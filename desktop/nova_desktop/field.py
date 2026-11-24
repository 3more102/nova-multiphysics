"""Read the nodal field produced by nova-run and display selectable grid planes."""
from array import array
from dataclasses import dataclass
import math

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QImage, QLinearGradient, QPainter
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QSlider, QVBoxLayout, QWidget


@dataclass
class NodalField:
    axes: list
    values: array

    def value(self, index):
        x, y, z = index
        return self.values[x + len(self.axes[0]) * (y + len(self.axes[1]) * z)]

    @classmethod
    def read(cls, path, divisions):
        shape = tuple(n + 1 for n in divisions)
        cells = 6 * math.prod(divisions)
        if min(divisions) < 1 or cells > 196608:
            raise ValueError("Unsupported result dimensions.")
        with path.open(encoding="ascii") as src:
            def line():
                data = src.readline(512)
                if not data or len(data) == 512:
                    raise ValueError("Truncated or invalid solver field record.")
                return data.strip()

            def expect(value):
                if line() != value:
                    raise ValueError(f"Expected solver field record {value!r}.")

            expect("# vtk DataFile Version 3.0")
            line()  # Dataset title.
            expect("ASCII")
            expect("DATASET UNSTRUCTURED_GRID")
            count = math.prod(shape)
            expect(f"POINTS {count} double")
            axes = [array("d"), array("d"), array("d")]
            for node in range(count):
                coords = tuple(map(float, line().split()))
                if len(coords) != 3 or not all(map(math.isfinite, coords)):
                    raise ValueError("Invalid mesh coordinate.")
                index = (node % shape[0], node // shape[0] % shape[1], node // (shape[0] * shape[1]))
                for axis, i, coordinate in zip(axes, index, coords):
                    if i == len(axis):
                        if axis and coordinate <= axis[-1]:
                            raise ValueError("A structured mesh axis is not increasing.")
                        axis.append(coordinate)
                    elif axis[i] != coordinate:
                        raise ValueError("Result coordinates disagree with the structured mesh.")
            expect(f"CELLS {cells} {5 * cells}")
            referenced = bytearray(count)
            # Canonical corner IDs use x + 2*y + 4*z. These are the six
            # body-diagonal tetrahedra emitted by make_box_mesh_3d().
            expected_tetrahedra = {
                (0, 1, 3, 7): 0,
                (0, 2, 3, 7): 1,
                (0, 2, 6, 7): 2,
                (0, 4, 6, 7): 3,
                (0, 4, 5, 7): 4,
                (0, 1, 5, 7): 5,
            }
            cube_masks = bytearray(math.prod(divisions))
            for _ in range(cells):
                words = line().split()
                if len(words) != 5 or words[0] != "4":
                    raise ValueError("Invalid tetrahedral connectivity record.")
                try:
                    nodes = tuple(int(word, 10) for word in words[1:])
                except ValueError as error:
                    raise ValueError("Invalid tetrahedral node index.") from error
                if len(set(nodes)) != 4 or any(node < 0 or node >= count for node in nodes):
                    raise ValueError("Tetrahedron contains duplicate or out-of-range node indices.")
                grid = tuple(
                    (node % shape[0], node // shape[0] % shape[1],
                     node // (shape[0] * shape[1]))
                    for node in nodes
                )
                lower = tuple(min(point[axis] for point in grid) for axis in range(3))
                if any(max(point[axis] for point in grid) - lower[axis] != 1
                       for axis in range(3)):
                    raise ValueError("Tetrahedron does not occupy one structured mesh cell.")
                local_tetrahedron = tuple(sorted(
                    (point[0] - lower[0]) + 2 * (point[1] - lower[1]) +
                    4 * (point[2] - lower[2]) for point in grid
                ))
                pattern = expected_tetrahedra.get(local_tetrahedron)
                if pattern is None:
                    raise ValueError("Tetrahedron does not match NOVA's six-tetrahedron cell split.")
                cube = lower[0] + divisions[0] * (lower[1] + divisions[1] * lower[2])
                bit = 1 << pattern
                if cube_masks[cube] & bit:
                    raise ValueError("Duplicate tetrahedron in structured mesh cell.")
                cube_masks[cube] |= bit
                for node in nodes:
                    referenced[node] = 1
            if not all(referenced) or any(mask != 0b111111 for mask in cube_masks):
                raise ValueError("Tetrahedral connectivity does not cover the native structured mesh.")
            expect(f"CELL_TYPES {cells}")
            for _ in range(cells):
                expect("10")
            expect(f"POINT_DATA {count}")
            expect("SCALARS temperature double 1")
            expect("LOOKUP_TABLE default")
            values = array("d", (float(line()) for _ in range(count)))
            if not all(math.isfinite(t) and t >= 0 for t in values):
                raise ValueError("Invalid absolute temperature in result.")
        return cls(axes, values)


STOPS = ((24, 47, 86), (30, 105, 141), (43, 158, 151), (159, 203, 106), (253, 224, 86))


def colors():
    palette = []
    for index in range(256):
        position = index / 255 * (len(STOPS) - 1)
        left = min(int(position), len(STOPS) - 2)
        f = position - left
        palette.append(QColor(*(round(a + f * (b - a)) for a, b in zip(STOPS[left], STOPS[left + 1]))))
    return palette


class SliceCanvas(QWidget):
    probed = Signal(str)

    def __init__(self):
        super().__init__()
        self.field = None
        self.palette_colors = colors()
        self.image = QImage()
        self.normal = 2
        self.level = 0
        self.selected = None
        self.setMinimumSize(260, 230)
        self.setAccessibleName("Temperature field. Click a node to probe kelvin and position.")

    def set_field(self, field):
        self.field = field
        self.selected = None
        if field:
            self.low, self.high = min(field.values), max(field.values)
        self.update()

    def set_plane(self, normal, level):
        self.normal, self.level = normal, level
        self.selected = None
        if self.field:
            self.u, self.v = [axis for axis in range(3) if axis != normal]
            self.w = len(self.field.axes[self.u])
            self.h = len(self.field.axes[self.v])
            self.image = QImage(self.w, self.h, QImage.Format.Format_RGB32)
            index = [0, 0, 0]
            index[normal] = level
            span = self.high - self.low
            for j in range(self.h):
                index[self.v] = j
                for i in range(self.w):
                    index[self.u] = i
                    fraction = (self.field.value(index) - self.low) / span if span else 0.5
                    self.image.setPixelColor(i, self.h - 1 - j, self.palette_colors[round(255 * fraction)])
        self.update()

    def plot_rect(self):
        dx, dy = (self.field.axes[a][-1] - self.field.axes[a][0] for a in (self.u, self.v))
        scale = min((self.width() - 90) / dx, (self.height() - 110) / dy)
        w, h = dx * scale, dy * scale
        return QRectF((self.width() - w) / 2, (self.height() - h) / 2 - 12, w, h)

    def paintEvent(self, _event):
        painter = QPainter(self)
        painter.fillRect(self.rect(), self.palette().base())
        painter.setPen(self.palette().text().color())
        if not self.field or self.image.isNull():
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Solve the study to inspect temperature.")
            return
        rect = self.plot_rect()
        painter.drawImage(rect, self.image, QRectF(0.5, 0.5, self.w - 1, self.h - 1))
        painter.drawRect(rect)
        painter.drawText(QPointF(rect.left(), rect.bottom() + 20),
                         f"{'XYZ'[self.u]}   0 – {self.field.axes[self.u][-1]:.6g} m")
        painter.drawText(QPointF(rect.left(), rect.top() - 9),
                         f"{'XYZ'[self.v]}   0 – {self.field.axes[self.v][-1]:.6g} m")
        legend = QRectF(45, self.height() - 28, self.width() - 90, 8)
        gradient = QLinearGradient(legend.topLeft(), legend.topRight())
        for i, rgb in enumerate(STOPS):
            gradient.setColorAt(i / (len(STOPS) - 1), QColor(*rgb))
        painter.fillRect(legend, gradient)
        painter.drawText(QRectF(legend.left(), legend.bottom(), legend.width(), 17),
                         Qt.AlignmentFlag.AlignLeft, f"{self.low:.7g} K")
        painter.drawText(QRectF(legend.left(), legend.bottom(), legend.width(), 17),
                         Qt.AlignmentFlag.AlignRight, f"{self.high:.7g} K")
        if self.selected:
            i, j = self.selected
            point = QPointF(rect.left() + rect.width() * i / (self.w - 1),
                            rect.bottom() - rect.height() * j / (self.h - 1))
            painter.setPen(Qt.GlobalColor.white)
            painter.drawEllipse(point, 5, 5)
            painter.setPen(Qt.GlobalColor.black)
            painter.drawEllipse(point, 6, 6)

    def mousePressEvent(self, event):
        if not self.field or event.button() != Qt.MouseButton.LeftButton:
            return
        rect = self.plot_rect()
        point = event.position()
        if not rect.contains(point):
            return
        i = round((point.x() - rect.left()) * (self.w - 1) / rect.width())
        j = round((rect.bottom() - point.y()) * (self.h - 1) / rect.height())
        index = [0, 0, 0]
        index[self.normal], index[self.u], index[self.v] = self.level, i, j
        xyz = ", ".join(f"{self.field.axes[axis][node]:.7g}" for axis, node in enumerate(index))
        self.selected = i, j
        self.probed.emit(f"({xyz}) m    {self.field.value(index):.10g} K")
        self.update()


class FieldView(QWidget):
    def __init__(self):
        super().__init__()
        self.canvas = SliceCanvas()
        self.plane = QComboBox()
        self.plane.addItems(["YZ", "XZ", "XY"])
        self.plane.setCurrentIndex(2)
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.position = QLabel()
        self.probe = QLabel()
        controls = QHBoxLayout()
        controls.addWidget(QLabel("Slice"))
        controls.addWidget(self.plane)
        controls.addWidget(self.slider, 1)
        controls.addWidget(self.position)
        layout = QVBoxLayout(self)
        layout.addLayout(controls)
        layout.addWidget(self.canvas, 1)
        layout.addWidget(self.probe)
        self.canvas.probed.connect(self.probe.setText)
        self.plane.currentIndexChanged.connect(self._plane_changed)
        self.slider.valueChanged.connect(self._slice_changed)
        self.set_field(None)

    def set_field(self, field):
        self.canvas.set_field(field)
        self.plane.setEnabled(field is not None)
        self.slider.setEnabled(field is not None)
        self._plane_changed()

    def _plane_changed(self):
        field = self.canvas.field
        if field:
            count = len(field.axes[self.plane.currentIndex()])
            self.slider.blockSignals(True)
            self.slider.setRange(0, count - 1)
            self.slider.setValue((count - 1) // 2)
            self.slider.blockSignals(False)
        self._slice_changed()

    def _slice_changed(self):
        axis, index = self.plane.currentIndex(), self.slider.value()
        self.canvas.set_plane(axis, index)
        field = self.canvas.field
        self.position.setText(f"{'XYZ'[axis]} = {field.axes[axis][index]:.7g} m" if field else "")
        self.probe.setText("Click a node to read temperature." if field else "")
