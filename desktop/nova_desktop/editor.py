"""Editors for the records accepted by the schema-1 C++ model codec."""
from PySide6.QtCore import Signal, Qt
from PySide6.QtWidgets import (
    QComboBox, QFormLayout, QGridLayout, QLabel, QLineEdit, QStackedWidget,
    QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

SIDES = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")
BOUNDARY_TYPES = {"insulated": "Insulated", "temperature": "Fixed temperature (K)",
                  "outward_flux": "Outward flux (W/m²)", "convection": "Convection (W/m²·K)"}
PAGES = {
    "Project": [("project", "Project ID", "untitled")],
    "Geometry": [("geometry", "Box ID", "box"), ("lx", "Length X (m)", "0.5"),
                 ("ly", "Length Y (m)", "0.5"), ("lz", "Length Z (m)", "0.5")],
    "Mesh": [("mesh", "Mesh ID", "mesh"), ("nx", "X subdivisions", "8"),
             ("ny", "Y subdivisions", "3"), ("nz", "Z subdivisions", "2")],
    "Material": [("material", "Material ID", "solid"), ("k", "Conductivity (W/m·K)", "2")],
    "Boundaries": [],
    "Study": [("study", "Study ID", "heat"), ("source", "Generation (W/m³)", "0"),
              ("iterations", "Maximum iterations", "2000"), ("relative", "Relative tolerance", "1e-13"),
              ("absolute", "Absolute tolerance", "1e-14")],
}


class ModelEditor(QWidget):
    changed = Signal()

    def __init__(self):
        super().__init__()
        self.fields = {}
        self.sides = {}
        self.tree = QTreeWidget()
        self.tree.setHeaderLabel("Model builder")
        self.stack = QStackedWidget()
        layout = QVBoxLayout(self)
        layout.addWidget(self.tree, 1)
        layout.addWidget(self.stack, 2)
        for index, (page, entries) in enumerate(PAGES.items()):
            item = QTreeWidgetItem([page])
            item.setData(0, Qt.ItemDataRole.UserRole, index)
            self.tree.addTopLevelItem(item)
            pane = QWidget()
            form = QFormLayout(pane)
            for key, label, default in entries:
                field = QLineEdit(default)
                field.setMaxLength(100)
                field.textChanged.connect(self.changed)
                self.fields[key] = field
                form.addRow(label, field)
            if page == "Boundaries":
                table = QGridLayout()
                form.addRow(table)
                for column, label in enumerate(("Side", "Condition", "Value", "Ambient (K)")):
                    table.addWidget(QLabel(label), 0, column)
                for row, side in enumerate(SIDES, start=1):
                    kind, value, ambient = QComboBox(), QLineEdit(), QLineEdit()
                    for key, label in BOUNDARY_TYPES.items():
                        kind.addItem(label, key)
                    value.setMaxLength(100)
                    ambient.setMaxLength(100)
                    value.setAccessibleName(side + " boundary value")
                    ambient.setAccessibleName(side + " ambient kelvin")
                    self.sides[side] = kind, value, ambient
                    for column, widget in enumerate((QLabel(side), kind, value, ambient)):
                        table.addWidget(widget, row, column)
                    kind.currentIndexChanged.connect(lambda _index, s=side: self._boundary_changed(s))
                    value.textChanged.connect(self.changed)
                    ambient.textChanged.connect(self.changed)
            if page == "Study":
                self.preconditioner = QComboBox()
                self.preconditioner.addItems(["jacobi", "ic0", "none"])
                self.preconditioner.currentIndexChanged.connect(lambda _: self.changed.emit())
                form.addRow("CG preconditioner", self.preconditioner)
            self.stack.addWidget(pane)
        self.tree.currentItemChanged.connect(
            lambda item, _: self.stack.setCurrentIndex(item.data(0, Qt.ItemDataRole.UserRole)) if item else None)
        self.tree.setCurrentItem(self.tree.topLevelItem(0))
        self.reset_model()

    def _boundary_changed(self, side):
        kind, value, ambient = self.sides[side]
        value.setEnabled(kind.currentData() != "insulated")
        ambient.setEnabled(kind.currentData() == "convection")
        self.changed.emit()

    def set_boundary(self, side, kind, value="", ambient=""):
        combo, first, second = self.sides[side]
        combo.setCurrentIndex(combo.findData(kind))
        first.setText(value)
        second.setText(ambient)
        self._boundary_changed(side)

    def reset_model(self):
        for entries in PAGES.values():
            for key, _, default in entries:
                self.fields[key].setText(default)
        for side in SIDES:
            self.set_boundary(side, "insulated")
        self.set_boundary("xmin", "temperature", "400")
        self.set_boundary("xmax", "convection", "12", "300")
        self.preconditioner.setCurrentText("jacobi")

    def text(self):
        # Only single tokens are interpolated. Physical and reference validation
        # remains the responsibility of nova-run for both save and solve.
        def token(widget):
            value = widget.text().strip()
            if not value or not value.isascii() or "#" in value or any(c.isspace() for c in value):
                raise ValueError("Each active setting needs a single ASCII value.")
            return value

        f = {key: token(widget) for key, widget in self.fields.items()}
        lines = ["NOVA_PROJECT 1", f"project {f['project']}",
                 f"geometry {f['geometry']} box {f['lx']} {f['ly']} {f['lz']}",
                 f"mesh {f['mesh']} {f['geometry']} structured_tet {f['nx']} {f['ny']} {f['nz']}",
                 f"material {f['material']} isotropic {f['k']}",
                 f"study {f['study']} steady_heat {f['mesh']} {f['material']}",
                 f"source {f['study']} {f['source']}",
                 f"solver {f['study']} cg {self.preconditioner.currentText()} "
                 f"{f['iterations']} {f['relative']} {f['absolute']}"]
        for side, (combo, value, ambient) in self.sides.items():
            kind = combo.currentData()
            if kind != "insulated":
                line = f"boundary {f['study']} {f['geometry']} {side} {kind} {token(value)}"
                lines.append(line + (" " + token(ambient) if kind == "convection" else ""))
        return "\n".join([*lines, "end", ""])

    def load_canonical(self, text):
        # Input comes exclusively from the C++ runner's canonical 'show' response.
        rows = [line.split() for line in text.splitlines()]
        if not rows or rows[0] != ["NOVA_PROJECT", "1"] or rows[-1] != ["end"]:
            raise ValueError("The solver returned an unsupported project format.")
        project, geometry, mesh, material, study, source, solver = rows[1:8]
        values = {"project": project[1], "geometry": geometry[1], "mesh": mesh[1], "material": material[1],
                  "k": material[3], "study": study[1], "source": source[2], "iterations": solver[4],
                  "relative": solver[5], "absolute": solver[6]}
        values.update(zip(("lx", "ly", "lz"), geometry[3:6], strict=True))
        values.update(zip(("nx", "ny", "nz"), mesh[4:7], strict=True))
        selections = {row[3]: row[4:] for row in rows[8:-1]}
        for key, value in values.items():
            self.fields[key].setText(value)
        self.preconditioner.setCurrentText(solver[3])
        for side in SIDES:
            self.set_boundary(side, *selections.get(side, ["insulated"]))
