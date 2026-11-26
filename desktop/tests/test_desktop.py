"""Exercise the installed Qt desktop through real widgets and the installed solver."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEventLoop, QPoint, Qt, QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog, QToolBar

from nova_desktop.app import Window
from nova_desktop.field import NodalField


class DesktopWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])
        cls.app.setQuitOnLastWindowClosed(False)

    def setUp(self):
        self.files = tempfile.TemporaryDirectory()
        self.directory = Path(self.files.name)
        self.win = Window(RUNNER)
        self.win.show()
        self.app.processEvents()

    def tearDown(self):
        self.win.dirty = False
        self.win.close()
        self.files.cleanup()

    def click_action(self, action):
        button = self.win.findChild(QToolBar).widgetForAction(action)
        QTest.mouseClick(button, Qt.MouseButton.LeftButton)

    def wait_for_study(self, start):
        loop = QEventLoop()
        outcomes = []

        def done(ok):
            outcomes.append(ok)
            loop.quit()

        self.win.study_finished.connect(done)
        timer = QTimer()
        timer.setSingleShot(True)
        timer.timeout.connect(loop.quit)
        timer.start(45000)
        try:
            start()
            if not outcomes:
                loop.exec()
            self.assertEqual(len(outcomes), 1, "Timed out waiting for the solver")
            return outcomes[0]
        finally:
            timer.stop()
            self.win.study_finished.disconnect(done)

    def test_create_save_reopen_solve_probe_and_export(self):
        project = self.directory / "قياس thermal model.nova"
        exported = self.directory / "exported.vtk"
        field = self.win.editor.fields["project"]
        field.selectAll()
        QTest.keyClicks(field, "thermal-wall")
        with patch.object(QFileDialog, "getSaveFileName", return_value=(str(project), "")):
            self.click_action(self.win.save_action)
        self.assertTrue(project.is_file(), self.win.console.toPlainText())
        self.assertFalse(self.win.dirty)
        self.click_action(self.win.new_action)
        self.assertEqual(field.text(), "untitled")
        with patch.object(QFileDialog, "getOpenFileName", return_value=(str(project), "")):
            self.click_action(self.win.open_action)
        self.assertEqual(field.text(), "thermal-wall")
        self.assertTrue(self.wait_for_study(lambda: self.click_action(self.win.run_action)),
                        self.win.console.toPlainText())
        self.assertIn("Convective power 75 W", self.win.summary.text())
        view = self.win.view
        self.assertEqual(len(view.canvas.field.values), 108)
        for normal in range(3):
            view.plane.setCurrentIndex(normal)
            view.slider.setValue(view.slider.maximum())
            index = [0, 0, 0]
            index[normal] = view.slider.value()
            # Click the center of an actual node sample in the displayed axes.
            index[view.canvas.u] = (view.canvas.w - 1) // 2
            index[view.canvas.v] = (view.canvas.h - 1) // 2
            rect = view.canvas.plot_rect()
            point = QPoint(round(rect.left() + rect.width() * index[view.canvas.u] / (view.canvas.w - 1)),
                           round(rect.bottom() - rect.height() * index[view.canvas.v] / (view.canvas.h - 1)))
            QTest.mouseClick(view.canvas, Qt.MouseButton.LeftButton, pos=point)
            actual = float(view.probe.text().split()[-2])
            x = view.canvas.field.axes[0][index[0]]
            self.assertAlmostEqual(actual, 400 - 150 * x, delta=2e-7)
        with patch.object(QFileDialog, "getSaveFileName", return_value=(str(exported), "")):
            self.win.export_action.trigger()
        self.assertEqual(exported.read_bytes(), self.win.vtk.read_bytes())
        self.win.editor.fields["nx"].setText("10")
        self.assertFalse(self.win.export_action.isEnabled())
        self.assertIsNone(self.win.view.canvas.field)
        self.assertTrue(self.win.dirty)

    def test_invalid_save_solve_and_open_preserve_project(self):
        path = self.directory / "saved.nova"
        self.assertTrue(self.win.save_to(path))
        expected = path.read_bytes()
        self.win.editor.fields["k"].setText("-1")
        self.assertFalse(self.win.save_to(path))
        self.assertEqual(expected, path.read_bytes())
        self.assertFalse(self.wait_for_study(self.win.start_solve))
        self.assertIn("conductivity", self.win.console.toPlainText().lower())
        self.assertFalse(self.win.export_action.isEnabled())
        damaged = self.directory / "damaged.nova"
        damaged.write_text("NOVA_PROJECT 9000\n", encoding="ascii")
        self.assertFalse(self.win.open_path(damaged))
        self.assertEqual(self.win.editor.fields["k"].text(), "-1")
        self.assertEqual(self.win.path, path)

    def test_schema_two_anisotropic_editor_round_trip_and_solve(self):
        editor = self.win.editor
        editor.fields["project"].setText("anisotropic-wall")
        editor.fields["geometry"].setText("wall")
        editor.fields["lx"].setText("0.5")
        editor.fields["ly"].setText("0.25")
        editor.fields["lz"].setText("0.125")
        editor.fields["mesh"].setText("tets")
        editor.fields["nx"].setText("6")
        editor.fields["ny"].setText("2")
        editor.fields["nz"].setText("1")
        editor.fields["material"].setText("benchmark-solid")
        editor.fields["study"].setText("steady")
        editor.fields["source"].setText("0")
        editor.material_mode.setCurrentIndex(editor.material_mode.findData("anisotropic"))
        for key, value in {"kxx": "2", "kyy": "3", "kzz": "4",
                           "kxy": "0.25", "kxz": "-0.15", "kyz": "0.10"}.items():
            editor.fields[key].setText(value)
        editor.set_boundary("xmin", "temperature", "400")
        editor.set_boundary("xmax", "convection", "12", "300")
        editor.set_boundary("ymin", "outward_flux", "-37.5")
        editor.set_boundary("ymax", "outward_flux", "37.5")
        editor.set_boundary("zmin", "outward_flux", "22.5")
        editor.set_boundary("zmax", "outward_flux", "-22.5")

        project = self.directory / "anisotropic.nova"
        self.assertTrue(self.win.save_to(project), self.win.console.toPlainText())
        saved = project.read_text(encoding="ascii")
        self.assertTrue(saved.startswith("NOVA_PROJECT 2\n"))
        material = next(line.split() for line in saved.splitlines()
                        if line.startswith("material "))
        self.assertEqual(material[:3], ["material", "benchmark-solid", "anisotropic"])
        self.assertEqual(len(material), 9)

        self.win.new_project()
        self.assertTrue(self.win.open_path(project), self.win.console.toPlainText())
        self.assertEqual(editor.material_mode.currentData(), "anisotropic")
        for key, value in {"kxx": 2.0, "kyy": 3.0, "kzz": 4.0,
                           "kxy": 0.25, "kxz": -0.15, "kyz": 0.10}.items():
            self.assertAlmostEqual(float(editor.fields[key].text()), value, delta=1e-14)

        self.assertTrue(self.wait_for_study(self.win.start_solve), self.win.console.toPlainText())
        self.assertIn("Convective power 9.375 W", self.win.summary.text())
        field = self.win.view.canvas.field
        self.assertEqual(len(field.values), 42)
        self.assertAlmostEqual(min(field.values), 325.0, delta=2e-8)
        self.assertAlmostEqual(max(field.values), 400.0, delta=2e-8)

    def test_cancel_then_solve_a_new_snapshot(self):
        for key in ("nx", "ny", "nz"):
            self.win.editor.fields[key].setText("32")

        def start_and_cancel():
            self.win.start_solve()
            self.assertTrue(self.win.job.waitForStarted(5000))
            self.assertFalse(self.win.editor.isEnabled())
            self.click_action(self.win.cancel_action)

        self.assertFalse(self.wait_for_study(start_and_cancel))
        self.assertEqual(self.win.statusBar().currentMessage(), "Solve cancelled")
        self.assertFalse(self.win.export_action.isEnabled())
        for key in ("nx", "ny", "nz"):
            self.win.editor.fields[key].setText("2")
        self.assertTrue(self.wait_for_study(self.win.start_solve))
        self.assertEqual(len(self.win.view.canvas.field.values), 27)

    def test_missing_solver_and_damaged_result_are_reported(self):
        self.win.runner = str(self.directory / "absent-solver")
        self.assertFalse(self.wait_for_study(self.win.start_solve))
        self.assertTrue(self.win.run_action.isEnabled())
        self.win.runner = RUNNER
        self.assertTrue(self.wait_for_study(self.win.start_solve))
        corrupted = self.directory / "truncated.vtk"
        corrupted.write_bytes(self.win.vtk.read_bytes()[:150])
        with self.assertRaises(ValueError):
            NodalField.read(corrupted, (8, 3, 2))

        records = self.win.vtk.read_text(encoding="ascii").splitlines()
        cells = 6 * 8 * 3 * 2
        first_cell = records.index(f"CELLS {cells} {5 * cells}") + 1
        for name, bad_record in (
            ("duplicate-node", "4 0 0 1 2"),
            ("out-of-range-node", "4 0 1 2 108"),
            ("nonlocal-tetrahedron", "4 0 1 2 107"),
            # Four corners with x+y=1 are coplanar but span the cell in x, y and z.
            ("invalid-local-tetrahedron", "4 1 9 37 45"),
        ):
            damaged_records = records.copy()
            damaged_records[first_cell] = bad_record
            damaged = self.directory / f"{name}.vtk"
            damaged.write_text("\n".join(damaged_records) + "\n", encoding="ascii")
            with self.subTest(connectivity=name), self.assertRaises(ValueError):
                NodalField.read(damaged, (8, 3, 2))

        duplicate_records = records.copy()
        duplicate_records[first_cell + 1] = duplicate_records[first_cell]
        duplicate = self.directory / "duplicate-tetrahedron.vtk"
        duplicate.write_text("\n".join(duplicate_records) + "\n", encoding="ascii")
        with self.assertRaisesRegex(ValueError, "Duplicate tetrahedron"):
            NodalField.read(duplicate, (8, 3, 2))

    def test_round_trip_all_settings_and_boundary_types(self):
        editor = self.win.editor
        for name, value in {"mesh": "mesh-two", "geometry": "new-box", "material": "sample",
                            "study": "run-two", "lx": "0.3", "nz": "3", "source": "0.125"}.items():
            editor.fields[name].setText(value)
        editor.set_boundary("ymin", "outward_flux", "-12.5")
        editor.set_boundary("zmax", "convection", "2", "290.25")
        editor.preconditioner.setCurrentText("ic0")
        path = self.directory / "settings.nova"
        self.assertTrue(self.win.save_to(path))
        original = path.read_bytes()
        self.win.new_project()
        self.assertTrue(self.win.open_path(path))
        self.assertEqual(editor.sides["ymin"][0].currentData(), "outward_flux")
        self.assertEqual(editor.sides["zmax"][2].text(), "290.25")
        self.assertEqual(editor.preconditioner.currentText(), "ic0")
        self.assertTrue(self.win.save_to(path))
        self.assertEqual(path.read_bytes(), original)
        self.assertTrue(self.wait_for_study(self.win.start_solve))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python test_desktop.py INSTALLED_SOLVER_BIN_DIR")
    name = "nova-run.exe" if os.name == "nt" else "nova-run"
    RUNNER = str((Path(sys.argv.pop(1)) / name).resolve(strict=True))
    unittest.main(verbosity=2)
