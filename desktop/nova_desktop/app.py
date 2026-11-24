"""Qt desktop front end for the versioned NOVA thermal solver process."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from PySide6.QtCore import QIODevice, QProcess, QSaveFile, QSettings, Qt, Signal
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication, QFileDialog, QLabel, QMainWindow, QMessageBox, QPlainTextEdit,
    QProgressBar, QSplitter, QToolBar, QVBoxLayout, QWidget,
)

from .editor import ModelEditor
from .field import FieldView, NodalField


def find_runner(explicit):
    if explicit:
        return str(Path(explicit).resolve())
    saved = QSettings().value("solver/path", "")
    name = "nova-run.exe" if os.name == "nt" else "nova-run"
    for candidate in (Path(sys.executable).with_name(name), Path(sys.argv[0]).resolve().with_name(name)):
        if candidate.is_file():
            return str(candidate)
    return saved if saved and Path(saved).is_file() else shutil.which(name)


class Window(QMainWindow):
    study_finished = Signal(bool)

    def __init__(self, runner=None):
        super().__init__()
        self.runner = find_runner(runner)
        self.workspace = tempfile.TemporaryDirectory(prefix="nova-thermal-")
        self.snapshot = Path(self.workspace.name) / "snapshot.nova"
        self.vtk = Path(self.workspace.name) / "solution.vtk"
        self.path = None
        self.dirty = False
        self.loading = False
        self.job = None
        self.cancelled = False
        self.editor = ModelEditor()
        self.view = FieldView()
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.console = QPlainTextEdit()
        self.console.setReadOnly(True)
        self.console.setAccessibleName("Solver diagnostics")
        self.console.setMaximumBlockCount(80)
        self.console.setMaximumHeight(95)
        result = QWidget()
        layout = QVBoxLayout(result)
        layout.addWidget(self.view, 1)
        layout.addWidget(self.summary)
        layout.addWidget(self.console)
        split = QSplitter(Qt.Orientation.Horizontal)
        split.addWidget(self.editor)
        split.addWidget(result)
        split.setSizes([410, 650])
        self.setCentralWidget(split)
        self.resize(1080, 690)
        self.activity = QProgressBar()
        self.activity.setRange(0, 0)
        self.activity.setMaximumWidth(110)
        self.statusBar().addPermanentWidget(self.activity)
        self.activity.hide()

        files = self.menuBar().addMenu("&File")
        study = self.menuBar().addMenu("&Study")
        toolbar = QToolBar("Project and study")
        toolbar.setMovable(False)
        self.addToolBar(toolbar)

        def action(menu, title, callback, shortcut=None, on_toolbar=False):
            item = QAction(title, self)
            if shortcut:
                item.setShortcut(QKeySequence(shortcut))
            item.triggered.connect(callback)
            menu.addAction(item)
            if on_toolbar:
                toolbar.addAction(item)
            return item

        self.new_action = action(files, "New", self.new_project, QKeySequence.StandardKey.New, True)
        self.open_action = action(files, "Open…", self.open_dialog, QKeySequence.StandardKey.Open, True)
        self.save_action = action(files, "Save", self.save, QKeySequence.StandardKey.Save, True)
        self.save_as_action = action(files, "Save as…", self.save_as, QKeySequence.StandardKey.SaveAs)
        self.export_action = action(files, "Export result VTK…", self.export_dialog)
        toolbar.addSeparator()
        self.run_action = action(study, "Solve", self.start_solve, "F5", True)
        self.cancel_action = action(study, "Cancel solve", self.cancel, "Esc", True)
        self.runner_action = action(study, "Locate solver…", self.choose_runner)
        self.editor.changed.connect(self._edited)
        self._controls()
        self.statusBar().showMessage("Ready" if self.runner else "Choose the installed nova-run under Study → Locate solver.")

    def _controls(self):
        idle = self.job is None
        self.editor.setEnabled(idle)
        for action in (self.new_action, self.open_action, self.save_action, self.save_as_action,
                       self.run_action, self.runner_action):
            action.setEnabled(idle)
        self.cancel_action.setEnabled(not idle and not self.cancelled)
        self.export_action.setEnabled(idle and self.view.canvas.field is not None)
        self.activity.setVisible(not idle)
        name = self.path.name if self.path else "Untitled"
        self.setWindowTitle(f"{name}{' *' if self.dirty else ''} — NOVA Thermal")

    def _clear_result(self):
        self.view.set_field(None)
        self.summary.clear()

    def _edited(self):
        if not self.loading:
            self.dirty = True
            self._clear_result()
            self._controls()
            self.statusBar().showMessage("Model edited")

    def _error(self, error):
        message = str(error).strip()
        self.console.appendPlainText(message)
        self.statusBar().showMessage(message)

    def choose_runner(self):
        path, _ = QFileDialog.getOpenFileName(self, "Locate the installed nova-run")
        if path:
            self.runner = path
            QSettings().setValue("solver/path", path)
            self.statusBar().showMessage("Solver selected")
        return bool(path)

    def _runner_ready(self):
        return bool(self.runner) or self.choose_runner()

    def _invoke(self, *args):
        result = subprocess.run([self.runner, *map(str, args)], capture_output=True,
                                encoding="utf-8", errors="replace", timeout=15,
                                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        if result.returncode:
            raise ValueError(result.stderr.strip() or f"Solver exited with code {result.returncode}.")
        return result.stdout

    def _discard_or_save(self):
        if not self.dirty:
            return True
        buttons = QMessageBox.StandardButton
        choice = QMessageBox.question(self, "Unsaved project", "Save the current project?",
                                      buttons.Save | buttons.Discard | buttons.Cancel, buttons.Save)
        return choice == buttons.Discard or (choice == buttons.Save and self.save())

    def new_project(self):
        if not self._discard_or_save():
            return
        self.loading = True
        try:
            self.editor.reset_model()
        finally:
            self.loading = False
        self.path = None
        self.dirty = False
        self._clear_result()
        self.console.clear()
        self._controls()
        self.statusBar().showMessage("New thermal project")

    def open_dialog(self):
        if self._discard_or_save():
            path, _ = QFileDialog.getOpenFileName(self, "Open thermal project", "", "NOVA projects (*.nova)")
            if path:
                self.open_path(Path(path))

    def open_path(self, path):
        if not self._runner_ready():
            return False
        try:
            canonical = self._invoke("show", path)
            self.loading = True
            self.editor.load_canonical(canonical)
        except (ValueError, IndexError, OSError, subprocess.TimeoutExpired) as error:
            self._error(error)
            return False
        finally:
            self.loading = False
        self.path = path
        self.dirty = False
        self._clear_result()
        self.console.clear()
        self._controls()
        self.statusBar().showMessage("Opened project")
        return True

    def save(self):
        return self.save_to(self.path) if self.path else self.save_as()

    def save_as(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save thermal project", str(self.path or "untitled.nova"),
                                               "NOVA projects (*.nova)")
        if not path:
            return False
        return self.save_to(Path(path if path.endswith(".nova") else path + ".nova"))

    def save_to(self, path):
        if not self._runner_ready():
            return False
        try:
            self.snapshot.write_text(self.editor.text(), encoding="ascii")
            self._invoke("save-as", self.snapshot, path)
        except (ValueError, OSError, subprocess.TimeoutExpired) as error:
            self._error(error)
            return False
        self.path = path
        self.dirty = False
        self._controls()
        self.statusBar().showMessage("Saved project")
        return True

    def start_solve(self):
        if self.job is not None or not self._runner_ready():
            return
        try:
            self.snapshot.write_text(self.editor.text(), encoding="ascii")
        except (OSError, ValueError) as error:
            self._error(error)
            return
        self._clear_result()
        self.console.clear()
        self.cancelled = False
        job = QProcess(self)
        self.job = job
        job.setProgram(self.runner)
        job.setArguments(["solve", str(self.snapshot), "--vtk", str(self.vtk)])
        job.finished.connect(self._finished)
        job.errorOccurred.connect(self._process_error)
        self._controls()
        self.statusBar().showMessage("Solving stationary heat conduction…")
        job.start()

    def cancel(self):
        if self.job is not None:
            self.cancelled = True
            self._controls()
            self.statusBar().showMessage("Cancelling solve…")
            self.job.kill()

    def _process_error(self, error):
        if error == QProcess.ProcessError.FailedToStart and self.job is not None:
            self._error(self.job.errorString())
            self._release_job(False)

    def _finished(self, code, exit_status):
        if self.job is None:
            return
        ok = False
        try:
            if self.cancelled:
                self.statusBar().showMessage("Solve cancelled")
            elif code != 0 or exit_status != QProcess.ExitStatus.NormalExit:
                diagnostic = bytes(self.job.readAllStandardError()).decode("utf-8", "replace").strip()
                raise ValueError(diagnostic or f"Solver exited with code {code}.")
            else:
                output = bytes(self.job.readAllStandardOutput()).decode("ascii").splitlines()
                if not output or output[0] != "NOVA_RESULT 1":
                    raise ValueError("Unsupported solver response.")
                data = dict(row.split() for row in output[1:])
                if (data["project"] != self.editor.fields["project"].text().strip() or
                        data["study"] != self.editor.fields["study"].text().strip()):
                    raise ValueError("Result belongs to a different project or study.")
                divisions = tuple(int(self.editor.fields[key].text()) for key in ("nx", "ny", "nz"))
                field = NodalField.read(self.vtk, divisions)
                if len(field.values) != int(data["nodes"]):
                    raise ValueError("Result and solver summary node counts differ.")
                summary = (f"{data['nodes']} nodes · {data['tetrahedra']} tetrahedra · {data['iterations']} CG iterations\n"
                           f"Temperature {float(data['min_K']):.8g} – {float(data['max_K']):.8g} K    "
                           f"Convective power {float(data['convection_out_W']):.8g} W    "
                           f"Heat balance defect {float(data['balance_defect_W']):.4g} W")
                self.view.set_field(field)
                self.summary.setText(summary)
                self.console.setPlainText("\n".join(output))
                self.statusBar().showMessage("Study complete")
                ok = True
        except (OSError, ValueError, KeyError, IndexError) as error:
            self._error(error)
        self._release_job(ok)

    def _release_job(self, ok):
        self.job.deleteLater()
        self.job = None
        self._controls()
        self.study_finished.emit(ok)

    def export_dialog(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export result", "temperature.vtk", "VTK result (*.vtk)")
        if path:
            self.export_to(Path(path if path.endswith(".vtk") else path + ".vtk"))

    def export_to(self, path):
        if self.view.canvas.field is None:
            return False
        try:
            if path.is_symlink() or (path.exists() and not path.is_file()):
                raise ValueError("Result destination must be a regular file.")
            if self.path and self.path.exists() and path.exists() and path.samefile(self.path):
                raise ValueError("Result destination is the open project.")
            target = QSaveFile(str(path))
            if not target.open(QIODevice.OpenModeFlag.WriteOnly):
                raise OSError(target.errorString())
            with self.vtk.open("rb") as src:
                while chunk := src.read(65536):
                    if target.write(chunk) != len(chunk):
                        raise OSError(target.errorString())
            if not target.commit():
                raise OSError(target.errorString())
        except (ValueError, OSError) as error:
            self._error(error)
            return False
        self.statusBar().showMessage("Exported VTK result")
        return True

    def closeEvent(self, event):
        if self.job is not None:
            job = self.job
            self.cancel()
            job.waitForFinished(3000)
            if self.job is not None:
                event.ignore()
                return
        if not self._discard_or_save():
            event.ignore()
            return
        self.workspace.cleanup()
        event.accept()


def main():
    parser = argparse.ArgumentParser(description="NOVA Thermal desktop project editor")
    parser.add_argument("project", nargs="?", type=Path)
    parser.add_argument("--runner", help="Path to the installed nova-run executable")
    args = parser.parse_args()
    app = QApplication(sys.argv[:1])
    app.setOrganizationName("NOVA Multiphysics")
    app.setApplicationName("NOVA Thermal")
    window = Window(args.runner)
    window.show()
    if args.project:
        window.open_path(args.project)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
