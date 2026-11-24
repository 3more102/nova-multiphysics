# NOVA Thermal desktop

The desktop edits native schema-1 thermal box projects. It uses Qt 6 for the controls and the installed `nova-run` process for project validation, canonical persistence, mesh generation, solving and VTK creation.

## Development installation

Build and install the NOVA C++ SDK from this revision. In a Python 3.10+ virtual environment, install the optional desktop client. The regression configuration uses Python 3.12 and PySide6 Essentials 6.8.3.

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
cmake --install build --config Release --prefix ./nova-sdk
python -m pip install ./desktop "PySide6-Essentials==6.8.3"
nova-desktop
~~~

The desktop launcher can also be opened from the environment's executable directory. In **Study → Locate solver**, select `nova-sdk/bin/nova-run` on Linux or `nova-sdk/bin/nova-run.exe` on Windows. The client remembers this path. For automation, pass `--runner PATH` and an optional project path to the launcher.

## Create, solve and inspect a study

1. Choose **New**. Enter a project ID under **Project**, box dimensions in metres under **Geometry**, and grid subdivisions under **Mesh**.
2. Set the isotropic conductivity under **Material**. In **Boundaries**, choose a condition per box side and enter the displayed SI values. **Study** controls heat generation, iteration budget, tolerances and preconditioner.
3. Choose **Save** to keep a native project. **Open** restores all supported settings by asking the solver to read the file. Saving uses the C++ canonical writer.
4. Choose **Solve** (F5). The desktop submits an immutable working copy and displays progress while the process runs. **Cancel solve** (Esc) stops the job. Errors appear in the diagnostics panel. Editing the model clears the previous field.
5. Choose XY, XZ or YZ and use the slider to select a nodal plane. Click a colored sample to read its actual coordinates in metres and nodal temperature in kelvin. The scale covers the entire computed field. The display uses nearest nodal samples. Cell heat flux and generation are available through **File → Export result VTK**.
6. For the shipped wall, the cold face is 325 K and the convective heat rate is 75 W. The status panel also reports the mesh size, CG iteration count and signed heat balance. [Model equations and supported scope](THERMAL_PROJECT.md).

**File formats and scope.** The client reads the canonical response of `nova-run show PROJECT` and accepts solver response protocol `NOVA_RESULT 1`. It reads temperatures from the solver's private VTK file and checks its coordinates against the structured grid. Save, open and solve use argument arrays with no shell. The project file contains no executable commands or external output paths. The viewer supports structured box slices, with the full VTK available for external visualization.

The desktop is a separately installed development client. Free-form CAD editing, arbitrary mesh viewing, resumable jobs, recovery after a machine failure, signed installers and user acceptance testing remain open. Qt dependency provenance, license choice and distribution obligations require the owner's review before bundling. [Qt for Python commercial use](https://doc.qt.io/qtforpython-6/commercial/index.html).

## Installed-client verification

~~~sh
python desktop/tests/test_desktop.py ./nova-sdk/bin
~~~

The tests run Qt with the offscreen platform plugin. They exercise actual widget actions, native save/open (including a Unicode path), a separate solver process, all three slice directions, analytical temperature probes, VTK export, invalid input, cancellation, a subsequent solve and failure to launch the solver. They also verify that edits invalidate the result and that a failed save preserves the original project. Hardware-rendered visualization and interactive user usability are separate qualification work.
