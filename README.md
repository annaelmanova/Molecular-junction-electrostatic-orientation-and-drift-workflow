# Molecular junction electrostatic orientation and drift workflow

Reproducible GUI and single-notebook workflow for one-dimensional molecular-junction electrostatics, isotropic force/drift analysis, and exponential-field anisotropic orientation calculations. Defaults reproduce the Y6/ZnO/SAM/PM6:Y6 system.

The molecule name and the left-layer, interlayer, and active-layer names are editable in the GUI or `CFG` cell. They propagate to device-strip labels, run metadata, and report text. Molecular figures are not bundled assets: enable them explicitly and provide an XYZ geometry plus a matching OpenDX electrostatic-potential file. The notebook then plots the input geometry and selected calculated minimum-energy orientations.

## What the workflow produces

- Exponential and parabolic potential/field comparisons.
- Analytical electric-field derivatives.
- Isotropic dipolar, polarizability, quadrupolar, and total forces.
- Isotropic drift-length sweeps for user-selected diffusion coefficients and times.
- Full-sphere molecular orientation-energy minimization and three distinct minima.
- Exponential-only anisotropic forces and drift.
- PNG and vector-PDF figures with fixed publication panel dimensions.
- CSV data, a configured notebook, a Word report, and a PDF report.

All numerical plots are regenerated from the GUI inputs. Anisotropic calculations intentionally use only the exponential electrostatic model.

## Repository contents

| Path | Purpose |
|---|---|
| `Molecular_Junction_GUI.py` | Desktop GUI for inputs, QC text, execution, and output generation. |
| `Molecular_Junction_Workflow.ipynb` | Self-contained calculation notebook. |
| `build_molecular_junction_notebook.py` | Rebuilds the notebook from source cells and embedded assets. |
| `build_molecular_junction_report.py` | Creates the structured Word report from a completed run. |
| `repair_molecular_junction_run.py` | Recreates a configured notebook from a saved run configuration. |
| `assets/` | System-independent workflow and junction schematics. |
| `examples/sample_user_inputs.json` | Example parameter configuration. |

## Installation

Python 3.10 or newer is recommended.

```powershell
git clone <repository-url>
cd Molecular_Junction_Workflow
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Linux/macOS, activate the environment with `source .venv/bin/activate`.

## Running the GUI

```powershell
python Molecular_Junction_GUI.py
```

The GUI provides tabs for device parameters, electrostatic sweeps, diffusion/time values, the ground-state dipole, polarizability tensor, optional quadrupole tensor, QC output, and instructions.

1. Enter or paste the QC information.
2. Review the molecular dipole and tensors in the same Cartesian frame as the XYZ geometry.
3. Set the voltage, interlayer-step, layer decay-length, diffusion, and time sweeps.
4. Select **Validate inputs**.
5. Select **Create and run notebook**.
6. Use **Open output folder** after completion.

The GUI normally uses the Python interpreter that launched it. To select another interpreter, define `MOLECULAR_WORKFLOW_PYTHON` before starting the GUI.

```powershell
$env:MOLECULAR_WORKFLOW_PYTHON = "C:\path\to\python.exe"
python Molecular_Junction_GUI.py
```

## Running the notebook directly

Open `Molecular_Junction_Workflow.ipynb`, edit the `CFG` cell, and run all cells from top to bottom. The notebook creates its output directory automatically.

## Output layout

Each GUI run is stored under `Molecular_Junction_GUI_runs/<timestamp>/`:

```text
Molecular_Junction_GUI_runs/<timestamp>/
├── Configured_Molecular_Junction_Workflow.ipynb
├── user_inputs.json
├── QC_input.txt
└── outputs/
    ├── Molecular_Junction_Report.docx
    ├── Molecular_Junction_Report.pdf
    ├── UPDATED_PUBLICATION_PDF_PLOTS/
    ├── figures/
    └── data/
```

`UPDATED_PUBLICATION_PDF_PLOTS` contains the clean vector-PDF figure set. The `figures` directory contains both PNG and PDF versions. Numerical tables and orientation minima are written to `data` as CSV files.

## Important model conventions

- Dipole, polarizability tensor, quadrupole tensor, and XYZ geometry must use the same Cartesian coordinate system.
- The dipole is the ground-state dipole moment, not a transition dipole moment.
- Isotropic quantities use `|mu|`, `Tr(alpha)/3`, and `Tr(Q)/3`.
- Orientation minimization sweeps the supplied orientation sphere at every spatial coordinate.
- Anisotropic forces and drift are evaluated only with the exponential electrostatic closure.
- Drift lengths are local Einstein-relation displacement scales, not full morphology trajectories.

## Word and PDF export

The Word report requires `python-docx`. On Windows, the GUI exports PDF through Microsoft Word automation. If Word is unavailable, the DOCX and all figure PDFs are still produced; export the DOCX manually or add a LibreOffice-based conversion step.

## Rebuilding the notebook

After editing `build_molecular_junction_notebook.py` or files in `assets/`:

```powershell
python build_molecular_junction_notebook.py
```

## Reproducibility

Archive the configured notebook, `user_inputs.json`, QC input, XYZ geometry, CSV data, and report together. Record the Python/package versions and the orientation-sphere sample count used for publication calculations.

## License and citation

No license has been selected in this package. Add an appropriate `LICENSE` file before public release. Update `CITATION.cff` with the author names, repository URL, and release DOI when available.
