# FCC-ee W-pair polarization analysis

This repository contains Jupyter notebooks for angular, kinematic, likelihood, and TMVA studies of polarized W-pair samples. The analysis workflows are provided as **Jupyter notebooks only**; there are no standalone analysis scripts. The shared `analysis_config.py` module centralizes paths and analysis settings used by the notebooks.

## Repository contents

- `study_angles.ipynb` — W production and decay angles, lepton kinematics, and signal-versus-full-background comparisons. It also creates the ROOT templates used by the conventional likelihood fit.
- `likelihood_fit.ipynb` — one-, two-, and three-parameter likelihood fits and likelihood scans using the leptonic templates.
- `tmva_analysis.ipynb` — leptonic signal/background BDT training, test-sample plots, BDT templates, and a likelihood fit.
- `hadronic_analysis.ipynb` — hadronic-channel observables, BDT training, and likelihood analysis.
- `analysis_config.py` — shared input/output paths, sample filenames, tree names, cross-sections, luminosity, and binning. This is a configuration/helper module imported by the notebooks, not a standalone analysis workflow.
- `requirements.txt` — Python packages used by the notebooks.
- `data/README.md` — instructions for placing input ROOT files.

## Requirements

You need a working ROOT installation with PyROOT available in the Python interpreter used by Jupyter. The TMVA notebooks additionally require ROOT to be built with TMVA support. ROOT is an external prerequisite and is intentionally not installed through `requirements.txt`, because its installation depends on the operating system and Python environment.

Install the Python packages with:

```bash
python -m pip install -r requirements.txt
```

Before starting Jupyter, activate the environment in which this succeeds:

```bash
python -c "import ROOT; print(ROOT.gROOT.GetVersion())"
```

For TMVA, you can also check:

```bash
python -c "import ROOT; print(ROOT.TMVA.Tools.Instance())"
```

## Input files and configuration

Input ROOT files are not included in the repository. By default, notebooks look for input files in the `data/` directory next to `analysis_config.py`. Put your ROOT files there, or set `FCC_INPUT_DIR` to another directory.

The expected default filenames are listed in `SAMPLE_FILES` and `HADRONIC_SAMPLE_FILES` in `analysis_config.py`. Edit those mappings if your local files have different names. Absolute file paths are also supported, but using the shared input directory is recommended.

The shared configuration also exposes the tree names and physics settings. For example:

```bash
export FCC_INPUT_DIR=/path/to/my/root/files
export FCC_OUTPUT_DIR=/path/to/my/analysis/output
export FCC_EVENTS_TREE=Events
export FCC_BDT_TREE=TreeS
```

The input trees must provide the branches used by the relevant notebook. The leptonic and hadronic analyses use LHE particle information, including `nLHEPart`, `LHEPart_pdgId`, `LHEPart_pt`, `LHEPart_eta`, `LHEPart_phi`, and `LHEPart_mass`. A different tree or branch schema may require adapting the analysis code.

## Running the notebooks

Run Jupyter from the repository directory so that `analysis_config.py` can be imported and TMVA's generated `dataset/` and `dataset_adr/` folders remain inside the project:

```bash
jupyter lab
```

Open the notebooks and run their cells from top to bottom. The recommended order is:

1. `study_angles.ipynb`
2. `likelihood_fit.ipynb`
3. `tmva_analysis.ipynb`
4. `hadronic_analysis.ipynb`

`likelihood_fit.ipynb` expects the leptonic template file produced by `study_angles.ipynb`. The leptonic and hadronic TMVA notebooks require their corresponding event-level ROOT input files. The TMVA workflows are otherwise independent of the conventional likelihood workflow.

## Output files

Generated ROOT files and plots are written under `outputs/` by default. TMVA weight XML files are stored in the relative `dataset/weights/` and `dataset_adr/weights/` folders because TMVA manages those directories itself. These generated folders are ignored by Git.

Set `FCC_OUTPUT_DIR` to move regular ROOT and plot outputs elsewhere. TMVA dataset folders remain relative to the repository working directory.

## Physics configuration

Default luminosities and cross-sections reproduce the numerical setup encoded in the source analysis. They are configuration inputs, not universal constants. Update the dictionaries and environment variables in `analysis_config.py` to match the generator samples and physics scenario being studied.

The likelihood workflows use Asimov templates where specified in the source notebooks. They therefore estimate expected sensitivity rather than fitting observed collision data.
