# FCC-ee W-pair polarization analysis

This repository contains Jupyter notebooks for generator-level studies of polarized \(W^+W^-\) production and decay at FCC-ee. It covers leptonic and hadronic final states, kinematic distributions, signal-versus-background comparisons, binned likelihood fits, and TMVA Boosted Decision Trees (BDTs).

All analysis workflows are run from Jupyter notebooks. The shared `analysis_config.py` module contains the input-file mapping, physics settings, tree names, and output locations used by the notebooks.

## Repository contents

| File | Purpose |
|---|---|
| `study_angles.ipynb` | Leptonic production and decay angles, lepton kinematics, and comparisons of signal polarization components with the full dilepton background. Produces templates for the conventional likelihood fit. |
| `likelihood_fit.ipynb` | One-, two-, and three-parameter binned likelihood fits and one- and two-dimensional likelihood scans using the leptonic templates. |
| `tmva_analysis.ipynb` | Leptonic signal/background BDT training, test-sample distributions, independent holdout templates, and a likelihood fit to the BDT score. |
| `hadronic_analysis.ipynb` | Hadronic event reconstruction, kinematic templates, likelihood fits and scans, BDT training, and a likelihood fit to the BDT score. |
| `analysis_config.py` | Shared configuration for input and output directories, sample filenames, tree names, cross-sections, luminosity, and histogram binning. |
| `requirements.txt` | Python packages used by the notebooks. |
| `data/README.md` | Instructions for providing the input ROOT files. |

## Requirements

The analysis requires:

- Python and Jupyter.
- ROOT with PyROOT available in the Python environment used by the Jupyter kernel.
- ROOT built with TMVA support for `tmva_analysis.ipynb` and the BDT section of `hadronic_analysis.ipynb`.
- The Python packages listed in `requirements.txt`.

ROOT is an external prerequisite and is not installed through `requirements.txt`. Its installation depends on the operating system and Python environment. Install the Python packages in the same environment that provides PyROOT:

```bash
python -m pip install -r requirements.txt
```

Check that the active Python interpreter can import ROOT:

```bash
python -c "import ROOT; print(ROOT.gROOT.GetVersion())"
```

For notebooks that use TMVA, also check that TMVA is available:

```bash
python -c "import ROOT; ROOT.TMVA.Tools.Instance(); print('TMVA is available')"
```

Start Jupyter from this environment. In Jupyter, select the kernel that uses the same Python interpreter; otherwise the notebook may report `ModuleNotFoundError: No module named 'ROOT'` even when ROOT is installed elsewhere on the machine.

## Input ROOT files

Input ROOT files are not included in this repository. Obtain or generate compatible samples separately and place them in `data/`, or configure another input directory with `FCC_INPUT_DIR`.

The default filenames are defined in `SAMPLE_FILES` and `HADRONIC_SAMPLE_FILES` in `analysis_config.py`. They are examples for the supplied analysis configuration; replace them with the filenames of your own datasets when needed.

### Leptonic and production samples

| Configuration key | Default filename | Used for |
|---|---|---|
| `production_00` | `both_longitudinal.root` | Longitudinal production-angle sample |
| `production_tt` | `both_transverse.root` | Transverse production-angle sample |
| `production_mixed_1` | `mixed.root` | First mixed production sample (0T) |
| `production_mixed_2` | `mixed2.root` | Second mixed production sample (T0) |
| `decay_00` | `decay_00.root` | Longitudinal leptonic decay sample |
| `decay_tt` | `decay_TT.root` | Transverse leptonic decay sample |
| `decay_0t` | `decay_0T.root` | Mixed 0T leptonic decay sample |
| `decay_t0` | `decay_T0.root` | Mixed T0 leptonic decay sample |
| `background` | `bckgrnd_dollaro.root` | Leptonic background sample |

The production-angle samples are used by `study_angles.ipynb`. The four decay samples and the background are used by the leptonic angular, likelihood-template, and BDT workflows. `likelihood_fit.ipynb` reads the template file created by `study_angles.ipynb` and does not need to reopen the original event samples.

### Hadronic samples

| Configuration key | Default filename |
|---|---|
| `00` | `ww_hadronic_00.root` |
| `TT` | `ww_hadronic_TT.root` |
| `0T` | `ww_hadronic_0T.root` |
| `T0` | `ww_hadronic_T0.root` |
| `background` | `bckgrnd_had.root` |

These five files are used by `hadronic_analysis.ipynb`.

### Changing filenames and directories

Edit the mappings in `analysis_config.py` if your files have different names. Relative filenames are resolved under `INPUT_DIR`; absolute paths are also accepted. By default, `INPUT_DIR` is the repository's `data/` directory.

You can instead set environment variables before starting Jupyter:

```bash
export FCC_INPUT_DIR=/path/to/my/root/files
export FCC_OUTPUT_DIR=/path/to/my/analysis/output
```

On Windows PowerShell, the equivalent commands are:

```powershell
$env:FCC_INPUT_DIR = "C:\path\to\my\root\files"
$env:FCC_OUTPUT_DIR = "C:\path\to\my\analysis\output"
```

The shared configuration also allows you to change tree names and common physics settings. Relevant variables include `FCC_EVENTS_TREE`, `FCC_BDT_TREE`, `FCC_LUMI_PB`, `FCC_COLLISION_ENERGY_GEV`, and the cross-section variables documented in `analysis_config.py`. You can edit the Python configuration file directly instead of using environment variables.

## ROOT tree and branch structure

The input event trees must contain the particle-level information accessed by the notebooks. The standard event tree is named `Events` and can be changed with `FCC_EVENTS_TREE`.

The leptonic and hadronic event selections use these branches:

- `nLHEPart`
- `LHEPart_pdgId`
- `LHEPart_pt`
- `LHEPart_eta`
- `LHEPart_phi`
- `LHEPart_mass`

The leptonic TMVA workflow first looks for the tree named by `FCC_BDT_TREE` (default: `TreeS`). If it is not present, it searches the file for a `TTree` containing the `nLHEPart` branch. The selected tree must expose the LHE particle branches listed above.

If your files use a different tree or branch schema, adapt the tree configuration or the relevant branch access in the notebook. Changing only the filenames is sufficient only when the input ROOT files have a compatible structure and contain the required physical samples.

## Running the notebooks

1. Configure the input filenames and physics settings in `analysis_config.py`.
2. Place the required ROOT files in `data/` or set `FCC_INPUT_DIR`.
3. Activate the Python environment containing ROOT, PyROOT, TMVA (where required), and the packages in `requirements.txt`.
4. Start Jupyter from the repository directory:

   ```bash
   jupyter lab
   ```

5. Open the required notebook and use **Run All** (or execute its cells from top to bottom).

The notebooks can be run as follows:

| Notebook | Inputs / prerequisites | What it produces |
|---|---|---|
| `study_angles.ipynb` | Four production samples, four leptonic decay samples, and the leptonic background | Angular and kinematic plots; `outputs/histos.root`, used by the conventional leptonic likelihood fit |
| `likelihood_fit.ipynb` | `outputs/histos.root` from `study_angles.ipynb` | Fit summaries, likelihood scans, and plots |
| `tmva_analysis.ipynb` | Four leptonic decay samples and the leptonic background | BDT training output, test-sample plots, BDT-score templates, and a BDT likelihood fit |
| `hadronic_analysis.ipynb` | Four hadronic signal samples and one hadronic background sample | Hadronic distributions, `outputs/histos_hadronic.root`, likelihood results, and BDT outputs |

`likelihood_fit.ipynb` depends on the output of `study_angles.ipynb`. The leptonic BDT and hadronic workflows use their own event-level input samples and can be run independently of the conventional leptonic likelihood workflow.

## Analysis details

### Leptonic angular and kinematic study

`study_angles.ipynb` studies the absolute W production-angle cosine, charge-separated W decay angles, leading- and subleading-lepton transverse momentum, and dilepton invariant mass. It also compares the signal polarization components with a background sample retaining events with exactly two charged electrons and/or muons. The full-background mode includes the available dilepton flavour combinations; it is not restricted to the `e+ mu-` state.

### Conventional likelihood fit

`likelihood_fit.ipynb` fits the longitudinal, transverse, and mixed signal-strength parameters (`mu_00`, `mu_tt`, and `mu_mix`) using the leptonic templates. It performs one-, two-, and three-parameter fits and likelihood scans for the available observables. The fit uses an Asimov dataset formed from the nominal signal and background templates, so the results describe expected sensitivity rather than a fit to observed collision data.

### Leptonic BDT

`tmva_analysis.ipynb` trains a binary signal/background BDT using the leading-lepton transverse momentum, dilepton invariant mass, and observable `cos(theta*)`. The event split reserves 30% of each input file for independent holdout templates; the split is assigned before the event selection and is reproducible for a given filename and configured seed. The remaining 70% is supplied to TMVA for its internal training/test split.

The `SOLO_E_MU` setting in the notebook controls the background flavour selection. Its default value is `False`, which retains all supported two-lepton electron/muon final states. Set it to `True` to restrict the background to the `e+ mu-` state.

### Hadronic analysis

`hadronic_analysis.ipynb` selects events with four final-state quarks/gluons and reconstructs two W candidates by testing the three possible pairings of the four partons. The pairing is chosen by its proximity to the configured W mass. The main observables are the absolute decay-angle cosine and the transverse momentum of the leading parton.

For the BDT workflow, events are split deterministically by their entry index: 70% are used for training and 30% are reserved as holdout events for the likelihood templates. The hadronic notebook also performs the conventional likelihood fits and scans using the hadronic templates.

## Output files

Regular ROOT outputs and plots are written to `outputs/` by default. The main template files are:

- `outputs/histos.root` — leptonic templates created by `study_angles.ipynb`.
- `outputs/histos_hadronic.root` — hadronic templates created by `hadronic_analysis.ipynb`.

Plots are saved under `outputs/plots/` and likelihood-scan plots under `outputs/likelihood_plots/`. TMVA also creates intermediate ROOT files under `outputs/`. The TMVA weight XML files are stored in `dataset/weights/` and `dataset_adr/weights/`, relative to the repository directory, because TMVA manages these working directories.

Set `FCC_OUTPUT_DIR` to place the regular ROOT files and plots in another directory. The TMVA `dataset/` and `dataset_adr/` working directories remain relative to the repository directory. Generated files and local ROOT datasets are excluded by `.gitignore`; do not commit large input samples or generated products unless you intend to publish them.

## Physics configuration and normalization

The default integrated luminosity, collision energy, cross-sections, histogram ranges, and binning are configured in `analysis_config.py`. These values describe the supplied analysis setup and are not universal constants. Update them to match the generator samples and physics scenario being studied. In particular, luminosity normalization uses the configured cross-section in pb and the generated event count of each input sample.

The two mixed polarization samples, 0T and T0, are handled as separate contributions where appropriate and then combined for the mixed component. Check that the configured cross-sections correspond to the conventions used for your generated samples before interpreting expected event yields.
