# FCCee-Polarization

This repository contains Python scripts for the polarization analysis of the \\(W^+W^-\\) signal in the FCC-ee framework. It includes angular and kinematic studies, signal and background template extraction, binned likelihood fits, and a binary Boosted Decision Tree (BDT) analysis.

## Logical steps

1. **Angular and kinematic studies (`study_angles.py`)**: studies the production angle \\(\\cos\\theta_W\\), decay-angle distributions, leading-lepton transverse momentum, and dilepton invariant mass. It also compares the signal polarization components with the background.
2. **Template extraction (`extract_histograms.py`)**: builds the signal and background templates, normalizes them to the configured integrated luminosity, and saves them in a ROOT file.
3. **Statistical likelihood fit (`likelihood_fit.py`)**: performs binned Poisson likelihood fits with `iminuit` to estimate the signal-strength parameters and their uncertainties and correlations.
4. **BDT analysis (`bdt_analysis.py`)**: trains a binary signal/background classifier with TMVA, builds BDT-score templates from a holdout sample, and performs a likelihood fit.

## Repository structure

- `study_angles.py` — angular and kinematic studies.
- `extract_histograms.py` — extraction of templates for the conventional likelihood fit.
- `likelihood_fit.py` — one-, two-, and three-parameter fits for the kinematic observables.
- `bdt_analysis.py` — BDT training and likelihood fit using the BDT score.
- `analysis_config.py` — input sample mapping, ROOT tree names, and configurable settings for the angular and kinematic study.
- `requirements.txt` — Python dependencies.

## Input data and configuration

Input ROOT datasets are not included. You do **not** need to use the same files or file names as the original analysis: edit `SAMPLE_FILES` in `analysis_config.py` to map each sample role to your own file name or absolute path. Relative paths are resolved from `FCC_ANALYSIS_DIR` (default: the current directory). Generated files are written to `FCC_OUTPUT_DIR` (default: the current directory).

The same configuration file also exposes the settings used by `study_angles.py`: `CROSS_SECTIONS_PB`, `LUMI_PB`, `N_BINS`, `COLLISION_ENERGY_GEV`, and `OBSERVABLE_RANGES`. Adapt these values to the physical setup and datasets being analyzed. The cross-sections are required for luminosity normalization; they are configuration inputs rather than values embedded in the analysis workflow.

The input samples must represent the required physical components and provide the event-level information used by the analysis. Tree names can be configured through `FCC_EVENTS_TREE` and `FCC_BDT_TREE`. The ROOT trees must contain the particle branches referenced in the scripts; changing file names does not require changing the analysis, but a different data schema may require adapting the branch access.

## Prerequisites

The analysis requires ROOT with PyROOT. The BDT workflow additionally requires TMVA enabled in ROOT. Python dependencies are listed in `requirements.txt` and can be installed with:

```bash
python -m pip install -r requirements.txt
```

## Running the analysis

After configuring the sample paths and setting up ROOT, run the required steps:

```bash
python study_angles.py
python extract_histograms.py
python likelihood_fit.py
python bdt_analysis.py
```

`study_angles.py` runs the full angular and kinematic study and also writes the templates used by the conventional fit. `extract_histograms.py` can be run independently when only the templates are needed. Run `likelihood_fit.py` after creating the templates. The BDT workflow is independent and creates its own intermediate files and trained weights.

## Normalization and statistical model

The angular and kinematic study uses the integrated luminosity and cross-sections configured in `analysis_config.py`. The likelihood fits use an Asimov dataset formed from the nominal signal and background templates, so they estimate expected statistical sensitivity rather than fitting observed collision data.
