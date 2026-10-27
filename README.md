# FCCee-Polarization

This repository contains the Python scripts used for the FCC-ee study of the polarization components of the \(W^+W^-\) signal. The analysis considers the longitudinal (`00`), transverse (`TT`), and mixed (`0T` and `T0`) samples, together with a background sample. It includes kinematic studies, template extraction, binned likelihood fits, and a binary signal/background Boosted Decision Tree (BDT).

## Analysis workflow

1. **Angular and kinematic studies (`study_angles.py`)**  
   Studies the production angle \(\cos\theta_W\) in generated events and the decay-angle distributions for the signal polarization components. It also compares the generator-level and observable definitions of \(\cos\theta^*\), studies the leading-lepton transverse momentum \(p_T\) and dilepton invariant mass \(m_{\ell\ell}\), and compares signal yields with the background. Histograms are normalized to the expected luminosity and plots and yield summaries are produced.

2. **Template extraction (`extract_histograms.py`)**  
   Builds the signal and background templates for \(\cos\theta^*\), leading-lepton \(p_T\), and \(m_{\ell\ell}\). The histograms are scaled to the expected integrated luminosity and saved to `histos.root` for the conventional likelihood fit. This script can be run independently when only the fit templates are needed.

3. **Statistical likelihood fits (`likelihood_fit.py`)**  
   Reads the templates from `histos.root` and performs binned Poisson likelihood fits with `iminuit`. For each of the three observables, it runs three one-parameter fits, three two-parameter fits, and one three-parameter fit. The fitted parameters are the signal-strength factors \(\mu_{00}\), \(\mu_{TT}\), and \(\mu_{\mathrm{mix}}\). The script reports fit validity, parameter uncertainties, covariance quality, and correlations, and compares the uncertainties from the three-dimensional fits.

4. **BDT classification and likelihood fit (`bdt_analysis.py`)**  
   Prepares ROOT trees from the signal and background samples, splits events into a TMVA training/test sample and an independent holdout sample, and trains a binary BDT using the leading-lepton \(p_T\), dilepton invariant mass, and observable \(\cos\theta^*\). The holdout sample is used to build the BDT-score templates, which are then used in a binned likelihood fit to extract the signal-strength factors.

## Repository structure

- **`study_angles.py`** — angular and kinematic studies, signal/background comparisons, yield and efficiency summaries, and template output.
- **`extract_histograms.py`** — standalone extraction of the templates used by the conventional likelihood fit.
- **`likelihood_fit.py`** — one-, two-, and three-parameter likelihood fits for the three kinematic observables, including correlation and uncertainty comparisons.
- **`bdt_analysis.py`** — TMVA BDT training, holdout-template construction, and likelihood fit using the BDT score.

## Input samples

The scripts use the following ROOT files:

- `both_longitudinal.root`
- `both_transverse.root`
- `mixed.root`
- `mixed2.root`
- `decay_00.root`
- `decay_TT.root`
- `decay_0T.root`
- `decay_T0.root`
- `bckgrnd_dollaro.root`

The angular studies and conventional template extraction expect an `Events` tree. For the BDT preparation, the code first looks for a `TreeS` tree and otherwise searches for a `TTree` containing the `nLHEPart` branch. The BDT normalization also reads the `Events` tree to obtain the number of generated events.

## Configuration and prerequisites

The analysis requires ROOT with PyROOT and TMVA enabled, as well as Python, NumPy, and `iminuit` (version 2 or later). `pandas` and IPython are optional and are used to display fit-result tables in notebook environments.

Install the Python packages with:

```bash
python -m pip install "iminuit>=2.0" numpy
```

If table display in Jupyter is desired, install the optional packages as well:

```bash
python -m pip install pandas ipython
```

The scripts currently contain the ROOT library path `/home/davide/root/lib`. The input paths for `study_angles.py` and `extract_histograms.py` are set by `BASE_DIR = "/home/davide/FCCee-TreeProducer"`. Update these paths in the scripts if ROOT or the input files are stored elsewhere. `bdt_analysis.py` opens its input files by name relative to the current working directory, so run it from a directory containing the ROOT samples or adjust the file paths accordingly.

## Running the analysis

Run the commands from the repository directory after configuring ROOT and the input-file paths.

### 1. Angular and kinematic studies

```bash
python study_angles.py
```

This runs the studies in the first analysis script and writes `histos.root` with the templates used by the conventional fit.

To produce only the conventional templates, without running the full set of studies and plots, use:

```bash
python extract_histograms.py
```

Both scripts write to `histos.root` in the current working directory.

### 2. Conventional likelihood fits

After creating `histos.root`, run:

```bash
python likelihood_fit.py
```

The script reads the signal and background templates for `cos_theta_star`, `lep_pt`, and `mll`, then runs the seven fit configurations for each observable.

### 3. BDT analysis

Run the BDT workflow from a directory containing the input ROOT files:

```bash
python bdt_analysis.py
```

The script creates prepared ROOT trees, the TMVA output file `TMVA_Binario_Output.root`, and the BDT weights at `dataset/weights/AnalisiTesi_BDT.weights.xml`. It then creates the holdout BDT templates and performs the likelihood fit.

## Normalization and statistical model

The scripts use an integrated luminosity of \(19.2\,\mathrm{ab}^{-1}\) and the cross-sections configured in the source code:

- `00`: \(0.01276\,\mathrm{pb}\)
- `TT`: \(0.05131\,\mathrm{pb}\)
- `0T`: \(0.01975\,\mathrm{pb}\)
- `T0`: \(0.01975\,\mathrm{pb}\)
- Background: \(0.01676\,\mathrm{pb}\)

The likelihood fits use an Asimov dataset formed from the sum of the nominal signal and background templates. They therefore estimate the expected statistical sensitivity rather than performing a fit to observed collision data. In the conventional fits, signal-strength parameters have bounds \([0,5]\), and parameters not included in a given fit are fixed to one.
