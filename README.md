# FCCee-Polarization

This repository contains Python scripts for the polarization analysis of $W^+W^-$ signals ($00$, $TT$, and mixed states).

## Logical steps:
1. **Preliminary Angular Studies (`study_angles.py`)**: 
   - Analysis of the production angle ($\cos\theta_W$) for non-decayed events.
   - Investigation of the decay angle ($\cos\theta^*$) distributions for pure longitudinal, transverse, and mixed polarization states using decayed events (without background)[cite: 1, 3].
3. **Template Extraction & Luminosity Scaling (`extract_histograms.py`)**: Normalization of MC samples to the expected luminosity ($19.2 \text{ ab}^{-1}$) including backgrounds, generating the final $\cos\theta^*$ templates.
4. **Statistical Likelihood Fit (`likelihood_fit.py`)**: Implementation of a binned Poisson profile likelihood to extract polarization fractions ($\mu$) using `iminuit`[cite: 2].

## Repository Structure
- **`study_angles.py`**: Analyzes pure signal production and decay angles ($\theta_W$ and $\theta^*$).
- **`extract_histograms.py`**: Extracts normalized $\cos\theta^*$ templates including backgrounds and saves them to `histos.root`.
- **`likelihood_fit.py`**: Performs binned Poisson maximum likelihood fits and correlation checks.
- **`requirements.txt`**: List of Python dependencies.

## Prerequisites
Ensure you have ROOT (with PyROOT) installed, then install the required Python packages:
```bash
pip install -r requirements.txt --break-system-packages
