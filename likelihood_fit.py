import numpy as np
import ROOT
from iminuit import Minuit

# --- 1. OPEN INPUT ROOT FILE AND EXTRACT TEMPLATES ---
f_input = ROOT.TFile.Open("histos.root", "READ")
if not f_input or f_input.IsZombie():
    raise FileNotFoundError("Could not open 'histos.root'. Run extraction script first.")

h_bkg   = f_input.Get("h_cos_theta_star_bkg")
h_00    = f_input.Get("h_cos_theta_star_00")
h_tt    = f_input.Get("h_cos_theta_star_tt")
h_mixed = f_input.Get("h_cos_theta_star_mixed")

n_bins = h_bkg.GetNbinsX()

# Convert ROOT histograms into NumPy arrays for high-performance computation
bkg_template   = np.array([h_bkg.GetBinContent(i) for i in range(1, n_bins + 1)])
sig00_template = np.array([h_00.GetBinContent(i) for i in range(1, n_bins + 1)])
sigtt_template = np.array([h_tt.GetBinContent(i) for i in range(1, n_bins + 1)])
sigmix_template= np.array([h_mixed.GetBinContent(i) for i in range(1, n_bins + 1)])

# --- 2. DEFINE DATA (Asimov dataset baseline) ---
# Sum of Standard Model components / Monte Carlo expectation
data = bkg_template + sig00_template + sigtt_template + sigmix_template

# --- 3. DEFINE POISSON NEGATIVE LOG-LIKELIHOOD (-2 ln L) ---
def negative_log_likelihood(mu_00, mu_tt, mu_mix):
    # Expected bin-by-bin yield model
    nu = bkg_template + mu_00 * sig00_template + mu_tt * sigtt_template + mu_mix * sigmix_template
    
    # Mathematical protection against log(0) or zero division
    nu = np.where(nu <= 0, 1e-5, nu)
    
    # Binned Poisson NLL formulation (-2 ln L)
    term_data = np.where(data > 0, data * np.log(data / nu), 0.0)
    nll = 2.0 * np.sum(nu - data + term_data)
    
    return nll

# --- 4. MINUIT INITIALIZATION & CONFIGURATION ---
m = Minuit(negative_log_likelihood, mu_00=1.0, mu_tt=1.0, mu_mix=1.0)

# Set physical parameter boundaries
m.limits["mu_00"] = (0.0, 5.0)
m.limits["mu_tt"] = (0.0, 5.0)
m.limits["mu_mix"] = (0.0, 5.0)

# Set error definition for -2lnL (errordef = 1.0 yields 68% CL intervals)
m.errordef = 1.0  

# --- 5. RUN MINIMIZATION (FIT) ---
m.migrad()  # Execute MIGRAD minimizer
m.hesse()   # Compute covariance matrix for uncertainties

# --- 6. PRINT RESULTS AND DIAGNOSTICS ---
print("--- 3D FIT RESULTS ---")
print(m.values)
print(m.errors)
print("Is the fit valid?", m.valid)
print("\nCorrelation Matrix:")
print(m.covariance.correlation())