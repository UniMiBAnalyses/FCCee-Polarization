import os
import numpy as np
import ROOT

# --- CONFIGURATION & GLOBAL SETTINGS ---
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)

# Path configuration (can be adapted or set to local directory)
BASE_DIR = os.getenv("FCC_ANALYSIS_DIR", ".")
BINS = 40
LUMI_PB = 19.2e6  # 19.2 ab^-1 expressed in pb^-1[cite: 1]

# Cross-sections (in pb)
SIGMA_00 = 0.01276
SIGMA_TT = 0.05131
SIGMA_0T = 0.01975
SIGMA_BKG = 0.01676

def open_tree(filename, tree_name="Events"):
    """Open a ROOT file and return the file and tree objects."""
    path = os.path.join(BASE_DIR, filename)
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise FileNotFoundError(f"Cannot open ROOT file: {path}")
    tree = f.Get(tree_name)
    if not tree:
        raise KeyError(f"Tree '{tree_name}' not found in {path}")
    return f, tree

def get_leptons(event):
    """Return charged leptons (e, mu) as TLorentzVectors."""
    leptons = []
    for p in range(event.nLHEPart):
        if abs(event.LHEPart_pdgId[p]) in (11, 13) and event.LHEPart_pt[p] > 0.1:
            v = ROOT.TLorentzVector()
            v.SetPtEtaPhiM(
                event.LHEPart_pt[p],
                event.LHEPart_eta[p],
                event.LHEPart_phi[p],
                event.LHEPart_mass[p],
            )
            leptons.append(v)
    return leptons

def fill_lepton_variable(tree, hist, variable):
    """Fill histogram with the decay lepton variable (e.g. cos_theta_star)."""
    if not tree:
        return
    for event in tree:
        leptons = get_leptons(event)
        if len(leptons) != 2:
            continue
        if leptons[0].Pt() < 10 or leptons[1].Pt() < 10:
            continue
            
        if variable == "cos_theta_star":
            v_ll = leptons[0] + leptons[1]
            l1_rest = ROOT.TLorentzVector(leptons[0])
            l1_rest.Boost(-v_ll.BoostVector())
            hist.Fill(l1_rest.CosTheta())

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    # Open input files (assuming they are placed in the working directory)
    _, tree_dec_00 = open_tree("decay_00.root")
    _, tree_dec_tt = open_tree("decay_TT.root")
    _, tree_dec_0t = open_tree("decay_0T.root")
    _, tree_dec_t0 = open_tree("decay_T0.root")
    _, tree_bkg    = open_tree("bckgrnd_dollaro.root")

    # Define histograms
    h_00    = ROOT.TH1F("h_cos_theta_star_00", ";cos(#theta^{*});Events / bin", BINS, -1.0, 1.0)
    h_tt    = ROOT.TH1F("h_cos_theta_star_tt", ";cos(#theta^{*});Events / bin", BINS, -1.0, 1.0)
    h_mixed = ROOT.TH1F("h_cos_theta_star_mixed", ";cos(#theta^{*});Events / bin", BINS, -1.0, 1.0)
    h_bkg   = ROOT.TH1F("h_cos_theta_star_bkg", ";cos(#theta^{*});Events / bin", BINS, -1.0, 1.0)

    dataset_map = [
        (tree_dec_00, h_00, SIGMA_00),
        (tree_dec_tt, h_tt, SIGMA_TT),
        ([tree_dec_0t, tree_dec_t0], h_mixed, SIGMA_0T * 2),
        (tree_bkg, h_bkg, SIGMA_BKG),
    ]

    for source, hist, sig in dataset_map:
        trees = source if isinstance(source, list) else [source]
        n_gen_total = sum(t.GetEntries() for t in trees)
        
        for t in trees:
            fill_lepton_variable(t, hist, "cos_theta_star")
            
        if n_gen_total > 0:
            hist.Scale((sig * LUMI_PB) / n_gen_total)

    # Save histograms to output ROOT file for the fit stage
    f_out = ROOT.TFile.Open("histos.root", "RECREATE")
    h_00.Write()
    h_tt.Write()
    h_mixed.Write()
    h_bkg.Write()
    f_out.Close()
    
    print("Histograms successfully extracted and saved to 'histos.root'!")