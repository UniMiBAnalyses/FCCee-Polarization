import os
import ROOT

# --- CONFIGURATION & GLOBAL SETTINGS ---
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)

BASE_DIR = os.getenv("FCC_ANALYSIS_DIR", ".")
BINS = 40
LUMI_PB = 19.2e6  # 19.2 ab^-1 in pb^-1[cite: 1]

# Cross-sections (in pb)
SIGMA_00 = 0.01276
SIGMA_TT = 0.05131
SIGMA_0T = 0.01975

def open_tree(filename, tree_name="Events"):
    """Open a ROOT file from BASE_DIR and return the tree."""
    path = os.path.join(BASE_DIR, filename)
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise FileNotFoundError(f"Cannot open ROOT file: {path}")
    return f, f.Get(tree_name)

def get_leptons(event):
    """Return the two charged leptons as TLorentzVectors."""
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

def fill_production_angle(tree, hist):
    """Fill histogram with cos(theta_W) of the W flight direction (Pz / |P|)."""
    if not tree:
        return
    for event in tree:
        for p in range(event.nLHEPart):
            pdg = event.LHEPart_pdgId[p]
            if abs(pdg) != 24:
                continue
            v = ROOT.TLorentzVector()
            v.SetPtEtaPhiM(
                event.LHEPart_pt[p],
                event.LHEPart_eta[p],
                event.LHEPart_phi[p],
                event.LHEPart_mass[p],
            )
            p_tot = v.Vect().Mag()
            if p_tot > 0:
                hist.Fill(v.Pz() / p_tot)

def fill_decay_angle(tree, hist):
    """Fill histogram with cos(theta*) of the decay lepton in W rest frame."""
    if not tree:
        return
    for event in tree:
        leptons = get_leptons(event)
        if len(leptons) != 2:
            continue
        if leptons[0].Pt() < 10 or leptons[1].Pt() < 10:
            continue
        v_ll = leptons[0] + leptons[1]
        l1_rest = ROOT.TLorentzVector(leptons[0])
        l1_rest.Boost(-v_ll.BoostVector())
        hist.Fill(l1_rest.CosTheta())

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    # Open trees for production & decay studies
    _, tree_l = open_tree("both_longitudinal.root")
    _, tree_t = open_tree("both_transverse.root")
    _, tree_m1 = open_tree("mixed.root")
    _, tree_m2 = open_tree("mixed2.root")

    _, tree_dec_00 = open_tree("decay_00.root")
    _, tree_dec_tt = open_tree("decay_TT.root")
    _, tree_dec_0t = open_tree("decay_0T.root")
    _, tree_dec_t0 = open_tree("decay_T0.root")

    print("Files successfully loaded. Ready for angular distribution studies!")