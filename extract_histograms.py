# ============================================================
# Estrazione dei template per i fit di likelihood.
# Versione script del blocco "STUDIO DEL BACKGROUND (1)" del primo notebook.
# La selezione e la normalizzazione sono mantenute come nel notebook.
# ============================================================
import ROOT

from analysis_config import SAMPLE_FILES, EVENTS_TREE_NAME, resolve_input, output_path, TEMPLATE_FILE

# --- CONFIGURAZIONE E IMPOSTAZIONI GLOBALI ---
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetEndErrorSize(0)

N_BINS = 80
LUMI_PB = 19.2e6

# Sezioni d'urto in pb
SIGMA_00 = 0.01276
SIGMA_TT = 0.05131
SIGMA_0T = 0.01975
SIGMA_MIXED = 2 * SIGMA_0T
SIGMA_BKG = 0.01676

Y_TITLE = "Events / bin"
TITLE_COS_THETA_STAR = "cos(#theta^{*})"
TITLE_LEP_PT = "Lepton p_{T} [GeV]"
TITLE_MLL = "M_{ll} [GeV]"


def open_tree(filename, tree_name=EVENTS_TREE_NAME):
    """Apre un file ROOT da BASE_DIR e ritorna il file e l'albero."""
    path = resolve_input(filename)
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise FileNotFoundError(f"Impossibile aprire il file ROOT: {path}")
    tree = f.Get(tree_name)
    if not tree:
        raise KeyError(f"Albero '{tree_name}' non trovato in {path}")
    return f, tree


def make_hist(name, x_title, y_title, x_min, x_max, sumw2=True):
    """Crea un TH1F con N_BINS bin, come nel primo notebook."""
    hist = ROOT.TH1F(name, f";{x_title};{y_title}", N_BINS, x_min, x_max)
    if sumw2:
        hist.Sumw2()
    return hist


def get_leptons(event):
    """Ritorna i leptoni carichi (e, mu) dell'evento come TLorentzVector."""
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
    """Riempie l'istogramma con una entry per evento usando i due leptoni.

    Selezione: esattamente due leptoni carichi e pT > 10 GeV per entrambi.
    La definizione di cos_theta_star è quella osservabile del notebook.
    """
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

        elif variable == "lep_pt":
            hist.Fill(max(leptons[0].Pt(), leptons[1].Pt()))

        elif variable == "mll":
            v_ll = leptons[0] + leptons[1]
            hist.Fill(v_ll.M())


def event_weight(sigma_pb, n_generated):
    """Peso nominale per evento: sigma * luminosità / N_generati."""
    return sigma_pb * LUMI_PB / n_generated


def scale_to_luminosity(hist, sigma_pb, n_generated):
    """Normalizza l'istogramma alla luminosità integrata."""
    if n_generated > 0:
        hist.Scale(event_weight(sigma_pb, n_generated))
    else:
        print(
            f"Attenzione: nessun evento generato per {hist.GetName()}, "
            "normalizzazione non eseguita."
        )


def fill_and_scale(trees, hist, variable, sigma_pb):
    """Riempie tutti gli alberi e normalizza usando la somma degli eventi generati.

    Per il campione misto il notebook applica SIGMA_MIXED al totale dei due file,
    con N_generati = N(0T) + N(T0). Questo comportamento è mantenuto.
    """
    n_generated = sum(tree.GetEntries() for tree in trees)
    for tree in trees:
        fill_lepton_variable(tree, hist, variable)
    scale_to_luminosity(hist, sigma_pb, n_generated)


# --- Apertura degli alberi ---
f_dec_00, tree_dec_00 = open_tree(SAMPLE_FILES["decay_00"])
f_dec_tt, tree_dec_tt = open_tree(SAMPLE_FILES["decay_tt"])
f_dec_0t, tree_dec_0t = open_tree(SAMPLE_FILES["decay_0t"])
f_dec_t0, tree_dec_t0 = open_tree(SAMPLE_FILES["decay_t0"])
f_bkg, tree_bkg = open_tree(SAMPLE_FILES["background"])

# --- Estrazione delle tre osservabili per segnale e background ---
# sumw2=False mantiene esattamente la scelta del blocco di estrazione del notebook.

def crea_template(variable, x_min, x_max, x_title, name):
    h_00 = make_hist(f"h_{name}_00", x_title, Y_TITLE, x_min, x_max, sumw2=False)
    h_tt = make_hist(f"h_{name}_tt", x_title, Y_TITLE, x_min, x_max, sumw2=False)
    h_mixed = make_hist(f"h_{name}_mixed", x_title, Y_TITLE, x_min, x_max, sumw2=False)
    h_bkg = make_hist(f"h_{name}_bkg", x_title, Y_TITLE, x_min, x_max, sumw2=False)

    fill_and_scale([tree_dec_00], h_00, variable, SIGMA_00)
    fill_and_scale([tree_dec_tt], h_tt, variable, SIGMA_TT)
    fill_and_scale([tree_dec_0t, tree_dec_t0], h_mixed, variable, SIGMA_MIXED)
    fill_and_scale([tree_bkg], h_bkg, variable, SIGMA_BKG)

    return {"00": h_00, "TT": h_tt, "mixed": h_mixed, "bkg": h_bkg}


res_comp_cos = crea_template(
    "cos_theta_star", -1.0, 1.0, TITLE_COS_THETA_STAR, "cos_theta_star"
)
res_comp_pt = crea_template("lep_pt", 0.0, 80.0, TITLE_LEP_PT, "lep_pt")
res_comp_mll = crea_template("mll", 0.0, 140.0, TITLE_MLL, "mll")

# --- Salvataggio degli istogrammi su file ROOT con i nomi usati dal fit ---
f_out = ROOT.TFile.Open(TEMPLATE_FILE, "RECREATE")
if not f_out or f_out.IsZombie():
    raise OSError(f"Impossibile creare {TEMPLATE_FILE!r}")
f_out.cd()

for name, result in [
    ("cos_theta_star", res_comp_cos),
    ("lep_pt", res_comp_pt),
    ("mll", res_comp_mll),
]:
    result["00"].Write(f"h_{name}_00")
    result["TT"].Write(f"h_{name}_tt")
    result["mixed"].Write(f"h_{name}_mixed")
    result["bkg"].Write(f"h_{name}_bkg")

f_out.Close()
print(f"Tutti gli istogrammi (cos_theta_star, lep_pt, mll) sono stati salvati in {TEMPLATE_FILE!r}.")
