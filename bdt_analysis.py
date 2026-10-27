# %% 1. SETUP E FUNZIONI DI SUPPORTO
import array
import random
import zlib
import numpy as np
from iminuit import Minuit

import ROOT

from analysis_config import (
    SAMPLE_FILES, EVENTS_TREE_NAME, BDT_TREE_NAME, resolve_input, output_path, OUTPUT_DIR
)

ROOT.TMVA.Tools.Instance()
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetEndErrorSize(0)

# ---------------- CONFIGURAZIONE ----------------
TRAIN_FRACTION = 0.70
SPLIT_SEED = 20261026
BINS_BDT = 50
LUMI_PB = 19.2e6  # 19.2 ab^-1 in pb^-1

SIGMA_LL = 0.01276
SIGMA_TT = 0.05131
SIGMA_0T = 0.01975
SIGMA_T0 = 0.01975
SIGMA_BKG = 0.01676

COLORS = {
    "LL": ROOT.kGreen + 2,
    "TT": ROOT.kOrange + 1,
    "LT": ROOT.kBlue + 1,
    "Bkg": ROOT.kRed + 1,
}

# Frazione holdout di ciascun file originale, calcolata prima dei tagli.
holdout_fraction_by_file = {}


def get_leptons(event):
    """Estrae i leptoni carichi e li restituisce come TLorentzVector."""
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


def open_event_tree(filename):
    """Apre il file originale e trova TreeS (o un TTree con nLHEPart)."""
    input_path = resolve_input(filename)
    f_in = ROOT.TFile.Open(input_path, "READ")
    if not f_in or f_in.IsZombie():
        raise OSError(f"Impossibile aprire il file ROOT: {input_path}")

    tree_in = f_in.Get(BDT_TREE_NAME)
    if not tree_in:
        for key in f_in.GetListOfKeys():
            obj = f_in.Get(key.GetName())
            if obj and obj.InheritsFrom("TTree") and obj.GetBranch("nLHEPart"):
                tree_in = obj
                break
    if not tree_in:
        f_in.Close()
        raise KeyError(f"Nessun TTree con branch nLHEPart trovato in {input_path}")
    return f_in, tree_in


def crea_tree_per_tmva(input_files, output_filename):
    """Crea un tree con le tre variabili e il flag is_train.

    is_train=1: eventi ammessi al training/test interno di TMVA.
    is_train=0: eventi riservati ai template indipendenti della likelihood.
    La selezione casuale è riproducibile e viene assegnata prima dei tagli.
    """
    if isinstance(input_files, str):
        input_files = [input_files]

    f_out = ROOT.TFile.Open(output_filename, "RECREATE")
    if not f_out or f_out.IsZombie():
        raise OSError(f"Impossibile creare il file ROOT: {output_filename}")

    tree_out = ROOT.TTree("T", "Tree per TMVA binario")
    lep_pt_f = array.array('f', [0.0])
    mll_f = array.array('f', [0.0])
    cos_star_f = array.array('f', [0.0])
    is_train_f = array.array('i', [0])

    tree_out.Branch("lep_pt", lep_pt_f, "lep_pt/F")
    tree_out.Branch("mll", mll_f, "mll/F")
    tree_out.Branch("cos_theta_star", cos_star_f, "cos_theta_star/F")
    tree_out.Branch("is_train", is_train_f, "is_train/I")

    try:
        for filename in input_files:
            f_in, tree_in = open_event_tree(filename)
            n_all = int(tree_in.GetEntries())
            if n_all <= 0:
                f_in.Close()
                raise ValueError(f"Il tree di input è vuoto: {filename}")

            # Seed stabile per file, così lo stesso evento riceve lo stesso flag
            # quando il file viene usato nel campione segnale globale e nei template.
            seed_file = SPLIT_SEED + zlib.crc32(filename.encode("utf-8"))
            rng = random.Random(seed_file)
            n_holdout = 0

            try:
                for event in tree_in:
                    is_train = 1 if rng.random() < TRAIN_FRACTION else 0
                    if is_train == 0:
                        n_holdout += 1

                    leptons = get_leptons(event)
                    if len(leptons) != 2:
                        continue
                    if leptons[0].Pt() < 10 or leptons[1].Pt() < 10:
                        continue

                    lep_pt_f[0] = max(leptons[0].Pt(), leptons[1].Pt())
                    v_ll = leptons[0] + leptons[1]
                    mll_f[0] = v_ll.M()
                    l1_rest = ROOT.TLorentzVector(leptons[0])
                    l1_rest.Boost(-v_ll.BoostVector())
                    cos_star_f[0] = l1_rest.CosTheta()
                    is_train_f[0] = is_train
                    tree_out.Fill()
            finally:
                f_in.Close()

            holdout_fraction_by_file[filename] = n_holdout / n_all
            print(
                f"{filename}: {n_all} eventi sorgente, "
                f"frazione holdout = {holdout_fraction_by_file[filename]:.3f}"
            )
    finally:
        f_out.cd()
        tree_out.Write()
        f_out.Close()

    print(f"Creato {output_filename}")


def n_eventi_generati(filename):
    """Numero di eventi generati, letto dal tree Events del file originale."""
    input_path = resolve_input(filename)
    f = ROOT.TFile.Open(input_path, "READ")
    if not f or f.IsZombie():
        raise OSError(f"Impossibile aprire {input_path}")
    t = f.Get(EVENTS_TREE_NAME)
    if not t:
        f.Close()
        raise KeyError(f"Tree {EVENTS_TREE_NAME!r} non trovato in {input_path}")
    n_gen = int(t.GetEntries())
    f.Close()
    if n_gen <= 0:
        raise ValueError(f"N_gen non positivo per {filename}")
    return n_gen


def make_canvas(name, title, width=900, height=650):
    c = ROOT.TCanvas(name, title, width, height)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    return c


def draw_fcc_label():
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(61)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.15, 0.92, "FCC-ee")
    latex.SetTextFont(52)
    latex.SetTextSize(0.035)
    latex.DrawLatex(0.28, 0.92, "Simulation")
    latex.SetTextFont(42)
    latex.SetTextSize(0.030)
    latex.DrawLatex(0.45, 0.92, "#sqrt{s} = 164 GeV, L = 19.2 ab^{-1}")


def style_hist(hist, color):
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.SetFillColorAlpha(color, 0.12)


def make_legend(x1, y1, x2, y2):
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    return leg


# %% 2. CREA GLI ALBERI ROOT
file_segnale = [SAMPLE_FILES["decay_00"], SAMPLE_FILES["decay_0t"], SAMPLE_FILES["decay_t0"], SAMPLE_FILES["decay_tt"]]

# Campioni usati per il training binario
crea_tree_per_tmva(file_segnale, output_path("tmva_signal.root"))
crea_tree_per_tmva(SAMPLE_FILES["background"], output_path("tmva_background.root"))

# Campioni separati, destinati esclusivamente ai template fisici della likelihood
crea_tree_per_tmva(SAMPLE_FILES["decay_00"], output_path("tmva_signal_00.root"))
crea_tree_per_tmva(SAMPLE_FILES["decay_tt"], output_path("tmva_signal_TT.root"))
crea_tree_per_tmva(SAMPLE_FILES["decay_0t"], output_path("tmva_signal_0T.root"))
crea_tree_per_tmva(SAMPLE_FILES["decay_t0"], output_path("tmva_signal_T0.root"))

print("\nPreparazione dei campioni completata.")


# %% 3. TRAINING BDT SEGNALE/FONDO
f_signal = ROOT.TFile.Open(output_path("tmva_signal.root"), "READ")
f_bkg = ROOT.TFile.Open(output_path("tmva_background.root"), "READ")
if not f_signal or f_signal.IsZombie() or not f_bkg or f_bkg.IsZombie():
    raise OSError("Impossibile aprire i file ROOT per il training.")

signal_tree = f_signal.Get("T")
bkg_tree = f_bkg.Get("T")
if not signal_tree or not bkg_tree:
    raise KeyError("Tree 'T' non trovato nei file di training.")

output_file = ROOT.TFile.Open(output_path("TMVA_Binario_Output.root"), "RECREATE")
weights_file = output_path("dataset/weights/AnalisiTesi_BDT.weights.xml")
factory = ROOT.TMVA.Factory(
    "AnalisiTesi",
    output_file,
    "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification"
)
dataloader = ROOT.TMVA.DataLoader(str(OUTPUT_DIR / "dataset"))

dataloader.AddVariable("lep_pt", "F")
dataloader.AddVariable("mll", "F")
dataloader.AddVariable("cos_theta_star", "F")
dataloader.AddSignalTree(signal_tree, 1.0)
dataloader.AddBackgroundTree(bkg_tree, 1.0)

# Solo gli eventi is_train=1 entrano nella divisione training/test interna di TMVA.
cut_training = ROOT.TCut("is_train==1")
dataloader.PrepareTrainingAndTestTree(
    cut_training,
    cut_training,
    "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:!V"
)

factory.BookMethod(
    dataloader,
    ROOT.TMVA.Types.kBDT,
    "BDT",
    "!H:!V:NTrees=100:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost"
)

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
output_file.Close()
f_signal.Close()
f_bkg.Close()

print("\nTraining BDT binario completato.")
print(f"Pesi TMVA: {weights_file}")


# %% 4. TEMPLATE DEL BDT E FIT DI LIKELIHOOD
reader = ROOT.TMVA.Reader("!Color:!Silent")
v_pt = array.array('f', [0.0])
v_mll = array.array('f', [0.0])
v_cos = array.array('f', [0.0])
reader.AddVariable("lep_pt", v_pt)
reader.AddVariable("mll", v_mll)
reader.AddVariable("cos_theta_star", v_cos)
reader.BookMVA("BDT", weights_file)


def get_bdt_histogram(contributions, hist_name):
    """Costruisce un template pesato; ogni contributo si normalizza separatamente.

    contributions: lista di tuple (file_preparato, file_originale, sezione_d'urto).
    Legge solo is_train=0, riservato ai template della likelihood.
    """
    hist = ROOT.TH1F(hist_name, ";BDT score;Eventi attesi / bin", BINS_BDT, -1.0, 1.0)
    hist.Sumw2()

    for i_contrib, (prepared_file, original_file, sigma) in enumerate(contributions):
        frac = holdout_fraction_by_file.get(original_file)
        if frac is None or frac <= 0:
            raise ValueError(f"Frazione holdout non disponibile per {original_file}")

        n_gen = n_eventi_generati(original_file)
        weight = sigma * LUMI_PB / (n_gen * frac)
        h_part = ROOT.TH1F(
            f"{hist_name}_part{i_contrib}",
            ";BDT score;Eventi attesi / bin",
            BINS_BDT, -1.0, 1.0
        )
        h_part.Sumw2()

        f_in = ROOT.TFile.Open(prepared_file, "READ")
        if not f_in or f_in.IsZombie():
            raise OSError(f"Impossibile aprire {prepared_file}")
        tree = f_in.Get("T")
        if not tree:
            f_in.Close()
            raise KeyError(f"Tree 'T' non trovato in {prepared_file}")

        is_train = array.array('i', [0])
        tree.SetBranchAddress("lep_pt", v_pt)
        tree.SetBranchAddress("mll", v_mll)
        tree.SetBranchAddress("cos_theta_star", v_cos)
        tree.SetBranchAddress("is_train", is_train)

        n_holdout_selected = 0
        for i in range(tree.GetEntries()):
            tree.GetEntry(i)
            if is_train[0] != 0:
                continue
            score = reader.EvaluateMVA("BDT")
            h_part.Fill(float(score))
            n_holdout_selected += 1

        f_in.Close()
        h_part.Scale(weight)
        hist.Add(h_part)
        print(
            f"{original_file}: N_gen={n_gen}, holdout={frac:.3f}, "
            f"eventi selezionati holdout={n_holdout_selected}"
        )

    hist.SetDirectory(0)
    return hist


# LL e TT: un singolo campione per componente
h_LL = get_bdt_histogram(
    [(output_path("tmva_signal_00.root"), SAMPLE_FILES["decay_00"], SIGMA_LL)], "h_bdt_LL"
)
h_TT = get_bdt_histogram(
    [(output_path("tmva_signal_TT.root"), SAMPLE_FILES["decay_tt"], SIGMA_TT)], "h_bdt_TT"
)

# LT: normalizzazione fisica separata dei due contributi misti
h_LT = get_bdt_histogram(
    [
        (output_path("tmva_signal_0T.root"), SAMPLE_FILES["decay_0t"], SIGMA_0T),
        (output_path("tmva_signal_T0.root"), SAMPLE_FILES["decay_t0"], SIGMA_T0),
    ],
    "h_bdt_LT"
)
h_bkg = get_bdt_histogram(
    [(output_path("tmva_background.root"), SAMPLE_FILES["background"], SIGMA_BKG)], "h_bdt_bkg"
)

# Template per la likelihood: tutti gli istogrammi hanno lo stesso binning.
N_BINS = h_bkg.GetNbinsX()
template_bkg = np.array([h_bkg.GetBinContent(i) for i in range(1, N_BINS + 1)], dtype=float)
template_LL = np.array([h_LL.GetBinContent(i) for i in range(1, N_BINS + 1)], dtype=float)
template_TT = np.array([h_TT.GetBinContent(i) for i in range(1, N_BINS + 1)], dtype=float)
template_LT = np.array([h_LT.GetBinContent(i) for i in range(1, N_BINS + 1)], dtype=float)

# Dati Asimov: somma esatta dei template nominali.
data = template_bkg + template_LL + template_TT + template_LT


def negative_log_likelihood(mu_00, mu_tt, mu_mix):
    expected = template_bkg + mu_00 * template_LL + mu_tt * template_TT + mu_mix * template_LT
    mask = data > 0
    if np.any(expected[mask] <= 0):
        return 1e100
    term = np.zeros_like(data)
    term[mask] = data[mask] * np.log(data[mask] / expected[mask])
    return float(2.0 * np.sum(expected - data + term))


fit_bdt = Minuit(negative_log_likelihood, mu_00=1.0, mu_tt=1.0, mu_mix=1.0)
fit_bdt.errordef = 1.0
for par in ("mu_00", "mu_tt", "mu_mix"):
    fit_bdt.limits[par] = (0.0, 5.0)
fit_bdt.migrad()
fit_bdt.hesse()

print("\n========== FIT LIKELIHOOD — BDT BINARIO ==========")
print(f"μ_LL = {fit_bdt.values['mu_00']:.5f} ± {fit_bdt.errors['mu_00']:.5f}")
print(f"μ_TT = {fit_bdt.values['mu_tt']:.5f} ± {fit_bdt.errors['mu_tt']:.5f}")
print(f"μ_LT = {fit_bdt.values['mu_mix']:.5f} ± {fit_bdt.errors['mu_mix']:.5f}")
print(f"Fit valido: {'sì' if fit_bdt.valid else 'no'}")
print(f"Covarianza accurata: {'sì' if getattr(fit_bdt.fmin, 'has_accurate_covar', False) else 'no'}")
print(f"EDM: {fit_bdt.fmin.edm:.3g}")
print("Correlazioni:")
print(f"  ρ(μ_LL, μ_TT) = {fit_bdt.covariance.correlation()['mu_00', 'mu_tt']:+.4f}")
print(f"  ρ(μ_LL, μ_LT) = {fit_bdt.covariance.correlation()['mu_00', 'mu_mix']:+.4f}")
print(f"  ρ(μ_TT, μ_LT) = {fit_bdt.covariance.correlation()['mu_tt', 'mu_mix']:+.4f}")
