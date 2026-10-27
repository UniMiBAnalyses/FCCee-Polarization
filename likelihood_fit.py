# %% SETUP E FUNZIONI GRAFICHE
import ROOT

from analysis_config import TEMPLATE_FILE
import numpy as np
from iminuit import Minuit

# --- Stile dei grafici e funzioni di supporto, come nel secondo notebook ---
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)       # no stats box on the plots
ROOT.gStyle.SetPadTickX(1)      # ticks also on the top axis
ROOT.gStyle.SetPadTickY(1)      # ticks also on the right axis
ROOT.gStyle.SetEndErrorSize(0)  # error bars without end caps

# colours of the samples (same as in semplificato.ipynb)
COLORS = {"00": ROOT.kGreen + 2, "TT": ROOT.kOrange + 1, "mixed": ROOT.kBlue + 1, "bkg": ROOT.kRed + 1}


def make_canvas(name, title, width=900, height=650):
    """Create a canvas with the margins used throughout the notebooks."""
    c = ROOT.TCanvas(name, title, width, height)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    return c


def draw_fcc_label():
    """Draw the 'FCC-ee Simulation' label with sqrt(s) and luminosity on the current pad."""
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


def style_hist(hist, color, line_style=1, fill=True):
    """Set the line colour and width of a histogram (plus a light fill of the same colour)."""
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
    hist.SetLineStyle(line_style)
    if fill:
        hist.SetFillColorAlpha(color, 0.15)


def make_legend(x1, y1, x2, y2, header=None):
    """Return a legend without border and background, with an optional header line."""
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    if header:
        leg.SetHeader(header)
    return leg


# %% FIT 1D / 2D / 3D SU TUTTE LE OSSERVABILI
from itertools import combinations

PARAMS = ("mu_00", "mu_tt", "mu_mix")
VARIABILI = ("cos_theta_star", "lep_pt", "mll")

# Nomi fisici usati nelle tabelle, mantenendo i nomi interni del codice.
NOMI_PARAMETRI = {
    "mu_00": "μ_LL",
    "mu_tt": "μ_TT",
    "mu_mix": "μ_mix",
}

NOMI_VARIABILI = {
    "cos_theta_star": r"cos(θ*)",
    "lep_pt": r"pT leptone",
    "mll": r"massa dileptonica mll",
}

try:
    import pandas as pd
    from IPython import get_ipython
    from IPython.display import display, Markdown
    IN_NOTEBOOK = get_ipython() is not None
except ImportError:
    pd = None
    IN_NOTEBOOK = False


def carica_template(var_name, root_file=TEMPLATE_FILE):
    """Legge i template ROOT e restituisce bkg, LL(00), TT e LT(mixed) come array NumPy."""
    f_input = ROOT.TFile.Open(root_file, "READ")
    if not f_input or f_input.IsZombie():
        raise OSError(f"Non riesco ad aprire il file ROOT: {root_file}")

    try:
        nomi = {
            "bkg": f"h_{var_name}_bkg",
            "00": f"h_{var_name}_00",
            "tt": f"h_{var_name}_tt",
            "mixed": f"h_{var_name}_mixed",
        }
        hist = {}
        for chiave, nome in nomi.items():
            h = f_input.Get(nome)
            if not h:
                raise KeyError(f"Istogramma '{nome}' non trovato in {root_file}")
            hist[chiave] = h

        n_bins = hist["bkg"].GetNbinsX()
        for chiave, h in hist.items():
            if h.GetNbinsX() != n_bins:
                raise ValueError(
                    f"Numero di bin non coerente per '{nomi[chiave]}' "
                    f"({h.GetNbinsX()} invece di {n_bins})"
                )

        arrays = {
            chiave: np.array(
                [h.GetBinContent(i) for i in range(1, n_bins + 1)], dtype=float
            )
            for chiave, h in hist.items()
        }
    finally:
        f_input.Close()

    return arrays


def crea_nll(template):
    """Crea la devianza di Poisson (-2 Delta log L) per i template forniti."""
    bkg = template["bkg"]
    sig00 = template["00"]
    sigtt = template["tt"]
    sigmix = template["mixed"]

    # Dati Asimov: in assenza di fluttuazioni, il punto atteso è mu=(1,1,1).
    data = bkg + sig00 + sigtt + sigmix
    mask_data = data > 0

    def negative_log_likelihood(mu_00, mu_tt, mu_mix):
        nu = bkg + mu_00 * sig00 + mu_tt * sigtt + mu_mix * sigmix
        if np.any((nu <= 0) & mask_data):
            return 1e100

        term_data = np.zeros_like(data)
        term_data[mask_data] = data[mask_data] * np.log(data[mask_data] / nu[mask_data])
        return float(2.0 * np.sum(nu - data + term_data))

    return negative_log_likelihood


def esegui_fit_likelihood(var_name, parametri_liberi=PARAMS,
                          root_file=TEMPLATE_FILE, stampa=False):
    """Esegue il fit richiesto; i parametri non liberi sono fissati a 1."""
    parametri_liberi = tuple(parametri_liberi)
    if not parametri_liberi:
        raise ValueError("Deve esserci almeno un parametro libero.")
    if any(p not in PARAMS for p in parametri_liberi):
        raise ValueError(f"Parametri non riconosciuti: {parametri_liberi}")

    nll = crea_nll(carica_template(var_name, root_file))
    m = Minuit(nll, mu_00=1.0, mu_tt=1.0, mu_mix=1.0)
    m.errordef = 1.0
    for p in PARAMS:
        m.limits[p] = (0.0, 5.0)
        m.fixed[p] = p not in parametri_liberi
        if p not in parametri_liberi:
            m.values[p] = 1.0

    m.migrad()
    m.hesse()

    # La stampa dettagliata resta disponibile solo se richiesta esplicitamente.
    if stampa:
        print(f"Fit {len(parametri_liberi)}D — {var_name} — liberi: {', '.join(parametri_liberi)}")
        for p in PARAMS:
            if p in parametri_liberi:
                print(f"  {NOMI_PARAMETRI[p]} = {m.values[p]:.5f} ± {m.errors[p]:.5f}")
            else:
                print(f"  {NOMI_PARAMETRI[p]} = 1 (fissato)")
        print(f"  valido={m.valid}; covarianza accurata={getattr(m.fmin, 'has_accurate_covar', None)}; EDM={m.fmin.edm:.3g}")

    return m


def leggi_correlazione(fit, p1, p2, parametri_liberi):
    """Estrae una correlazione in modo compatibile con diverse versioni di iminuit."""
    corr = fit.covariance.correlation()
    try:
        return float(corr[p1, p2])
    except Exception:
        arr = np.asarray(corr, dtype=float)
        ordine = list(PARAMS) if arr.shape[0] == len(PARAMS) else list(parametri_liberi)
        return float(arr[ordine.index(p1), ordine.index(p2)])


def stampa_tabella_testo(titolo, righe, colonne):
    """Stampa una tabella allineata quando il codice è eseguito fuori da Jupyter."""
    print(f"\n{titolo}")
    if not righe:
        print("(nessun risultato)")
        return
    larghezze = {
        col: max(len(str(col)), *(len(str(r.get(col, ""))) for r in righe))
        for col in colonne
    }
    intestazione = " | ".join(str(col).ljust(larghezze[col]) for col in colonne)
    separatore = "-+-".join("-" * larghezze[col] for col in colonne)
    print(intestazione)
    print(separatore)
    for riga in righe:
        print(" | ".join(str(riga.get(col, "")).ljust(larghezze[col]) for col in colonne))


# Un fit 1D per ciascun parametro, tre fit 2D per ciascuna coppia e un fit 3D,
# ripetuti per tutte e tre le osservabili. I singoli fit non stampano blocchi ripetitivi.
fit_results = {}
righe_riepilogo = []
righe_correlazioni = []

for var_name in VARIABILI:
    combinazioni_fit = [
        combo
        for n_liberi in (1, 2, 3)
        for combo in combinations(PARAMS, n_liberi)
    ]

    for parametri_liberi in combinazioni_fit:
        fit = esegui_fit_likelihood(var_name, parametri_liberi, stampa=False)
        fit_results[(var_name, tuple(parametri_liberi))] = fit

        dimensione = len(parametri_liberi)
        riga = {
            "Fit": f"{dimensione}D",
            "Parametri liberi": ", ".join(NOMI_PARAMETRI[p] for p in parametri_liberi),
            "μ_LL": (f"{fit.values['mu_00']:.5f} ± {fit.errors['mu_00']:.5f}"
                     if "mu_00" in parametri_liberi else "1 (fisso)"),
            "μ_TT": (f"{fit.values['mu_tt']:.5f} ± {fit.errors['mu_tt']:.5f}"
                     if "mu_tt" in parametri_liberi else "1 (fisso)"),
            "μ_LT": (f"{fit.values['mu_mix']:.5f} ± {fit.errors['mu_mix']:.5f}"
                     if "mu_mix" in parametri_liberi else "1 (fisso)"),
            "Valido": "OK" if fit.valid else "NO",
            "Cov. accurata": "sì" if getattr(fit.fmin, "has_accurate_covar", False) else "no",
            "EDM": f"{fit.fmin.edm:.2g}",
            "_variabile": var_name,
        }
        righe_riepilogo.append(riga)

        if dimensione >= 2 and fit.covariance is not None:
            for p1, p2 in combinations(parametri_liberi, 2):
                rho = leggi_correlazione(fit, p1, p2, parametri_liberi)
                righe_correlazioni.append({
                    "Fit": f"{dimensione}D",
                    "Coppia": f"ρ({NOMI_PARAMETRI[p1]}, {NOMI_PARAMETRI[p2]})",
                    "ρ": f"{rho:+.4f}",
                    "_variabile": var_name,
                })

# Manteniamo i nomi usati dal grafico comparativo già presente sotto.
fit_cos = fit_results[("cos_theta_star", PARAMS)]
fit_pt  = fit_results[("lep_pt", PARAMS)]
fit_mll = fit_results[("mll", PARAMS)]

COLONNE_FIT = ["Fit", "Parametri liberi", "μ_LL", "μ_TT", "μ_LT", "Valido", "Cov. accurata", "EDM"]
COLONNE_CORR = ["Fit", "Coppia", "ρ"]

# Output raggruppato per dimensione del fit: prima tutti i 1D, poi i 2D e infine i 3D.
# Le correlazioni vengono mostrate dopo, raggruppate per 2D e 3D.
COLONNE_FIT_ORDINATE = [
    "Variabile", "Parametri liberi", "μ_LL", "μ_TT", "μ_LT",
    "Valido", "Cov. accurata", "EDM"
]
COLONNE_CORR_ORDINATE = ["Variabile", "Coppia", "ρ"]


def righe_per_display(righe, nome_colonna_variabile="Variabile"):
    """Prepara le righe rendendo leggibile il nome dell'osservabile."""
    output = []
    for r in righe:
        item = {k: v for k, v in r.items() if k != "_variabile" and k != "Fit"}
        item[nome_colonna_variabile] = NOMI_VARIABILI[r["_variabile"]]
        output.append(item)
    return output


if IN_NOTEBOOK and pd is not None:
    for dimensione in (1, 2, 3):
        titolo = f"Fit {dimensione}D"
        display(Markdown(f"## {titolo}"))
        subset = [r for r in righe_riepilogo if r["Fit"] == f"{dimensione}D"]
        rows_display = righe_per_display(subset)
        display(pd.DataFrame(rows_display, columns=COLONNE_FIT_ORDINATE))

    display(Markdown("## Correlazioni"))
    for dimensione in (2, 3):
        display(Markdown(f"### Correlazioni dei fit {dimensione}D"))
        subset = [r for r in righe_correlazioni if r["Fit"] == f"{dimensione}D"]
        rows_display = righe_per_display(subset)
        display(pd.DataFrame(rows_display, columns=COLONNE_CORR_ORDINATE))
else:
    for dimensione in (1, 2, 3):
        subset = [r for r in righe_riepilogo if r["Fit"] == f"{dimensione}D"]
        stampa_tabella_testo(
            f"Fit {dimensione}D",
            righe_per_display(subset),
            COLONNE_FIT_ORDINATE,
        )

    print("\nCORRELAZIONI")
    for dimensione in (2, 3):
        subset = [r for r in righe_correlazioni if r["Fit"] == f"{dimensione}D"]
        stampa_tabella_testo(
            f"Correlazioni dei fit {dimensione}D",
            righe_per_display(subset),
            COLONNE_CORR_ORDINATE,
        )

# DataFrame disponibili se in seguito vuoi salvarli o filtrarli:
# riepilogo_fit e riepilogo_correlazioni contengono una riga per fit/correlazione.
riepilogo_fit = pd.DataFrame(righe_riepilogo) if pd is not None else righe_riepilogo
riepilogo_correlazioni = pd.DataFrame(righe_correlazioni) if pd is not None else righe_correlazioni


# %% PLOT COMPARATIVO DELLE INCERTEZZE
# --- Confronto delle incertezze sui parametri mu nei tre fit ---
# (usa gli oggetti fit_cos, fit_pt, fit_mll restituiti dalla cella precedente)

parametri = ["mu_00", "mu_tt", "mu_mix"]
etichette = ["00", "TT", "Mixed"]
fit_variabili = [
    (fit_cos, ROOT.kViolet + 1, "cos #theta^{*}"),
    (fit_pt,  ROOT.kTeal + 3,   "Lepton p_{T}"),
    (fit_mll, ROOT.kGray + 2,   "M_{ll}"),
]

massimo = max(fit.errors[par] for fit, _, _ in fit_variabili for par in parametri)

c_err = make_canvas("c_err", "Incertezze sui parametri mu")
leg_err = make_legend(0.55, 0.72, 0.88, 0.88)
hist_err = []  # tengo da parte gli istogrammi, così non spariscono

for k, (fit, colore, nome) in enumerate(fit_variabili):
    h = ROOT.TH1F(f"h_err_{k}", ";Polarisation component;#sigma(#mu) from the fit",
                  len(parametri), 0, len(parametri))
    for i, par in enumerate(parametri):
        h.SetBinContent(i + 1, fit.errors[par])
        h.GetXaxis().SetBinLabel(i + 1, etichette[i])

    style_hist(h, colore, fill=False)
    h.SetFillColor(colore)             # barre piene
    h.SetBarWidth(0.28)                # tre barre affiancate dentro ogni bin
    h.SetBarOffset(0.06 + 0.30 * k)
    h.SetMinimum(0.0)
    h.SetMaximum(massimo * 1.6)        # spazio per la legenda
    h.GetXaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleOffset(1.4)
    h.Draw("B" if k == 0 else "B SAME")

    leg_err.AddEntry(h, nome, "f")
    hist_err.append(h)

leg_err.Draw()
draw_fcc_label()
c_err.Draw()
