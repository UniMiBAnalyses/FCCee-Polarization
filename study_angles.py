# ============================================================
# SETUP: libraries, global settings, constants, and functions
# Script version of the first notebook; the analysis logic is preserved.
# ============================================================
import math
import ROOT

from analysis_config import (
    SAMPLE_FILES,
    CROSS_SECTIONS_PB,
    OBSERVABLE_RANGES,
    EVENTS_TREE_NAME,
    N_BINS,
    LUMI_PB,
    COLLISION_ENERGY_GEV,
    resolve_input,
    TEMPLATE_FILE,
)

# --- ROOT global settings ---
ROOT.TH1.AddDirectory(False)    # histograms are not attached to open files
ROOT.gStyle.SetOptStat(0)       # disable the statistics box
ROOT.gStyle.SetPadTickX(1)      # show ticks on the top axis
ROOT.gStyle.SetPadTickY(1)      # show ticks on the right axis
ROOT.gStyle.SetEndErrorSize(0)  # error bars without end caps


# --- Plot styling ---
# Use a consistent colour for each sample across all plots.
COLOR_00 = ROOT.kGreen + 2
COLOR_TT = ROOT.kOrange + 1
COLOR_MIXED = ROOT.kBlue + 1
COLOR_SIGNAL = ROOT.kBlue + 1   # total signal; it is not plotted together with the mixed component
COLOR_BKG = ROOT.kRed + 1

# True and observable variables use colours distinct from the sample colours.
COLOR_TRUE = ROOT.kBlack
COLOR_OBS = ROOT.kMagenta + 1

# Legend labels
LABEL_00 = "Longitudinal (00)"
LABEL_TT = "Transverse (TT)"
LABEL_MIXED = "Mixed (0T + T0)"
LABEL_SIGNAL = "Signal (00 + TT + 0T + T0)"
LABEL_BKG = "Background"

# Axis titles
Y_TITLE = "Events / bin"
TITLE_COS_THETA_W = "cos(#theta_{W})"
TITLE_COS_THETA_STAR = "cos(#theta^{*})"
TITLE_LEP_PT = "Lepton p_{T} [GeV]"
TITLE_MLL = "M_{ll} [GeV]"

FILL_ALPHA = 0.15          # histogram fill transparency
LEGEND_ROW_HEIGHT = 0.06   # legend row height in pad coordinates


# --- File and histogram utilities ---

def open_tree(filename, tree_name=EVENTS_TREE_NAME):
    """Open a ROOT file and return the file handle and tree.

    Keep the file handle alive while the tree is in use.
    """
    path = resolve_input(filename)
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise FileNotFoundError(f"Could not open ROOT file: {path}")
    tree = f.Get(tree_name)
    if not tree:
        raise KeyError(f"Tree '{tree_name}' not found in {path}")
    return f, tree


def make_hist(name, x_title, y_title, x_min, x_max, sumw2=True):
    """Create a histogram with N_BINS bins.

    sumw2=True enables weighted statistical errors, allowing them to be propagated
    correctly when the histogram is scaled or added to another histogram.
    """
    hist = ROOT.TH1F(name, f";{x_title};{y_title}", N_BINS, x_min, x_max)
    if sumw2:
        hist.Sumw2()
    return hist


# --- Histogram filling functions ---

def fill_prod(tree, hist, w_pdg=None):
    """Fill `hist` with cos(theta) of the W flight direction (Pz / |P|).

    w_pdg = None -> W+ and W- [default]
    w_pdg = 24   -> W+ only
    w_pdg = -24  -> W- only
    """
    if not tree:
        return
    for event in tree:
        for p in range(event.nLHEPart):
            pdg = event.LHEPart_pdgId[p]
            if abs(pdg) != 24:
                continue
            if w_pdg is not None and pdg != w_pdg:
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


def fill_cos_theta_star(tree, hist):
    """Fill `hist` with the generator-level (true) cos(theta*) using the W boson.

    It is the angle between the positron momentum in the W+ rest frame and the
    W+ flight direction in the laboratory frame.
    """
    if not tree:
        return
    for event in tree:
        w_vec, lep_vec = None, None

        for p in range(event.nLHEPart):
            pdg = event.LHEPart_pdgId[p]
            if pdg not in (24, -11):
                continue

            v = ROOT.TLorentzVector()
            v.SetPtEtaPhiM(
                event.LHEPart_pt[p],
                event.LHEPart_eta[p],
                event.LHEPart_phi[p],
                event.LHEPart_mass[p],
            )

            if pdg == 24:    # W+ boson
                w_vec = v
            else:            # pdg == -11: positron from the W+ decay
                lep_vec = v

        if w_vec is not None and lep_vec is not None:
            lep_rest = ROOT.TLorentzVector(lep_vec)
            lep_rest.Boost(-w_vec.BoostVector())    # boost to the W rest frame

            p_lep_rest = lep_rest.Vect()
            p_w_lab = w_vec.Vect()

            if p_lep_rest.Mag() > 0 and p_w_lab.Mag() > 0:
                cos_theta = p_lep_rest.Dot(p_w_lab) / (
                    p_lep_rest.Mag() * p_w_lab.Mag()
                )
                hist.Fill(cos_theta)


def get_leptons(event):
    """Return the event's charged leptons (electrons and muons) as TLorentzVectors.

    Incoming beams (pT = 0) are skipped.
    """
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
    """Fill `hist` with one entry per event, calculated from the two leptons.

    variable = "cos_theta_star" -> observable decay-angle cosine
    variable = "lep_pt"         -> leading-lepton pT
    variable = "mll"            -> dilepton invariant mass

    Selection: exactly two charged leptons, both with pT > 10 GeV. The same
    selection is applied to signal and background samples, as in the analysis.
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
            # angle of the first lepton with respect to the beam axis,
            # in the dilepton rest frame
            v_ll = leptons[0] + leptons[1]
            l1_rest = ROOT.TLorentzVector(leptons[0])
            l1_rest.Boost(-v_ll.BoostVector())
            hist.Fill(l1_rest.CosTheta())

        elif variable == "lep_pt":
            hist.Fill(max(leptons[0].Pt(), leptons[1].Pt()))

        elif variable == "mll":
            v_ll = leptons[0] + leptons[1]    # sum of the two four-vectors
            hist.Fill(v_ll.M())


# --- Luminosity normalization functions ---

def event_weight(sigma_pb, n_generated):
    """Weight per generated event: sigma * L / N_generated (expected events per generated event)."""
    return sigma_pb * LUMI_PB / n_generated


def scale_to_luminosity(hist, sigma_pb, n_generated):
    """Normalize `hist` to the expected yield at the integrated luminosity."""
    if n_generated > 0:
        hist.Scale(event_weight(sigma_pb, n_generated))
    else:
        print(f"Warning: no generated events for {hist.GetName()}, normalization skipped.")


def fill_and_scale(trees, hist, variable, sigma_pb):
    """Fill `hist` using all trees in the list and normalize it to the integrated luminosity.

    N_generated is the sum of the entries in all trees in the list.
    """
    n_generated = sum(tree.GetEntries() for tree in trees)
    for tree in trees:
        fill_lepton_variable(tree, hist, variable)
    scale_to_luminosity(hist, sigma_pb, n_generated)


# --- Plotting functions ---

def make_canvas(name, title, width=900, height=650):
    """Create a canvas with the margins used throughout the analysis."""
    c = ROOT.TCanvas(name, title, width, height)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    return c


def draw_fcc_label():
    """Draw the FCC-ee Simulation label, collision energy, and luminosity on the current pad."""
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
    latex.DrawLatex(0.45, 0.92, f"#sqrt{{s}} = {COLLISION_ENERGY_GEV:g} GeV, L = {LUMI_PB / 1e6:g} ab^{{-1}}")


def style_hist(hist, color):
    """Standard style: thick coloured line and transparent fill."""
    hist.SetLineColor(color)
    hist.SetLineStyle(1)
    hist.SetLineWidth(2)
    hist.SetFillColorAlpha(color, FILL_ALPHA)
    hist.SetMarkerColor(color)    # markers, if used, match the line colour


def make_legend(n_rows, header=None):
    """Create a legend without a border or background, in the upper-right area.

    n_rows is the number of rows, including the header.
    """
    legend = ROOT.TLegend(0.25, 0.88 - LEGEND_ROW_HEIGHT * n_rows, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    if header:
        legend.SetHeader(header)
    return legend


def y_axis_max(max_value, n_rows):
    """Set the y-axis maximum, leaving space for a legend with n_rows rows."""
    if max_value <= 0:
        return 1.0
    return max_value * (1.0 + 0.13 * n_rows)


def plot_overlay(name, title, entries, y_title=None, draw_option="HIST",
                 int_format=".3e", header=None):
    """Draw overlaid histograms using the standard analysis style.

    entries:     list of (histogram, colour, label) tuples
    y_title:     y-axis title (the histogram title is kept if None)
    draw_option: "HIST" or "HIST E" (with statistical error bars)
    int_format:  format of the integral displayed in the legend
    header:      optional legend header

    Return (canvas, legend); keep them referenced, otherwise ROOT may delete them.
    """
    canvas = make_canvas(name, title)
    canvas.cd()

    for hist, color, _ in entries:
        style_hist(hist, color)

    # y axis: start at zero and leave space above for the legend
    n_rows = len(entries) + (1 if header else 0)
    max_value = max(hist.GetMaximum() for hist, _, _ in entries)
    first_hist = entries[0][0]
    first_hist.GetYaxis().SetRangeUser(0.0, y_axis_max(max_value, n_rows))
    first_hist.GetYaxis().SetTitleOffset(1.4)
    if y_title:
        first_hist.GetYaxis().SetTitle(y_title)

    for i, (hist, _, _) in enumerate(entries):
        hist.Draw(draw_option if i == 0 else draw_option + " SAME")

    legend = make_legend(n_rows, header)
    for hist, _, label in entries:
        legend.AddEntry(hist, f"{label} [Int: {format(hist.Integral(), int_format)}]", "l")
    legend.Draw()

    draw_fcc_label()
    canvas.Draw()
    return canvas, legend


def plot_stack(name, title, entries, x_title, y_title=Y_TITLE):
    """Draw stacked histograms with a statistical uncertainty band on the total.

    entries: list of (histogram, colour, label) tuples, ordered bottom to top.

    Return (canvas, legend, stack, total); keep references to all four objects
    in the caller. The stack and total are created inside this function; if they are not returned,
    Python may delete them when the function exits, leaving an empty canvas.
    """
    canvas = make_canvas(name, title)
    canvas.cd()

    stack = ROOT.THStack(f"{name}_stack", f";{x_title};{y_title}")
    for hist, color, _ in entries:
        style_hist(hist, color)
        stack.Add(hist)

    # Total histogram, used only to draw the uncertainty band on top of the stack
    total = entries[0][0].Clone(f"{name}_total")
    for hist, _, _ in entries[1:]:
        total.Add(hist)
    total.SetFillStyle(3004)
    total.SetFillColor(ROOT.kGray + 2)
    total.SetLineColor(ROOT.kGray + 2)
    total.SetMarkerSize(0)

    stack.Draw("HIST")
    stack.GetXaxis().SetTitle(x_title)
    stack.GetYaxis().SetTitle(y_title)
    stack.GetYaxis().SetTitleOffset(1.4)

    # y axis: start at zero and leave space above for the legend
    n_rows = len(entries) + 1    # one row per sample plus the uncertainty band
    stack.SetMinimum(0)
    stack.SetMaximum(y_axis_max(stack.GetMaximum(), n_rows))
    total.Draw("E2 SAME")

    legend = make_legend(n_rows)
    for hist, _, label in reversed(entries):    # legend entries are ordered from top to bottom
        legend.AddEntry(hist, f"{label} [Int: {hist.Integral():.3e}]", "lf")
    legend.AddEntry(total, "MC Stat. Uncertainty", "f")
    legend.Draw()

    draw_fcc_label()
    canvas.Draw()
    return canvas, legend, stack, total


def plot_signal_components(variable, x_min, x_max, x_title, name):
    """Compare the 00, TT, and mixed signal components for a dilepton observable.

    Histograms are normalized to the integrated luminosity and include statistical
    uncertainties (Sumw2). The 0T and T0 channels are normalized separately and added.
    """
    h_00 = make_hist(f"h_{name}_dec_00", x_title, Y_TITLE, x_min, x_max)
    h_tt = make_hist(f"h_{name}_dec_tt", x_title, Y_TITLE, x_min, x_max)
    h_0t = make_hist(f"h_{name}_dec_0t", x_title, Y_TITLE, x_min, x_max)
    h_t0 = make_hist(f"h_{name}_dec_t0", x_title, Y_TITLE, x_min, x_max)

    fill_and_scale([tree_dec_00], h_00, variable, CROSS_SECTIONS_PB["00"])
    fill_and_scale([tree_dec_tt], h_tt, variable, CROSS_SECTIONS_PB["TT"])
    fill_and_scale([tree_dec_0t], h_0t, variable, CROSS_SECTIONS_PB["0T"])
    fill_and_scale([tree_dec_t0], h_t0, variable, CROSS_SECTIONS_PB["0T"])

    # mixed = 0T + T0 (uncertainties are added in quadrature)
    h_mixed = h_0t.Clone(f"h_{name}_dec_mixed")
    h_mixed.Add(h_t0)

    canvas, legend = plot_overlay(
        f"c_decay_{name}", f"Decay distribution - {name}",
        [(h_00, COLOR_00, LABEL_00),
         (h_tt, COLOR_TT, LABEL_TT),
         (h_mixed, COLOR_MIXED, LABEL_MIXED)],
        draw_option="HIST E",
    )

    return {"canvas": canvas, "legend": legend,
            "00": h_00, "TT": h_tt, "mixed": h_mixed}


def plot_components_vs_bkg(variable, x_min, x_max, x_title, name):
    """Compare each signal component with the background.

    Histograms are normalized to the integrated luminosity. Sumw2 is disabled, so
    statistical uncertainties are not shown. The mixed component uses the combined
    number of generated events from the 0T and T0 trees, as in the original analysis.
    """
    h_00 = make_hist(f"h_{name}_00", x_title, Y_TITLE, x_min, x_max, sumw2=False)
    h_tt = make_hist(f"h_{name}_tt", x_title, Y_TITLE, x_min, x_max, sumw2=False)
    h_mixed = make_hist(f"h_{name}_mixed", x_title, Y_TITLE, x_min, x_max, sumw2=False)
    h_bkg = make_hist(f"h_{name}_bkg", x_title, Y_TITLE, x_min, x_max, sumw2=False)

    fill_and_scale([tree_dec_00], h_00, variable, CROSS_SECTIONS_PB["00"])
    fill_and_scale([tree_dec_tt], h_tt, variable, CROSS_SECTIONS_PB["TT"])
    fill_and_scale([tree_dec_0t, tree_dec_t0], h_mixed, variable, (2 * CROSS_SECTIONS_PB["0T"]))
    fill_and_scale([tree_bkg], h_bkg, variable, CROSS_SECTIONS_PB["background"])

    canvas, legend = plot_overlay(
        f"c_comp_{name}", f"Signal vs background - {variable}",
        [(h_00, COLOR_00, LABEL_00),
         (h_tt, COLOR_TT, LABEL_TT),
         (h_mixed, COLOR_MIXED, LABEL_MIXED),
         (h_bkg, COLOR_BKG, LABEL_BKG)],
    )

    signal = h_00.Integral() + h_tt.Integral() + h_mixed.Integral()
    print_signal_background(variable, signal, h_bkg.Integral())

    return {"canvas": canvas, "legend": legend,
            "00": h_00, "TT": h_tt, "mixed": h_mixed, "bkg": h_bkg}


def plot_signal_vs_bkg(variable, x_min, x_max, x_title, name):
    """Compare the total signal (00 + TT + 0T + T0) with the background.

    Requires the 00, TT, 0T, T0, and background trees.
    Histograms are normalized to the integrated luminosity and include statistical errors (Sumw2).
    The 0T and T0 channels are normalized separately.

    Return a dictionary containing the canvas, legend, and histograms
    ("signal", "00", "TT", "mixed", "bkg").
    """
    h_00 = make_hist(f"h_{name}_00", x_title, Y_TITLE, x_min, x_max)
    h_tt = make_hist(f"h_{name}_tt", x_title, Y_TITLE, x_min, x_max)
    h_0t = make_hist(f"h_{name}_0t", x_title, Y_TITLE, x_min, x_max)
    h_t0 = make_hist(f"h_{name}_t0", x_title, Y_TITLE, x_min, x_max)
    h_bkg = make_hist(f"h_{name}_bkg", x_title, Y_TITLE, x_min, x_max)

    fill_and_scale([tree_dec_00], h_00, variable, CROSS_SECTIONS_PB["00"])
    fill_and_scale([tree_dec_tt], h_tt, variable, CROSS_SECTIONS_PB["TT"])
    fill_and_scale([tree_dec_0t], h_0t, variable, CROSS_SECTIONS_PB["0T"])
    fill_and_scale([tree_dec_t0], h_t0, variable, CROSS_SECTIONS_PB["0T"])
    fill_and_scale([tree_bkg], h_bkg, variable, CROSS_SECTIONS_PB["background"])

    # mixed = 0T + T0; signal = 00 + TT + mixed (uncertainties are added in quadrature)
    h_mixed = h_0t.Clone(f"h_{name}_mixed")
    h_mixed.Add(h_t0)

    h_signal = h_00.Clone(f"h_{name}_signal")
    h_signal.Add(h_tt)
    h_signal.Add(h_mixed)

    canvas, legend = plot_overlay(
        f"c_sig_{name}", f"Signal vs background - {variable}",
        [(h_signal, COLOR_SIGNAL, LABEL_SIGNAL),
         (h_bkg, COLOR_BKG, LABEL_BKG)],
        draw_option="HIST E",
    )

    print_signal_background(variable, h_signal.Integral(), h_bkg.Integral())

    return {"canvas": canvas, "legend": legend, "signal": h_signal,
            "00": h_00, "TT": h_tt, "mixed": h_mixed, "bkg": h_bkg}


# --- Reporting and tables ---

def print_signal_background(variable, signal, background):
    """Print expected signal and background yields and their ratio."""
    ratio = background / signal if signal > 0 else float("nan")
    print(f"{variable}: signal = {signal:.3e}, background = {background:.3e}, background/signal = {ratio:.4f}")


def error_on_integral(hist):
    """Calculate the integral uncertainty by adding bin uncertainties in quadrature."""
    return math.sqrt(
        sum(hist.GetBinError(i) ** 2 for i in range(1, hist.GetNbinsX() + 1))
    )


def print_yield_table(variable_name, result):
    """Print expected pre-selection and post-selection yields and efficiency per sample.

    `result` is the dictionary returned by `plot_signal_vs_bkg()`.
    """
    rows = [
        ("00", result["00"], CROSS_SECTIONS_PB["00"]),
        ("TT", result["TT"], CROSS_SECTIONS_PB["TT"]),
        ("Mixed (0T + T0)", result["mixed"], (2 * CROSS_SECTIONS_PB["0T"])),
        ("Background", result["bkg"], CROSS_SECTIONS_PB["background"]),
    ]

    print(f"\n\n{'=' * 100}")
    print(f"EFFICIENCY ANALYSIS: {variable_name}")
    print("=" * 100)
    print(
        f"{'Sample':<18}{'sigma [pb]':>13}{'Produced':>16}"
        f"{'After selection':>18}{'MC uncertainty':>15}{'Efficiency':>14}"
    )
    print("-" * 100)

    for name, hist, sigma in rows:
        n_produced = sigma * LUMI_PB    # expected events before selection
        n_selected = hist.Integral()    # expected events after selection, within the histogram range
        error = error_on_integral(hist)
        efficiency = n_selected / n_produced if n_produced > 0 else float("nan")
        print(
            f"{name:<18}{sigma:>13.5f}{n_produced:>16.3e}"
            f"{n_selected:>18.3e}{error:>15.2e}{efficiency:>13.2%}"
        )
    print("-" * 100)


# ============================================================
# STUDY OF THE PRODUCTION ANGLE cos(theta_W)
# 1. Generated events (not normalized)
# 2. Arbitrary units (unit area)
# 3. Normalized to the integrated luminosity
# ============================================================

# --- Open input trees ---
f_prod_00, tree_prod_00 = open_tree(SAMPLE_FILES["production_00"])
f_prod_tt, tree_prod_tt = open_tree(SAMPLE_FILES["production_tt"])
f_prod_m1, tree_prod_m1 = open_tree(SAMPLE_FILES["production_mixed_1"])
f_prod_m2, tree_prod_m2 = open_tree(SAMPLE_FILES["production_mixed_2"])

# --- Raw histograms (generated events) ---
h_prod_00_raw = make_hist("h_prod_00_raw", TITLE_COS_THETA_W, "Entries (Gen)", *OBSERVABLE_RANGES["cos_theta_w"], sumw2=False)
h_prod_tt_raw = make_hist("h_prod_tt_raw", TITLE_COS_THETA_W, "Entries (Gen)", *OBSERVABLE_RANGES["cos_theta_w"], sumw2=False)
h_prod_mixed_raw = make_hist("h_prod_mixed_raw", TITLE_COS_THETA_W, "Entries (Gen)", *OBSERVABLE_RANGES["cos_theta_w"], sumw2=False)

fill_prod(tree_prod_00, h_prod_00_raw)
fill_prod(tree_prod_tt, h_prod_tt_raw)
# the mixed production sample is the sum of the two configured input files
fill_prod(tree_prod_m1, h_prod_mixed_raw)
fill_prod(tree_prod_m2, h_prod_mixed_raw)

# Make copies for the other plots before modifying any histogram.
h_prod_00_au = h_prod_00_raw.Clone("h_prod_00_au")
h_prod_tt_au = h_prod_tt_raw.Clone("h_prod_tt_au")
h_prod_mixed_au = h_prod_mixed_raw.Clone("h_prod_mixed_au")

h_prod_00_lum = h_prod_00_raw.Clone("h_prod_00_lum")
h_prod_tt_lum = h_prod_tt_raw.Clone("h_prod_tt_lum")
h_prod_mixed_lum = h_prod_mixed_raw.Clone("h_prod_mixed_lum")

# --- Plot 1: generated events ---
c_prod_raw, leg_prod_raw = plot_overlay(
    "c_prod_raw", "Production cos(theta_W) - Raw",
    [(h_prod_00_raw, COLOR_00, LABEL_00),
     (h_prod_tt_raw, COLOR_TT, LABEL_TT),
     (h_prod_mixed_raw, COLOR_MIXED, LABEL_MIXED)],
    int_format=".0f",
)

# --- Plot 2: arbitrary units ---
for hist in (h_prod_00_au, h_prod_tt_au, h_prod_mixed_au):
    if hist.Integral() > 0:
        hist.Scale(1.0 / hist.Integral())

c_prod_au, leg_prod_au = plot_overlay(
    "c_prod_au", "Production cos(theta_W) - a.u.",
    [(h_prod_00_au, COLOR_00, LABEL_00),
     (h_prod_tt_au, COLOR_TT, LABEL_TT),
     (h_prod_mixed_au, COLOR_MIXED, LABEL_MIXED)],
    y_title="a.u.",
    int_format=".2f",
)

# --- Plot 3: normalized to the integrated luminosity ---
n_gen_00 = tree_prod_00.GetEntries()
n_gen_tt = tree_prod_tt.GetEntries()
n_gen_mixed = tree_prod_m1.GetEntries() + tree_prod_m2.GetEntries()    # sum of the two mixed input files

scale_to_luminosity(h_prod_00_lum, CROSS_SECTIONS_PB["00"], n_gen_00)
scale_to_luminosity(h_prod_tt_lum, CROSS_SECTIONS_PB["TT"], n_gen_tt)
scale_to_luminosity(h_prod_mixed_lum, (2 * CROSS_SECTIONS_PB["0T"]), n_gen_mixed)

c_prod_lum, leg_prod_lum = plot_overlay(
    "c_prod_lum", "Production cos(theta_W) - Luminosity",
    [(h_prod_00_lum, COLOR_00, LABEL_00),
     (h_prod_tt_lum, COLOR_TT, LABEL_TT),
     (h_prod_mixed_lum, COLOR_MIXED, LABEL_MIXED)],
    y_title=Y_TITLE,
)


# ============================================================
# STUDY OF THE DECAY ANGLE cos(theta*)
# 1. Distributions normalized to the integrated luminosity (with uncertainties)
# 2. Stacked distribution (with uncertainties)
# ============================================================

# --- Open input trees ---
f_dec_00, tree_dec_00 = open_tree(SAMPLE_FILES["decay_00"])
f_dec_tt, tree_dec_tt = open_tree(SAMPLE_FILES["decay_tt"])
f_dec_0t, tree_dec_0t = open_tree(SAMPLE_FILES["decay_0t"])
f_dec_t0, tree_dec_t0 = open_tree(SAMPLE_FILES["decay_t0"])

# --- Raw histograms (generated events) ---
h_dec_00_raw = make_hist("h_dec_00_raw", TITLE_COS_THETA_STAR, Y_TITLE, *OBSERVABLE_RANGES["cos_theta_star"])
h_dec_tt_raw = make_hist("h_dec_tt_raw", TITLE_COS_THETA_STAR, Y_TITLE, *OBSERVABLE_RANGES["cos_theta_star"])
h_dec_0t_raw = make_hist("h_dec_0t_raw", TITLE_COS_THETA_STAR, Y_TITLE, *OBSERVABLE_RANGES["cos_theta_star"])
h_dec_t0_raw = make_hist("h_dec_t0_raw", TITLE_COS_THETA_STAR, Y_TITLE, *OBSERVABLE_RANGES["cos_theta_star"])

fill_cos_theta_star(tree_dec_00, h_dec_00_raw)
fill_cos_theta_star(tree_dec_tt, h_dec_tt_raw)
fill_cos_theta_star(tree_dec_0t, h_dec_0t_raw)
fill_cos_theta_star(tree_dec_t0, h_dec_t0_raw)

# --- Plot 1: normalized to the integrated luminosity ---
h_dec_00_lum = h_dec_00_raw.Clone("h_dec_00_lum")
h_dec_tt_lum = h_dec_tt_raw.Clone("h_dec_tt_lum")
h_dec_0t_lum = h_dec_0t_raw.Clone("h_dec_0t_lum")
h_dec_t0_lum = h_dec_t0_raw.Clone("h_dec_t0_lum")

# Errors are rescaled automatically thanks to Sumw2
scale_to_luminosity(h_dec_00_lum, CROSS_SECTIONS_PB["00"], tree_dec_00.GetEntries())
scale_to_luminosity(h_dec_tt_lum, CROSS_SECTIONS_PB["TT"], tree_dec_tt.GetEntries())
scale_to_luminosity(h_dec_0t_lum, CROSS_SECTIONS_PB["0T"], tree_dec_0t.GetEntries())
scale_to_luminosity(h_dec_t0_lum, CROSS_SECTIONS_PB["0T"], tree_dec_t0.GetEntries())

# Mixed component: sum of the expected yields of both channels (uncertainties add in quadrature)
h_dec_mixed_lum = h_dec_0t_lum.Clone("h_dec_mixed_lum")
h_dec_mixed_lum.Add(h_dec_t0_lum)

c_decay_lum, leg_decay_lum = plot_overlay(
    "c_decay_lum", "Lepton Decay Angle cos(theta*) - Luminosity",
    [(h_dec_00_lum, COLOR_00, LABEL_00),
     (h_dec_tt_lum, COLOR_TT, LABEL_TT),
     (h_dec_mixed_lum, COLOR_MIXED, LABEL_MIXED)],
    draw_option="HIST E",
)

# --- Plot 2: stacked distribution ---
h_dec_00_stack = h_dec_00_lum.Clone("h_dec_00_stack")
h_dec_mixed_stack = h_dec_mixed_lum.Clone("h_dec_mixed_stack")
h_dec_tt_stack = h_dec_tt_lum.Clone("h_dec_tt_stack")

# Order from bottom to top
c_stack, leg_stack, stack_decay, h_total_stack = plot_stack(
    "c_stack", "Lepton Decay Angle cos(theta*) - Stacked",
    [(h_dec_00_stack, COLOR_00, LABEL_00),
     (h_dec_mixed_stack, COLOR_MIXED, LABEL_MIXED),
     (h_dec_tt_stack, COLOR_TT, LABEL_TT)],
    TITLE_COS_THETA_STAR,
)

# ============================================================
# STUDY OF OTHER OBSERVABLES: leading-lepton pT and invariant mass
# Signal split into 00, TT, and mixed components, normalized to luminosity (with uncertainties)
# ============================================================

# Requires the decay trees opened in the previous section
res_dec_pt = plot_signal_components("lep_pt", *OBSERVABLE_RANGES["lep_pt"], TITLE_LEP_PT, "lep_pt")
res_dec_mll = plot_signal_components("mll", *OBSERVABLE_RANGES["mll_signal"], TITLE_MLL, "mll")

# ============================================================
# TRUE vs OBSERVABLE cos(theta*)
# One plot per sample, normalized to the integrated luminosity
# ============================================================

# (sample name, trees, label, cross-section in pb)
samples = [
    ("00", [tree_dec_00], LABEL_00, CROSS_SECTIONS_PB["00"]),
    ("TT", [tree_dec_tt], LABEL_TT, CROSS_SECTIONS_PB["TT"]),
    ("mixed", [tree_dec_0t, tree_dec_t0], LABEL_MIXED, (2 * CROSS_SECTIONS_PB["0T"])),
]

plots_true_obs = []    # keep canvases and histograms alive so ROOT does not delete them

for name, trees, label, sigma in samples:
    h_true = make_hist(f"h_true_{name}", TITLE_COS_THETA_STAR, Y_TITLE, *OBSERVABLE_RANGES["cos_theta_star"])
    h_obs = make_hist(f"h_obs_{name}", TITLE_COS_THETA_STAR, Y_TITLE, *OBSERVABLE_RANGES["cos_theta_star"])

    # number of generated events in this sample
    n_generated = sum(tree.GetEntries() for tree in trees)
    if n_generated <= 0:
        print(f"{name}: no generated events; sample skipped.")
        continue

    for tree in trees:
        fill_cos_theta_star(tree, h_true)                    # uses the W boson from the generator
        fill_lepton_variable(tree, h_obs, "cos_theta_star")  # uses only the two leptons

    # luminosity normalization: sigma * L / N_generated
    weight = event_weight(sigma, n_generated)
    print(f"{name}: event weight = {weight:.2f}")
    h_true.Scale(weight)
    h_obs.Scale(weight)

    # the sample label is used as the legend header
    canvas, legend = plot_overlay(
        f"c_true_obs_{name}", f"cos theta* true vs observable - {name}",
        [(h_true, COLOR_TRUE, "True (W from generator)"),
         (h_obs, COLOR_OBS, "Observable (two leptons)")],
        draw_option="HIST E",
        header=label,
    )
    plots_true_obs.append((canvas, h_true, h_obs, legend))

    integral_true = h_true.Integral()
    integral_obs = h_obs.Integral()
    ratio = integral_obs / integral_true if integral_true > 0 else float("nan")
    print(f"{name}: expected true yield = {integral_true:.3e}, observable yield = {integral_obs:.3e}, ratio = {ratio:.1%}")

# ============================================================
# BACKGROUND STUDY (1): signal components vs background
# One plot per observable, yield/error tables, and template output
# ============================================================

# --- Open the background tree ---
f_bkg, tree_bkg = open_tree(SAMPLE_FILES["background"])

# --- Plots: cos(theta*), leading-lepton pT, and dilepton invariant mass ---
res_comp_cos = plot_components_vs_bkg("cos_theta_star", *OBSERVABLE_RANGES["cos_theta_star"], TITLE_COS_THETA_STAR, "cos_theta_star")
res_comp_pt = plot_components_vs_bkg("lep_pt", *OBSERVABLE_RANGES["lep_pt"], TITLE_LEP_PT, "lep_pt")
res_comp_mll = plot_components_vs_bkg("mll", *OBSERVABLE_RANGES["mll_components"], TITLE_MLL, "mll")

# --- Yield, statistical uncertainty, and efficiency tables ---
print_yield_table("cos(theta*)", res_comp_cos)
print_yield_table("pT leading", res_comp_pt)
print_yield_table("mll", res_comp_mll)

# --- Save histograms to a ROOT file ---
f_out = ROOT.TFile.Open(TEMPLATE_FILE, "RECREATE")
f_out.cd()

for name, result in [("cos_theta_star", res_comp_cos), ("lep_pt", res_comp_pt), ("mll", res_comp_mll)]:
    result["00"].Write(f"h_{name}_00")
    result["TT"].Write(f"h_{name}_tt")
    result["mixed"].Write(f"h_{name}_mixed")
    result["bkg"].Write(f"h_{name}_bkg")

f_out.Close()
print(f"\nAll histograms (cos_theta_star, lep_pt, mll) were saved to '{TEMPLATE_FILE}'.")

# ============================================================
# BACKGROUND STUDY (2): total signal vs background
# One plot per observable, with yield and efficiency tables
# ============================================================

# --- Plots: cos(theta*), leading-lepton pT, and dilepton invariant mass ---
res_sig_cos = plot_signal_vs_bkg("cos_theta_star", *OBSERVABLE_RANGES["cos_theta_star"], TITLE_COS_THETA_STAR, "cos_theta_star")
res_sig_pt = plot_signal_vs_bkg("lep_pt", *OBSERVABLE_RANGES["lep_pt"], TITLE_LEP_PT, "lep_pt")
res_sig_mll = plot_signal_vs_bkg("mll", *OBSERVABLE_RANGES["mll_signal"], TITLE_MLL, "mll")

# --- Yield and efficiency tables ---
print_yield_table("cos(theta*)", res_sig_cos)
print_yield_table("pT leading", res_sig_pt)
print_yield_table("mll", res_sig_mll)

print("\nAnalysis completed for cos(theta*), pT, and mll.")
