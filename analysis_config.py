"""Shared configuration for the analysis notebooks and scripts.

Input ROOT files are searched in ``INPUT_DIR`` (default: ``<project>/data``).
Generated files are written to ``OUTPUT_DIR`` (default: ``<project>/outputs``).

You can change the directories either by editing the values below or by setting:

    FCC_INPUT_DIR=/path/to/root/files
    FCC_OUTPUT_DIR=/path/to/output

The input file names are analysis defaults. Replace them in the mappings below
when your own ROOT files use different names. No ROOT installation path is
hard-coded: activate a Python environment in which ``import ROOT`` works.
"""

import os
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent


def _directory_from_env(variable: str, default: Path) -> Path:
    value = os.environ.get(variable)
    path = Path(value).expanduser() if value else default
    if not path.is_absolute():
        path = PROJECT_DIR / path
    return path.resolve()


INPUT_DIR = _directory_from_env("FCC_INPUT_DIR", PROJECT_DIR / "data")
OUTPUT_DIR = _directory_from_env("FCC_OUTPUT_DIR", PROJECT_DIR / "outputs")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

EVENTS_TREE_NAME = os.environ.get("FCC_EVENTS_TREE", "Events")
BDT_TREE_NAME = os.environ.get("FCC_BDT_TREE", "TreeS")

# Common settings
N_BINS = int(os.environ.get("FCC_N_BINS", "80"))
LUMI_PB = float(os.environ.get("FCC_LUMI_PB", "19200000"))
COLLISION_ENERGY_GEV = float(os.environ.get("FCC_COLLISION_ENERGY_GEV", "164"))
BINS_BDT = int(os.environ.get("FCC_BINS_BDT", "80"))
REBIN_BDT = int(os.environ.get("FCC_REBIN_BDT", "2"))
W_MASS_GEV = float(os.environ.get("FCC_W_MASS_GEV", "80.4"))

# Leptonic-channel cross-sections in pb.
CROSS_SECTIONS_PB = {
    "00": float(os.environ.get("FCC_SIGMA_00_PB", "0.01276")),
    "TT": float(os.environ.get("FCC_SIGMA_TT_PB", "0.05131")),
    "0T": float(os.environ.get("FCC_SIGMA_0T_PB", "0.01975")),
    "background": float(os.environ.get("FCC_SIGMA_BKG_PB", "0.01676")),
}

# Production-level W-pair cross-sections before W decays, in pb.
PRODUCTION_CROSS_SECTIONS_PB = {
    "00": float(os.environ.get("FCC_SIGMA_PROD_00_PB", "1.292")),
    "TT": float(os.environ.get("FCC_SIGMA_PROD_TT_PB", "4.878")),
    "0T": float(os.environ.get("FCC_SIGMA_PROD_0T_PB", "1.977")),
    "T0": float(os.environ.get("FCC_SIGMA_PROD_T0_PB", "1.977")),
}

# Hadronic-channel cross-sections in pb.
HADRONIC_CROSS_SECTIONS_PB = {
    "00": float(os.environ.get("FCC_HAD_SIGMA_00_PB", "0.4592")),
    "TT": float(os.environ.get("FCC_HAD_SIGMA_TT_PB", "1.845")),
    "0T": float(os.environ.get("FCC_HAD_SIGMA_0T_PB", "0.7111")),
    "T0": float(os.environ.get("FCC_HAD_SIGMA_T0_PB", "0.7137")),
    "background": float(os.environ.get("FCC_HAD_SIGMA_BKG_PB", "0.3005")),
}

# Input files for the leptonic analysis. Edit this mapping for another dataset.
SAMPLE_FILES = {
    "production_00": "both_longitudinal.root",
    "production_tt": "both_transverse.root",
    "production_mixed_1": "mixed.root",
    "production_mixed_2": "mixed2.root",
    "decay_00": "decay_00.root",
    "decay_tt": "decay_TT.root",
    "decay_0t": "decay_0T.root",
    "decay_t0": "decay_T0.root",
    "background": "bckgrnd_dollaro.root",
}

LEPTONIC_SAMPLE_FILES = {
    "decay_00": SAMPLE_FILES["decay_00"],
    "decay_tt": SAMPLE_FILES["decay_tt"],
    "decay_0t": SAMPLE_FILES["decay_0t"],
    "decay_t0": SAMPLE_FILES["decay_t0"],
    "background": SAMPLE_FILES["background"],
}

# Input files for the hadronic analysis.
HADRONIC_SAMPLE_FILES = {
    "00": "ww_hadronic_00.root",
    "TT": "ww_hadronic_TT.root",
    "0T": "ww_hadronic_0T.root",
    "T0": "ww_hadronic_T0.root",
    "background": "bckgrnd_had.root",
}

OBSERVABLE_RANGES = {
    "cos_theta_w_abs": (0.0, 1.0),
    "cos_theta_w": (-1.0, 1.0),
    "cos_theta_star": (-1.0, 1.0),
    "lep_pt": (0.0, 80.0),
    "subleading_lepton_pt": (0.0, 60.0),
    "mll": (0.0, 140.0),
    "mll_components": (0.0, 140.0),
    "mll_signal": (0.0, 130.0),
}


def resolve_input(filename) -> str:
    """Resolve an input file; absolute paths are accepted unchanged."""
    path = Path(filename).expanduser()
    if not path.is_absolute():
        path = INPUT_DIR / path
    return str(path.resolve())


def sample_path(role: str) -> str:
    """Return the resolved path for a leptonic sample role."""
    try:
        return resolve_input(SAMPLE_FILES[role])
    except KeyError as exc:
        raise KeyError(f"Unknown leptonic sample role {role!r}.") from exc


def hadronic_sample_path(role: str) -> str:
    """Return the resolved path for a hadronic sample role."""
    try:
        return resolve_input(HADRONIC_SAMPLE_FILES[role])
    except KeyError as exc:
        raise KeyError(f"Unknown hadronic sample role {role!r}.") from exc


def output_path(filename) -> str:
    """Resolve an output file under OUTPUT_DIR and create its parent folder."""
    path = Path(filename).expanduser()
    if not path.is_absolute():
        path = OUTPUT_DIR / path
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    return str(path)


def output_dir(dirname) -> str:
    """Create and return an output directory under OUTPUT_DIR."""
    path = Path(dirname).expanduser()
    if not path.is_absolute():
        path = OUTPUT_DIR / path
    path = path.resolve()
    path.mkdir(parents=True, exist_ok=True)
    return str(path)


TEMPLATE_FILE = output_path("histos.root")
HADRONIC_TEMPLATE_FILE = output_path("histos_hadronic.root")
PLOTS_DIR = Path(output_dir("plots"))
