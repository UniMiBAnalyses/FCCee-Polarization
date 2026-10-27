"""User configuration for input ROOT files and output locations.

Edit SAMPLE_FILES to point to datasets with the required physical meaning and
compatible ROOT tree/branch structure. Relative paths are resolved from
FCC_ANALYSIS_DIR; absolute paths are also accepted.

Environment variables:
  FCC_ANALYSIS_DIR  Directory used for relative input paths (default: current directory)
  FCC_OUTPUT_DIR    Directory for generated files (default: current directory)
  FCC_EVENTS_TREE   Tree containing generated events (default: Events)
  FCC_BDT_TREE      Preferred tree for BDT input files (default: TreeS)
"""

import os
from pathlib import Path

BASE_DIR = Path(os.environ.get("FCC_ANALYSIS_DIR", ".")).expanduser().resolve()
OUTPUT_DIR = Path(os.environ.get("FCC_OUTPUT_DIR", ".")).expanduser().resolve()
EVENTS_TREE_NAME = os.environ.get("FCC_EVENTS_TREE", "Events")
BDT_TREE_NAME = os.environ.get("FCC_BDT_TREE", "TreeS")

# Shared analysis settings. Adapt these values to the dataset and analysis setup.
N_BINS = int(os.environ.get("FCC_N_BINS", "80"))
LUMI_PB = float(os.environ.get("FCC_LUMI_PB", "19200000"))
COLLISION_ENERGY_GEV = float(os.environ.get("FCC_COLLISION_ENERGY_GEV", "164"))

# Cross-sections in pb, keyed by physics component rather than file name.
# As in the original analysis, both mixed channels use sigma_0T and their total is 2 * sigma_0T.
CROSS_SECTIONS_PB = {
    "00": 0.01276,
    "TT": 0.05131,
    "0T": 0.01975,
    "background": 0.01676,
}

# Histogram ranges can be adapted without editing the analysis workflow.
# The two mll ranges remain distinct because the original analysis used them
# in different plot groups.
OBSERVABLE_RANGES = {
    "cos_theta_w": (-1.0, 1.0),
    "cos_theta_star": (-1.0, 1.0),
    "lep_pt": (0.0, 80.0),
    "mll_components": (0.0, 140.0),
    "mll_signal": (0.0, 130.0),
}

# Map analysis roles to your own file names or absolute file paths.
# The values below preserve the defaults of the original analysis; they can be
# replaced without editing the analysis scripts.
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

TEMPLATE_FILENAME = "histos.root"


def resolve_input(filename):
    """Resolve a relative input path against BASE_DIR; preserve absolute paths."""
    path = Path(filename).expanduser()
    if not path.is_absolute():
        path = BASE_DIR / path
    return str(path)


def sample_path(role):
    """Return the resolved path for a sample role defined in SAMPLE_FILES."""
    try:
        filename = SAMPLE_FILES[role]
    except KeyError as exc:
        raise KeyError(f"Unknown sample role {role!r}. Configure it in SAMPLE_FILES.") from exc
    return resolve_input(filename)


def output_path(filename):
    """Resolve an output file under OUTPUT_DIR and create its parent directory."""
    path = Path(filename).expanduser()
    if not path.is_absolute():
        path = OUTPUT_DIR / path
    path.parent.mkdir(parents=True, exist_ok=True)
    return str(path)


TEMPLATE_FILE = output_path(TEMPLATE_FILENAME)
