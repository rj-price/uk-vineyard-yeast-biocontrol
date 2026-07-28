"""Shared manuscript figure style — one visual system across all main figures.
SVG (editable) + PNG (300 dpi preview) output; DejaVu Sans; recessive axes.
Palette rules follow the dataviz method: categorical hues assigned in fixed order
(species via _palette), sequential = one hue light->dark, pathogens = their own
fixed colourblind-safe triplet (validated, distinct from the species set)."""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ---- global rc ----
def apply():
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",      # keep text as text in the SVG (editable)
        "pdf.fonttype": 42,
        "figure.dpi": 120,
    })

# ---- recurring entities ----
# Pathogens are a recurring entity across figures -> one fixed colour each.
# Validated triplet (ported dataviz checks, light mode): worst adjacent CVD dE
# 13.4 (>=8), worst normal-vision dE 21.6 (>=15), all in the OKLCH lightness band.
# Pathogens are always direct-labelled, so cross-panel overlap with the species
# Okabe-Ito hues is acceptable.
PATHO_ORDER = ["B. cinerea", "F. graminearum", "F. oxysporum"]
PATHO_COL = {"B. cinerea": "#764E9F",       # purple
             "F. graminearum": "#CC9911",   # gold
             "F. oxysporum": "#009988"}      # teal

def patho_it(p):
    return f"$\\it{{{p.replace(' ', chr(92)+' ')}}}$"

# sequential ramp for %MI (single hue, light -> dark teal-blue)
MI_CMAP = LinearSegmentedColormap.from_list(
    "mi_seq", ["#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"])

INK = "#222222"; MUTED = "#666666"; GRID = "#dddddd"

def panel_letter(ax, letter, x=-0.06, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12, fontweight="bold",
            va="bottom", ha="right")

def save(fig, path_noext):
    fig.savefig(f"{path_noext}.svg", bbox_inches="tight")
    fig.savefig(f"{path_noext}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
