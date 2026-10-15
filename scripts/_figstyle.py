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
    _draw_label_chips(fig)
    fig.savefig(f"{path_noext}.svg", bbox_inches="tight")
    fig.savefig(f"{path_noext}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


_PENDING_CHIPS = []


def label_chips(ax, ys, colors, w_pt=4.5, gap_pt=2.5):
    """Queue a colour chip just LEFT of each y-tick label. Positions are measured from the
    rendered labels at save time, so chips never overlap names whatever the panel width."""
    _PENDING_CHIPS.append((ax, list(ys), list(colors), w_pt, gap_pt))


def _draw_label_chips(fig):
    from matplotlib.patches import Rectangle
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for ax, ys, colors, w_pt, gap_pt in [c for c in _PENDING_CHIPS if c[0].figure is fig]:
        labs = {round(t.get_position()[1], 6): t for t in ax.get_yticklabels() if t.get_text()}
        px = fig.dpi / 72.0
        inv = ax.transAxes.inverted()
        # one aligned column per panel: left of the widest label
        hits = [labs.get(round(float(y), 6)) for y in ys]
        x0s = [t.get_window_extent(r).x0 for t in hits if t is not None]
        if not x0s:
            continue
        x0_disp = min(x0s) - gap_pt * px
        for y, col, t in zip(ys, colors, hits):
            if t is None:
                continue
            xa1 = inv.transform((x0_disp, 0))[0]
            xa0 = inv.transform((x0_disp - w_pt * px, 0))[0]
            ax.add_patch(Rectangle((xa0, y - 0.42), xa1 - xa0, 0.84, color=col,
                                   transform=ax.get_yaxis_transform(), clip_on=False, lw=0))
    _PENDING_CHIPS[:] = [c for c in _PENDING_CHIPS if c[0].figure is not fig]
