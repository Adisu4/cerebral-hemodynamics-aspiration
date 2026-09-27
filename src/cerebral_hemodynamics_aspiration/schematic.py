"""Draw the topology implemented in cerebral_hemodynamics_aspiration.model.

This diagram introduces no numerical model changes. Resistive connections are
unarrowed because their flow direction follows the pressure difference. Only
the one-way CSF terms and prescribed return have arrowheads.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse, FancyArrowPatch

SCHEMATIC_STYLE = {
    "font.family": "Arial", "font.size": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "savefig.facecolor": "white",
}

SCHEMATIC_CAPTION = 'Figure 1. Model structure for the supine simulations. (A) Intracranial circulation, CSF exchange, and autoregulation. (B) Extracranial venous network. Solid boxes identify dynamic states, ovals algebraic pressures, and dashed boxes fixed pressure boundaries. Each bilateral jugular box represents two separate pressure states. The venous sinus shown in both panels is the same compartment. Bold outlines identify the four aspiration sites; an equal flow is returned to the lower SVC. For bilateral jugular aspiration, the total flow is divided equally between sides. Resistance and conductance labels correspond to the Supplementary Material. Autoregulation responds to cerebral arterial inflow and adjusts Cpa and Rpa. Dashed connections in panel A indicate constitutive coupling rather than flow.'


def box(ax, name, xy, text, *, width=1.75, height=.60,
        kind="state", aspiration=False, fontsize=7.7):
    x, y = xy
    style = {"facecolor": "white", "edgecolor": "black",
             "linewidth": 1.75 if aspiration else .85, "zorder": 4}
    if kind == "algebraic":
        patch = Ellipse(xy, width, height, **style)
    else:
        patch = Rectangle((x-width/2, y-height/2), width, height,
                          linestyle=(0, (3, 2)) if kind == "fixed" else "-",
                          **style)
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
            linespacing=1.18, zorder=5)


def wire(ax, points, label=None, label_at=None, *, arrow=False,
         dashed=False, fontsize=7.2, rotation=0):
    style = (0, (3, 2)) if dashed else "-"
    if arrow:
        ax.plot(*zip(*points[:-1]), color="black", lw=.85,
                linestyle=style, zorder=1)
        ax.add_patch(FancyArrowPatch(points[-2], points[-1],
                    arrowstyle="-|>", mutation_scale=8, linewidth=.85,
                    color="black", linestyle=style, zorder=2))
    else:
        ax.plot(*zip(*points), color="black", lw=.85,
                linestyle=style, zorder=1)
    if label:
        ax.text(*label_at, label, ha="center", va="center", fontsize=fontsize,
                rotation=rotation, zorder=3,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.1})


def panel(ax, letter, title, ylim):
    ax.set_xlim(-.2, 12.15)
    ax.set_ylim(*ylim)
    ax.axis("off")
    ax.text(0, 1.025, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="bottom")
    ax.text(.055, 1.025, title, transform=ax.transAxes,
            fontsize=9.5, va="bottom")


def _draw(output: Path) -> None:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.1, 7.85))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.5, 5.35], hspace=.21)
    a, b = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    panel(a, "A", "Intracranial circulation and CSF", (-.15, 5.20))
    panel(b, "B", "Extracranial venous network", (-.25, 7.0))

    # A. The exact inflow resistance is Rla + Rpa/2; there is no independent Rc.
    box(a, "Pa_a", (.85, 4.35), "Pa = 100\nmmHg", width=1.50, kind="fixed")
    box(a, "Ppa", (3.30, 4.35), "Pial arteries\nPpa")
    box(a, "Pc", (5.80, 4.35), "Capillary\nPc", kind="algebraic")
    box(a, "Pv", (8.30, 4.35), "Cerebral veins\nPv", aspiration=True)
    box(a, "Pvs_a", (10.80, 4.35), "Venous sinus\nPvs", aspiration=True)
    wire(a, [(1.60,4.35),(2.425,4.35)], "Rla + Rpa/2", (2.01,4.82), fontsize=6.8)
    wire(a, [(4.175,4.35),(4.925,4.35)], "Rpa/2", (4.55,4.78))
    wire(a, [(6.675,4.35),(7.425,4.35)], "Rpv", (7.05,4.78))
    wire(a, [(9.175,4.35),(9.925,4.35)], "Rvs", (9.55,4.78))
    box(a, "Pic", (5.80, 2.15), "Intracranial pressure\nPic", width=2.10)
    box(a, "xaut", (1.45, 2.15), "Autoregulation\nxaut", width=2.15)
    # One-way CSF production and absorption as implemented by max(...,0).
    wire(a, [(5.80,4.05),(5.80,2.45)], "Qf (Rf)", (5.80,3.28), arrow=True)
    wire(a, [(6.85,2.15),(10.80,2.15),(10.80,4.05)],
         "Q0 (R0)", (8.65,2.15), arrow=True)
    # Constitutive pressure coupling, not an additional flow path.
    wire(a, [(5.0,2.45),(3.30,3.20),(3.30,4.05)],
         "transmural\npressure", (3.86,3.21), dashed=True)
    wire(a, [(6.65,2.45),(8.30,3.20),(8.30,4.05)],
         "transmural\npressure", (7.57,3.21), dashed=True)
    wire(a, [(1.45,2.45),(1.45,3.30),(2.43,3.30),(2.85,4.05)],
         "Cpa, Rpa", (1.47,3.57), dashed=True)
    a.text(5.9, .95,
           "Rvs depends on Pv, Pic and Pvs; Pv and Pic are separate states.",
           ha="center", va="center", fontsize=7.5)
    a.text(5.9, .34,
           "Solid lines: vascular connections     Arrows: CSF flow     Dashed lines: constitutive coupling",
           ha="center", va="center", fontsize=7.0)

    # B. Right/left jugular states are grouped only for diagram economy.
    box(b, "Pvs_b", (4.45, 6.48), "Venous sinus\nPvs (from A)", aspiration=True)
    box(b, "J3", (1.45, 5.04), "Jugular J3\nPjr3, Pjl3", aspiration=True)
    box(b, "C3", (4.45, 5.04), "Collateral C3\nPc3")
    box(b, "Pa_b", (6.63, 5.04), "Pa = 100\nmmHg", width=1.50, kind="fixed")
    box(b, "VV", (8.65, 5.04), "Vertebral veins\nPvv")
    box(b, "Piv", (11.0, 5.04), "Piv\nfixed", width=1.48, kind="fixed")
    box(b, "J2", (1.45, 3.29), "Jugular J2\nPjr2, Pjl2", aspiration=True)
    box(b, "C2", (4.45, 3.29), "Collateral C2\nPc2")
    box(b, "Pcv_c", (6.63, 3.29), "Pcv = 5\nmmHg", width=1.50, kind="fixed")
    box(b, "Azy", (8.65, 3.29), "Azygos\nPazy")
    box(b, "Plv", (11.0, 3.29), "Plv\nfixed", width=1.48, kind="fixed")
    box(b, "Svc1", (1.45, 1.55), "Upper-SVC junction\nPsvc1", width=2.08, kind="algebraic", fontsize=7.1)
    box(b, "Svc", (4.45, 1.55), "Lower SVC\nPsvc")
    box(b, "Pcv_s", (4.45, .10), "Pcv = 5 mmHg", height=.48, kind="fixed")
    wire(b, [(3.575,6.48),(1.45,6.48),(1.45,5.34)],
         "Gjr3, Gjl3", (1.45,5.92))
    wire(b, [(4.45,6.18),(4.45,5.34)], "Gc3", (4.45,5.76))
    wire(b, [(5.325,6.48),(8.65,6.48),(8.65,5.34)],
         "Gvvr + Gvvl", (8.65,5.92))
    wire(b, [(2.325,5.04),(3.575,5.04)], "Gcjr3, Gcjl3", (2.95,5.42), fontsize=6.9)
    wire(b, [(5.325,5.04),(5.88,5.04)], "Gex", (5.63,5.44))
    wire(b, [(9.525,5.04),(10.26,5.04)], "Gvv2", (9.89,5.44))
    wire(b, [(1.45,4.74),(1.45,3.59)], "Gjr2, Gjl2", (1.45,4.165))
    wire(b, [(4.45,4.74),(4.45,3.59)], "Gc2", (4.45,4.165))
    wire(b, [(8.65,4.74),(8.65,3.59)], "Gazy1", (8.65,4.165))
    wire(b, [(2.325,3.29),(3.575,3.29)], "Gcjr2, Gcjl2", (2.95,3.68), fontsize=6.9)
    wire(b, [(5.325,3.29),(5.88,3.29)], "Gc1", (5.63,3.69))
    wire(b, [(9.525,3.29),(10.26,3.29)], "Glv", (9.89,3.69))
    wire(b, [(1.45,2.99),(1.45,1.85)], "Gjr1, Gjl1", (1.45,2.42))
    wire(b, [(2.49,1.55),(3.575,1.55)], "Gsvc1", (3.03,1.93))
    wire(b, [(8.65,2.99),(8.65,2.35),(4.45,2.35),(4.45,1.85)],
         "Gazy2", (6.37,2.35))
    wire(b, [(4.45,1.25),(4.45,.34)], "Gsvc2", (4.45,.79))
    wire(b, [(9.30,.95),(6.37,.95),(6.37,1.55),(5.325,1.55)], arrow=True)
    b.text(9.30,1.36,"Equal returned flow",ha="center",fontsize=7.7)
    b.text(9.30,.52,"Qreturn = Qasp",ha="center",fontsize=8.0)

    fig.text(.50,.036,
             "Solid boxes: dynamic states   Ovals: algebraic pressures   Dashed boxes: fixed pressures",
             ha="center", va="center", fontsize=7.0)
    fig.text(.50,.016,
             "Bold boxes: aspiration sites (one per run); total J3 or J2 aspiration is divided equally between sides.",
             ha="center", va="center", fontsize=7.0)
    fig.subplots_adjust(left=.027,right=.985,top=.953,bottom=.073)
    stem = output / "figure_1_model_schematic"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_model_schematic(output: Path) -> None:
    """Write vector PDF and 600-dpi PNG without modifying plotting defaults."""
    with plt.rc_context(SCHEMATIC_STYLE):
        _draw(output)
