"""Visualize saved hemodynamic simulations and sensitivity analyses."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
from matplotlib import font_manager
import matplotlib.pyplot as plt
import numpy as np

RESULTS = Path('.')
FIGURES = Path('figures')

ARIAL_FILES = [
    Path(r"C:\Windows\Fonts\arial.ttf"),
    Path(r"C:\Windows\Fonts\arialbd.ttf"),
    Path(r"C:\Windows\Fonts\ariali.ttf"),
    Path(r"C:\Windows\Fonts\arialbi.ttf"),
]
for font_file in ARIAL_FILES:
    if font_file.exists():
        font_manager.fontManager.addfont(font_file)

# Okabe-Ito colours: restrained, colour-vision-safe, and distinct in print.
# Line styles and markers remain redundant encodings for grayscale reproduction.
COLORS = {
    "Pv": "#0072B2",       # blue
    "Pv480": "#56B4E9",    # sky blue
    "Pvs": "#E69F00",      # orange
    "J3": "#009E73",       # bluish green
    "J2": "#CC79A7",       # reddish purple
    "ink": "#202124",
    "gray": "#666666",
    "light": "#E6E6E6",
}
MARKERS = {"Pv": "o", "Pvs": "s", "J3": "^", "J2": "D"}
LINESTYLES = {"Pv": "-", "Pvs": "--", "J3": "-.", "J2": ":"}
SITE_LABELS = {"Pv": "Cerebral veins (Pv)", "Pvs": "Venous sinus (Pvs)",
               "J3": "Bilateral J3", "J2": "Bilateral J2"}

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 10.0, "axes.labelsize": 11.0,
    "legend.fontsize": 9.5, "xtick.labelsize": 10.0, "ytick.labelsize": 10.0,
    "axes.linewidth": 0.8, "lines.linewidth": 1.8, "pdf.fonttype": 42,
    "ps.fonttype": 42, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "mathtext.fontset": "custom",
    "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})

def read_json(name: str) -> dict:
    return json.loads((RESULTS / f"{name}.json").read_text(encoding="utf-8"))


def read_rows() -> list[dict[str, str]]:
    with (RESULTS / "study_summary.csv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def terminal_state(name: str) -> np.ndarray:
    return np.asarray(read_json(name)["terminal_window_mean_state"], dtype=float)


def load_trajectory(name: str) -> tuple[np.ndarray, np.ndarray, list[str]]:
    with np.load(RESULTS / f"{name}.npz") as data:
        return data["t"].copy(), data["y"].copy(), data["state_names"].astype(str).tolist()


def finish_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3.5, width=0.8, color=COLORS["ink"])
    ax.grid(axis="y", color="#D9D9D9", linestyle=(0, (3, 3)), linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)


def save_figure(fig: plt.Figure, stem: str) -> None:
    pdf = FIGURES / f"{stem}.pdf"
    png = FIGURES / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight", metadata={"Creator": "Matplotlib; cerebral hemodynamics model"})
    fig.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(fig)


def make_timecourse() -> None:
    baseline_report = read_json("tbi_baseline")
    baseline_mean = float(baseline_report["summary"]["icp_mean_mmhg"])
    t_base, y_base, names = load_trajectory("tbi_baseline")
    pic_index = names.index("Pic")
    pre_mask = t_base >= t_base[-1] - 1200.0
    pre_time = (t_base[pre_mask] - t_base[-1]) / 60.0
    pre_icp = y_base[pic_index, pre_mask]

    fig, ax = plt.subplots(figsize=(7.1, 4.15))
    ax.plot(pre_time, pre_icp, color="#000000", linewidth=2.0)

    series = [
        ("dose_pv_240", "Pv, 240 mL/min", "Pv", "-", "o"),
        ("dose_pv_480", "Pv, 480 mL/min", "Pv480", "--", "s"),
        ("dose_pvs_240", "Pvs, 240 mL/min", "Pvs", "-.", "^")
    ]
    for name, label, site, linestyle, marker in series:
        t, y, state_names = load_trajectory(name)
        index = state_names.index("Pic")
        mark_every = max(1, len(t) // 9)
        ax.plot(
            t / 60.0,
            y[index],
            color=COLORS[site],
            linestyle=linestyle,
            marker=marker,
            markevery=mark_every,
            markersize=4.5,
            markerfacecolor=COLORS[site],
            markeredgecolor=COLORS[site],
            label=label,
        )

    ax.axvline(0, color="#777777", linestyle=(0, (4, 3)), linewidth=1.2)
    ax.axhline(baseline_mean, color="#9A9A9A", linestyle=(0, (1, 2)), linewidth=1.0)
    ax.text(
        -18,
        23.22,
        "Elevated baseline ICP: 22.70 mmHg",
        ha="left",
        va="center",
        fontsize=9.2,
        style="italic",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
    )
    ax.annotate(
        "Withdrawal onset",
        xy=(0, 21.35),
        xytext=(16, 21.35),
        ha="left",
        va="center",
        fontsize=9.0,
        color="#555555",
        arrowprops={"arrowstyle": "-|>", "color": "#777777", "lw": 0.9},
    )
    ax.set_xlim(-20, 185)
    ax.set_ylim(14.25, 23.55)
    ax.set_xticks([-20, 0, 30, 60, 90, 120, 150, 180])
    ax.set_xlabel("Simulation time relative to withdrawal onset (min)")
    ax.set_ylabel("Intracranial pressure (mmHg)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.23), ncol=3)
    finish_axes(ax)
    fig.subplots_adjust(bottom=0.28)
    save_figure(fig, "figure_03_icp_time_course")


def make_dose_response(rows: list[dict[str, str]]) -> None:
    fig, ax = plt.subplots(figsize=(7.1, 4.35))
    baseline = float(next(r for r in rows if r["family"] == "baseline")["baseline_icp_mmhg"])
    for site in ("Pv", "Pvs", "J3", "J2"):
        subset = sorted(
            (r for r in rows if r["family"] == "dose_response" and r["site"] == site),
            key=lambda r: float(r["flow_ml_min"]),
        )
        x = [0.0] + [float(r["flow_ml_min"]) for r in subset]
        y = [baseline] + [float(r["outcome_icp_mmhg"]) for r in subset]
        ax.plot(
            x,
            y,
            color=COLORS[site],
            marker=MARKERS[site],
            linestyle=LINESTYLES[site],
            markersize=5.2,
            markerfacecolor=COLORS[site],
            markeredgecolor=COLORS[site],
            markeredgewidth=1.1,
            label=SITE_LABELS[site],
        )
    ax.axhline(baseline, color="#8A8A8A", linestyle=(0, (1, 2)), linewidth=1.0)
    label_box = {"facecolor": "white", "edgecolor": "none", "pad": 1.2}
    ax.text(8, 22.93, "No withdrawal: 22.70 mmHg", ha="left", va="center", fontsize=9.0,
            color="#555555", bbox=label_box)
    ax.text(226, 19.18, "19.63", ha="center", va="top", fontsize=9.2, bbox=label_box)
    ax.text(466, 16.24, "16.49", ha="right", va="bottom", fontsize=9.2, bbox=label_box)
    ax.set_xlabel("Withdrawal rate (mL/min)")
    ax.set_ylabel("Final intracranial pressure (mmHg)")
    ax.set_xticks([0, 60, 120, 180, 240, 300, 360, 480])
    ax.set_xlim(-8, 492)
    ax.set_ylim(15.8, 23.25)
    ax.legend(frameon=False, loc="lower left")
    finish_axes(ax)
    fig.tight_layout()
    save_figure(fig, "figure_04_flow_response")


def make_pressure_coupling() -> None:
    baseline = terminal_state("tbi_baseline")
    sites = ("Pv", "Pvs", "J3", "J2")
    outcomes = {site: terminal_state(f"dose_{site.lower()}_240") for site in sites}
    panels = ((0, "ICP ($P_{ic}$)"), (2, "Cerebral veins ($P_v$)"),
              (3, "Venous sinus ($P_{vs}$)"))
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.15))
    x = np.arange(len(sites))
    for panel_index, (state_index, title) in enumerate(panels):
        values = [baseline[state_index] - outcomes[site][state_index] for site in sites]
        axes[panel_index].bar(
            x, values, color=[COLORS[site] for site in sites], width=.66,
            edgecolor=COLORS["ink"], linewidth=.7,
        )
        axes[panel_index].axhline(0, color=COLORS["ink"], linewidth=.7)
        axes[panel_index].set_xticks(x, sites)
        axes[panel_index].set_title(title)
        axes[panel_index].set_ylabel("Pressure reduction (mmHg)")
        axes[panel_index].text(
            -.18, 1.04, chr(ord("A") + panel_index), transform=axes[panel_index].transAxes,
            weight="bold", fontsize=10,
        )
        ymax = max(values)
        axes[panel_index].set_ylim(0, ymax * 1.18 if ymax > 0 else 1.0)
        for xi, value in zip(x, values):
            axes[panel_index].text(
                xi, value + ymax * .025, f"{value:.2f}", ha="center", va="bottom", fontsize=8,
            )
        finish_axes(axes[panel_index])
    fig.tight_layout(w_pad=1.6)
    save_figure(fig, "figure_05_pressure_coupling")


def make_sensitivity(rows: list[dict[str, str]]) -> None:
    definitions = [
        ("R0", "CSF outflow resistance (R0)", 1800.0),
        ("Gaut", "Autoregulatory gain (Gaut)", 0.30),
        ("kE", "Elastance coefficient (kE)", 0.077),
        ("A", "Jugular sigmoid width (A)", 0.80),
        ("J3_remaining_fraction", "Remaining J3 conductance", 0.15),
        ("Rvs1_multiplier", "Terminal venous resistance multiplier", 3.0),
        ("Gc3", "Collateral conductance (Gc3)", 21.43),
    ]
    nominal_delta = 3.0682335071234057
    fig, ax = plt.subplots(figsize=(7.1, 4.6))
    for position, (parameter, label, nominal_value) in enumerate(definitions):
        subset = [r for r in rows if r["family"] == "matched_sensitivity" and r["parameter"] == parameter]
        subset.sort(key=lambda r: float(r["parameter_value"]))
        parameter_values = np.array([float(r["parameter_value"]) for r in subset])
        outcomes = np.array([float(r["delta_icp_mmhg"]) for r in subset])
        low, high = float(np.min(outcomes)), float(np.max(outcomes))
        ax.hlines(position, low, high, color="#000000", linewidth=1.7, zorder=1)
        ax.vlines([low, high], position - 0.12, position + 0.12, color="#000000", linewidth=1.1, zorder=1)
        ax.scatter(outcomes, np.full_like(outcomes, position), s=28, facecolor="white", edgecolor="#000000", linewidth=1.0, zorder=2)
        nominal_index = int(np.argmin(np.abs(parameter_values - nominal_value)))
        ax.scatter(outcomes[nominal_index], position, s=85, marker="*", color="#000000", edgecolor="#000000", linewidth=0.6, zorder=3)
    ax.axvline(nominal_delta, color="#9A9A9A", linestyle=(0, (3, 3)), linewidth=1.0, zorder=0)
    ax.set_yticks(np.arange(len(definitions)), [d[1] for d in definitions])
    ax.invert_yaxis()
    ax.set_xlim(1.55, 5.15)
    ax.set_xlabel("ICP reduction at Pv 240 mL/min, ΔICP (mmHg)")
    ax.grid(axis="x", color=COLORS["light"], linewidth=0.65)
    ax.grid(axis="y", visible=False)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.scatter([], [], s=28, facecolor="white", edgecolor="#000000", label="Tested values")
    ax.scatter([], [], s=85, marker="*", color="#000000", label="Nominal value")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    save_figure(fig, "figure_06_parameter_sensitivity")


def make_source_benchmark() -> None:
    manifest = json.loads((RESULTS / "study_manifest.json").read_text(encoding="utf-8"))
    ref = manifest["reference_validation"]
    report = read_json("reference_baseline")
    flows = report["terminal_flows_ml_s"]

    flow_labels = ["Jugular\n(Qj3)", "Vertebral\n(Qvv)", "Collateral\n(Qc3)", "Cerebral\n(Q)"]
    flow_source = np.array([11.74, 0.79, 0.0, 12.5])
    flow_model = np.array([
        float(flows["Qjr3"]) + float(flows["Qjl3"]),
        float(flows["Qvvr"]) + float(flows["Qvvl"]),
        float(flows["Qc3"]),
        float(flows["Q_cerebral"]),
    ])
    pressure_labels = ["ICP", "Venous sinus\npressure"]
    pressure_source = np.array([float(ref["targets"]["Pic_mmhg"]), float(ref["targets"]["Pvs_mmhg"])])
    pressure_model = np.array([float(ref["reproduced"]["Pic_mmhg"]), float(ref["reproduced"]["Pvs_mmhg"])])

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.8), gridspec_kw={"width_ratios": [1.7, 1.0]})
    panels = [
        (axes[0], flow_labels, flow_source, flow_model, "Flow (mL/s)"),
        (axes[1], pressure_labels, pressure_source, pressure_model, "Pressure (mmHg)"),
    ]
    for panel_index, (ax, labels, source, model, ylabel) in enumerate(panels):
        x = np.arange(len(labels))
        ax.bar(x, source, width=0.64, color="black", edgecolor="black", linewidth=1.2,
               label="Published reference")
        ax.scatter(x, model, s=58, facecolor="white", edgecolor="black", linewidth=1.2,
                   marker="o", zorder=3,
                   label="Model reproduction")
        ax.set_xticks(x, labels)
        ax.set_ylabel(ylabel)
        upper = max(float(np.max(source)), float(np.max(model))) * 1.18
        lower = min(-0.35 if panel_index == 0 else 0.0, float(np.min(model)) * 1.25)
        ax.set_ylim(lower, upper)
        ax.text(-0.12, 1.03, chr(ord("A") + panel_index), transform=ax.transAxes,
                fontsize=11, weight="bold")
        finish_axes(ax)
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0.10, 1, 1), w_pad=2.0)
    save_figure(fig, "figure_02_reference_reproduction")



def make_all() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    rows = read_rows()
    make_timecourse()
    make_dose_response(rows)
    make_pressure_coupling()
    make_sensitivity(rows)
    make_source_benchmark()
    input_names = [
        "study_summary.csv", "study_manifest.json", "reference_baseline.json",
        "tbi_baseline.json", "tbi_baseline.npz", "dose_pv_240.json", "dose_pv_240.npz",
        "dose_pv_480.npz", "dose_pvs_240.json", "dose_pvs_240.npz",
        "dose_j3_240.json", "dose_j2_240.json",
    ]
    output_files = sorted(
        path for path in FIGURES.iterdir()
        if path.name.startswith("figure_") and path.suffix.lower() in {".pdf", ".png"}
    )
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    provenance = {
        "description": "Figures generated from saved deterministic model outputs",
        "generator": "src/cerebral_hemodynamics_aspiration/visualization.py",
        "generator_sha256": digest(Path(__file__).resolve()),
        "inputs": {name: digest(RESULTS / name) for name in input_names},
        "outputs": {path.name: digest(path) for path in output_files},
    }
    (FIGURES / "figure_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )


def cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default_data = Path(__file__).resolve().parents[2] / "reference_results" / "figure_data"
    parser.add_argument("--results-dir", type=Path, default=default_data,
                        help="Directory containing the saved summary, reports, and selected trajectories")
    parser.add_argument("--output-dir", type=Path, default=Path("figures"),
                        help="Destination for PDF and 600-dpi PNG files")
    args = parser.parse_args(argv)
    global RESULTS, FIGURES
    RESULTS = args.results_dir.resolve()
    FIGURES = args.output_dir.resolve()
    make_all()
    print(f"Wrote figures to {FIGURES}")
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
