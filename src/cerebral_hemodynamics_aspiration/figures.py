"""Generate manuscript figures and tables from saved study results.

This module reads saved simulations; it does not change or rerun the model.
Every nonzero-flow report must document equal return to the lower SVC.
"""

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


RATES = (60, 120, 180, 240, 300, 360, 480)
SITES = ("Pv", "Pvs", "J3", "J2")
SITE_NAMES = {
    "Pv": "Cerebral veins (Pv)",
    "Pvs": "Venous sinus (Pvs)",
    "J3": "Bilateral J3",
    "J2": "Bilateral J2",
}
COLORS = {
    "Pv": "#0072B2",
    "Pv480": "#56B4E9",
    "Pvs": "#E69F00",
    "J3": "#009E73",
    "J2": "#CC79A7",
    "ink": "#202124",
    "gray": "#666666",
}
MARKERS = {"Pv": "o", "Pvs": "s", "J3": "^", "J2": "D"}
LINESTYLES = {"Pv": "-", "Pvs": "--", "J3": "-.", "J2": ":"}
SENSITIVITY = (
    ("R0", "CSF outflow resistance (R0)", 1800.0, (1350.0, 1800.0, 2250.0)),
    ("Gaut", "Autoregulatory gain (Gaut)", 0.30, (0.15, 0.30, 0.60)),
    ("kE", "Elastance coefficient (kE)", 0.077, (0.0616, 0.077, 0.0924)),
    ("A", "Jugular sigmoid width (A)", 0.80, (0.40, 0.80, 1.60)),
    ("J3_remaining_fraction", "Remaining J3 conductance", 0.15, (0.075, 0.15, 0.30)),
    ("Rvs1_multiplier", "Terminal venous resistance multiplier", 3.0, (1.50, 3.00, 6.00)),
    ("Gc3", "Collateral conductance (Gc3)", 21.43, (0.0, 3.50, 10.715, 21.43, 42.86)),
)


def safe_token(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_style() -> None:
    for path in ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"):
        font = Path("C:/Windows/Fonts") / path
        if font.exists():
            font_manager.fontManager.addfont(font)
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 10, "axes.labelsize": 11,
        "legend.fontsize": 9, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "axes.linewidth": 0.8, "lines.linewidth": 1.8,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "mathtext.fontset": "custom",
        "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
    })


def finish_axes(ax: plt.Axes, *, grid_axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3.5, width=0.8,
                   color=COLORS["ink"])
    ax.grid(axis=grid_axis, color="#D9D9D9", linestyle=(0, (3, 3)),
            linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)


class Source:
    def __init__(self, reports: Path, trajectories: Path):
        self.reports = reports
        self.trajectories = trajectories
        self.used: set[Path] = set()

    def report(self, name: str, *, site: str | None = None, rate: int | None = None,
               intervention: bool = False) -> dict:
        path = self.reports / f"{name}.json"
        report = json.loads(path.read_text(encoding="utf-8"))
        self.used.add(path)
        if not report.get("converged"):
            raise ValueError(f"{path}: simulation did not converge")
        if intervention:
            if report.get("return_fraction") != 1.0:
                raise ValueError(f"{path}: intervention does not use full return")
            if report.get("return_site") != "lower_SVC_state":
                raise ValueError(f"{path}: return is not to lower SVC")
            if report.get("integrator", {}).get("method") != "BDF":
                raise ValueError(f"{path}: primary solver is not BDF")
            if report.get("ramp_s") != 60.0:
                raise ValueError(f"{path}: ramp duration changed")
        if site is not None and report.get("site") != site:
            raise ValueError(f"{path}: expected site {site}")
        if rate is not None and report.get("flow_ml_min") != float(rate):
            raise ValueError(f"{path}: expected rate {rate}")
        return report

    def trajectory(self, name: str) -> tuple[np.ndarray, np.ndarray, list[str]]:
        path = self.trajectories / f"{name}.npz"
        with np.load(path, allow_pickle=False) as data:
            t = data["t"].copy()
            y = data["y"].copy()
            names = data["state_names"].astype(str).tolist()
        self.used.add(path)
        return t, y, names


def icp(report: dict) -> float:
    return float(report["summary"]["icp_mean_mmhg"])


def state(report: dict, index: int) -> float:
    return float(report["terminal_window_mean_state"][index])


def save_figure(fig: plt.Figure, output: Path, stem: str) -> None:
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight",
                metadata={"Creator": "Matplotlib; SVC-return model results"})
    fig.savefig(output / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_reference_hemodynamics(source: Source, output: Path) -> None:
    """Reproduce the published supine reference flows and pressures."""
    report = source.report("reference_baseline")
    flows = report["terminal_flows_ml_s"]
    flow_source = np.array([11.74, 0.79, 0.0, 12.5])
    flow_model = np.array([
        float(flows["Qjr3"]) + float(flows["Qjl3"]),
        float(flows["Qvvr"]) + float(flows["Qvvl"]),
        float(flows["Qc3"]),
        float(flows["Q_cerebral"]),
    ])
    panels = (
        (["Jugular\n(Qj3)", "Vertebral\n(Qvv)", "Collateral\n(Qc3)", "Cerebral\n(Q)"],
         flow_source, flow_model, "Flow (mL/s)"),
        (["ICP", "Venous sinus\npressure"],
         np.array([9.44, 6.00]),
         np.array([icp(report), state(report, 3)]), "Pressure (mmHg)"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.8),
                             gridspec_kw={"width_ratios": [1.7, 1.0]})
    for index, (ax, labels, published, reproduced, ylabel) in enumerate(
        (axes[i], *panels[i]) for i in range(2)
    ):
        x = np.arange(len(labels))
        ax.bar(x, published, width=.64, color="black", edgecolor="black",
               linewidth=1.2, label="Published reference")
        ax.scatter(x, reproduced, s=58, facecolor="white", edgecolor="black",
                   linewidth=1.2, marker="o", zorder=3,
                   label="Model reproduction")
        ax.set_xticks(x, labels)
        ax.set_ylabel(ylabel)
        upper = max(float(np.max(published)), float(np.max(reproduced))) * 1.18
        lower = min(-.35 if index == 0 else 0., float(np.min(reproduced)) * 1.25)
        ax.set_ylim(lower, upper)
        ax.text(-.12, 1.03, chr(ord("A") + index), transform=ax.transAxes,
                fontsize=11, weight="bold")
        finish_axes(ax)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, -.02),
               ncol=2, frameon=False)
    fig.tight_layout(rect=(0, .10, 1, 1), w_pad=2)
    save_figure(fig, output, "figure_s1_reference_hemodynamics")


def plot_icp_response(source: Source, reports: dict, output: Path,
             baseline_icp: float) -> None:
    t_base, y_base, names = source.trajectory("tbi_baseline")
    index = names.index("Pic")
    mask = t_base >= t_base[-1] - 1200.0
    fig, ax = plt.subplots(figsize=(7.1, 4.15))
    ax.plot((t_base[mask] - t_base[-1]) / 60.0, y_base[index, mask],
            color="black", linewidth=2)

    series = (
        ("Pv", 240, "Pv", "-", "o"),
        ("Pv", 480, "Pv480", "--", "s"),
        ("Pvs", 240, "Pvs", "-.", "^"),
    )
    all_icp = [y_base[index, mask]]
    max_time = 0.0
    for site, rate, color_key, linestyle, marker in series:
        name = f"dose_{site.lower()}_{rate}"
        t, y, state_names = source.trajectory(name)
        if state_names.index("Pic") != index:
            raise ValueError(f"{name}: Pic state ordering changed")
        max_time = max(max_time, float(t[-1] / 60))
        all_icp.append(y[index])
        ax.plot(t / 60, y[index], color=COLORS[color_key],
                linestyle=linestyle, marker=marker,
                markevery=max(1, len(t) // 9), markersize=4.5,
                markerfacecolor=COLORS[color_key],
                markeredgecolor=COLORS[color_key],
                label=f"{site}, {rate} mL/min")

    extent = np.concatenate(all_icp)
    span = max(float(np.ptp(extent)), 1)
    ax.set_ylim(float(np.min(extent)) - 0.06 * span,
                float(np.max(extent)) + 0.09 * span)
    ax.set_xlim(-20, max_time + 3)
    ax.axvline(0, color="#777777", linestyle=(0, (4, 3)), linewidth=1.2)
    ax.axhline(baseline_icp, color="#9A9A9A", linestyle=(0, (1, 2)), linewidth=1)
    ax.text(-18, baseline_icp + 0.038 * span,
            f"Elevated baseline ICP: {baseline_icp:.2f} mmHg",
            ha="left", va="bottom", fontsize=9, style="italic",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5})
    ax.set_xlabel("Time relative to aspiration onset (min)")
    ax.set_ylabel("Intracranial pressure (mmHg)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.23),
              ncol=3)
    finish_axes(ax)
    fig.subplots_adjust(bottom=0.28)
    save_figure(fig, output, "figure_2_icp_response")


def plot_venous_pressures(reports: dict, baseline: dict, output: Path) -> None:
    panels = ((2, "Cerebral veins (Pv)"), (3, "Venous sinus (Pvs)"))
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.2))
    x = np.arange(len(SITES))
    for panel_index, (state_index, title) in enumerate(panels):
        ax = axes[panel_index]
        reductions = [state(baseline, state_index)
                      - state(reports[(site, 240)], state_index)
                      for site in SITES]
        bars = ax.bar(x, reductions, width=.66,
                      color=[COLORS[site] for site in SITES],
                      edgecolor=COLORS["ink"], linewidth=.7)
        ax.axhline(0, color=COLORS["ink"], linewidth=.7)
        ax.set_xticks(x, SITES)
        ax.set_title(title)
        ax.set_ylabel("Pressure reduction (mmHg)")
        ax.text(-.18, 1.04, chr(ord("A") + panel_index),
                transform=ax.transAxes, weight="bold", fontsize=10)
        bottom = min(0., min(reductions))
        top = max(0., max(reductions))
        span = max(top - bottom, .2)
        ax.set_ylim(bottom - .08 * span, top + .18 * span)
        for bar, value in zip(bars, reductions):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    value + (.025 * span if value >= 0 else -.025 * span),
                    f"{value:.2f}", ha="center",
                    va="bottom" if value >= 0 else "top", fontsize=8)
        finish_axes(ax)
    fig.tight_layout(w_pad=1.6)
    save_figure(fig, output, "figure_3_venous_pressure")


def plot_parameter_sensitivity(sensitivity: dict, nominal_delta: float,
             output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.1, 4.6))
    all_outcomes = []
    for position, (parameter, label, nominal, tested) in enumerate(SENSITIVITY):
        outcomes = np.array([sensitivity[(parameter, value)] for value in tested])
        all_outcomes.extend(outcomes)
        low, high = float(np.min(outcomes)), float(np.max(outcomes))
        ax.hlines(position, low, high, color="black", linewidth=1.7, zorder=1)
        ax.vlines([low, high], position - .12, position + .12,
                  color="black", linewidth=1.1, zorder=1)
        ax.scatter(outcomes, np.full_like(outcomes, position), s=28,
                   facecolor="white", edgecolor="black", linewidth=1,
                   zorder=2)
        nominal_outcome = sensitivity[(parameter, nominal)]
        if abs(nominal_outcome - nominal_delta) > .001:
            raise ValueError(f"{parameter}: nominal sensitivity result does not match primary")
        ax.scatter(nominal_outcome, position, s=85, marker="*",
                   color="black", edgecolor="black", linewidth=.6,
                   zorder=3)
    ax.axvline(nominal_delta, color="#9A9A9A", linestyle=(0, (3, 3)),
               linewidth=1, zorder=0)
    ax.set_yticks(np.arange(len(SENSITIVITY)), [d[1] for d in SENSITIVITY])
    ax.invert_yaxis()
    xmin, xmax = min(all_outcomes), max(all_outcomes)
    span = max(xmax - xmin, .2)
    ax.set_xlim(xmin - .08 * span, xmax + .08 * span)
    ax.set_xlabel("ICP reduction, ΔICP (mmHg)")
    finish_axes(ax, grid_axis="x")
    ax.grid(axis="y", visible=False)
    ax.scatter([], [], s=28, facecolor="white", edgecolor="black",
               label="Tested values")
    ax.scatter([], [], s=85, marker="*", color="black",
               label="Nominal value")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    save_figure(fig, output, "figure_s2_parameter_sensitivity")


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_convergence_table(reports: dict, baseline: dict, output: Path) -> None:
    """Report attained convergence diagnostics, not just acceptance limits."""
    fields = ["site", "rate_ml_min", "duration_s", "icp_sd_mmhg",
              "max_pressure_drift_mmhg_min", "max_pressure_derivative_mmhg_min"]
    cases = [("Baseline", 0, baseline)] + [
        (site, rate, reports[(site, rate)]) for site in SITES for rate in RATES
    ]
    rows = []
    lines = ["Table S9. Convergence of the post-traumatic baseline and primary aspiration simulations.", "",
             "| Site | Rate (mL/min) | Duration (s) | ICP SD (10⁻⁵ mmHg) | Maximum absolute drift (10⁻⁵ mmHg/min) | Maximum absolute derivative (10⁻⁵ mmHg/min) |",
             "|---|---:|---:|---:|---:|---:|"]
    for site, rate, report in cases:
        stats = report["summary"]
        if stats["window_s"] != 120.0 or not all(
            check["passed"] for check in report["convergence_checks"][-2:]
        ):
            raise ValueError(f"{report['label']}: final-window convergence is inconsistent")
        if report["domain"]["branches_seen"] != [0]:
            raise ValueError(f"{report['label']}: resistance branch changed")
        values = [stats["icp_std_mmhg"], stats["max_pressure_drift_mmhg_min"],
                  stats["max_pressure_residual_mmhg_min"]]
        rows.append(dict(zip(fields, [site, rate, int(report["duration_s"]),
                                      *[f"{value:.12e}" for value in values]])))
        lines.append(f"| {site} | {rate} | {report['duration_s']:.0f} | "
                     + " | ".join(f"{value * 1e5:.3f}" for value in values) + " |")
    lines += ["", "Baseline duration is the equilibration time before aspiration. Intervention durations are measured from aspiration onset and include the 60-s ramp. ICP SD is temporal standard deviation over the final 120 s, not biological variability. Drift is the largest absolute fitted slope among the 13 pressure states over that window; derivative is the largest absolute pressure derivative at the endpoint. All runs met all three convergence criteria in two consecutive 600-s segments. Pv > Pvs and Pv > ICP held throughout every primary aspiration run."]
    write_csv(output / "table_s9_convergence.csv", fields, rows)
    (output / "table_s9_convergence.md").write_bytes(("\n".join(lines) + "\n").encode("utf-8"))


def plot_terminal_resistance(source: Source, output: Path) -> None:
    """Show the saved pressure and resistance trajectories without resimulation."""
    series = (("Pv", 240, "Pv", "-", "o"),
              ("Pv", 480, "Pv480", "--", "s"),
              ("Pvs", 240, "Pvs", "-.", "^"))
    panels = (("Pic", "Intracranial pressure\n(mmHg)"),
              ("Pv", "Cerebral venous pressure\n(mmHg)"),
              ("Pvs", "Venous sinus pressure\n(mmHg)"),
              ("Rvs", "Terminal resistance\n(mmHg s/mL)"))

    def trajectory(label: str) -> tuple[np.ndarray, dict[str, np.ndarray]]:
        t, y, names = source.trajectory(label)
        with np.load(source.trajectories / f"{label}.npz", allow_pickle=False) as data:
            resistance = data["Rvs"].copy()
            branches = data["branch_code"].copy()
        if not np.all(branches == 0):
            raise ValueError(f"{label}: expected the pressure-dependent resistance branch")
        values = {name: y[names.index(name)] for name in ("Pic", "Pv", "Pvs")}
        values["Rvs"] = resistance
        if not np.all(values["Pv"] > values["Pvs"]):
            raise ValueError(f"{label}: branch codes disagree with the saved pressures")
        if not np.all(values["Pv"] > values["Pic"]):
            raise ValueError(f"{label}: nonpositive cerebral venous transmural pressure")
        return t, values

    t_base, base = trajectory("tbi_baseline")
    mask = t_base >= t_base[-1] - 1200.0
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.8), sharex=True)
    for panel, (ax, (variable, ylabel)) in enumerate(zip(axes.flat, panels)):
        ax.plot((t_base[mask] - t_base[-1]) / 60, base[variable][mask],
                color="black", linewidth=2)
        for site, rate, color, style, marker in series:
            t, values = trajectory(f"dose_{site.lower()}_{rate}")
            ax.plot(t / 60, values[variable], color=COLORS[color],
                    linestyle=style, marker=marker, markevery=max(1, len(t) // 9),
                    markersize=3.3, label=f"{site}, {rate} mL/min")
        ax.axvline(0, color="#777777", linestyle=(0, (4, 3)), linewidth=1)
        ax.set_ylabel(ylabel)
        ax.set_xlim(-20, 183)
        ax.set_xticks([0, 60, 120, 180])
        ax.text(-.17, 1.04, chr(ord("A") + panel), transform=ax.transAxes,
                weight="bold", fontsize=10)
        finish_axes(ax)
    for ax in axes[1]:
        ax.set_xlabel("Time from aspiration onset (min)", fontsize=10)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="lower center",
               bbox_to_anchor=(.52, .015), ncol=3)
    fig.subplots_adjust(left=.11, right=.99, bottom=.16, top=.96,
                        hspace=.25, wspace=.38)
    save_figure(fig, output, "figure_s3_terminal_resistance")


def write_primary_table(reports: dict, baseline_icp: float, output: Path) -> None:
    fields = ["site", "delta_icp_240_mmhg", "delta_icp_480_mmhg"]
    rows = []
    md = ["Table 2. ICP reduction at the nominal and maximum aspiration flow rates.", "",
          "| Aspiration site | ΔICP at 240 mL/min (mmHg) | ΔICP at 480 mL/min (mmHg) |",
          "|---|---:|---:|"]
    for site in SITES:
        values = [baseline_icp - icp(reports[(site, rate)]) for rate in (240, 480)]
        rows.append(dict(zip(fields, [site, *[f"{v:.9f}" for v in values]])))
        md.append(f"| {SITE_NAMES[site]} | {values[0]:.3f} | {values[1]:.3f} |")
    md += ["", f"ΔICP is the no-aspiration post-traumatic baseline ICP ({baseline_icp:.3f} mmHg) minus the mean ICP over the final 120 s after convergence. Values are deterministic model outputs."]
    write_csv(output / "table_2_icp_reduction.csv", fields, rows)
    (output / "table_2_icp_reduction.md").write_bytes(("\n".join(md) + "\n").encode("utf-8"))


def write_aspiration_table(reports: dict, baseline_icp: float, output: Path) -> None:
    fields = ["rate_ml_min", "site", "final_icp_mmhg", "delta_icp_mmhg"]
    rows = []
    md = ["Table S7. Final ICP and ICP reduction across aspiration rates and sites.", "",
          "| Rate (mL/min) | Aspiration site | Final ICP (mmHg) | ΔICP (mmHg) |",
          "|---:|---|---:|---:|"]
    for rate in RATES:
        for site in SITES:
            final = icp(reports[(site, rate)])
            delta = baseline_icp - final
            rows.append(dict(zip(fields, [rate, site, f"{final:.9f}", f"{delta:.9f}"])))
            md.append(f"| {rate} | {SITE_NAMES[site]} | {final:.3f} | {delta:.3f} |")
    md += ["", f"The no-aspiration post-traumatic baseline ICP was {baseline_icp:.3f} mmHg. Final ICP is the mean over the final 120 s after convergence; ΔICP equals baseline minus final ICP. Each aspiration flow rate was simulated separately."]
    write_csv(output / "table_s7_aspiration_response.csv", fields, rows)
    (output / "table_s7_aspiration_response.md").write_bytes(("\n".join(md) + "\n").encode("utf-8"))


def write_comparison_table(source: Source, primary: dict, baseline: dict,
                  output: Path) -> None:
    comparisons = (
        ("Post-traumatic primary analysis", "tbi_baseline", "dose"),
        ("Physiological baseline state", "reference_baseline", "reference"),
        ("Fixed terminal resistance", "fixed_rvs_baseline", "fixed_rvs"),
    )
    fields = ["comparison", "baseline_icp_mmhg"] + [f"{site}_delta_icp_mmhg" for site in SITES]
    rows = []
    md = ["Table S8. Model comparisons of ICP reduction at 240 mL/min.", "",
          "| Comparison | Cerebral vein ΔICP | Venous sinus ΔICP | J3 ΔICP | J2 ΔICP |",
          "|---|---:|---:|---:|---:|"]
    for description, baseline_label, prefix in comparisons:
        matching_baseline = baseline if baseline_label == "tbi_baseline" else source.report(baseline_label)
        base_icp = icp(matching_baseline)
        row = {"comparison": description, "baseline_icp_mmhg": f"{base_icp:.9f}"}
        text = [description]
        for site in SITES:
            intervention = (primary[(site, 240)] if prefix == "dose" else
                            source.report(f"{prefix}_{site.lower()}_240",
                                          site=site, rate=240, intervention=True))
            delta = base_icp - icp(intervention)
            row[f"{site}_delta_icp_mmhg"] = f"{delta:.9f}"
            text.append(f"{delta:.3f}")
        rows.append(row)
        md.append("| " + " | ".join(text) + " |")
    md += ["", "Values are ΔICP relative to the corresponding baseline without aspiration. The fixed-resistance comparison has a different baseline equilibrium. It tests sensitivity to the terminal-resistance law and does not simulate a transition into flow limitation."]
    write_csv(output / "table_s8_model_comparison.csv", fields, rows)
    (output / "table_s8_model_comparison.md").write_bytes(
        ("\n".join(md) + "\n").encode("utf-8"))


def read_sensitivity(source: Source) -> dict[tuple[str, float], float]:
    result = {}
    for parameter, _label, _nominal, tested in SENSITIVITY:
        for value in tested:
            prefix = f"sens_{parameter.lower()}_{safe_token(value)}"
            baseline = source.report(prefix + "_baseline")
            intervention = source.report(prefix + "_pv240",
                                         site="Pv", rate=240, intervention=True)
            result[(parameter, value)] = icp(baseline) - icp(intervention)
    if len(result) != 23:
        raise ValueError("Expected 23 matched sensitivity settings")
    return result


def cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    reference_reports = root / "data" / "reports"
    reference_trajectories = root / "data" / "trajectories"
    parser.add_argument("--results-dir", type=Path, default=reference_reports,
                        help="Directory of converged JSON run reports")
    parser.add_argument("--trajectories-dir", type=Path, default=None,
                        help="NPZ trajectories; defaults to data/trajectories for reference reports, or the run directory")
    parser.add_argument("--output-dir", type=Path, default=root / "figures",
                        help="Figure output directory")
    parser.add_argument("--tables-dir", type=Path, default=root / "tables",
                        help="Table output directory")
    args = parser.parse_args(argv)
    reports_dir = args.results_dir.resolve()
    trajectories_dir = (args.trajectories_dir.resolve() if args.trajectories_dir
                        else (reference_trajectories if reports_dir == reference_reports.resolve()
                              else reports_dir))
    figure_dir, table_dir = args.output_dir.resolve(), args.tables_dir.resolve()
    if not reports_dir.is_dir() or not trajectories_dir.is_dir():
        raise FileNotFoundError("The reports and trajectories directories are required")
    for directory in (figure_dir, table_dir, figure_dir / "supplementary", table_dir / "supplementary"):
        directory.mkdir(parents=True, exist_ok=True)
    set_style()
    source = Source(reports_dir, trajectories_dir)
    baseline = source.report("tbi_baseline")
    baseline_icp = icp(baseline)
    primary = {
        (site, rate): source.report(f"dose_{site.lower()}_{rate}", site=site,
                                   rate=rate, intervention=True)
        for site in SITES for rate in RATES
    }
    sensitivity = read_sensitivity(source)
    write_primary_table(primary, baseline_icp, table_dir)
    write_aspiration_table(primary, baseline_icp, table_dir / "supplementary")
    write_comparison_table(source, primary, baseline, table_dir / "supplementary")
    write_convergence_table(primary, baseline, table_dir / "supplementary")
    plot_icp_response(source, primary, figure_dir, baseline_icp)
    plot_venous_pressures(primary, baseline, figure_dir)
    plot_reference_hemodynamics(source, figure_dir / "supplementary")
    nominal_delta = baseline_icp - icp(primary[("Pv", 240)])
    plot_parameter_sensitivity(sensitivity, nominal_delta, figure_dir / "supplementary")
    plot_terminal_resistance(source, figure_dir / "supplementary")

    captions = (
        "Figure 2. ICP response after onset of venous aspiration. "
        f"The black trace shows the final 20 min of the elevated no-aspiration baseline ({baseline_icp:.3f} mmHg). "
        "Aspiration flow increases linearly over 60 s; intervention curves show cerebral-vein aspiration at 240 and 480 mL/min and venous-sinus aspiration at 240 mL/min.\n\n"
        "Figure 3. Venous pressure reductions after convergence at 240 mL/min. "
        "Panels show the reductions from the matched no-aspiration baseline in cerebral venous pressure (A) and venous sinus pressure (B).\n\n"
        "Figure S1. Reproduction of selected published supine hemodynamic values. "
        "Black bars show the published reference values; open circles show the model reproduction. This is a source-reference benchmark, not independent validation of the aspiration intervention.\n\n"
        "Figure S2. One-at-a-time sensitivity of ICP reduction during cerebral-vein aspiration at 240 mL/min. "
        "Open circles show tested parameter settings; horizontal segments span their deterministic outputs, and stars mark nominal settings. "
        f"The dashed line marks the nominal {nominal_delta:.3f}-mmHg reduction. The ranges are not statistical uncertainty intervals.\n\n"
        "Figure S3. Pressures and terminal cerebral-vein resistance before and during aspiration. "
        "Panels show ICP (A), cerebral venous pressure Pv (B), venous sinus pressure Pvs (C), and terminal resistance Rvs (D). "
        "The black traces show the final 20 min of the elevated baseline. Curves show cerebral-vein aspiration at 240 and 480 mL/min and venous-sinus aspiration at 240 mL/min, ending when convergence criteria were met. "
        "Pv and Pvs in the legend identify aspiration sites. The vertical line marks aspiration onset; flow increases over 60 s. "
        "The pressure-dependent branch (Pv > Pvs, with Pv > ICP) remained active throughout the displayed trajectories and all 28 primary aspiration runs. No transition to the other resistance branch occurred.\n"
    )
    (figure_dir / "captions.md").write_bytes(captions.encode("utf-8"))
    def files_in(directory: Path) -> dict[str, str]:
        return {p.relative_to(directory).as_posix(): sha256(p)
                for p in sorted(directory.rglob("*"))
                if p.is_file() and p.name not in {"provenance.json", "SHA256SUMS.txt"}}
    def source_label(directory: Path) -> str:
        try:
            return directory.relative_to(root).as_posix()
        except ValueError:
            return str(directory)
    provenance = {
        "generator": "src/cerebral_hemodynamics_aspiration/figures.py",
        "generator_sha256": sha256(Path(__file__).resolve()),
        "model_run_boundary_condition": "return_fraction=1.0; lower_SVC_state",
        "reports_dir": source_label(reports_dir),
        "trajectories_dir": source_label(trajectories_dir),
        "baseline_icp_mmhg": baseline_icp,
        "nominal_pv240_delta_icp_mmhg": nominal_delta,
        "inputs_sha256": {p.name: sha256(p) for p in sorted(source.used)},
        "figures_sha256": files_in(figure_dir),
        "tables_sha256": files_in(table_dir),
    }
    (figure_dir / "provenance.json").write_bytes((json.dumps(provenance, indent=2) + "\n").encode("utf-8"))
    for directory in (figure_dir, table_dir):
        paths = [p for p in sorted(directory.rglob("*")) if p.is_file() and p.name != "SHA256SUMS.txt"]
        manifest = "\n".join(f"{sha256(p)}  {p.relative_to(directory).as_posix()}" for p in paths) + "\n"
        (directory / "SHA256SUMS.txt").write_bytes(manifest.encode("utf-8"))
    print(f"Wrote Figures 2–3 and S1–S3 to {figure_dir}; Tables 2 and S7–S9 to {table_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
