"""Build manuscript figures and tables from converged full-return run reports.

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


def figure_2(source: Source, output: Path) -> None:
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
    save_figure(fig, output, "figure_02_reference_reproduction")


def figure_3(source: Source, reports: dict, output: Path,
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
    ax.set_xlabel("Time relative to extraction onset (min)")
    ax.set_ylabel("Intracranial pressure (mmHg)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.23),
              ncol=3)
    finish_axes(ax)
    fig.subplots_adjust(bottom=0.28)
    save_figure(fig, output, "figure_03_icp_time_course_svc_return")


def figure_4(reports: dict, output: Path, baseline_icp: float) -> None:
    fig, ax = plt.subplots(figsize=(7.1, 4.35))
    values = [baseline_icp]
    for site in SITES:
        subset = [reports[(site, rate)] for rate in RATES]
        y = [baseline_icp] + [icp(report) for report in subset]
        values.extend(y)
        ax.plot([0, *RATES], y, color=COLORS[site], marker=MARKERS[site],
                linestyle=LINESTYLES[site], markersize=5.2,
                markerfacecolor=COLORS[site], markeredgecolor=COLORS[site],
                markeredgewidth=1.1, label=SITE_NAMES[site])
    ax.axhline(baseline_icp, color="#8A8A8A", linestyle=(0, (1, 2)),
               linewidth=1, label=f"No extraction ({baseline_icp:.2f} mmHg)")
    span = max(max(values) - min(values), 1)
    ax.set_ylim(min(values) - .08 * span, max(values) + .08 * span)
    ax.set_xlim(-8, 492)
    ax.set_xticks([0, *RATES])
    ax.set_xlabel("Extraction and SVC-return rate (mL/min)")
    ax.set_ylabel("Final intracranial pressure (mmHg)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.23),
              ncol=3)
    finish_axes(ax)
    fig.subplots_adjust(bottom=0.3)
    save_figure(fig, output, "figure_04_flow_response_svc_return")


def figure_5(reports: dict, baseline: dict, output: Path) -> None:
    panels = ((0, "ICP (Pic)"), (2, "Cerebral veins (Pv)"),
              (3, "Venous sinus (Pvs)"))
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.15))
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
    save_figure(fig, output, "figure_05_pressure_coupling_svc_return")


def figure_6(sensitivity: dict, nominal_delta: float,
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
    save_figure(fig, output, "figure_06_parameter_sensitivity_svc_return")


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_table_2(reports: dict, baseline_icp: float, output: Path) -> None:
    fields = ["rate_ml_min"]
    for site in SITES:
        fields.extend((f"{site}_final_icp_mmhg", f"{site}_delta_icp_mmhg"))
    rows = []
    md = ["Table 2. Final ICP and ICP reduction for complete lower-SVC return.",
          f"Matched no-extraction TBI baseline ICP: {baseline_icp:.6f} mmHg.", "",
          "| Rate (mL/min) | Cerebral vein final ICP | ΔICP | Venous sinus final ICP | ΔICP | J3 final ICP | ΔICP | J2 final ICP | ΔICP |",
          "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for rate in RATES:
        row = {"rate_ml_min": rate}
        text = [str(rate)]
        for site in SITES:
            final = icp(reports[(site, rate)])
            delta = baseline_icp - final
            row[f"{site}_final_icp_mmhg"] = f"{final:.9f}"
            row[f"{site}_delta_icp_mmhg"] = f"{delta:.9f}"
            text.extend((f"{final:.3f}", f"{delta:.3f}"))
        rows.append(row)
        md.append("| " + " | ".join(text) + " |")
    md += ["", "Values are deterministic final 120-s mean ICP results. The extraction rate is returned in full to the lower-SVC state. ΔICP equals the matched no-extraction baseline minus final ICP."]
    write_csv(output / "table_02_dose_response_svc_return.csv", fields, rows)
    (output / "table_02_dose_response_svc_return.md").write_bytes(
        ("\n".join(md) + "\n").encode("utf-8"))


def write_table_3(source: Source, primary: dict, baseline: dict,
                  output: Path) -> None:
    comparisons = (
        ("Post-traumatic primary analysis", "tbi_baseline", "dose"),
        ("Physiological baseline state", "reference_baseline", "reference"),
        ("Fixed terminal resistance", "fixed_rvs_baseline", "fixed_rvs"),
    )
    fields = ["comparison", "baseline_icp_mmhg"] + [f"{site}_delta_icp_mmhg" for site in SITES]
    rows = []
    md = ["Table 3. Model comparisons of ICP reduction at 240 mL/min with complete lower-SVC return.", "",
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
    md += ["", "Values are ΔICP relative to each comparison's matched no-extraction baseline. All aspiration cases use complete return to the lower-SVC state. The fixed-resistance comparison has a different baseline equilibrium."]
    write_csv(output / "table_03_model_comparisons_svc_return.csv", fields, rows)
    (output / "table_03_model_comparisons_svc_return.md").write_bytes(
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
    reference_reports = root / "reference_results" / "reports"
    reference_trajectories = root / "reference_results" / "figure_data"
    parser.add_argument("--results-dir", type=Path, default=reference_reports,
                        help="Directory of converged JSON run reports")
    parser.add_argument("--trajectories-dir", type=Path, default=None,
                        help="Directory of NPZ trajectories (defaults to the saved reference trajectories or the results directory)")
    parser.add_argument("--output-dir", type=Path, default=root / "figures",
                        help="Destination for figures and Tables 2-3")
    args = parser.parse_args(argv)
    reports_dir = args.results_dir.resolve()
    trajectories_dir = (args.trajectories_dir.resolve() if args.trajectories_dir
                        else (reference_trajectories if reports_dir == reference_reports.resolve()
                              else reports_dir))
    output_dir = args.output_dir.resolve()
    if not reports_dir.is_dir() or not trajectories_dir.is_dir():
        raise FileNotFoundError("The reports and trajectories directories are required")
    output_dir.mkdir(parents=True, exist_ok=True)
    set_style()
    source = Source(reports_dir, trajectories_dir)
    figure_2(source, output_dir)
    baseline = source.report("tbi_baseline")
    baseline_icp = icp(baseline)
    primary = {
        (site, rate): source.report(f"dose_{site.lower()}_{rate}",
                                    site=site, rate=rate, intervention=True)
        for site in SITES for rate in RATES
    }
    sensitivity = read_sensitivity(source)
    write_table_2(primary, baseline_icp, output_dir)
    write_table_3(source, primary, baseline, output_dir)
    figure_3(source, primary, output_dir, baseline_icp)
    figure_4(primary, output_dir, baseline_icp)
    figure_5(primary, baseline, output_dir)
    nominal_delta = baseline_icp - icp(primary[("Pv", 240)])
    figure_6(sensitivity, nominal_delta, output_dir)

    captions = (
        "Figure 3. ICP response after onset of prescribed venous extraction with complete lower-SVC return. "
        f"The black pre-intervention trace shows the final 20 min of the no-extraction post-traumatic baseline ({baseline_icp:.3f} mmHg). "
        "Extraction and return flows ramp linearly over the first 60 s from t = 0; intervention curves show Pv at 240 and 480 mL/min and Pvs at 240 mL/min.\n\n"
        "Figure 4. Final ICP across extraction locations and rates with complete lower-SVC return. "
        f"The dotted horizontal line marks the no-extraction baseline of {baseline_icp:.2f} mmHg.\n\n"
        "Figure 5. Terminal pressure reductions at 240 mL/min with complete lower-SVC return, relative to the matched no-extraction post-traumatic baseline. "
        "Panels show ICP (Pic), cerebral venous pressure (Pv), and venous sinus pressure (Pvs).\n\n"
        "Figure 6. One-at-a-time sensitivity of ICP reduction during 240 mL/min cerebral-vein extraction with complete lower-SVC return. "
        "Open circles show tested values; horizontal segments span each parameter's predicted reductions; "
        f"stars mark nominal settings; the dashed line marks the nominal {nominal_delta:.3f}-mmHg reduction.\n"
    )
    (output_dir / "figure_captions_svc_return.md").write_bytes(captions.encode("utf-8"))
    inputs = {path.name: sha256(path) for path in sorted(source.used)}
    outputs = {path.name: sha256(path) for path in sorted(output_dir.iterdir())
               if path.is_file() and path.name not in
               {"artifact_provenance.json", "SHA256SUMS.txt"}}
    try:
        results_label = reports_dir.relative_to(root).as_posix()
    except ValueError:
        results_label = str(reports_dir)
    try:
        trajectory_label = trajectories_dir.relative_to(root).as_posix()
    except ValueError:
        trajectory_label = str(trajectories_dir)
    provenance = {
        "generator": "src/cerebral_hemodynamics_aspiration/visualization.py",
        "generator_sha256": sha256(Path(__file__).resolve()),
        "model_run_boundary_condition": "return_fraction=1.0; lower_SVC_state",
        "reports_dir": results_label,
        "trajectories_dir": trajectory_label,
        "baseline_icp_mmhg": baseline_icp,
        "nominal_pv240_delta_icp_mmhg": nominal_delta,
        "inputs_sha256": inputs,
        "outputs_sha256": outputs,
    }
    (output_dir / "artifact_provenance.json").write_bytes(
        (json.dumps(provenance, indent=2) + "\n").encode("utf-8"))
    print(f"Wrote Figures 2-6 and Tables 2-3 to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
