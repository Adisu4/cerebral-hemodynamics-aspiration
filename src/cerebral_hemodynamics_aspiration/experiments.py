"""Run the predefined baseline, aspiration, and sensitivity experiments.

Every aspiration result is paired with a converged baseline generated with
the same parameter set.  The script is resumable; pass ``--force`` to replace
existing files.
"""

from __future__ import annotations

import argparse
import copy
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from .parameters import make_reference_parameters, make_tbi_parameters
from .simulation import (
    DEFAULT_RESULTS_DIR,
    IntegratorConfig,
    integrate,
    load_report,
    source_hashes,
)


DOSE_RATES = (60.0, 120.0, 180.0, 240.0, 300.0, 360.0, 480.0)
SITES = ("Pv", "Pvs", "J3", "J2")


def icp(report: dict) -> float:
    return float(report["summary"]["icp_mean_mmhg"])


def safe_token(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def maybe_run(force: bool, output_dir: Path, *args, **kwargs) -> dict:
    label = kwargs["label"]
    cached_pair = (
        (output_dir / f"{label}.json").exists()
        and (output_dir / f"{label}.npz").exists()
    )
    report = integrate(
        *args,
        output_dir=output_dir,
        overwrite=force,
        **kwargs,
    )
    if not force and cached_pair:
        print(f"{label}: reused", flush=True)
    return report


def sensitivity_parameter(name: str, value: float):
    p = make_tbi_parameters()
    if name == "R0":
        p.R0 = value
    elif name == "Gaut":
        p.Gaut = value
    elif name == "kE":
        p.kE = value
    elif name == "A":
        p.A = value
    elif name == "J3_remaining_fraction":
        p.kjr3 = make_reference_parameters().kjr3 * value
        p.kjl3 = make_reference_parameters().kjl3 * value
    elif name == "Rvs1_multiplier":
        p.Rvs1 = make_reference_parameters().Rvs1 * value
    elif name == "Gc3":
        p.Gc3 = value
    else:
        raise KeyError(name)
    return p


def main(force: bool = False, output_dir: Path = DEFAULT_RESULTS_DIR) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []

    reference_parameters = make_reference_parameters()
    reference = maybe_run(
        force, output_dir, reference_parameters, reference_parameters.initial_state(),
        label="reference_baseline",
    )
    reference_flows = reference["terminal_flows_ml_s"]
    validation = {
        "targets": {
            "Pic_mmhg": 9.44,
            "Pvs_mmhg": 6.00,
            "Qvv_total_ml_s": 0.79,
        },
        "reproduced": {
            "Pic_mmhg": icp(reference),
            "Pvs_mmhg": float(reference["terminal_window_mean_state"][3]),
            "Qvv_total_ml_s": float(reference_flows["Qvvr"]) + float(reference_flows["Qvvl"]),
        },
        "source": "Gadda et al. 2015, Tables 5-6",
    }
    validation["absolute_error"] = {
        key: validation["reproduced"][key] - value
        for key, value in validation["targets"].items()
    }

    reference_state = np.asarray(reference["terminal_window_mean_state"], dtype=float)
    reference_site_reports: dict[str, dict] = {}
    for site in SITES:
        label = f"reference_{site.lower()}_240"
        result = maybe_run(
            force, output_dir, copy.deepcopy(reference_parameters), reference_state.copy(),
            label=label, site=site, flow_ml_min=240.0, return_fraction=1.0,
        )
        reference_site_reports[site] = result
        rows.append({
            "family": "reference_site_comparator", "label": label, "site": site,
            "flow_ml_min": 240.0, "parameter": "reference_state", "parameter_value": "",
            "baseline_icp_mmhg": icp(reference), "outcome_icp_mmhg": icp(result),
            "delta_icp_mmhg": icp(reference) - icp(result),
        })

    tbi_parameters = make_tbi_parameters()
    tbi_baseline = maybe_run(
        force, output_dir, tbi_parameters, tbi_parameters.initial_state(),
        label="tbi_baseline",
    )
    tbi_state = np.asarray(tbi_baseline["terminal_window_mean_state"], dtype=float)
    baseline_icp = icp(tbi_baseline)
    rows.append({
        "family": "baseline", "label": "tbi_baseline", "site": "none",
        "flow_ml_min": 0.0, "parameter": "nominal", "parameter_value": "",
        "baseline_icp_mmhg": baseline_icp, "outcome_icp_mmhg": baseline_icp,
        "delta_icp_mmhg": 0.0,
    })

    control = maybe_run(
        force, output_dir, copy.deepcopy(tbi_parameters), tbi_state.copy(),
        label="control", site="Pv", flow_ml_min=0.0,
    )
    rows.append({
        "family": "control", "label": "control", "site": "none",
        "flow_ml_min": 0.0, "parameter": "nominal", "parameter_value": "",
        "baseline_icp_mmhg": baseline_icp, "outcome_icp_mmhg": icp(control),
        "delta_icp_mmhg": baseline_icp - icp(control),
    })

    dose_reports: dict[tuple[str, float], dict] = {}
    for site in SITES:
        for rate in DOSE_RATES:
            label = f"dose_{site.lower()}_{safe_token(rate)}"
            report = maybe_run(
                force, output_dir, copy.deepcopy(tbi_parameters), tbi_state.copy(),
                label=label, site=site, flow_ml_min=rate, return_fraction=1.0,
            )
            dose_reports[(site, rate)] = report
            rows.append({
                "family": "dose_response", "label": label, "site": site,
                "flow_ml_min": rate, "parameter": "nominal", "parameter_value": "",
                "baseline_icp_mmhg": baseline_icp, "outcome_icp_mmhg": icp(report),
                "delta_icp_mmhg": baseline_icp - icp(report),
            })

    fixed_parameters = make_tbi_parameters()
    fixed_parameters.terminal_resistance_mode = "fixed_Rvs1"
    fixed_baseline = maybe_run(
        force, output_dir, fixed_parameters, fixed_parameters.initial_state(),
        label="fixed_rvs_baseline",
    )
    fixed_state = np.asarray(fixed_baseline["terminal_window_mean_state"], dtype=float)
    rows.append({
        "family": "fixed_resistance_comparator", "label": "fixed_rvs_baseline", "site": "none",
        "flow_ml_min": 0.0, "parameter": "terminal_resistance_mode", "parameter_value": "fixed_Rvs1",
        "baseline_icp_mmhg": icp(fixed_baseline), "outcome_icp_mmhg": icp(fixed_baseline),
        "delta_icp_mmhg": 0.0,
    })
    fixed_site_reports: dict[str, dict] = {}
    for site in SITES:
        label = f"fixed_rvs_{site.lower()}_240"
        result = maybe_run(
            force, output_dir, copy.deepcopy(fixed_parameters), fixed_state.copy(),
            label=label, site=site, flow_ml_min=240.0, return_fraction=1.0,
        )
        fixed_site_reports[site] = result
        rows.append({
            "family": "fixed_resistance_comparator", "label": label, "site": site,
            "flow_ml_min": 240.0, "parameter": "terminal_resistance_mode", "parameter_value": "fixed_Rvs1",
            "baseline_icp_mmhg": icp(fixed_baseline), "outcome_icp_mmhg": icp(result),
            "delta_icp_mmhg": icp(fixed_baseline) - icp(result),
        })

    sensitivity_grid = {
        "R0": (1350.0, 1800.0, 2250.0),
        "Gaut": (0.15, 0.30, 0.60),
        "kE": (0.0616, 0.0770, 0.0924),
        "A": (0.40, 0.80, 1.60),
        "J3_remaining_fraction": (0.075, 0.150, 0.300),
        "Rvs1_multiplier": (1.50, 3.00, 6.00),
        "Gc3": (0.0, 3.50, 10.715, 21.43, 42.86),
    }
    sensitivity_pairs: list[dict] = []
    for parameter_name, values in sensitivity_grid.items():
        for value in values:
            token = safe_token(value)
            p = sensitivity_parameter(parameter_name, value)
            base_label = f"sens_{parameter_name.lower()}_{token}_baseline"
            baseline = maybe_run(
                force, output_dir, copy.deepcopy(p), p.initial_state(), label=base_label,
            )
            matched_state = np.asarray(baseline["terminal_window_mean_state"], dtype=float)
            aspiration_label = f"sens_{parameter_name.lower()}_{token}_pv240"
            aspiration = maybe_run(
                force, output_dir, copy.deepcopy(p), matched_state,
                label=aspiration_label, site="Pv", flow_ml_min=240.0,
                return_fraction=1.0,
            )
            pair = {
                "parameter": parameter_name,
                "parameter_value": value,
                "baseline_label": base_label,
                "aspiration_label": aspiration_label,
                "baseline_icp_mmhg": icp(baseline),
                "outcome_icp_mmhg": icp(aspiration),
                "delta_icp_mmhg": icp(baseline) - icp(aspiration),
                "baseline_domain": baseline["domain"],
                "aspiration_domain": aspiration["domain"],
            }
            sensitivity_pairs.append(pair)
            rows.append({
                "family": "matched_sensitivity", "label": aspiration_label,
                "site": "Pv", "flow_ml_min": 240.0,
                "parameter": parameter_name, "parameter_value": value,
                "baseline_icp_mmhg": pair["baseline_icp_mmhg"],
                "outcome_icp_mmhg": pair["outcome_icp_mmhg"],
                "delta_icp_mmhg": pair["delta_icp_mmhg"],
            })

    strict_config = IntegratorConfig(
        method="Radau", rtol=5e-9, atol=5e-11, max_step_s=1.0,
        save_step_s=1.0, chunk_s=600.0, max_duration_s=18000.0,
    )
    solver_pairs: list[dict] = []
    for site, rate in (("Pv", 240.0), ("Pv", 480.0), ("Pvs", 240.0)):
        label = f"verify_radau_{site.lower()}_{safe_token(rate)}"
        report = maybe_run(
            force, output_dir, copy.deepcopy(tbi_parameters), tbi_state.copy(),
            label=label, site=site, flow_ml_min=rate, return_fraction=1.0,
            config=strict_config,
        )
        primary = dose_reports[(site, rate)]
        solver_pairs.append({
            "site": site,
            "flow_ml_min": rate,
            "bdf_icp_mmhg": icp(primary),
            "radau_icp_mmhg": icp(report),
            "absolute_difference_mmhg": abs(icp(primary) - icp(report)),
        })

    for row in rows:
        report = load_report(row["label"], output_dir)
        row.update({
            "converged": report["converged"],
            "duration_s": report["duration_s"],
            "icp_std_mmhg": report["summary"]["icp_std_mmhg"],
            "icp_drift_mmhg_min": report["summary"]["icp_drift_mmhg_min"],
            "max_pressure_drift_mmhg_min": report["summary"]["max_pressure_drift_mmhg_min"],
            "max_pressure_residual_mmhg_min": report["summary"]["max_pressure_residual_mmhg_min"],
            "branches_seen": ";".join(str(x) for x in report["domain"]["branches_seen"]),
            "minimum_Pv_minus_Pic_mmhg": report["domain"]["minimum_Pv_minus_Pic_mmhg"],
        })

    csv_path = output_dir / "study_summary.csv"
    fieldnames = list(rows[0].keys())
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    key_results = {
        "baseline_icp_mmhg": baseline_icp,
        "pv_240_delta_icp_mmhg": baseline_icp - icp(dose_reports[("Pv", 240.0)]),
        "pv_480_delta_icp_mmhg": baseline_icp - icp(dose_reports[("Pv", 480.0)]),
        "pvs_240_delta_icp_mmhg": baseline_icp - icp(dose_reports[("Pvs", 240.0)]),
        "j3_240_delta_icp_mmhg": baseline_icp - icp(dose_reports[("J3", 240.0)]),
        "j2_240_delta_icp_mmhg": baseline_icp - icp(dose_reports[("J2", 240.0)]),
    }
    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_scope": "supine, computational, TBI-like parameter state",
        "primary_aspiration_boundary": "equal flow returned to lower_SVC_state (return_fraction=1.0)",
        "primary_endpoint": "baseline ICP minus terminal 120-s mean ICP",
        "key_results": key_results,
        "reference_validation": validation,
        "reference_site_comparator": {
            site: {
                "baseline_icp_mmhg": icp(reference),
                "outcome_icp_mmhg": icp(result),
                "delta_icp_mmhg": icp(reference) - icp(result),
            }
            for site, result in reference_site_reports.items()
        },
        "fixed_resistance_comparator": {
            "baseline_icp_mmhg": icp(fixed_baseline),
            "sites": {
                site: {
                    "outcome_icp_mmhg": icp(result),
                    "delta_icp_mmhg": icp(fixed_baseline) - icp(result),
                }
                for site, result in fixed_site_reports.items()
            },
        },
        "solver_verification": solver_pairs,
        "sensitivity_grid": sensitivity_grid,
        "sensitivity_pairs": sensitivity_pairs,
        "all_report_count": len([
            path for path in output_dir.glob("*.json")
            if path.name not in {"study_manifest.json", "verification_report.json"}
        ]),
        "all_converged": all(bool(row["converged"]) for row in rows),
        "source_hashes_sha256": source_hashes(),
        "summary_csv": csv_path.name,
    }
    (output_dir / "study_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(key_results, indent=2), flush=True)


def cli() -> None:
    """Run the predefined experiments from the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force", action="store_true",
        help="Replace existing result files after recomputing each simulation",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_RESULTS_DIR,
        help="Destination for simulation reports and trajectories",
    )
    args = parser.parse_args()
    main(force=args.force, output_dir=args.output_dir)


if __name__ == "__main__":
    cli()
