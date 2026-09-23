"""Check saved run reports against the modular equations and summaries."""

from __future__ import annotations

import csv
from dataclasses import fields
import json
from pathlib import Path

import numpy as np

from cerebral_hemodynamics_aspiration import model
from cerebral_hemodynamics_aspiration.parameters import ModelParameters
from cerebral_hemodynamics_aspiration.simulation import source_hashes


ROOT = Path(__file__).resolve().parents[1]
REPORT_DIR = ROOT / "reference_results" / "reports"


def _parameters(report: dict) -> ModelParameters:
    names = {item.name for item in fields(ModelParameters)}
    return ModelParameters(
        **{name: value for name, value in report["parameters"].items() if name in names}
    )


def _reports() -> list[dict]:
    return [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(REPORT_DIR.glob("*.json"))
    ]


def test_all_89_reports_are_converged_and_equation_consistent() -> None:
    reports = _reports()
    assert len(reports) == 89
    current_hashes = source_hashes()

    for report in reports:
        assert report["converged"] is True
        assert report["source_hashes_sha256"] == current_hashes
        if float(report["flow_ml_min"]) > 0:
            assert report["return_fraction"] == 1.0
            assert report["return_site"] == "lower_SVC_state"
        state = np.asarray(report["final_state"], dtype=float)
        flow_ml_s = float(report["flow_ml_min"]) / 60.0
        derivatives, _ = model.evaluate(
            state,
            _parameters(report),
            model.Aspiration(report["site"], flow_ml_s),
            float(report["return_fraction"]) * flow_ml_s,
        )
        residual = float(np.max(np.abs(derivatives[: model.PRESSURE_STATE_COUNT])) * 60.0)
        expected = float(report["summary"]["max_pressure_residual_mmhg_min"])
        assert np.isclose(residual, expected, rtol=1e-9, atol=1e-12), report["label"]
        assert int(report["domain"]["invalid_open_branch_samples"]) == 0


def test_summary_csv_matches_all_referenced_reports() -> None:
    reports = {report["label"]: report for report in _reports()}
    with (ROOT / "reference_results" / "study_summary.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 62
    for row in rows:
        report = reports[row["label"]]
        observed = float(report["summary"]["icp_mean_mmhg"])
        assert np.isclose(observed, float(row["outcome_icp_mmhg"]), atol=1e-12)


def test_key_results_match_saved_manifest() -> None:
    study = json.loads(
        (ROOT / "reference_results" / "study_manifest.json").read_text(encoding="utf-8")
    )
    expected = json.loads(
        (ROOT / "reference_results" / "key_results.json").read_text(encoding="utf-8")
    )
    for key, value in study["key_results"].items():
        assert value == expected[key]
