"""Independent steady-state cross-checks for the nominal model results."""

from __future__ import annotations

from dataclasses import fields
import json
from pathlib import Path

import numpy as np
from scipy.optimize import root

from cerebral_hemodynamics_aspiration import model
from cerebral_hemodynamics_aspiration.parameters import ModelParameters


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_ICP = {
    "tbi_baseline": 22.698884822524803,
    "dose_pv_240": 19.666231102126197,
    "dose_pv_480": 16.566100463002616,
}


def _load(name: str) -> tuple[dict, ModelParameters]:
    report = json.loads(
        (ROOT / "data" / "reports" / f"{name}.json").read_text(
            encoding="utf-8"
        )
    )
    field_names = {field.name for field in fields(ModelParameters)}
    parameters = ModelParameters(
        **{
            key: value
            for key, value in report["parameters"].items()
            if key in field_names
        }
    )
    return report, parameters


def test_nominal_local_equilibria_match_independent_audit() -> None:
    solved: dict[str, float] = {}
    for name, expected_icp in EXPECTED_ICP.items():
        report, parameters = _load(name)
        flow_ml_s = float(report["flow_ml_min"]) / 60.0
        aspiration = model.Aspiration(report["site"], flow_ml_s)
        returned = float(report["return_fraction"]) * flow_ml_s

        result = root(
            lambda state: model.evaluate(state, parameters, aspiration, returned)[0],
            np.asarray(report["terminal_window_mean_state"], dtype=float),
        )
        assert result.success, result.message
        assert np.max(np.abs(result.fun)) < 1e-10
        assert np.isclose(float(result.x[0]), expected_icp, rtol=0.0, atol=1e-9)
        solved[name] = float(result.x[0])

    assert np.isclose(
        solved["tbi_baseline"] - solved["dose_pv_240"],
        3.0326537202788195,
        atol=1e-9,
    )
    assert np.isclose(
        solved["tbi_baseline"] - solved["dose_pv_480"],
        6.1327843594024,
        atol=1e-9,
    )
