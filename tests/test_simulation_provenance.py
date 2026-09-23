"""Checks for modular-run provenance and cache discrimination."""

from __future__ import annotations

import copy

from cerebral_hemodynamics_aspiration.parameters import make_tbi_parameters
from cerebral_hemodynamics_aspiration.simulation import (
    IntegratorConfig,
    input_fingerprint,
    source_hashes,
)


def test_modular_source_hashes_cover_the_numerical_pipeline() -> None:
    names = set(source_hashes())
    assert names == {
        "src/cerebral_hemodynamics_aspiration/model.py",
        "src/cerebral_hemodynamics_aspiration/parameters.py",
        "src/cerebral_hemodynamics_aspiration/simulation.py",
        "src/cerebral_hemodynamics_aspiration/experiments.py",
    }


def test_input_fingerprint_changes_with_numerical_intent() -> None:
    parameters = make_tbi_parameters()
    state = parameters.initial_state()
    base = input_fingerprint(parameters, state, label="case", flow_ml_min=240.0)
    assert base == input_fingerprint(
        parameters, state, label="case", flow_ml_min=240.0
    )

    changed_parameters = copy.deepcopy(parameters)
    changed_parameters.R0 += 1.0
    assert base != input_fingerprint(
        changed_parameters, state, label="case", flow_ml_min=240.0
    )
    assert base != input_fingerprint(
        parameters, state, label="case", flow_ml_min=480.0
    )
    assert base != input_fingerprint(
        parameters,
        state,
        label="case",
        flow_ml_min=240.0,
        config=IntegratorConfig(rtol=1e-7),
    )
