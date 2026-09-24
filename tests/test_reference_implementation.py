"""Regression checks against the archived result-producing implementation."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from cerebral_hemodynamics_aspiration import model
from cerebral_hemodynamics_aspiration.parameters import (
    ModelParameters,
    make_reference_parameters,
    make_tbi_parameters,
)


ROOT = Path(__file__).resolve().parents[1]
REFERENCE_SOURCE = ROOT / "tests" / "reference"


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


source_parameters = _load_module(
    "publication_parameters", REFERENCE_SOURCE / "parameters.py"
)
source_model = _load_module("publication_model", REFERENCE_SOURCE / "model.py")
source_runner = _load_module("run_publication", REFERENCE_SOURCE / "simulation.py")


def test_archived_source_hashes() -> None:
    manifest = json.loads((REFERENCE_SOURCE / "source_manifest.json").read_text(encoding="utf-8"))
    for filename, expected in manifest["files"].items():
        actual = hashlib.sha256((REFERENCE_SOURCE / filename).read_bytes()).hexdigest()
        assert actual == expected


def test_parameter_sets_match_archived_source_exactly() -> None:
    assert asdict(make_reference_parameters()) == asdict(
        source_parameters.make_reference_parameters()
    )
    assert asdict(make_tbi_parameters()) == asdict(
        source_parameters.make_tbi_parameters()
    )


def test_model_evaluation_matches_archived_source_exactly() -> None:
    report_names = ("tbi_baseline", "dose_pv_240", "dose_pv_480", "dose_pvs_240")
    reference = make_reference_parameters()
    source_reference = source_parameters.make_reference_parameters()
    tbi = make_tbi_parameters()
    source_tbi = source_parameters.make_tbi_parameters()
    cases = [
        (reference.initial_state(), reference, source_reference, "Pv", 0.0, 0.0),
        (tbi.initial_state(), tbi, source_tbi, "Pv", 0.0, 0.0),
    ]
    for name in report_names:
        report = json.loads(
            (ROOT / "data" / "reports" / f"{name}.json").read_text(
                encoding="utf-8"
            )
        )
        field_names = set(ModelParameters.__dataclass_fields__)
        values = {
            key: value for key, value in report["parameters"].items() if key in field_names
        }
        cases.append(
            (
                np.asarray(report["terminal_window_mean_state"], dtype=float),
                ModelParameters(**values),
                source_parameters.PublicationParameters(**values),
                report["site"],
                float(report["flow_ml_min"]) / 60.0,
                float(report["return_fraction"]),
            )
        )

    for state, package_parameters, archived_parameters, site, flow_ml_s, return_fraction in cases:
        package_dy, package_aux = model.evaluate(
            state,
            package_parameters,
            model.Aspiration(site, flow_ml_s),
            return_fraction * flow_ml_s,
        )
        source_dy, source_aux = source_model.evaluate(
            state,
            archived_parameters,
            source_model.Aspiration(site, flow_ml_s),
            return_fraction * flow_ml_s,
        )
        np.testing.assert_array_equal(package_dy, source_dy)
        assert package_aux.keys() == source_aux.keys()
        for key in package_aux:
            if isinstance(package_aux[key], str):
                assert package_aux[key] == source_aux[key]
            else:
                assert package_aux[key] == source_aux[key]


def test_short_trajectory_matches_archived_source_exactly(tmp_path: Path) -> None:
    from cerebral_hemodynamics_aspiration.simulation import (
        IntegratorConfig,
        integrate,
        runtime_environment,
    )

    shared = {
        "method": "BDF",
        "rtol": 1e-8,
        "atol": 1e-10,
        "max_step_s": 2.0,
        "save_step_s": 1.0,
        "chunk_s": 600.0,
        "max_duration_s": 600.0,
        "averaging_window_s": 120.0,
        "icp_std_limit_mmhg": 1e9,
        "pressure_drift_limit_mmhg_min": 1e9,
        "pressure_residual_limit_mmhg_min": 1e9,
        "consecutive_passes": 1,
    }
    package_parameters = make_tbi_parameters()
    archived_parameters = source_parameters.make_tbi_parameters()
    package_config = IntegratorConfig(**shared)
    package_report = integrate(
        package_parameters,
        package_parameters.initial_state(),
        label="short_package",
        site="Pv",
        flow_ml_min=240.0,
        config=package_config,
        output_dir=tmp_path / "package",
    )
    source_report = source_runner.integrate(
        archived_parameters,
        archived_parameters.initial_state(),
        label="short_archived",
        site="Pv",
        flow_ml_min=240.0,
        config=source_runner.IntegratorConfig(**shared),
        output_dir=tmp_path / "archived",
    )
    assert package_report["summary"] == source_report["summary"]
    assert package_report["environment"] == runtime_environment()
    reused = integrate(
        package_parameters,
        package_parameters.initial_state(),
        label="short_package",
        site="Pv",
        flow_ml_min=240.0,
        config=package_config,
        output_dir=tmp_path / "package",
        overwrite=False,
    )
    assert reused["input_fingerprint_sha256"] == package_report["input_fingerprint_sha256"]

    report_path = tmp_path / "package" / "short_package.json"
    cached_report = json.loads(report_path.read_text(encoding="utf-8"))
    nonconverged = dict(cached_report)
    nonconverged["converged"] = False
    report_path.write_text(json.dumps(nonconverged), encoding="utf-8")
    with pytest.raises(RuntimeError, match="cached run did not converge"):
        integrate(
            package_parameters,
            package_parameters.initial_state(),
            label="short_package",
            site="Pv",
            flow_ml_min=240.0,
            config=package_config,
            output_dir=tmp_path / "package",
            overwrite=False,
        )
    report_path.write_text(json.dumps(cached_report), encoding="utf-8")

    with pytest.raises(RuntimeError, match="does not match current inputs/source"):
        integrate(
            package_parameters,
            package_parameters.initial_state(),
            label="short_package",
            site="Pv",
            flow_ml_min=480.0,
            config=package_config,
            output_dir=tmp_path / "package",
            overwrite=False,
        )
    with np.load(tmp_path / "package" / "short_package.npz") as package_data:
        with np.load(tmp_path / "archived" / "short_archived.npz") as source_data:
            assert package_data.files == source_data.files
            for name in package_data.files:
                np.testing.assert_array_equal(package_data[name], source_data[name])

    trajectory = tmp_path / "package" / "short_package.npz"
    trajectory.write_bytes(trajectory.read_bytes() + b"tampered")
    with pytest.raises(RuntimeError, match="cached trajectory hash mismatch"):
        integrate(
            package_parameters,
            package_parameters.initial_state(),
            label="short_package",
            site="Pv",
            flow_ml_min=240.0,
            config=package_config,
            output_dir=tmp_path / "package",
            overwrite=False,
        )
