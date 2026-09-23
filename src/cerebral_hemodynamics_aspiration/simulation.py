"""Numerical integration, convergence assessment, and result serialization."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import platform
from pathlib import Path
import time

import numpy as np
import scipy
from scipy.integrate import solve_ivp

from . import model
from .parameters import ModelParameters


PACKAGE_ROOT = Path(__file__).resolve().parent
DEFAULT_RESULTS_DIR = Path("results") / "simulations"


@dataclass(frozen=True)
class IntegratorConfig:
    method: str = "BDF"
    rtol: float = 1e-8
    atol: float = 1e-10
    max_step_s: float = 2.0
    save_step_s: float = 1.0
    chunk_s: float = 600.0
    max_duration_s: float = 18000.0
    averaging_window_s: float = 120.0
    icp_std_limit_mmhg: float = 5e-4
    pressure_drift_limit_mmhg_min: float = 1e-4
    pressure_residual_limit_mmhg_min: float = 1e-4
    consecutive_passes: int = 2


def _jsonable(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def runtime_environment() -> dict[str, str]:
    """Return solver-relevant runtime identifiers recorded with each run."""
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
    }


def source_hashes() -> dict[str, str]:
    files = (
        PACKAGE_ROOT / "model.py",
        PACKAGE_ROOT / "parameters.py",
        PACKAGE_ROOT / "simulation.py",
        PACKAGE_ROOT / "experiments.py",
    )
    return {
        f"src/cerebral_hemodynamics_aspiration/{path.name}": hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        for path in files if path.exists()
    }


def input_fingerprint(
    parameters: ModelParameters,
    initial_state: np.ndarray,
    *,
    label: str,
    site: str = "Pv",
    flow_ml_min: float = 0.0,
    return_fraction: float = 0.0,
    ramp_s: float = 60.0,
    config: IntegratorConfig | None = None,
) -> str:
    """Return a stable hash of every numerical input and current source file."""
    config = config or IntegratorConfig()
    payload = {
        "label": label,
        "site": site,
        "flow_ml_min": float(flow_ml_min),
        "return_fraction": float(return_fraction),
        "ramp_s": float(ramp_s),
        "initial_state": np.asarray(initial_state, dtype=float).tolist(),
        "parameters": parameters.serializable(),
        "integrator": asdict(config),
        "source_hashes_sha256": source_hashes(),
        "environment": runtime_environment(),
    }
    encoded = json.dumps(
        _jsonable(payload), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _window_diagnostics(
    t: np.ndarray,
    y: np.ndarray,
    rhs,
    averaging_window_s: float,
) -> dict:
    cutoff = max(float(t[0]), float(t[-1]) - averaging_window_s)
    mask = t >= cutoff
    if np.count_nonzero(mask) < 3:
        raise RuntimeError("Not enough saved samples for convergence diagnostics")
    window_t = t[mask]
    window_y = y[:, mask]
    slopes = np.array([
        np.polyfit(window_t, row, 1)[0] * 60.0 for row in window_y[:model.PRESSURE_STATE_COUNT]
    ])
    terminal_rhs = np.asarray(rhs(float(t[-1]), y[:, -1]), dtype=float)
    return {
        "window_start_s": float(window_t[0]),
        "window_end_s": float(window_t[-1]),
        "window_s": float(window_t[-1] - window_t[0]),
        "mean_state": np.mean(window_y, axis=1),
        "std_state": np.std(window_y, axis=1),
        "pressure_drift_mmhg_min": slopes,
        "icp_mean_mmhg": float(np.mean(window_y[0])),
        "icp_std_mmhg": float(np.std(window_y[0])),
        "icp_drift_mmhg_min": float(slopes[0]),
        "max_pressure_drift_mmhg_min": float(np.max(np.abs(slopes))),
        "max_pressure_residual_mmhg_min": float(
            np.max(np.abs(terminal_rhs[:model.PRESSURE_STATE_COUNT])) * 60.0
        ),
    }


def _convergence_passed(stats: dict, config: IntegratorConfig) -> bool:
    return bool(
        stats["icp_std_mmhg"] < config.icp_std_limit_mmhg
        and stats["max_pressure_drift_mmhg_min"] < config.pressure_drift_limit_mmhg_min
        and stats["max_pressure_residual_mmhg_min"] < config.pressure_residual_limit_mmhg_min
    )


def integrate(
    parameters: ModelParameters,
    initial_state: np.ndarray,
    *,
    label: str,
    site: str = "Pv",
    flow_ml_min: float = 0.0,
    return_fraction: float = 0.0,
    ramp_s: float = 60.0,
    config: IntegratorConfig | None = None,
    output_dir: Path | str = DEFAULT_RESULTS_DIR,
    overwrite: bool = True,
) -> dict:
    """Integrate until strict equilibrium criteria pass twice in succession."""
    config = config or IntegratorConfig()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / f"{label}.json"
    npz_path = output_dir / f"{label}.npz"
    if site not in model.VALID_ASPIRATION_SITES:
        raise ValueError(f"site must be one of {model.VALID_ASPIRATION_SITES}")
    if flow_ml_min < 0.0 or not 0.0 <= return_fraction <= 1.0:
        raise ValueError("flow must be nonnegative and return_fraction must be in [0, 1]")
    state = np.asarray(initial_state, dtype=float).copy()
    if state.shape != (len(model.STATE_NAMES),):
        raise ValueError(f"initial_state must have shape ({len(model.STATE_NAMES)},)")

    fingerprint = input_fingerprint(
        parameters,
        state,
        label=label,
        site=site,
        flow_ml_min=flow_ml_min,
        return_fraction=return_fraction,
        ramp_s=ramp_s,
        config=config,
    )
    if not overwrite:
        if json_path.exists() != npz_path.exists():
            raise RuntimeError(f"{label}: cached output pair is incomplete")
        if json_path.exists():
            cached = json.loads(json_path.read_text(encoding="utf-8"))
            if cached.get("input_fingerprint_sha256") != fingerprint:
                raise RuntimeError(
                    f"{label}: cached output does not match current inputs/source"
                )
            if not cached.get("converged"):
                raise RuntimeError(f"{label}: cached run did not converge")
            if int(cached.get("domain", {}).get("invalid_open_branch_samples", -1)):
                raise RuntimeError(f"{label}: cached run has an invalid model domain")
            expected_trajectory_hash = cached.get("trajectory_sha256")
            if (
                not expected_trajectory_hash
                or _sha256_file(npz_path) != expected_trajectory_hash
            ):
                raise RuntimeError(f"{label}: cached trajectory hash mismatch")
            return cached

    target_flow_ml_s = flow_ml_min / 60.0

    def delivered_flow(t: float) -> float:
        if ramp_s <= 0.0:
            return target_flow_ml_s
        return target_flow_ml_s * min(max(t / ramp_s, 0.0), 1.0)

    aspiration_callback = lambda t, y: model.Aspiration(site, delivered_flow(t))
    return_callback = lambda t, y: return_fraction * delivered_flow(t)
    ode = model.rhs(parameters, aspiration_callback, return_callback)

    start_wall = time.monotonic()
    current_time = 0.0
    stable_passes = 0
    time_parts: list[np.ndarray] = []
    state_parts: list[np.ndarray] = []
    checks: list[dict] = []

    while current_time < config.max_duration_s:
        end_time = min(current_time + config.chunk_s, config.max_duration_s)
        samples = np.arange(current_time, end_time + config.save_step_s * 0.5, config.save_step_s)
        if samples[-1] > end_time:
            samples[-1] = end_time
        elif samples[-1] < end_time:
            samples = np.append(samples, end_time)
        solution = solve_ivp(
            ode,
            (current_time, end_time),
            state,
            method=config.method,
            rtol=config.rtol,
            atol=config.atol,
            max_step=config.max_step_s,
            t_eval=samples,
        )
        if not solution.success:
            raise RuntimeError(f"{label}: integrator failed: {solution.message}")
        chunk_t = solution.t
        chunk_y = solution.y
        stats = _window_diagnostics(chunk_t, chunk_y, ode, config.averaging_window_s)
        passed = _convergence_passed(stats, config)
        stable_passes = stable_passes + 1 if passed else 0
        stats["passed"] = passed
        stats["consecutive_passes"] = stable_passes
        checks.append(_jsonable(stats))
        if time_parts:
            chunk_t = chunk_t[1:]
            chunk_y = chunk_y[:, 1:]
        time_parts.append(chunk_t)
        state_parts.append(chunk_y)
        current_time = end_time
        state = solution.y[:, -1]
        if stable_passes >= config.consecutive_passes:
            break

    t = np.concatenate(time_parts)
    y = np.concatenate(state_parts, axis=1)
    final_stats = _window_diagnostics(t, y, ode, config.averaging_window_s)
    final_stats["passed"] = _convergence_passed(final_stats, config)
    mean_state = np.asarray(final_stats["mean_state"], dtype=float)
    final_aspiration = model.Aspiration(site, target_flow_ml_s)
    final_flows = model.flow_snapshot(
        mean_state, parameters, final_aspiration, return_fraction * target_flow_ml_s
    )

    open_mask = y[2] > y[3]
    if parameters.terminal_resistance_mode == "fixed_Rvs1":
        invalid_open_mask = np.zeros_like(open_mask, dtype=bool)
        branch_code = np.full(t.shape, 2, dtype=np.int8)
        Rvs = np.full(t.shape, parameters.Rvs1, dtype=float)
    else:
        invalid_open_mask = open_mask & (y[2] <= y[0])
        branch_code = np.where(open_mask, 0, 1).astype(np.int8)
        Rvs = np.where(
            open_mask,
            parameters.Rvs1 * (y[2] - y[3]) / (y[2] - y[0]),
            parameters.Rvs1,
        )
    Psvc1 = np.empty_like(t)
    for index in range(t.size):
        conductances = model.jugular_conductances(y[:, index], parameters)
        Psvc1[index] = model.superior_vena_cava_junction(y[:, index], conductances, parameters)
    delivered = np.array([delivered_flow(value) * 60.0 for value in t])

    np.savez_compressed(
        npz_path,
        t=t,
        y=y,
        state_names=np.asarray(model.STATE_NAMES),
        flow_ml_min=delivered,
        branch_code=branch_code,
        Rvs=Rvs,
        Psvc1=Psvc1,
    )

    report = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "label": label,
        "site": site,
        "flow_ml_min": float(flow_ml_min),
        "ramp_s": float(ramp_s),
        "return_fraction": float(return_fraction),
        "return_site": "lower_SVC_state" if return_fraction else None,
        "converged": stable_passes >= config.consecutive_passes,
        "duration_s": float(current_time),
        "integrator": asdict(config),
        "convergence_checks": checks,
        "summary": _jsonable(final_stats),
        "final_state": _jsonable(state),
        "terminal_window_mean_state": _jsonable(mean_state),
        "terminal_flows_ml_s": _jsonable(final_flows),
        "parameters": parameters.serializable(),
        "domain": {
            "branches_seen": sorted(set(branch_code.tolist())),
            "branch_code_definition": {
                "0": "open_starling", "1": "reverse_or_equal", "2": "fixed_Rvs1_comparator"
            },
            "invalid_open_branch_samples": int(np.count_nonzero(invalid_open_mask)),
            "minimum_Pv_minus_Pic_mmhg": float(np.min(y[2] - y[0])),
            "minimum_Cvi_denominator_pressure_mmhg": float(
                np.min(y[2] - y[0] - parameters.Pv1)
            ),
            "minimum_Ppa_minus_Pic_mmhg": float(np.min(y[1] - y[0])),
        },
        "environment": runtime_environment(),
        "source_hashes_sha256": source_hashes(),
        "input_fingerprint_sha256": fingerprint,
        "trajectory_file": npz_path.name,
        "trajectory_sha256": _sha256_file(npz_path),
        "wall_time_s": time.monotonic() - start_wall,
    }
    json_path.write_text(json.dumps(_jsonable(report), indent=2), encoding="utf-8")
    print(
        f"{label}: converged={report['converged']} duration={current_time:.0f}s "
        f"ICP={final_stats['icp_mean_mmhg']:.6f} "
        f"drift={final_stats['icp_drift_mmhg_min']:.3g} mmHg/min",
        flush=True,
    )
    if not report["converged"]:
        raise RuntimeError(f"{label}: failed convergence criteria by {current_time:.0f} s")
    if report["domain"]["invalid_open_branch_samples"]:
        raise RuntimeError(f"{label}: Eq. 14 domain was violated")
    return report


def load_report(label: str, output_dir: Path | str = DEFAULT_RESULTS_DIR) -> dict:
    return json.loads((Path(output_dir) / f"{label}.json").read_text(encoding="utf-8"))


__all__ = [
    "DEFAULT_RESULTS_DIR", "IntegratorConfig", "input_fingerprint", "integrate",
    "load_report", "runtime_environment", "source_hashes",
]
