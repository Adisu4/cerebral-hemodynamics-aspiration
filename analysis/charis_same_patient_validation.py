#!/usr/bin/env python3
"""CHARIS same-patient ABP→model ICP validation pilot.

Purpose
-------
Use one actual TBI patient's simultaneous arterial blood pressure (ABP) and
intracranial pressure (ICP) recording to test whether the maintained cerebral
hemodynamics model produces a realistic cardiac ICP pulse when driven by that
patient's measured ABP.

Final waveform data are raw clinical samples. No smoothing, filtering, Fourier
fitting, amplitude scaling, or waveform normalization is applied to the ABP or
measured ICP used for validation. Peak detection and phase resampling are used
only to select a stable 7-beat window objectively from the measured record.

For visual comparison, one additional model trace is saved after a *constant
vertical shift* to the measured ICP mean. This preserves model pulse amplitude
and morphology exactly; it only removes the expected offset from using a
generic TBI parameter set rather than patient-specific model calibration.

After the validation window is selected using measured-data quality only, a
virtual Pv aspiration test is run with the same raw ABP waveform repeated
periodically, 240 mL/min aspiration, 60-s ramp, and equal lower-SVC return.
The maintained publication model is not modified.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from cerebral_hemodynamics_aspiration import model
from cerebral_hemodynamics_aspiration.parameters import make_tbi_parameters


N_BEATS = 7
COARSE_WINDOW_S = 30.0
COARSE_STEP_S = 60.0
REFINE_STEP_S = 5.0
ASPIRATION_ML_MIN = 240.0
RAMP_S = 60.0
ASPIRATION_DURATION_S = 180.0
ICP_VALIDATION_RANGE_MMHG = (20.0, 25.0)


def parse_header(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    first = lines[0].split()
    record = first[0]
    n_sig = int(first[1])
    fs = float(first[2])
    n_samples = int(first[3])
    signals = []
    pat = re.compile(r"(?P<gain>[0-9.]+)\((?P<baseline>-?\d+)\)/(?P<unit>\S+)")
    for i in range(n_sig):
        parts = lines[i + 1].split()
        m = pat.fullmatch(parts[2])
        if not m:
            raise ValueError(f"Could not parse gain/baseline from {parts[2]!r}")
        signals.append(
            {
                "index": i,
                "gain": float(m.group("gain")),
                "baseline": int(m.group("baseline")),
                "unit": m.group("unit"),
                "name": parts[-1],
            }
        )
    comments = [line[2:] for line in lines if line.startswith("# ")]
    return {
        "record": record,
        "n_sig": n_sig,
        "fs": fs,
        "n_samples": n_samples,
        "signals": signals,
        "comments": comments,
    }


def phys(segment: np.ndarray, spec: dict) -> np.ndarray:
    return (segment[:, spec["index"]].astype(float) - spec["baseline"]) / spec["gain"]


def phase_resample(x: np.ndarray, n: int = 128) -> np.ndarray:
    old = np.linspace(0.0, 1.0, len(x), endpoint=True)
    new = np.linspace(0.0, 1.0, n, endpoint=False)
    return np.interp(new, old, x)


def median_beat_correlation(signal: np.ndarray, peaks: np.ndarray) -> float:
    beats = []
    for i in range(len(peaks) - 1):
        beat = signal[peaks[i] : peaks[i + 1] + 1]
        if len(beat) < 10:
            continue
        beats.append(phase_resample(beat))
    if len(beats) < N_BEATS:
        return float("nan")
    B = np.vstack(beats)
    template = np.mean(B, axis=0)
    corr = []
    for b in B:
        if np.std(b) > 0 and np.std(template) > 0:
            corr.append(float(np.corrcoef(b, template)[0, 1]))
    return float(np.median(corr)) if corr else float("nan")


def slope_per_minute(x: np.ndarray, fs: float) -> float:
    t = np.arange(len(x), dtype=float) / fs
    return float(np.polyfit(t, x, 1)[0] * 60.0)


def evaluate_window(abp: np.ndarray, icp: np.ndarray, fs: float, global_start: int) -> dict | None:
    if len(abp) < int(COARSE_WINDOW_S * fs) or len(icp) != len(abp):
        return None
    if not (np.all(np.isfinite(abp)) and np.all(np.isfinite(icp))):
        return None

    # Broad physiologic/artifact gates only. Selection is not based on model output.
    if not (45.0 <= np.mean(abp) <= 155.0):
        return None
    if np.quantile(abp, 0.005) < 20.0 or np.quantile(abp, 0.995) > 230.0:
        return None
    if not (ICP_VALIDATION_RANGE_MMHG[0] <= np.mean(icp) <= ICP_VALIDATION_RANGE_MMHG[1]):
        return None
    if np.quantile(icp, 0.001) < -20.0 or np.quantile(icp, 0.999) > 80.0:
        return None

    peaks, _ = find_peaks(
        abp,
        distance=int(0.45 * fs),
        prominence=12.0,
    )
    if len(peaks) < N_BEATS + 1:
        return None

    best = None
    for i in range(len(peaks) - N_BEATS):
        p = peaks[i : i + N_BEATS + 1]
        rr = np.diff(p) / fs
        mean_rr = float(np.mean(rr))
        if mean_rr <= 0:
            continue
        hr = 60.0 / mean_rr
        rr_cv = float(np.std(rr) / mean_rr)
        if not (45.0 <= hr <= 105.0) or rr_cv > 0.06:
            continue

        start, stop = int(p[0]), int(p[-1])
        a = abp[start:stop]
        c = icp[start:stop]
        if len(a) < int(N_BEATS * 0.45 * fs):
            continue

        abp_pp = float(np.quantile(a, 0.995) - np.quantile(a, 0.005))
        icp_pp = float(np.quantile(c, 0.995) - np.quantile(c, 0.005))
        if not (20.0 <= abp_pp <= 110.0):
            continue
        if not (0.25 <= icp_pp <= 15.0):
            continue

        local_peaks = p - start
        abp_corr = median_beat_correlation(a, local_peaks)
        icp_corr = median_beat_correlation(c, local_peaks)
        if not np.isfinite(abp_corr) or not np.isfinite(icp_corr):
            continue
        if abp_corr < 0.90 or icp_corr < 0.65:
            continue

        abp_drift = slope_per_minute(a, fs)
        icp_drift = slope_per_minute(c, fs)

        # Objective signal-quality score. No model quantity appears here.
        score = (
            12.0 * rr_cv
            + 2.0 * (1.0 - abp_corr)
            + 2.0 * (1.0 - icp_corr)
            + 0.015 * abs(abp_drift)
            + 0.05 * abs(icp_drift)
        )
        row = {
            "score": float(score),
            "global_start_sample": int(global_start + start),
            "global_stop_sample": int(global_start + stop),
            "local_peak_samples": p.tolist(),
            "hr_bpm": float(hr),
            "rr_cv": rr_cv,
            "abp_mean_mmhg": float(np.mean(a)),
            "abp_pulse_range_mmhg": abp_pp,
            "icp_mean_mmhg": float(np.mean(c)),
            "icp_pulse_range_mmhg": icp_pp,
            "abp_median_beat_corr": abp_corr,
            "icp_median_beat_corr": icp_corr,
            "abp_drift_mmhg_min": float(abp_drift),
            "icp_drift_mmhg_min": float(icp_drift),
        }
        if best is None or row["score"] < best["score"]:
            best = row
    return best


def choose_segment(data: np.memmap, specs: dict, fs: float, n_samples: int) -> dict:
    abp_spec, icp_spec = specs["ABP"], specs["ICP"]
    win = int(COARSE_WINDOW_S * fs)
    step = int(COARSE_STEP_S * fs)
    coarse = []

    for start in range(0, n_samples - win, step):
        seg = data[start : start + win]
        abp = phys(seg, abp_spec)
        icp = phys(seg, icp_spec)
        result = evaluate_window(abp, icp, fs, start)
        if result is not None:
            coarse.append(result)

    if not coarse:
        raise RuntimeError("No coarse CHARIS windows passed QC.")

    coarse.sort(key=lambda x: x["score"])
    top = coarse[: min(15, len(coarse))]

    refined = []
    refine_step = int(REFINE_STEP_S * fs)
    radius = int(COARSE_STEP_S * fs)
    seen = set()
    for seed in top:
        center = seed["global_start_sample"]
        lo = max(0, center - radius)
        hi = min(n_samples - win, center + radius)
        for start in range(lo, hi + 1, refine_step):
            if start in seen:
                continue
            seen.add(start)
            seg = data[start : start + win]
            result = evaluate_window(phys(seg, abp_spec), phys(seg, icp_spec), fs, start)
            if result is not None:
                refined.append(result)

    candidates = refined if refined else coarse
    candidates.sort(key=lambda x: x["score"])
    best = candidates[0]
    best["n_coarse_candidates"] = len(coarse)
    best["n_refined_candidates"] = len(refined)
    return best


class RawPeriodicPressure:
    def __init__(self, samples: np.ndarray, fs: float):
        self.samples = np.asarray(samples, float)
        self.fs = float(fs)
        self.period_s = len(self.samples) / self.fs

    def __call__(self, t: float) -> float:
        pos = (float(t) % self.period_s) * self.fs
        i = int(math.floor(pos)) % len(self.samples)
        frac = pos - math.floor(pos)
        j = (i + 1) % len(self.samples)
        return float((1.0 - frac) * self.samples[i] + frac * self.samples[j])


def make_rhs(paw: RawPeriodicPressure, target_flow_ml_min: float, ramp_s: float, return_fraction: float):
    p = make_tbi_parameters()
    p_forced = copy.copy(p)
    target = target_flow_ml_min / 60.0

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        p_forced.Pa = paw(t)
        delivered = target
        if ramp_s > 0:
            delivered *= min(max(float(t) / ramp_s, 0.0), 1.0)
        aspiration = model.Aspiration("Pv", delivered)
        return model.evaluate(y, p_forced, aspiration, return_fraction * delivered)[0]

    return rhs


def equilibrate_at_mean_pressure(initial: np.ndarray, mean_pa_mmhg: float) -> np.ndarray:
    """Settle the slow model states at the patient's measured mean ABP first.

    This avoids asking a short cardiac-period convergence loop to also absorb
    hours-scale mean-state adaptation. The raw pulsatile waveform is introduced
    only after this constant-mean equilibration.
    """
    p = make_tbi_parameters()
    p.Pa = float(mean_pa_mmhg)
    aspiration = model.Aspiration("Pv", 0.0)

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        return model.evaluate(y, p, aspiration, 0.0)[0]

    sol = solve_ivp(
        rhs,
        (0.0, 12000.0),
        np.asarray(initial, float),
        method="BDF",
        rtol=1e-8,
        atol=1e-10,
        max_step=2.0,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1]


def converge_periodic(initial: np.ndarray, paw: RawPeriodicPressure) -> tuple[np.ndarray, dict]:
    """Converge the *cardiac waveform shape* under raw repeated ABP.

    Absolute mean state adaptation is much slower than the cardiac cycle in this
    model. For waveform validation we therefore require convergence of the
    mean-centered ICP pulse shape and pulse amplitude, while recording any
    residual block-to-block mean drift explicitly.
    """
    rhs = make_rhs(paw, 0.0, 0.0, 0.0)
    state = np.asarray(initial, float).copy()
    previous_centered = None
    previous_amp = None
    previous_mean = None
    stable = 0
    diagnostics = {}

    for block in range(1, 81):
        a = (block - 1) * paw.period_s
        b = block * paw.period_s
        te = a + np.arange(len(paw.samples), dtype=float) / paw.fs
        sol = solve_ivp(
            rhs,
            (a, b),
            state,
            method="BDF",
            rtol=1e-8,
            atol=1e-10,
            max_step=min(0.01, 0.5 / paw.fs),
            t_eval=te,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        state = sol.y[:, -1]
        wave = sol.y[0]
        mean = float(np.mean(wave))
        centered = wave - mean
        amp = float(np.ptp(wave))

        if previous_centered is not None and len(previous_centered) == len(centered):
            shape_diff = float(np.max(np.abs(centered - previous_centered)))
            amp_diff = float(abs(amp - previous_amp))
            mean_drift = float(mean - previous_mean)
            stable = stable + 1 if shape_diff < 0.001 and amp_diff < 0.001 else 0
        else:
            shape_diff = float("inf")
            amp_diff = float("inf")
            mean_drift = float("inf")
            stable = 0

        diagnostics = {
            "blocks_run": block,
            "shape_converged": bool(stable >= 3),
            "final_centered_wave_max_diff_mmhg": shape_diff,
            "final_pulse_amplitude_diff_mmhg": amp_diff,
            "final_block_mean_drift_mmhg": mean_drift,
            "final_block_mean_icp_mmhg": mean,
            "final_block_peak_to_peak_mmhg": amp,
        }
        previous_centered = centered.copy()
        previous_amp = amp
        previous_mean = mean
        if stable >= 3:
            break

    return state, diagnostics

def simulate_one_block(initial: np.ndarray, paw: RawPeriodicPressure) -> tuple[np.ndarray, np.ndarray]:
    rhs = make_rhs(paw, 0.0, 0.0, 0.0)
    t = np.arange(len(paw.samples), dtype=float) / paw.fs
    sol = solve_ivp(
        rhs,
        (0.0, paw.period_s),
        np.asarray(initial, float),
        method="BDF",
        rtol=1e-8,
        atol=1e-10,
        max_step=min(0.01, 0.5 / paw.fs),
        t_eval=t,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def correlation_metrics(measured: np.ndarray, predicted: np.ndarray, fs: float) -> dict:
    n = min(len(measured), len(predicted))
    m = np.asarray(measured[:n], float)
    p = np.asarray(predicted[:n], float)
    mc = m - np.mean(m)
    pc = p - np.mean(p)
    direct = float(np.corrcoef(mc, pc)[0, 1])

    max_lag = int(0.50 * fs)
    best = (-np.inf, 0)
    for lag in range(-max_lag, max_lag + 1):
        if lag < 0:
            a, b = mc[-lag:], pc[: n + lag]
        elif lag > 0:
            a, b = mc[: n - lag], pc[lag:]
        else:
            a, b = mc, pc
        if len(a) < 20 or np.std(a) == 0 or np.std(b) == 0:
            continue
        r = float(np.corrcoef(a, b)[0, 1])
        if r > best[0]:
            best = (r, lag)

    aligned = p + (np.mean(m) - np.mean(p))
    return {
        "measured_mean_icp_mmhg": float(np.mean(m)),
        "predicted_mean_icp_mmhg": float(np.mean(p)),
        "constant_vertical_shift_for_overlay_mmhg": float(np.mean(m) - np.mean(p)),
        "measured_peak_to_peak_mmhg": float(np.ptp(m)),
        "predicted_peak_to_peak_mmhg": float(np.ptp(p)),
        "pulse_amplitude_ratio_predicted_over_measured": float(np.ptp(p) / np.ptp(m)),
        "direct_zero_lag_correlation_mean_centered": direct,
        "best_lag_correlation_mean_centered": float(best[0]),
        "best_lag_samples": int(best[1]),
        "best_lag_seconds": float(best[1] / fs),
        "mean_aligned_rmse_mmhg": float(np.sqrt(np.mean((m - aligned) ** 2))),
    }


def run_aspiration(initial: np.ndarray, paw: RawPeriodicPressure) -> tuple[np.ndarray, np.ndarray]:
    rhs = make_rhs(paw, ASPIRATION_ML_MIN, RAMP_S, 1.0)
    n_out = int(round(ASPIRATION_DURATION_S * paw.fs)) + 1
    t = np.linspace(0.0, ASPIRATION_DURATION_S, n_out)
    sol = solve_ivp(
        rhs,
        (0.0, ASPIRATION_DURATION_S),
        np.asarray(initial, float),
        method="BDF",
        rtol=1e-8,
        atol=1e-10,
        max_step=min(0.01, 0.5 / paw.fs),
        t_eval=t,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def extract_7beat_window(t: np.ndarray, icp: np.ndarray, paw: RawPeriodicPressure, start_s: float) -> dict:
    stop = start_s + paw.period_s
    mask = (t >= start_s) & (t < stop)
    return {
        "start_s": float(start_s),
        "stop_s": float(stop),
        "mean_icp_mmhg": float(np.mean(icp[mask])),
        "min_icp_mmhg": float(np.min(icp[mask])),
        "max_icp_mmhg": float(np.max(icp[mask])),
        "peak_to_peak_mmhg": float(np.ptp(icp[mask])),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dat", type=Path, required=True)
    parser.add_argument("--hea", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    header = parse_header(args.hea)
    fs = header["fs"]
    n_sig = header["n_sig"]
    n_samples = header["n_samples"]
    specs = {s["name"]: s for s in header["signals"]}
    if "ABP" not in specs or "ICP" not in specs:
        raise RuntimeError("Record lacks ABP or ICP.")

    raw = np.memmap(args.dat, dtype="<i2", mode="r")
    expected = n_samples * n_sig
    if raw.size != expected:
        raise RuntimeError(f"DAT sample count mismatch: {raw.size} vs {expected}")
    data = raw.reshape(n_samples, n_sig)

    selected = choose_segment(data, specs, fs, n_samples)
    s0, s1 = selected["global_start_sample"], selected["global_stop_sample"]
    seg = data[s0:s1]
    abp = phys(seg, specs["ABP"])
    icp_measured = phys(seg, specs["ICP"])
    time = np.arange(len(abp), dtype=float) / fs

    # Selection begins/ends on detected systolic peaks; final data remain raw.
    # The maintained TBI model is defined at Pa=100 mmHg. For waveform
    # validation, preserve this patient's measured beat-to-beat ABP exactly and
    # apply only a constant offset so the cycle-block mean equals 100 mmHg.
    # No pulse amplitude or shape rescaling is performed.
    abp_model_input = abp + (100.0 - float(np.mean(abp)))
    paw = RawPeriodicPressure(abp_model_input, fs)

    static_report = json.loads(Path("reference_results/reports/tbi_baseline.json").read_text())
    static_state = np.asarray(static_report["terminal_window_mean_state"], float)
    print("Selected measured-data window:", json.dumps(selected, indent=2), flush=True)
    periodic_state, periodic_diagnostics = converge_periodic(static_state, paw)
    print("Periodic-shape diagnostics:", json.dumps(periodic_diagnostics, indent=2), flush=True)
    tp, yp = simulate_one_block(periodic_state, paw)
    icp_pred = yp[0]

    metrics = correlation_metrics(icp_measured, icp_pred, fs)
    vertical_shift = metrics["constant_vertical_shift_for_overlay_mmhg"]
    icp_pred_mean_aligned = icp_pred + vertical_shift

    ta, ya = run_aspiration(periodic_state, paw)
    icp_asp = ya[0]

    # Seven-beat windows: baseline and two post-onset windows.
    last_start = max(0.0, ASPIRATION_DURATION_S - paw.period_s)
    aspiration_windows = {
        "baseline": {
            "mean_icp_mmhg": float(np.mean(icp_pred)),
            "min_icp_mmhg": float(np.min(icp_pred)),
            "max_icp_mmhg": float(np.max(icp_pred)),
            "peak_to_peak_mmhg": float(np.ptp(icp_pred)),
        },
        "around_60s": extract_7beat_window(ta, icp_asp, paw, 60.0),
        "around_120s": extract_7beat_window(ta, icp_asp, paw, 120.0),
        "final": extract_7beat_window(ta, icp_asp, paw, last_start),
    }

    # Save exact measured and model samples used for all figures.
    np.savez_compressed(
        out / "charis1_validation_and_aspiration.npz",
        validation_time_s=time,
        measured_abp_mmhg=abp,
        model_input_abp_mmhg=abp_model_input,
        measured_icp_mmhg=icp_measured,
        predicted_icp_mmhg=icp_pred,
        predicted_icp_mean_aligned_mmhg=icp_pred_mean_aligned,
        aspiration_time_s=ta,
        aspiration_icp_mmhg=icp_asp,
    )

    with (out / "charis1_selected_7beats.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["time_s", "measured_abp_mmhg", "model_input_abp_mmhg", "measured_icp_mmhg", "predicted_icp_mmhg", "predicted_icp_mean_aligned_mmhg"])
        for row in zip(time, abp, abp_model_input, icp_measured, icp_pred, icp_pred_mean_aligned):
            w.writerow([float(x) for x in row])

    summary = {
        "source": {
            "database": "CHARIS database, PhysioNet v1.0.0",
            "doi": "10.13026/C24G6F",
            "record": header["record"],
            "record_metadata": header["comments"],
            "sampling_hz": fs,
            "selection_rule": "eligibility required measured mean ICP 20-25 mmHg (matching the model's intracranial-hypertension operating range); within eligible windows, selection used measured ABP/ICP signal quality only and no model output",
            "validation_icp_mean_range_mmhg": list(ICP_VALIDATION_RANGE_MMHG),
            "final_waveform_processing": "measured ABP and ICP remain raw; no smoothing/filtering/ensemble averaging/Fourier fit",
            "model_input_abp": "raw measured ABP plus one constant offset to make its 7-beat mean 100 mmHg; pulse amplitude and shape unchanged",
            "model_input_abp_constant_shift_mmhg": float(100.0 - np.mean(abp)),
            "model_input_interpolation": "piecewise linear between raw 50-Hz ABP samples",
        },
        "selection": selected,
        "periodic_shape_convergence": periodic_diagnostics,
        "validation": metrics,
        "aspiration": {
            "site": "Pv",
            "flow_ml_min": ASPIRATION_ML_MIN,
            "ramp_s": RAMP_S,
            "equal_lower_svc_return": True,
            "duration_s": ASPIRATION_DURATION_S,
            "windows": aspiration_windows,
            "global_min_icp_mmhg": float(np.min(icp_asp)),
            "time_of_global_min_s": float(ta[np.argmin(icp_asp)]),
        },
    }
    (out / "charis1_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
