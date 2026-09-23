#!/usr/bin/env python3
"""Raw single-patient arterial-pressure pulsatility + aspiration transient pilot.

Uses one high-quality CONTRAST record (CON634) selected previously for stable
rhythm, repeatable beats, visible notch, and low pressure drift.

Important: the arterial waveform used to drive the model is NOT smoothed,
ensemble-averaged, or Fourier-fitted. Four consecutive measured Pa cycles are
taken directly from the raw 100-Hz recording, shifted only by a constant so
their block mean equals 100 mmHg, then linearly interpolated at the solver's
query times and repeated periodically. Linear interpolation is not filtering.

The maintained publication model is not modified.
"""

from __future__ import annotations

import argparse, csv, io, json, math, tarfile, copy
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from cerebral_hemodynamics_aspiration import model
from cerebral_hemodynamics_aspiration.parameters import make_tbi_parameters

FS = 100.0
PATIENT = "CON634.csv"
REST_WINDOW_S = (0.0, 15.0)
N_CYCLES = 4
ASPIRATION_ML_MIN = 240.0
RAMP_S = 60.0
TRANSIENT_S = 600.0


def read_patient(archive: Path) -> tuple[np.ndarray, np.ndarray]:
    with tarfile.open(archive, "r:gz") as tar:
        fh = tar.extractfile(PATIENT)
        if fh is None:
            raise FileNotFoundError(PATIENT)
        reader = csv.reader(io.TextIOWrapper(fh, encoding="utf-8", newline=""))
        header = next(reader)
        pa_i = header.index("Pa")
        ecg_i = header.index("ECG")
        pa, ecg = [], []
        stop_n = int(REST_WINDOW_S[1] * FS)
        for row in reader:
            if len(pa) >= stop_n:
                break
            pa.append(float(row[pa_i]) / 10.0)
            ecg.append(float(row[ecg_i]))
    a = int(REST_WINDOW_S[0] * FS)
    return np.asarray(pa[a:], float), np.asarray(ecg[a:], float)


def choose_four_cycles(pa: np.ndarray) -> tuple[np.ndarray, dict]:
    # Detection only. Values used below remain exactly the raw pressure samples.
    peaks, props = find_peaks(pa, distance=int(0.55 * FS), prominence=10.0)
    peaks = peaks[(peaks > int(0.8 * FS)) & (peaks < len(pa) - int(0.8 * FS))]
    if len(peaks) < N_CYCLES + 1:
        raise RuntimeError("Insufficient clean systolic peaks.")

    best = None
    for i in range(len(peaks) - N_CYCLES):
        p = peaks[i:i+N_CYCLES+1]
        rr = np.diff(p) / FS
        hr = 60.0 / np.mean(rr)
        if not 45 <= hr <= 100:
            continue
        block = pa[p[0]:p[-1]]  # 4 complete peak-to-peak cycles; no endpoint duplication
        if len(block) < 200:
            continue
        pp = np.percentile(block, 99.5) - np.percentile(block, 0.5)
        if not 25 <= pp <= 90:
            continue
        rr_cv = float(np.std(rr) / np.mean(rr))
        endpoint_jump = abs(float(block[-1] - block[0]))
        score = 30.0 * rr_cv + endpoint_jump / max(pp, 1.0)
        if best is None or score < best[0]:
            best = (score, p.copy(), block.copy(), rr.copy(), hr, pp, endpoint_jump)
    if best is None:
        raise RuntimeError("No suitable four-cycle block found.")

    _, peaks5, raw, rr, hr, pp, jump = best
    shifted = raw + (100.0 - float(np.mean(raw)))  # constant shift only
    meta = {
        "patient": PATIENT,
        "source_sampling_hz": FS,
        "peak_sample_indices": peaks5.tolist(),
        "rr_intervals_s": rr.tolist(),
        "heart_rate_bpm": float(hr),
        "raw_block_mean_mmhg": float(np.mean(raw)),
        "constant_shift_mmhg": float(100.0 - np.mean(raw)),
        "shifted_block_mean_mmhg": float(np.mean(shifted)),
        "raw_min_mmhg": float(np.min(raw)),
        "raw_max_mmhg": float(np.max(raw)),
        "raw_pulse_range_mmhg": float(np.ptp(raw)),
        "shifted_min_mmhg": float(np.min(shifted)),
        "shifted_max_mmhg": float(np.max(shifted)),
        "endpoint_jump_raw_mmhg": float(jump),
        "samples": int(len(raw)),
        "block_duration_s": float(len(raw) / FS),
    }
    return shifted, meta


class RawPeriodicPressure:
    def __init__(self, samples: np.ndarray, fs: float):
        self.samples = np.asarray(samples, float)
        self.fs = float(fs)
        self.dt = 1.0 / self.fs
        self.period = len(self.samples) / self.fs

    def __call__(self, t: float) -> float:
        u = float(t) % self.period
        pos = u * self.fs
        i = int(math.floor(pos))
        frac = pos - i
        j = (i + 1) % len(self.samples)
        return float((1.0 - frac) * self.samples[i] + frac * self.samples[j])


def rhs_factory(paw, flow_ml_min: float, ramp_s: float, return_fraction: float):
    p = make_tbi_parameters()
    p_forced = copy.copy(p)
    target = flow_ml_min / 60.0

    def rhs(t, y):
        p_forced.Pa = paw(t)
        delivered = target if ramp_s <= 0 else target * min(max(t / ramp_s, 0.0), 1.0)
        asp = model.Aspiration("Pv", delivered)
        return model.evaluate(y, p_forced, asp, return_fraction * delivered)[0]
    return rhs


def periodic_baseline(initial: np.ndarray, paw: RawPeriodicPressure) -> np.ndarray:
    rhs = rhs_factory(paw, 0.0, 0.0, 0.0)
    state = np.asarray(initial, float).copy()
    prev = None
    stable = 0
    for block in range(1, 121):
        a, b = (block-1)*paw.period, block*paw.period
        te = np.arange(a, b + 0.005, 0.01)
        te[-1] = b
        sol = solve_ivp(rhs, (a,b), state, method="BDF", rtol=1e-8, atol=1e-10,
                        max_step=0.01, t_eval=te)
        if not sol.success:
            raise RuntimeError(sol.message)
        state = sol.y[:,-1]
        wave = sol.y[0,:-1]
        if prev is not None:
            n = min(len(prev), len(wave))
            diff = float(np.max(np.abs(prev[:n] - wave[:n])))
            mean_diff = float(abs(np.mean(prev[:n]) - np.mean(wave[:n])))
            stable = stable + 1 if diff < 5e-4 and mean_diff < 1e-4 else 0
        prev = wave.copy()
        if stable >= 3:
            return state
    raise RuntimeError("Raw-waveform baseline did not reach block-periodic state.")


def run_transient(initial, paw, flow_ml_min, ramp_s, return_fraction):
    rhs = rhs_factory(paw, flow_ml_min, ramp_s, return_fraction)
    te = np.arange(0.0, TRANSIENT_S + 0.005, 0.01)
    te[-1] = TRANSIENT_S
    sol = solve_ivp(rhs, (0.0, TRANSIENT_S), np.asarray(initial,float),
                    method="BDF", rtol=1e-8, atol=1e-10, max_step=0.01, t_eval=te)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def block_stats(t, icp, paw, center):
    half = paw.period / 2
    start = max(0.0, center - half)
    stop = min(float(t[-1]), start + paw.period)
    mask = (t >= start) & (t < stop)
    x = icp[mask]
    return {
        "start_s": float(start), "stop_s": float(stop), "samples": int(mask.sum()),
        "mean_icp_mmhg": float(np.mean(x)),
        "min_icp_mmhg": float(np.min(x)),
        "max_icp_mmhg": float(np.max(x)),
        "pulse_amplitude_mmhg": float(np.ptp(x)),
    }


def domain_check(t, y, paw, flow_ml_min, ramp_s, ret):
    p = make_tbi_parameters(); pp = copy.copy(p); target=flow_ml_min/60.0
    invalid=0; branches=set()
    # Check every 1 s plus first/last four-cycle blocks to keep this inexpensive.
    idx = np.unique(np.concatenate([
        np.arange(0,len(t),100),
        np.arange(0,min(len(t),int(paw.period*FS*2)),1),
        np.arange(max(0,len(t)-int(paw.period*FS*2)),len(t),1)
    ]))
    for j in idx:
        pp.Pa=paw(float(t[j]))
        delivered=target*min(max(float(t[j])/ramp_s,0.0),1.0) if ramp_s>0 else target
        try:
            _,aux=model.evaluate(y[:,j],pp,model.Aspiration("Pv",delivered),ret*delivered)
            branches.add(str(aux["Rvs_branch"]))
        except model.ModelDomainError:
            invalid += 1
    return invalid, sorted(branches)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--output-dir",type=Path,required=True)
    a=ap.parse_args(); out=a.output_dir; out.mkdir(parents=True,exist_ok=True)

    pa,ecg = read_patient(a.archive)
    block, meta = choose_four_cycles(pa)
    paw = RawPeriodicPressure(block, FS)

    base_report=json.loads(Path("reference_results/reports/tbi_baseline.json").read_text())
    static_state=np.asarray(base_report["terminal_window_mean_state"],float)
    periodic_state=periodic_baseline(static_state,paw)

    t0,y0=run_transient(periodic_state,paw,0.0,0.0,0.0)
    t1,y1=run_transient(periodic_state,paw,ASPIRATION_ML_MIN,RAMP_S,1.0)

    invalid0,branches0=domain_check(t0,y0,paw,0.0,0.0,0.0)
    invalid1,branches1=domain_check(t1,y1,paw,ASPIRATION_ML_MIN,RAMP_S,1.0)

    # Save every solver output. No smoothing/downsampling.
    np.savez_compressed(out/"raw_patient_transient.npz",
        time_s=t0, control_icp=y0[0], aspiration_icp=y1[0],
        raw_shifted_pa=np.array([paw(x) for x in t0]),
        patient_block_pa=block)

    with (out/"raw_patient_4cycles.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["time_s","pa_mmhg"])
        for i,v in enumerate(block): w.writerow([i/FS,float(v)])

    times=[0.5*paw.period, 60.0, 120.0, 300.0, 599.0]
    windows={f"around_{int(c)}s":block_stats(t1,y1[0],paw,c) for c in times}
    cstart=block_stats(t0,y0[0],paw,0.5*paw.period)
    cend=block_stats(t0,y0[0],paw,599.0)

    summary={
      "input":meta,
      "method":{
        "waveform_processing":"raw measured Pa samples; constant mean shift to 100 mmHg only; no smoothing, no filtering, no ensemble averaging, no Fourier fit",
        "interpolation":"piecewise linear between original 100-Hz samples",
        "forcing":"selected four-cycle block repeated periodically",
        "solver":"BDF","rtol":1e-8,"atol":1e-10,"max_step_s":0.01,
        "aspiration_site":"Pv","aspiration_ml_min":ASPIRATION_ML_MIN,
        "equal_lower_svc_return":True,"aspiration_ramp_s":RAMP_S,
        "transient_duration_s":TRANSIENT_S
      },
      "control":{"initial_block":cstart,"final_block":cend,
                 "mean_shift_mmhg":cend["mean_icp_mmhg"]-cstart["mean_icp_mmhg"],
                 "invalid_domain_samples":invalid0,"rvs_branches":branches0},
      "aspiration":{"windows":windows,
                    "global_min_icp_mmhg":float(np.min(y1[0])),
                    "time_of_global_min_s":float(t1[np.argmin(y1[0])]),
                    "invalid_domain_samples":invalid1,"rvs_branches":branches1}
    }
    (out/"raw_patient_summary.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    main()
