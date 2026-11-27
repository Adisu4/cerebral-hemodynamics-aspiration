#!/usr/bin/env python3
"""CHARIS single-patient waveform validation and aspiration pilot.

Patient: charis1 (19M, SAH/SDH/TBI), chosen a priori because the header
explicitly identifies TBI and provides continuous ABP, ECG, and ICP.

No additional digital smoothing or filtering is applied to ABP or ICP.
Signal-quality metrics are used only to choose a clean 20-s segment. Five
consecutive raw cardiac cycles are then used for the comparison.

For model forcing, raw radial ABP is shifted by one constant so the selected
5-cycle mean equals the maintained model boundary mean of 100 mmHg. Pulse
amplitude and morphology are otherwise unchanged. Model output is not fitted
to the measured ICP waveform.
"""

from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from cerebral_hemodynamics_aspiration import model
from cerebral_hemodynamics_aspiration.parameters import make_tbi_parameters

FS = 50.0
GAIN_ABP, BASE_ABP = 85.953, -1332.0
GAIN_ECG, BASE_ECG = 6081.5219, -398.0
GAIN_ICP, BASE_ICP = 94.8784, -2316.0
NCHAN = 3
N_CYCLES = 5
TARGET_ICP_RANGE = (18.0, 26.0)
RAMP_S = 60.0
ASP_ML_MIN = 240.0
TRANSIENT_S = 100.0

def load_record(path: Path):
    raw = np.memmap(path, dtype="<i2", mode="r")
    n = (raw.size // NCHAN) * NCHAN
    x = np.asarray(raw[:n]).reshape(-1, NCHAN)
    abp = (x[:,0].astype(float) - BASE_ABP) / GAIN_ABP
    ecg = (x[:,1].astype(float) - BASE_ECG) / GAIN_ECG
    icp = (x[:,2].astype(float) - BASE_ICP) / GAIN_ICP
    return abp, ecg, icp

def window_quality(abp, icp, start, dur_s=20.0):
    n = int(dur_s*FS)
    a = abp[start:start+n]; p = icp[start:start+n]
    if len(a) < n: return None
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(p)): return None
    # Explicit artifact/physiology guardrails; not waveform fitting.
    if np.quantile(a,.005) < 40 or np.quantile(a,.995) > 220: return None
    if np.quantile(p,.005) < -5 or np.quantile(p,.995) > 60: return None
    ma, mp = float(np.mean(a)), float(np.mean(p))
    if not 70 <= ma <= 140: return None
    if not (TARGET_ICP_RANGE[0] <= mp <= TARGET_ICP_RANGE[1]): return None
    peaks,_ = find_peaks(a, distance=int(.45*FS), prominence=15.0)
    if len(peaks) < 14: return None
    rr = np.diff(peaks)/FS
    med = float(np.median(rr))
    hr = 60/med
    if not 50 <= hr <= 110: return None
    rr_cv = float(np.std(rr)/np.mean(rr))
    if rr_cv > .08: return None
    # Avoid obvious flat/flush traces.
    abp_pp = float(np.percentile(a,99)-np.percentile(a,1))
    icp_pp = float(np.percentile(p,99)-np.percentile(p,1))
    if not 20 <= abp_pp <= 100: return None
    if not .4 <= icp_pp <= 15: return None
    t=np.arange(n)/FS
    abp_drift=float(np.polyfit(t,a,1)[0]*60)
    icp_drift=float(np.polyfit(t,p,1)[0]*60)
    # We want a stationary validation window, not a clinical trend episode.
    if abs(abp_drift)>12 or abs(icp_drift)>4: return None
    score = 10*rr_cv + 0.03*abs(abp_drift) + 0.10*abs(icp_drift)
    return dict(score=score, mean_abp=ma, mean_icp=mp, hr=hr, rr_cv=rr_cv,
                abp_pp=abp_pp, icp_pp=icp_pp, abp_drift=abp_drift,
                icp_drift=icp_drift, peaks=peaks)

def select_clean_window(abp, icp):
    win=int(20*FS); step=int(30*FS)
    candidates=[]
    for start in range(0, len(abp)-win, step):
        q=window_quality(abp,icp,start)
        if q is not None:
            q["start"]=start
            candidates.append(q)
    if not candidates:
        raise RuntimeError("No clean 20-s segment in target ICP range.")
    candidates.sort(key=lambda x:x["score"])
    return candidates[0], len(candidates)

def select_five_cycles(abp, icp, q):
    start=q["start"]; n=int(20*FS)
    a=abp[start:start+n]; p=icp[start:start+n]
    peaks,_=find_peaks(a, distance=int(.45*FS), prominence=15.0)
    best=None
    for i in range(len(peaks)-N_CYCLES):
        pk=peaks[i:i+N_CYCLES+1]
        rr=np.diff(pk)/FS
        cv=float(np.std(rr)/np.mean(rr))
        if cv>.05: continue
        aa=a[pk[0]:pk[-1]]
        pp=p[pk[0]:pk[-1]]
        if len(aa)<100: continue
        # Reject blocks with endpoint discontinuity if repeated periodically.
        endpoint=abs(float(aa[-1]-aa[0]))
        pulse=max(float(np.ptp(aa)),1.0)
        score=20*cv + endpoint/pulse
        if best is None or score<best[0]:
            best=(score,pk.copy(),aa.copy(),pp.copy(),start+pk[0],start+pk[-1])
    if best is None:
        raise RuntimeError("No stable five-cycle block found.")
    _,pk,aa,pp,g0,g1=best
    return aa,pp,pk,g0,g1

class RawABP:
    def __init__(self,x,fs=FS):
        self.x=np.asarray(x,float); self.fs=fs; self.T=len(self.x)/fs
    def __call__(self,t):
        pos=((float(t)%self.T)*self.fs)
        i=int(math.floor(pos)); f=pos-i; j=(i+1)%len(self.x)
        return float((1-f)*self.x[i]+f*self.x[j])

def make_rhs(paw, target_ml_min=0.0, ramp_s=0.0, ret=0.0):
    p=make_tbi_parameters(); pp=copy.copy(p)
    target=target_ml_min/60.0
    def rhs(t,y):
        pp.Pa=paw(t)
        delivered=target if ramp_s<=0 else target*min(max(t/ramp_s,0.0),1.0)
        return model.evaluate(y,pp,model.Aspiration("Pv",delivered),ret*delivered)[0]
    return rhs

# local import kept here so header stays obvious
import copy

def periodic_state(initial,paw):
    rhs=make_rhs(paw)
    state=np.asarray(initial,float).copy(); prev=None; stable=0
    for k in range(120):
        a,b=k*paw.T,(k+1)*paw.T
        te=np.arange(a,b+1e-9,0.01)
        if te[-1] < b: te=np.append(te,b)
        sol=solve_ivp(rhs,(a,b),state,method="BDF",rtol=1e-8,atol=1e-10,max_step=.01,t_eval=te)
        if not sol.success: raise RuntimeError(sol.message)
        state=sol.y[:,-1]; w=sol.y[0,:-1]
        if prev is not None:
            n=min(len(prev),len(w)); d=float(np.max(np.abs(w[:n]-prev[:n])))
            md=float(abs(np.mean(w[:n])-np.mean(prev[:n])))
            stable=stable+1 if d<5e-4 and md<1e-4 else 0
        prev=w.copy()
        if stable>=3:return state
    raise RuntimeError("Baseline did not reach periodic state.")

def simulate_block(state,paw,flow=0.0,ramp=0.0,ret=0.0,duration=None):
    if duration is None: duration=paw.T
    rhs=make_rhs(paw,flow,ramp,ret)
    te=np.arange(0,duration+1e-9,.01)
    if te[-1] < duration: te=np.append(te,duration)
    sol=solve_ivp(rhs,(0,duration),np.asarray(state,float),method="BDF",rtol=1e-8,atol=1e-10,
                  max_step=.01,t_eval=te)
    if not sol.success: raise RuntimeError(sol.message)
    return sol

def resample_model_to_patient(sol_t,sol_icp,n):
    target=np.arange(n)/FS
    return np.interp(target,sol_t,sol_icp)

def beat_metrics(x,pk):
    vals=[]
    for i in range(len(pk)-1):
        b=x[pk[i]:pk[i+1]]
        vals.append(float(np.ptp(b)))
    return vals

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--dat",type=Path,required=True)
    ap.add_argument("--output-dir",type=Path,required=True)
    args=ap.parse_args(); out=args.output_dir; out.mkdir(parents=True,exist_ok=True)

    abp,ecg,icp=load_record(args.dat)
    q,ncand=select_clean_window(abp,icp)
    aa,meas_icp,pk,g0,g1=select_five_cycles(abp,icp,q)
    aa_shift=aa+(100.0-float(np.mean(aa)))
    paw=RawABP(aa_shift)

    base=json.loads(Path("reference_results/reports/tbi_baseline.json").read_text())
    state0=np.asarray(base["terminal_window_mean_state"],float)
    pstate=periodic_state(state0,paw)
    pred=simulate_block(pstate,paw)
    pred50=resample_model_to_patient(pred.t,pred.y[0],len(meas_icp))

    # Same-time comparison: no time shift and no amplitude fitting.
    m0=meas_icp-float(np.mean(meas_icp))
    p0=pred50-float(np.mean(pred50))
    corr=float(np.corrcoef(m0,p0)[0,1])
    rmse_centered=float(np.sqrt(np.mean((m0-p0)**2)))
    measured_amp=beat_metrics(meas_icp,pk-pk[0])
    predicted_amp=beat_metrics(pred50,pk-pk[0])

    # Diagnostic optimal lag only; it is NOT used to alter any plotted waveform.
    maxlag=int(.5*FS); best=(-2.0,0)
    for lag in range(-maxlag,maxlag+1):
        if lag<0: x,y=m0[-lag:],p0[:len(p0)+lag]
        elif lag>0: x,y=m0[:-lag],p0[lag:]
        else: x,y=m0,p0
        if len(x)>10:
            r=float(np.corrcoef(x,y)[0,1])
            if r>best[0]:best=(r,lag)

    trans=simulate_block(pstate,paw,ASP_ML_MIN,RAMP_S,1.0,TRANSIENT_S)

    # Save exact plotted source data at native patient sampling and solver sampling.
    patient_time=np.arange(len(aa))/FS
    np.savez_compressed(out/"charis1_validation_and_aspiration.npz",
        patient_time_s=patient_time, raw_abp_mmhg=aa, shifted_abp_mmhg=aa_shift,
        measured_icp_mmhg=meas_icp, predicted_icp_mmhg=pred50,
        model_block_time_s=pred.t, model_block_icp_mmhg=pred.y[0],
        aspiration_time_s=trans.t, aspiration_icp_mmhg=trans.y[0])

    summary={
      "record":"charis1",
      "patient_header":"19-year-old male; SAH/SDH/TBI; outcome Rehab",
      "sampling_hz":FS,
      "screening":{"target_mean_icp_mmhg":list(TARGET_ICP_RANGE),"eligible_20s_windows":ncand,
                   "selected_start_s":q["start"]/FS,"selected_mean_abp_mmhg":q["mean_abp"],
                   "selected_mean_icp_mmhg":q["mean_icp"],"selected_hr_bpm":q["hr"],
                   "selected_rr_cv":q["rr_cv"],"selected_abp_drift_mmhg_min":q["abp_drift"],
                   "selected_icp_drift_mmhg_min":q["icp_drift"]},
      "five_cycles":{"global_start_s":g0/FS,"global_stop_s":g1/FS,
                     "duration_s":len(aa)/FS,"raw_abp_mean_mmhg":float(np.mean(aa)),
                     "raw_abp_min_mmhg":float(np.min(aa)),"raw_abp_max_mmhg":float(np.max(aa)),
                     "model_input_constant_shift_mmhg":float(100-np.mean(aa)),
                     "measured_icp_mean_mmhg":float(np.mean(meas_icp)),
                     "predicted_icp_mean_mmhg":float(np.mean(pred50)),
                     "measured_beat_pulse_amplitudes_mmhg":measured_amp,
                     "predicted_beat_pulse_amplitudes_mmhg":predicted_amp,
                     "measured_median_pulse_amplitude_mmhg":float(np.median(measured_amp)),
                     "predicted_median_pulse_amplitude_mmhg":float(np.median(predicted_amp)),
                     "same_time_mean_centered_correlation":corr,
                     "same_time_mean_centered_rmse_mmhg":rmse_centered,
                     "diagnostic_optimal_lag_s":best[1]/FS,
                     "diagnostic_optimal_lag_correlation":best[0]},
      "aspiration":{"site":"Pv","flow_ml_min":ASP_ML_MIN,"equal_lower_svc_return":True,
                    "ramp_s":RAMP_S,"duration_s":TRANSIENT_S}
    }
    (out/"charis1_summary.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    main()
