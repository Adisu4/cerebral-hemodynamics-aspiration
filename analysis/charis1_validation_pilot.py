#!/usr/bin/env python3
"""CHARIS1 same-patient ABP->model ICP validation + aspiration pilot.

Selection is based ONLY on measured-signal quality and relevance to the model's
post-traumatic intracranial-hypertension state. Model output is not consulted
when selecting the segment.

Plotted/compared waveforms use raw acquired CHARIS samples. No smoothing,
filtering, Fourier fitting, or ensemble averaging is applied to ABP or measured
ICP. The CHARIS acquisition itself used the clinical monitor's 25-Hz filtered
outputs and sampled at 50 Hz, as documented by PhysioNet.

The selected block contains exactly seven consecutive cardiac cycles.
"""

from __future__ import annotations

import argparse, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from cerebral_hemodynamics_aspiration import model
from cerebral_hemodynamics_aspiration.parameters import make_tbi_parameters\nfrom cerebral_hemodynamics_aspiration.simulation import integrate

FS = 50.0
N_CYCLES = 7
FLOW_ML_MIN = 240.0
RAMP_S = 60.0

# charis1.hea:
# ABP gain=85.953, baseline=-1332; ICP gain=94.8784, baseline=-2316.
ABP_GAIN = 85.953
ABP_BASE = -1332.0
ICP_GAIN = 94.8784
ICP_BASE = -2316.0


def load_charis1(path: Path):
    raw = np.memmap(path, dtype="<i2", mode="r")
    if raw.size % 3:
        raise RuntimeError("WFDB data length is not divisible by 3 channels")
    x = raw.reshape(-1, 3)
    abp_d = np.asarray(x[:, 0], dtype=np.float64)
    ecg_d = np.asarray(x[:, 1], dtype=np.float64)
    icp_d = np.asarray(x[:, 2], dtype=np.float64)
    abp = (abp_d - ABP_BASE) / ABP_GAIN
    icp = (icp_d - ICP_BASE) / ICP_GAIN
    return abp, ecg_d, icp


def resample_cycle(x: np.ndarray, n=128):
    return np.interp(
        np.linspace(0, 1, n, endpoint=False),
        np.linspace(0, 1, len(x), endpoint=False),
        x,
    )


def beat_repeatability(signal: np.ndarray, peaks: np.ndarray):
    beats=[]
    for i in range(len(peaks)-1):
        b=signal[peaks[i]:peaks[i+1]]
        if len(b)>=15:
            beats.append(resample_cycle(b))
    if len(beats)<N_CYCLES:
        return -1.0, np.nan
    B=np.vstack(beats[:N_CYCLES])
    mean=B.mean(0)
    corrs=[]
    amps=[]
    for b in B:
        if np.std(b)>0 and np.std(mean)>0:
            corrs.append(np.corrcoef(b,mean)[0,1])
        amps.append(np.ptp(b))
    return float(np.median(corrs)), float(np.median(amps))


def evaluate_candidate(abp, icp, start, window_s=14.0):
    n=int(window_s*FS)
    a=int(start); b=min(len(abp),a+n)
    A=abp[a:b]; I=icp[a:b]
    if len(A)<n: return None
    if not (np.isfinite(A).all() and np.isfinite(I).all()): return None
    # Reject obvious disconnects/flushes/saturation.
    if np.min(A)<35 or np.max(A)>220: return None
    if np.min(I)<-5 or np.max(I)>60: return None
    mean_abp=float(np.mean(A)); mean_icp=float(np.mean(I))
    # Pre-specified relevance window: model is a TBI intracranial-hypertension state.
    if not (75<=mean_abp<=130): return None
    if not (18<=mean_icp<=30): return None

    peaks,_=find_peaks(A, distance=int(.42*FS), prominence=18.0)
    if len(peaks)<N_CYCLES+1: return None

    # Examine every possible seven-cycle block in this candidate window.
    best=None
    for j in range(len(peaks)-N_CYCLES):
        p=peaks[j:j+N_CYCLES+1]
        rr=np.diff(p)/FS
        hr=60/np.mean(rr)
        if not 45<=hr<=110: continue
        rr_cv=float(np.std(rr)/np.mean(rr))
        if rr_cv>.07: continue
        s0,s1=int(p[0]),int(p[-1])
        AA=A[s0:s1]; II=I[s0:s1]
        if len(AA)<N_CYCLES*20: continue
        abp_rep,abp_amp=beat_repeatability(A,p)
        icp_rep,icp_amp=beat_repeatability(I,p)
        if abp_rep<.94 or icp_rep<.65: continue
        if not 20<=abp_amp<=100: continue
        if not .3<=icp_amp<=12: continue
        endpoint_jump=abs(float(AA[-1]-AA[0]))
        # Quality score only; no model prediction is involved.
        score=(8*rr_cv + 3*(1-abp_rep) + 1.5*(1-icp_rep)
               + .02*endpoint_jump)
        rec=dict(
            global_start_sample=int(a+s0),
            global_stop_sample=int(a+s1),
            start_time_s=float((a+s0)/FS),
            stop_time_s=float((a+s1)/FS),
            cycles=N_CYCLES,
            hr_bpm=float(hr),
            rr_cv=rr_cv,
            mean_abp_mmhg=float(np.mean(AA)),
            mean_icp_mmhg=float(np.mean(II)),
            median_abp_pulse_mmhg=abp_amp,
            median_icp_pulse_mmhg=icp_amp,
            abp_repeatability=abp_rep,
            icp_repeatability=icp_rep,
            endpoint_abp_jump_mmhg=endpoint_jump,
            quality_score=float(score),
            local_peaks=(p-j*0).tolist(),
        )
        if best is None or score<best["quality_score"]:
            best=rec
    return best


def select_segment(abp,icp):
    # Scan the full record at 30-s spacing. Selection is frozen before model run.
    candidates=[]
    step=int(30*FS)
    for start in range(0,len(abp)-int(14*FS),step):
        rec=evaluate_candidate(abp,icp,start)
        if rec is not None:
            candidates.append(rec)
    if not candidates:
        raise RuntimeError("No CHARIS1 segment met the pre-specified quality/relevance criteria")
    candidates.sort(key=lambda r:r["quality_score"])
    return candidates[0], candidates


class RawBlock:
    def __init__(self,values,fs):
        self.v=np.asarray(values,float)
        self.fs=float(fs)
        self.T=len(self.v)/self.fs
    def __call__(self,t):
        pos=((float(t)%self.T)*self.fs)
        i=int(math.floor(pos))%len(self.v)
        frac=pos-math.floor(pos)
        j=(i+1)%len(self.v)
        return float((1-frac)*self.v[i]+frac*self.v[j])


def rhs_factory(paw, flow_ml_min=0.0, ramp_s=0.0, return_fraction=0.0):
    p=make_tbi_parameters(); pp=__import__("copy").copy(p)
    target=flow_ml_min/60.0
    def rhs(t,y):
        pp.Pa=paw(t)
        if target<=0: delivered=0.0
        elif ramp_s>0: delivered=target*min(max(t/ramp_s,0.0),1.0)
        else: delivered=target
        return model.evaluate(
            y, pp, model.Aspiration("Pv",delivered), return_fraction*delivered
        )[0]
    return rhs


def settle(initial,paw):
    state=np.asarray(initial,float).copy()
    rhs=rhs_factory(paw)
    prev=None; stable=0
    for k in range(1,121):
        a=(k-1)*paw.T; b=k*paw.T
        te=np.arange(a,b+0.0101,0.02); te[-1]=b
        sol=solve_ivp(rhs,(a,b),state,method="BDF",rtol=1e-8,atol=1e-10,
                      max_step=.02,t_eval=te)
        if not sol.success: raise RuntimeError(sol.message)
        state=sol.y[:,-1]; wave=sol.y[0,:-1]
        if prev is not None:
            n=min(len(prev),len(wave))
            wd=float(np.max(np.abs(wave[:n]-prev[:n])))
            md=float(abs(np.mean(wave[:n])-np.mean(prev[:n])))
            stable=stable+1 if wd<5e-4 and md<1e-4 else 0
        prev=wave.copy()
        if stable>=3: return state,k
    raise RuntimeError("Periodic baseline did not converge")


def simulate_block(initial,paw,flow=0.0,ramp=0.0,ret=0.0,duration=None):
    if duration is None: duration=paw.T
    rhs=rhs_factory(paw,flow,ramp,ret)
    te=np.arange(0,duration+0.0101,0.02); te[-1]=duration
    sol=solve_ivp(rhs,(0,duration),np.asarray(initial,float),method="BDF",
                  rtol=1e-8,atol=1e-10,max_step=.02,t_eval=te)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.t,sol.y


def interp_to_measured(t_model, y_model, n):
    tm=np.arange(n)/FS
    return np.interp(tm,t_model,y_model)


def metrics(pred,meas):
    pred=np.asarray(pred); meas=np.asarray(meas)
    corr=float(np.corrcoef(pred-np.mean(pred),meas-np.mean(meas))[0,1])
    shape_rmse=float(np.sqrt(np.mean(((pred-np.mean(pred))-(meas-np.mean(meas)))**2)))
    return dict(
      measured_mean_icp_mmhg=float(np.mean(meas)),
      predicted_mean_icp_mmhg=float(np.mean(pred)),
      mean_error_mmhg=float(np.mean(pred)-np.mean(meas)),
      measured_peak_to_peak_mmhg=float(np.ptp(meas)),
      predicted_peak_to_peak_mmhg=float(np.ptp(pred)),
      amplitude_ratio_predicted_to_measured=float(np.ptp(pred)/np.ptp(meas)),
      demeaned_waveform_correlation=corr,
      demeaned_rmse_mmhg=shape_rmse,
    )


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--dat",type=Path,required=True)
    ap.add_argument("--output-dir",type=Path,required=True)
    a=ap.parse_args(); out=a.output_dir; out.mkdir(parents=True,exist_ok=True)

    abp,ecg,icp=load_charis1(a.dat)
    chosen,cands=select_segment(abp,icp)
    s0=chosen["global_start_sample"]; s1=chosen["global_stop_sample"]
    A=abp[s0:s1].copy(); I=icp[s0:s1].copy()
    paw=RawBlock(A,FS)

    # Freeze selected raw segment before model run.
    with (out/"selected_raw_7cycles.csv").open("w") as f:
        f.write("time_s,abp_mmhg,measured_icp_mmhg\n")
        for i,(x,y) in enumerate(zip(A,I)):
            f.write(f"{i/FS:.5f},{x:.8f},{y:.8f}\n")

    baseline=json.loads(Path("reference_results/reports/tbi_baseline.json").read_text())
    initial=np.asarray(baseline["terminal_window_mean_state"],float)
    state,settle_blocks=settle(initial,paw)

    tb,yb=simulate_block(state,paw)
    pred=interp_to_measured(tb,yb[0],len(I))
    val=metrics(pred,I)

    # Aspiration transient on the same patient ABP.  Keep raw pulsatility.
    # 60-s ramp + 20 s after ramp so we can inspect seven-cycle windows.
    ta,ya=simulate_block(state,paw,FLOW_ML_MIN,RAMP_S,1.0,duration=80.0)
    # Seven-cycle windows at onset and immediately after full-ramp completion.
    win_n=len(A)
    onset_t=np.arange(win_n)/FS
    onset_pred=np.interp(onset_t,ta,ya[0])
    post_start=60.0
    post_t=post_start+np.arange(win_n)/FS
    post_pred=np.interp(post_t,ta,ya[0])
    post_mean=float(np.mean(post_pred))
    onset_mean=float(np.mean(onset_pred))

    # Domain check.
    p=make_tbi_parameters(); import copy; pp=copy.copy(p)
    invalid=0; branches=set(); target=FLOW_ML_MIN/60.0
    for j in range(0,len(ta),5):
        pp.Pa=paw(float(ta[j]))
        delivered=target*min(max(float(ta[j])/RAMP_S,0.0),1.0)
        try:
            _,aux=model.evaluate(ya[:,j],pp,model.Aspiration("Pv",delivered),delivered)
            branches.add(str(aux["Rvs_branch"]))
        except model.ModelDomainError:
            invalid+=1

    np.savez_compressed(
      out/"charis1_validation_and_aspiration.npz",
      time_7cycles_s=np.arange(len(A))/FS,
      measured_abp=A,
      measured_icp=I,
      predicted_icp=pred,
      aspiration_time_s=ta,
      aspiration_icp=ya[0],
      onset_7cycle_icp=onset_pred,
      post60_7cycle_icp=post_pred,
    )
    summary={
      "patient":{"record":"charis1","age":19,"sex":"M","diagnoses":"SAH/SDH/TBI"},
      "selection":{
        "rule":"best measured-signal-quality seven-cycle block among segments with mean ABP 75-130 mmHg and mean ICP 18-30 mmHg; model output not used for selection",
        "number_of_qualifying_candidates":len(cands),
        **{k:v for k,v in chosen.items() if k!="local_peaks"},
      },
      "validation":{
        "raw_waveforms":"no smoothing/filtering/Fourier/ensemble averaging applied in this analysis",
        "model_input":"absolute measured radial ABP, no mean normalization",
        "settling_blocks":settle_blocks,
        **val,
      },
      "aspiration":{
        "site":"Pv","flow_ml_min":FLOW_ML_MIN,"equal_lower_svc_return":True,
        "ramp_s":RAMP_S,
        "mean_icp_first_7cycles_mmhg":onset_mean,
        "mean_icp_7cycles_after_60s_mmhg":post_mean,
        "mean_drop_mmhg":float(onset_mean-post_mean),
        "onset_pulse_amplitude_mmhg":float(np.ptp(onset_pred)),
        "post60_pulse_amplitude_mmhg":float(np.ptp(post_pred)),
        "global_min_first80s_mmhg":float(np.min(ya[0])),
        "invalid_domain_samples":invalid,
        "rvs_branches":sorted(branches),
      }
    }
    (out/"summary.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    main()
