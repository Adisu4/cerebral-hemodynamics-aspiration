#!/usr/bin/env python3
from __future__ import annotations

import argparse, csv, io, json, math, re, tarfile, copy
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks, savgol_filter

from cerebral_hemodynamics_aspiration import model
from cerebral_hemodynamics_aspiration.parameters import make_tbi_parameters

FS=100.0
REST_S=60.0
NPHASE=256
NTRACE=12
NHARM=10
RX=re.compile(r"(CON\d{3})([abc]?)\.csv$", re.I)

def colidx(header,target):
    clean=[re.sub(r"[^a-z0-9]+","",h.lower()) for h in header]
    opts={"pa":{"pa","aorticpressure","aortic"},"ecg":{"ecg","ekg"}}[target]
    for i,x in enumerate(clean):
        if x in opts:return i
    if target=="pa":
        for i,x in enumerate(clean):
            if "aortic" in x or x.endswith("pa"): return i
    raise ValueError(header)

def sid(name):
    m=RX.search(Path(name).name)
    if not m: raise ValueError(name)
    return m.group(1).upper(),m.group(2).lower()

def read_rest(f):
    text=io.TextIOWrapper(f,encoding="utf-8",errors="replace",newline="")
    r=csv.reader(text); h=next(r); ip=colidx(h,"pa")
    vals=[]
    for row in r:
        if len(vals)>=int(FS*REST_S): break
        try: vals.append(float(row[ip])/10.0)
        except: pass
    return np.asarray(vals,float)

def resample(x,n=NPHASE):
    return np.interp(np.linspace(0,1,n,endpoint=False),
                     np.linspace(0,1,len(x),endpoint=True),x)

def notch_score(w):
    n=len(w); pp=np.ptp(w)
    if pp<=0:return 0.0,np.nan
    lo,hi=int(.12*n),int(.65*n)
    mins,_=find_peaks(-w[lo:hi],prominence=max(.25,.004*pp),distance=max(2,int(.03*n)))
    best=(0.0,np.nan)
    for rel in mins:
        m=lo+rel
        a,b=m+max(2,int(.02*n)),min(n-2,int(.85*n))
        if a>=b: continue
        maxs,_=find_peaks(w[a:b],prominence=max(.15,.0025*pp),distance=max(2,int(.025*n)))
        if maxs.size:
            s=a+maxs[0]
            sc=float((w[s]-w[m])/pp)
            if sc>best[0]: best=(sc,m/n)
    return best

def analyze(pa,name):
    """Choose the cleanest 15-s window within the initial 60-s resting record."""
    if len(pa)<4500 or np.mean(np.isfinite(pa))<.995:return None
    pa=np.asarray(pa,float)
    win=int(15*FS); step=int(5*FS)
    candidates=[]
    for a in range(0,max(1,len(pa)-win+1),step):
        x=pa[a:a+win]
        if len(x)<win or np.mean(np.isfinite(x))<.995: continue
        x=x[np.isfinite(x)]
        # Reject obvious catheter flushes/disconnections/spikes.
        if np.quantile(x,.005)<35 or np.quantile(x,.995)>220: continue
        s=savgol_filter(x,11,3)
        peaks,_=find_peaks(s,distance=int(.45*FS),prominence=8.0)
        if len(peaks)<9: continue
        rr=np.diff(peaks)/FS
        med=float(np.median(rr))
        if med<=0: continue
        hr=60/med
        # Robust rhythm variability; single detector errors should not reject an otherwise clean window.
        rr_mad=float(np.median(np.abs(rr-med)))
        robust_cv=(1.4826*rr_mad)/med
        if not 45<=hr<=100 or robust_cv>.06: continue
        beats=[]
        for i,r in enumerate(rr):
            if abs(r-med)/med>.12: continue
            b=s[peaks[i]:peaks[i+1]+1]
            if len(b)>=40: beats.append(resample(b))
        if len(beats)<8: continue
        B=np.vstack(beats); w=B.mean(0)
        corrs=[np.corrcoef(b,w)[0,1] for b in B if np.std(b)>0 and np.std(w)>0]
        corr=float(np.median(corrs)) if corrs else 0.0
        if corr<.93: continue
        mean=float(np.mean(x)); pp=float(np.ptp(w))
        if not 60<=mean<=150 or not 20<=pp<=100: continue
        t=np.arange(len(x))/FS; drift=float(np.polyfit(t,x,1)[0]*60)
        if abs(drift)>12: continue
        ns,np_=notch_score(w)
        # Prefer stable rhythm, repeatable beats, low drift, and a recognizable notch.
        score=8*robust_cv + 0.06*abs(drift) + 3*max(0,.98-corr) - 2*ns
        candidates.append(dict(subject=sid(name)[0],filename=Path(name).name,
            window_start_s=float(a/FS),window_end_s=float((a+win)/FS),
            hr_bpm=float(hr),rr_cv=float(robust_cv),mean_pa_mmhg=mean,
            pulse_pressure_mmhg=pp,median_beat_correlation=corr,
            drift_mmhg_min=drift,notch_score=float(ns),
            notch_phase=float(np_) if np.isfinite(np_) else None,
            n_beats=len(beats),wave=w,_score=float(score)))
    if not candidates:return None
    return min(candidates,key=lambda r:r["_score"])

def choose(recs,n=NTRACE):
    if not recs: raise RuntimeError("No traces passed QC")
    pool=[r for r in recs if r["notch_score"]>=.01]
    if len(pool)<5: pool=recs
    fields=["hr_bpm","mean_pa_mmhg","pulse_pressure_mmhg"]
    X=np.array([[r[f] for f in fields] for r in pool],float)
    med=np.median(X,0); mad=np.median(np.abs(X-med),0); mad[mad<1e-9]=1
    d=np.sqrt(np.sum(((X-med)/(1.4826*mad))**2,1))
    q=np.array([r["_score"] for r in pool])
    order=np.argsort(d+q)
    return [pool[i] for i in order[:min(n,len(pool))]]

def fourier(w,h=NHARM):
    p=np.arange(len(w))/len(w); cols=[np.ones_like(p)]
    for k in range(1,h+1): cols += [np.cos(2*np.pi*k*p),np.sin(2*np.pi*k*p)]
    A=np.column_stack(cols); beta,*_=np.linalg.lstsq(A,w,rcond=None)
    fit=A@beta
    coeff=[{"harmonic":0,"a_cos_mmhg":float(beta[0]),"b_sin_mmhg":0.0}]
    for k in range(1,h+1):
        coeff.append({"harmonic":k,"a_cos_mmhg":float(beta[2*k-1]),"b_sin_mmhg":float(beta[2*k])})
    rmse=float(np.sqrt(np.mean((fit-w)**2)))
    r2=float(1-np.sum((fit-w)**2)/np.sum((w-w.mean())**2))
    return fit,coeff,rmse,r2

class PaWave:
    def __init__(self,c,T): self.c,self.T=c,T
    def __call__(self,t):
        th=2*np.pi*((t%self.T)/self.T); v=self.c[0]["a_cos_mmhg"]
        for r in self.c[1:]:
            k=r["harmonic"]; v += r["a_cos_mmhg"]*np.cos(k*th)+r["b_sin_mmhg"]*np.sin(k*th)
        return float(v)

def periodic(initial,paw,flow=0.0,ret=0.0,maxcycles=120):
    p=make_tbi_parameters(); pp=copy.copy(p)
    asp=model.Aspiration("Pv",flow/60); rflow=ret*flow/60
    T=paw.T; state=np.array(initial,float); prev=None; stable=0
    last=None
    def rhs(t,y):
        pp.Pa=paw(t)
        return model.evaluate(y,pp,asp,rflow)[0]
    for cyc in range(1,maxcycles+1):
        a,b=(cyc-1)*T,cyc*T
        te=np.linspace(a,b,NPHASE+1)
        sol=solve_ivp(rhs,(a,b),state,method="BDF",rtol=1e-8,atol=1e-10,
                      max_step=min(.01,T/100),t_eval=te)
        if not sol.success: raise RuntimeError(sol.message)
        state=sol.y[:,-1]; icp=sol.y[0,:-1]
        if prev is not None:
            wd=float(np.max(np.abs(icp-prev))); md=float(abs(icp.mean()-prev.mean()))
            stable = stable+1 if wd<5e-4 and md<1e-4 else 0
        else: wd=md=float("inf")
        prev=icp.copy(); last=(sol,cyc,wd,md)
        if stable>=3: break
    if stable<3: raise RuntimeError("No periodic convergence")
    sol,cyc,wd,md=last; y=sol.y[:,:-1]; tt=sol.t[:-1]; icp=y[0]
    invalid=0; branches=set()
    for j,t in enumerate(tt):
        pp.Pa=paw(float(t))
        try:
            _,aux=model.evaluate(y[:,j],pp,asp,rflow); branches.add(str(aux["Rvs_branch"]))
        except model.ModelDomainError: invalid+=1
    pa=np.array([paw(float(t)) for t in tt])
    return dict(cycles=cyc,wave_diff_mmhg=wd,mean_diff_mmhg=md,
                mean_icp_mmhg=float(icp.mean()),min_icp_mmhg=float(icp.min()),
                max_icp_mmhg=float(icp.max()),icp_pulse_amplitude_mmhg=float(np.ptp(icp)),
                mean_pa_mmhg=float(pa.mean()),min_pa_mmhg=float(pa.min()),max_pa_mmhg=float(pa.max()),
                invalid_domain_samples=invalid,rvs_branches=sorted(branches),
                phase=((tt-tt[0])/T).tolist(),icp_wave=icp.tolist(),pa_wave=pa.tolist())

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--archive",type=Path,required=True); ap.add_argument("--output-dir",type=Path,required=True)
    a=ap.parse_args(); out=a.output_dir; out.mkdir(parents=True,exist_ok=True)
    recs=[]; seen=set()
    with tarfile.open(a.archive,"r|gz") as tar:
        for m in tar:
            if not m.isfile() or not m.name.lower().endswith(".csv"): continue
            try: s,suf=sid(m.name)
            except: continue
            if suf not in ("","a") or s in seen: continue
            f=tar.extractfile(m)
            if f is None: continue
            seen.add(s)
            try:r=analyze(read_rest(f),m.name)
            except Exception:r=None
            if r:recs.append(r)
    chosen=choose(recs)
    W=np.vstack([r["wave"] for r in chosen]); centered=W-W.mean(1,keepdims=True)
    rep=100+centered.mean(0); sd=centered.std(0,ddof=1)
    medhr=float(np.median([r["hr_bpm"] for r in chosen])); T=60/medhr
    fit,coeff,rmse,r2=fourier(rep); paw=PaWave(coeff,T)

    with (out/"selected_traces.csv").open("w",newline="") as f:
        fields=[k for k in chosen[0] if k not in ("wave","_score")]; w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for r in chosen:w.writerow({k:r[k] for k in fields})
    with (out/"fourier_coefficients.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["harmonic","a_cos_mmhg","b_sin_mmhg"]);w.writeheader();w.writerows(coeff)
    with (out/"representative_waveform.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["phase","time_s","representative_pa_mmhg","fourier_pa_mmhg","between_subject_sd_mmhg"])
        for i in range(NPHASE):w.writerow([i/NPHASE,i/NPHASE*T,rep[i],fit[i],sd[i]])

    br=json.loads(Path("reference_results/reports/tbi_baseline.json").read_text())
    pr=json.loads(Path("reference_results/reports/dose_pv_240.json").read_text())
    if float(pr.get("return_fraction",-1))!=1.0: raise RuntimeError("Pv240 reference is not full-return primary")
    bs=np.asarray(br["terminal_window_mean_state"],float); ps=np.asarray(pr["terminal_window_mean_state"],float)
    static_b=float(br["summary"]["icp_mean_mmhg"]); static_p=float(pr["summary"]["icp_mean_mmhg"])
    pb=periodic(bs,paw,0,0); pp=periodic(ps,paw,240,1)
    dstat=static_b-static_p; dp=pb["mean_icp_mmhg"]-pp["mean_icp_mmhg"]

    summary={
      "dataset":{"name":"CONTRAST","doi":"10.5061/dryad.f76nv","zenodo_record":"4992130","sampling_hz":FS,
                 "rest_window_s":REST_S,"subjects_screened":len(seen),"passed_qc":len(recs),
                 "selected_subjects":[r["subject"] for r in chosen]},
      "waveform":{"median_hr_bpm":medhr,"period_s":T,"mean_mmhg":float(rep.mean()),
                  "min_mmhg":float(rep.min()),"max_mmhg":float(rep.max()),"pulse_pressure_mmhg":float(np.ptp(rep)),
                  "fourier_harmonics":NHARM,"fourier_rmse_mmhg":rmse,"fourier_r2":r2},
      "static":{"baseline_icp_mmhg":static_b,"pv240_icp_mmhg":static_p,"delta_icp_mmhg":dstat},
      "pulsatile":{"baseline":{k:v for k,v in pb.items() if not k.endswith("_wave") and k!="phase"},
                   "pv240":{k:v for k,v in pp.items() if not k.endswith("_wave") and k!="phase"},
                   "delta_icp_cycle_mean_mmhg":dp,"difference_from_static_delta_mmhg":dp-dstat}
    }
    (out/"pilot_summary.json").write_text(json.dumps(summary,indent=2))
    np.savez_compressed(out/"pulsatile_cycles.npz",
        baseline_phase=np.array(pb["phase"]),baseline_icp=np.array(pb["icp_wave"]),baseline_pa=np.array(pb["pa_wave"]),
        pv240_phase=np.array(pp["phase"]),pv240_icp=np.array(pp["icp_wave"]),pv240_pa=np.array(pp["pa_wave"]))
    report=[
      "# CONTRAST pulsatility pilot","",
      f"- Screened subjects: {len(seen)}",f"- Passed QC: {len(recs)}",f"- Selected representative traces: {len(chosen)}",
      f"- Selected subjects: {', '.join(r['subject'] for r in chosen)}",
      f"- Representative HR: {medhr:.1f} bpm",
      f"- Representative Pa: {rep.min():.2f}-{rep.max():.2f} mmHg, mean {rep.mean():.2f}",
      f"- Fourier fit: {NHARM} harmonics, RMSE {rmse:.3f} mmHg, R2 {r2:.6f}","",
      "## ICP results",
      f"- Static baseline ICP: {static_b:.5f} mmHg",
      f"- Pulsatile cycle-mean baseline ICP: {pb['mean_icp_mmhg']:.5f} mmHg",
      f"- Static Pv240 full-return ICP: {static_p:.5f} mmHg",
      f"- Pulsatile cycle-mean Pv240 full-return ICP: {pp['mean_icp_mmhg']:.5f} mmHg",
      f"- Static delta ICP: {dstat:.5f} mmHg",
      f"- Pulsatile cycle-mean delta ICP: {dp:.5f} mmHg",
      f"- Change in delta ICP: {dp-dstat:+.5f} mmHg",
      f"- Baseline ICP pulse amplitude: {pb['icp_pulse_amplitude_mmhg']:.5f} mmHg",
      f"- Pv240 ICP pulse amplitude: {pp['icp_pulse_amplitude_mmhg']:.5f} mmHg",
      f"- Invalid-domain samples baseline/Pv240: {pb['invalid_domain_samples']}/{pp['invalid_domain_samples']}","",
      "Exploratory pilot only; the maintained publication model is unchanged."
    ]
    (out/"PILOT_REPORT.md").write_text("\n".join(report)+"\n")
    print("\n".join(report))

if __name__=="__main__": main()
