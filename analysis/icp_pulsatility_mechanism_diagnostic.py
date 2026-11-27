#!/usr/bin/env python3
"""Diagnose why the maintained Gadda-based model smooths cardiac ICP pulsatility.

This is an analysis-only script. It does not change the publication model.

Outputs:
1) Exact small-signal frequency response Pa -> Pic by numerical linearization
   around the maintained TBI equilibrium.
2) Controlled parameter/structure tests at cardiac frequency.
3) Comparison with the selected CHARIS1 seven-cycle measured ABP/ICP transfer.
"""

from __future__ import annotations
import copy, csv, io, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from cerebral_hemodynamics_aspiration import model
from cerebral_hemodynamics_aspiration.parameters import make_tbi_parameters

FS = 50.0
CHARIS_START = 390159
CHARIS_STOP = 390369


def load_charis1(dat: Path):
    raw = np.memmap(dat, dtype="<i2", mode="r")
    n = (raw.size // 3) * 3
    x = np.asarray(raw[:n]).reshape(-1,3)
    abp = (x[:,0].astype(float) - (-1332.0)) / 85.953
    icp = (x[:,2].astype(float) - (-2316.0)) / 94.8784
    return abp[CHARIS_START:CHARIS_STOP], icp[CHARIS_START:CHARIS_STOP]


def linearize(y0, p, eps_y=1e-5, eps_pa=1e-4):
    y0=np.asarray(y0,float)
    f0=model.evaluate(y0,p)[0]
    n=len(y0)
    A=np.zeros((n,n))
    for j in range(n):
        h=eps_y*max(1.0,abs(y0[j]))
        yp=y0.copy(); ym=y0.copy()
        yp[j]+=h; ym[j]-=h
        fp=model.evaluate(yp,p)[0]; fm=model.evaluate(ym,p)[0]
        A[:,j]=(fp-fm)/(2*h)
    pp=copy.copy(p); pm=copy.copy(p)
    pp.Pa += eps_pa; pm.Pa -= eps_pa
    B=(model.evaluate(y0,pp)[0]-model.evaluate(y0,pm)[0])/(2*eps_pa)
    return A,B


def frequency_response(A,B,freqs):
    C=np.zeros((1,A.shape[0])); C[0,0]=1.0
    I=np.eye(A.shape[0])
    rows=[]
    for f in freqs:
        H=(C @ np.linalg.solve(1j*2*np.pi*f*I-A,B[:,None]))[0,0]
        rows.append(dict(freq_hz=float(f), gain_mmhg_per_mmhg=float(abs(H)),
                         phase_deg=float(np.angle(H,deg=True))))
    return rows


def harmonics(x, fs=FS, f0=100/60, nh=6):
    x=np.asarray(x,float); n=len(x)
    X=np.fft.rfft(x-np.mean(x)); fr=np.fft.rfftfreq(n,1/fs)
    amp=2*np.abs(X)/n; phase=np.angle(X)
    out=[]
    for h in range(1,nh+1):
        i=int(np.argmin(np.abs(fr-h*f0)))
        out.append((h,float(fr[i]),float(amp[i]),float(phase[i])))
    return out


def periodic_output(abp, p, initial):
    # raw measured block, absolute values, repeated. Linear interpolation only.
    T=len(abp)/FS
    def pa(t):
        pos=(float(t)%T)*FS; i=int(math.floor(pos)); q=pos-i; j=(i+1)%len(abp)
        return float((1-q)*abp[i]+q*abp[j])
    pp=copy.copy(p)
    state=np.asarray(initial,float).copy(); prev=None; stable=0; last=None
    def rhs(t,y):
        pp.Pa=pa(t)
        return model.evaluate(y,pp)[0]
    for k in range(1,121):
        a,b=(k-1)*T,k*T
        te=np.arange(a,b+1e-12,.01)
        if te[-1] < b: te=np.append(te,b)
        sol=solve_ivp(rhs,(a,b),state,method="BDF",rtol=1e-8,atol=1e-10,max_step=.01,t_eval=te)
        if not sol.success: raise RuntimeError(sol.message)
        state=sol.y[:,-1]; w=sol.y[0,:-1]
        if prev is not None:
            n=min(len(w),len(prev)); d=np.max(np.abs(w[:n]-prev[:n])); md=abs(np.mean(w[:n])-np.mean(prev[:n]))
            stable=stable+1 if d<5e-4 and md<1e-4 else 0
        prev=w.copy(); last=sol
        if stable>=3: break
    if stable<3: raise RuntimeError("periodic convergence failed")
    # resample final model block to 50 Hz patient sample times
    tt=last.t-last.t[0]
    pic=np.interp(np.arange(len(abp))/FS,tt,last.y[0])
    return pic,state,k


def metrics(abp, measured, predicted):
    ha=harmonics(abp); hm=harmonics(measured); hp=harmonics(predicted)
    rows=[]
    for a,m,p in zip(ha,hm,hp):
        h,f,A,pha=a; _,_,M,phm=m; _,_,P,php=p
        rows.append(dict(harmonic=h,freq_hz=f,abp_amp=A,measured_icp_amp=M,predicted_icp_amp=P,
                         measured_gain=M/A if A else np.nan,predicted_gain=P/A if A else np.nan,
                         pred_over_measured=P/M if M else np.nan))
    mc=measured-np.mean(measured); pc=predicted-np.mean(predicted)
    return dict(mean_icp=float(np.mean(predicted)),peak_to_peak=float(np.ptp(predicted)),
                centered_corr=float(np.corrcoef(mc,pc)[0,1]),
                centered_rmse=float(np.sqrt(np.mean((mc-pc)**2)))),rows


def main():
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument("--charis-dat",type=Path,required=True); ap.add_argument("--output-dir",type=Path,required=True)
    args=ap.parse_args(); out=args.output_dir; out.mkdir(parents=True,exist_ok=True)

    ref=json.loads(Path("reference_results/reports/tbi_baseline.json").read_text())
    y0=np.asarray(ref["terminal_window_mean_state"],float)
    p=make_tbi_parameters()

    # Small-signal linearization at maintained TBI equilibrium (Pa=100).
    A,B=linearize(y0,p)
    eig=np.linalg.eigvals(A)
    freqs=np.unique(np.concatenate([np.logspace(-3,1.3,120),np.array([1.6666667,3.3333333,5,6.6666667,8.3333333,10])]))
    fr=frequency_response(A,B,freqs)
    with (out/"linear_frequency_response.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fr[0].keys()); w.writeheader(); w.writerows(fr)

    # Estimate direct arterial compartment time constants at equilibrium.
    aux=model.flow_snapshot(y0,p)
    Rup=p.Rla + float(aux["Rpa"])/2.0
    Rout=float(aux["Rpa"])/2.0
    Req=1.0/(1.0/Rup+1.0/Rout)
    tau_pial=float(aux["Cpa"])*Req
    Cic=float(aux["Cic"])
    Cvi=float(aux["Cvi"])

    abp,measured=load_charis1(args.charis_dat)

    variants=[]
    def add(name, mut):
        pv=make_tbi_parameters(); mut(pv); variants.append((name,pv))
    add("base_TBI",lambda q:None)
    add("kE_x0.5",lambda q:setattr(q,"kE",q.kE*.5))
    add("kE_x1.5",lambda q:setattr(q,"kE",q.kE*1.5))
    add("kE_x2",lambda q:setattr(q,"kE",q.kE*2))
    add("kE_x3",lambda q:setattr(q,"kE",q.kE*3))
    add("Rla_x0.5",lambda q:setattr(q,"Rla",q.Rla*.5))
    add("Rla_x2",lambda q:setattr(q,"Rla",q.Rla*2))
    add("arterial_compliance_x0.5",lambda q:(setattr(q,"Cpan",q.Cpan*.5),setattr(q,"Cpa1",q.Cpa1*.5),setattr(q,"Cpa2",q.Cpa2*.5)))
    add("arterial_compliance_x2",lambda q:(setattr(q,"Cpan",q.Cpan*2),setattr(q,"Cpa1",q.Cpa1*2),setattr(q,"Cpa2",q.Cpa2*2)))
    add("CSF_exchange_nearly_off",lambda q:(setattr(q,"Rf",q.Rf*1000),setattr(q,"R0",q.R0*1000)))
    add("autoreg_tau_x100",lambda q:setattr(q,"tau_aut",q.tau_aut*100))
    add("extracranial_Pa_path_off",lambda q:setattr(q,"Gex",0.0))

    summary_rows=[]; harmonic_rows=[]
    for name,pv in variants:
        try:
            pred,state,blocks=periodic_output(abp,pv,y0)
            mm,hh=metrics(abp,measured,pred)
            summary_rows.append(dict(variant=name,converged=True,blocks=blocks,**mm))
            for r in hh: harmonic_rows.append(dict(variant=name,**r))
        except Exception as e:
            summary_rows.append(dict(variant=name,converged=False,blocks="",mean_icp="",peak_to_peak="",centered_corr="",centered_rmse="",error=str(e)))

    with (out/"parameter_diagnostics.csv").open("w",newline="") as f:
        fields=sorted({k for r in summary_rows for k in r})
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(summary_rows)
    with (out/"harmonic_diagnostics.csv").open("w",newline="") as f:
        fields=list(harmonic_rows[0].keys()); w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(harmonic_rows)

    report={
      "equilibrium":{"Pic_mmhg":float(y0[0]),"Rpa":float(aux["Rpa"]),"Cpa_ml_per_mmhg":float(aux["Cpa"]),
                     "Cic_ml_per_mmhg":Cic,"Cvi_ml_per_mmhg":Cvi,
                     "arterial_Rup_mmhg_s_per_ml":Rup,"arterial_Rout_mmhg_s_per_ml":Rout,
                     "arterial_parallel_Req":Req,"pial_windkessel_tau_s":tau_pial},
      "linear_system":{"eigenvalues_real_imag":[[float(x.real),float(x.imag)] for x in eig]},
      "charis1":{"measured_mean_icp":float(np.mean(measured)),"measured_peak_to_peak":float(np.ptp(measured)),
                 "measured_abp_mean":float(np.mean(abp)),"measured_abp_peak_to_peak":float(np.ptp(abp))},
      "parameter_tests":summary_rows
    }
    (out/"diagnostic_report.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=="__main__":
    main()
