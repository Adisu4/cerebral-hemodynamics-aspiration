"""
Numerical solver for intracranial model

Integrates the model ODE system with adaptive time stepping and supports
stiff and non-stiff solvers.
"""

import numpy as np
from scipy.integrate import solve_ivp
from gadda_model.parameters import ModelParameters
from gadda_model.equations import gadda_ode_system_FIXED as gadda_ode_system

def get_flow_rate(protocol, t):
    if protocol is None or 'flow_schedule' not in protocol:
        return 0.0
    for segment in protocol['flow_schedule']:
        if segment['start'] <= t < segment['end']:
            frac = (t - segment['start']) / (segment['end'] - segment['start'])
            flow = segment['flow_start'] + frac * (segment['flow_end'] - segment['flow_start'])
            return flow
    return 0.0

class SimulationResult:
    def __init__(self, sol, params: ModelParameters):
        self.t = sol.t
        self.y = sol.y
        self.params = params
        self.success = sol.success
        self.message = sol.message
        
        # 15 State Unpacking
        self.Pic = self.y[0, :]
        self.Ppa = self.y[1, :]
        self.Pv = self.y[2, :]
        self.Pvs = self.y[3, :]
        self.Pjr3 = self.y[4, :]
        self.Pjl3 = self.y[5, :]
        self.Pjr2 = self.y[6, :]
        self.Pjl2 = self.y[7, :]
        self.Pc3 = self.y[8, :]
        self.Pc2 = self.y[9, :]
        self.Pvv = self.y[10, :]
        self.Pazy = self.y[11, :]
        self.Psvc = self.y[12, :]
        self.xaut = self.y[13, :]
        self.Cpa = self.y[14, :]
        
        self.ICP = self.Pic
        self.Pa = params.Pa * np.ones_like(self.t)
        self.CPP = self.Pa - self.ICP
    
    def get_final_state(self):
        return self.y[:, -1]
    
    def get_summary(self, window_s=60.0):
        if len(self.t) == 0: return {'ICP_mean': 0.0, 'success': False}
        idx = self.t >= (self.t[-1] - window_s)
        if not np.any(idx): idx = slice(None)
        
        return {
            'ICP_mean': float(np.mean(self.ICP[idx])),
            'ICP_std': float(np.std(self.ICP[idx])),
            'Pvs_mean': float(np.mean(self.Pvs[idx])),
            'CPP_mean': float(np.mean(self.CPP[idx])),
            'success': bool(self.success)
        }

def run_simulation(params, duration_s=900.0, t_start=0.0, y0=None, max_step=0.5, flow_rate_callback=None, method='Radau', rtol=None, atol=None):
    if y0 is None:
        y0 = params.get_initial_conditions()
    
    if rtol is None: rtol = 1e-6
    if atol is None: atol = 1e-8
    
    # Handle flow callback with absolute time
    def ode_func(t, y):

        return gadda_ode_system(t, y, params, flow_rate_callback=flow_rate_callback)
    
    try:
        # Time span must reflect absolute time for protocol timing to function correctly
        t_span = (t_start, t_start + duration_s)
        
        sol = solve_ivp(
            ode_func,
            t_span=t_span,
            y0=y0,
            method=method,
            dense_output=False,
            max_step=max_step, 
            rtol=rtol,
            atol=atol
        )
        return SimulationResult(sol, params)
    except Exception as e:
        print(f"    Solver {method} error: {e}")
        class FailedSol:
            t = np.array([t_start])
            y = np.zeros((15, 1))
            y[:,0] = y0
            success = False
            message = str(e)
        return SimulationResult(FailedSol(), params)

def run_full_protocol(params, stabilization_time_s=900.0, intervention_time_s=1500.0, max_stabilization_attempts=1):

    aspiration_backup = getattr(params, 'aspiration_protocol', None)
    params.aspiration_protocol = None
    

    baseline = run_simulation(params, duration_s=stabilization_time_s, t_start=0.0, method='RK45')
    if not baseline.success:
        print("  RK45 stabilization failed, retrying with LSODA...")
        baseline = run_simulation(params, duration_s=stabilization_time_s, t_start=0.0, method='LSODA')
    
    # Restore protocol
    params.aspiration_protocol = aspiration_backup
    
    # Define callback using the protocol object
    if aspiration_backup:
        # The callback receives ABSOLUTE time
        flow_cb = lambda t, pic: get_flow_rate(aspiration_backup, t)
    else:
        flow_cb = None
        
    y0 = baseline.get_final_state()
    t_intervention_start = baseline.t[-1]
    

    
    # Shift time for protocol timing
    def shifted_callback(t_abs, pic):
        t_rel = t_abs - t_intervention_start
        return flow_cb(t_rel, pic)

    intervention = run_simulation(
        params, 
        duration_s=intervention_time_s, 
        t_start=t_intervention_start,
        y0=y0, 
        flow_rate_callback=shifted_callback, 
        method='Radau' # Radau is safer for stiff intervention changes
    )
    
    return baseline, intervention

__all__ = ['SimulationResult', 'run_simulation', 'run_full_protocol']