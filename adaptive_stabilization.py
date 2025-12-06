"""
Adaptive Stabilization Algorithm

Replaces fixed stabilization times with convergence detection.
Implements chunked simulation with ICP drift monitoring for baseline equilibration.
"""

import numpy as np
from gadda_model.solver import run_simulation

def run_adaptive_stabilization(params,
                               max_time_s=3600,
                               check_interval_s=300,
                               window_s=120,
                               icp_tolerance_mmhg=0.1,
                               drift_tolerance_mmhg_per_min=0.02,
                               min_checks=2,
                               verbose=True):
    
    if params.R0 > 1000:
        solver_method = 'Radau'
        check_interval_s = min(check_interval_s * 3, 900)
    else:
        solver_method = 'RK45'

    y0 = params.get_initial_conditions()
    t_elapsed = 0.0
    consecutive_stable = 0

    while t_elapsed < max_time_s:
        chunk_duration = min(check_interval_s, max_time_s - t_elapsed)
        

        result = run_simulation(
            params,
            duration_s=chunk_duration,
            t_start=t_elapsed,
            y0=y0,
            method=solver_method,
            rtol=1e-4 if params.R0 > 1000 else 1e-6,
            atol=1e-6 if params.R0 > 1000 else 1e-8,
            max_step=2.0 if params.R0 > 1000 else 0.5
        )

        if not result.success:
            if verbose: print(f"Solver terminated at t={t_elapsed:.0f}s")
            return result

        t_elapsed = result.t[-1]
        y0 = result.get_final_state()


        idx = result.t >= (result.t[-1] - window_s)
        icp_window = result.ICP[idx]
        t_window = result.t[idx]

        if len(icp_window) < 5: continue

        icp_mean = np.mean(icp_window)
        icp_std = np.std(icp_window)
        
        if len(t_window) > 2:
            coeffs = np.polyfit(t_window - t_window[0], icp_window, 1)
            drift_rate_per_min = coeffs[0] * 60.0
        else:
            drift_rate_per_min = 0.0

        is_stable = (icp_std < icp_tolerance_mmhg and
                    abs(drift_rate_per_min) < drift_tolerance_mmhg_per_min)

        if is_stable:
            consecutive_stable += 1
        else:
            consecutive_stable = 0

        if verbose:
            print(f"t={t_elapsed:>6.0f}s: ICP={icp_mean:>6.2f}±{icp_std:>5.3f}mmHg, drift={drift_rate_per_min:>+7.4f}/min")

        if consecutive_stable >= min_checks:
            if verbose: print(f"Convergence achieved at t={t_elapsed:.0f}s")
            return result

    return result

def diagnose_equilibrium_time(test_case_id):
    """
    Diagnose equilibrium time for a test case.
    
    Args:
        test_case_id: ID of the test case to diagnose
    """
    from aspiration_study.test_cases import create_test_parameters, CLINICAL_TEST_CASES

    test = CLINICAL_TEST_CASES[test_case_id]
    params = create_test_parameters(test_case_id)

    print(f"\n{'='*80}")
    print(f"Equilibrium analysis: {test_case_id}")
    print(f"{'='*80}")
    print(f"Pathology: {test['name']}")
    print(f"Description: {test['description']}")
    print(f"Expected ICP: {test['expected_icp']}")
    print(f"\nParameters:")
    print(f"  R0={params.R0:.0f} mmHg·s/mL, Gaut={params.Gaut:.2f}, kE={params.kE:.3f}")

    estimated_Cic = 1.0 / (params.kE * 10.0)
    tau_csf = params.R0 * estimated_Cic
    tau_autoreg = params.tau_aut

    print(f"\nTime Constants:")
    print(f"  τ_autoreg={tau_autoreg:.0f}s, τ_CSF={tau_csf:.0f}s, 5×τ_CSF={5*tau_csf:.0f}s")
    print(f"  Recommended: {max(300, 5*tau_csf):.0f}s ({max(300, 5*tau_csf)/60:.1f}min)")
    print(f"{'='*80}\n")

    max_time = 14400 if params.R0 > 1000 else 7200
    result = run_adaptive_stabilization(
        params,
        max_time_s=max_time,
        check_interval_s=300,
        verbose=True
    )

    return result


def run_two_phase_aspiration_test(test_case_id, site='Pvs', flow_rate_mL_min=120.0, aspiration_duration_s=600.0):
    """
    Run two-phase aspiration test with forced timestep progression:
    1. Baseline stabilization with adaptive convergence detection
    2. Aspiration from stabilized baseline using BDF solver with t_eval
    
    Args:
        test_case_id: Test case ID (T0, T1, T2, T3)
        site: Aspiration site (Pvs, Pv, J3, J2, J1)
        flow_rate_mL_min: Target aspiration flow rate (mL/min), default 120
        aspiration_duration_s: Aspiration duration (seconds), default 600 (10 min)
    """
    from aspiration_study.test_cases import create_test_parameters, CLINICAL_TEST_CASES
    from gadda_model.aspiration import AspirationProtocol
    from scipy.integrate import solve_ivp
    import numpy as np
    
    print(f"\n{'='*80}")
    print(f"Two-phase aspiration analysis: {test_case_id}")
    print(f"{'='*80}")
    print(f"Pathology: {CLINICAL_TEST_CASES[test_case_id]['name']}")
    print(f"Site: {site}, Flow: {flow_rate_mL_min}mL/min, Duration: {aspiration_duration_s/60:.1f}min")
    
    params = create_test_parameters(test_case_id)
    
    print(f"\nPHASE 1: Baseline Stabilization (Adaptive Convergence)")
    
    max_time = 14400 if params.R0 > 1000 else 7200
    baseline = run_adaptive_stabilization(
        params,
        max_time_s=max_time,
        check_interval_s=300,
        verbose=True
    )
    
    if not baseline.success:
        print("Baseline simulation convergence failure")
        return None
    
    baseline_icp = baseline.get_summary()['ICP_mean']
    baseline_cpp = baseline.get_summary()['CPP_mean']
    baseline_pvs = baseline.get_summary()['Pvs_mean']
    
    print(f"Baseline stabilized: ICP={baseline_icp:.2f}mmHg, CPP={baseline_cpp:.2f}mmHg")
    print(f"  Baseline Pvs: {baseline_pvs:.2f} mmHg")
    print(f"  Pressure gradient (Pvs-ICP): {baseline_pvs - baseline_icp:.2f} mmHg")
    
    print(f"\nPHASE 2: Aspiration (Forced Timestep Progression)")
    
    ramp_duration_s = 60.0
    target_flow_mL_s = flow_rate_mL_min / 60.0
    
    # Smooth S-curve ramp
    def flow_callback(t, pic):
        if t < 0:
            return 0.0
        elif t < ramp_duration_s:
            frac = t / ramp_duration_s
            smooth_frac = 3*frac**2 - 2*frac**3  # S-curve
            return smooth_frac * target_flow_mL_s
        elif t < aspiration_duration_s:
            return target_flow_mL_s
        else:
            return 0.0
    
    protocol = AspirationProtocol.create_constant_flow(
        site=site,
        target_flow_mL_min=flow_rate_mL_min,
        ramp_duration_s=ramp_duration_s,
        total_duration_s=aspiration_duration_s
    )
    params.aspiration_protocol = protocol
    
    y0 = baseline.get_final_state()
    
    # Explicit time evaluation points for rapid ICP dynamics during aspiration
    num_eval_points = max(100, int(aspiration_duration_s / 5))
    t_eval = np.linspace(0, aspiration_duration_s, num_eval_points)
    
    print(f"\nSolver Configuration:")
    print(f"  Method: BDF (most stable for stiff systems)")
    print(f"  Forced evaluation points: {len(t_eval)}")
    print(f"  Max step: 5.0s")
    print(f"  Tolerances: rtol=1e-3, atol=1e-5 (relaxed for stability)")
    

    from gadda_model.equations import gadda_ode_system_FIXED
    
    def ode_with_callback(t, y):
        return gadda_ode_system_FIXED(t, y, params, flow_rate_callback=flow_callback)
    

    print(f"\nStarting intervention simulation...")
    print(f"Target: {flow_rate_mL_min}mL/min ({target_flow_mL_s:.2f}mL/s), Duration: {aspiration_duration_s/60:.0f}min")
    
    sol = solve_ivp(
        ode_with_callback,
        t_span=(0, aspiration_duration_s),
        y0=y0,
        method='BDF',
        t_eval=t_eval,
        max_step=5.0,
        rtol=1e-3,              # Relaxed tolerance
        atol=1e-5,              # Relaxed tolerance
        dense_output=False
    )
    
    if not sol.success:
        print(f"BDF solver failed: {sol.message}")
        print(f"  Attempting with Radau solver with tighter first_step...")
        
        sol = solve_ivp(
            ode_with_callback,
            t_span=(0, aspiration_duration_s),
            y0=y0,
            method='Radau',
            t_eval=t_eval,
            first_step=0.01,
            max_step=2.0,
            rtol=1e-4,
            atol=1e-6
        )
        
        if not sol.success:
            print(f"Solver convergence failure: Both BDF and Radau methods unsuccessful")
            return None
    
    print(f"Simulation completed successfully")
    print(f"  Timesteps: {len(sol.t)}")
    print(f"  Time range: {sol.t[0]:.1f}s to {sol.t[-1]:.1f}s")
    

    from gadda_model.solver import SimulationResult
    
    class MinimalSol:
        def __init__(self, t, y, success, message):
            self.t = t
            self.y = y
            self.success = success
            self.message = message
    
    minimal_sol = MinimalSol(sol.t, sol.y, sol.success, sol.message)
    intervention = SimulationResult(minimal_sol, params)
    

    steady_state_window_s = min(120.0, aspiration_duration_s * 0.5)
    intervention_summary = intervention.get_summary(window_s=steady_state_window_s)
    
    intervention_icp = intervention_summary['ICP_mean']
    intervention_cpp = intervention_summary['CPP_mean']
    
    reduction = baseline_icp - intervention_icp
    reduction_percent = (reduction / baseline_icp) * 100
    
    print(f"\nRESULTS:")
    print(f"  Baseline ICP:  {baseline_icp:.2f} mmHg")
    print(f"  Final ICP:     {intervention_icp:.2f} mmHg")
    print(f"  Reduction:     {reduction:.2f} mmHg ({reduction_percent:.1f}%)")
    print(f"  CPP:           {intervention_cpp:.2f} mmHg")
    
    therapeutic = reduction >= 5.0
    safe = intervention_cpp >= 50.0
    
    if therapeutic and safe:
        print(f"\nResult: Therapeutic reduction achieved with safe CPP")
    elif therapeutic:
        print(f"\nWarning: Therapeutic reduction achieved but CPP below safety threshold")
    elif reduction > 1.0:
        print(f"\nPartial response: {reduction:.2f} mmHg reduction (therapeutic threshold ≥5 mmHg)")
    else:
        print(f"\nInsufficient response: {reduction:.2f} mmHg reduction")
    
    return {
        'test_case_id': test_case_id,
        'site': site,
        'flow_rate_mL_min': flow_rate_mL_min,
        'baseline_icp': baseline_icp,
        'baseline_cpp': baseline_cpp,
        'intervention_icp': intervention_icp,
        'intervention_cpp': intervention_cpp,
        'reduction': reduction,
        'reduction_percent': reduction_percent,
        'therapeutic': therapeutic,
        'safe': safe,
        'baseline': baseline,
        'intervention': intervention
    }


if __name__ == '__main__':
    import sys
    
    print("\n" + "="*80)
    print("Aspiration simulation: Forced timestep with BDF solver")
    print("="*80)
    

    result = run_two_phase_aspiration_test(
        test_case_id='T1',
        site='Pv',
        flow_rate_mL_min=240.0,
        aspiration_duration_s=300.0  # 5 minutes
    )
    
    if result and result['reduction'] > 1.0:
        print(f"\n{'='*80}")
        print("Validation: aspiration protocol reduced ICP")
        print(f"{'='*80}")
        print(f"ICP reduction: {result['reduction']:.2f} mmHg")
        

        print(f"\n\nEvaluating elevated flow rates for T2 pathology...")
        for flow in [300, 360, 480]:
            print(f"\n{'='*80}")
            print(f"Flow rate: {flow} mL/min")
            print(f"{'='*80}")
            result2 = run_two_phase_aspiration_test('T2', 'Pv', flow, 300.0)
    else:
        print(f"\n{'='*80}")
        print("Response below threshold: consider parameter adjustments")
        print(f"{'='*80}")
        if result:
            print(f"Reduction: {result['reduction']:.2f} mmHg (confirmation threshold: >1 mmHg)")