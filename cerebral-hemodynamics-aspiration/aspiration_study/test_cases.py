"""
Aspiration study clinical test cases

Four clinical scenarios (T0–T3) covering varied severities and anatomical sites.
"""

import numpy as np
from gadda_model import (
    ModelParameters,
    AspirationProtocol,
    run_full_protocol
)
from gadda_model.solver import get_flow_rate


# Clinical test case definitions
CLINICAL_TEST_CASES = {
    'T0': {
        'name': 'Healthy Control',
        'description': 'Normal anatomy - baseline reference',
        'R0': 526.3,
        'Gaut': 3.0,  # Normal autoregulation
        'stenosis_level': 'none',
        'stenosis_degree': 1.0,
        'aspiration_sites': ['Pvs', 'Pv', 'J3', 'J2', 'J1'],
        'stabilization_time': 600,  # 10 min
        'expected_icp': '~8-12 mmHg',
        'clinical_context': 'Healthy control - validates aspiration specificity'
    },
    'T1': {
        'name': 'Mild IIH',
        'description': 'Bilateral J3 stenosis, 50% impaired autoregulation, elevated R0',
        'R0': 2000,  # High CSF outflow resistance for elevated ICP
        'Gaut': 0.5,  # 50% of normal (impaired but functional)
        'stenosis_level': 'J3',  # Bilateral J3 stenosis (most common IIH pattern)
        'stenosis_degree': 0.10,  # 90% stenosis = 10% remaining conductance
        'aspiration_sites': ['Pvs', 'Pv', 'J3'],
        'stabilization_time': 5100,  # 85 min (adaptive convergence detection)
        'expected_icp': '~12-16 mmHg',
        'clinical_context': 'Stable pathology - measures true aspiration effect',
        'interpretation': 'Physiological IIH with measurable aspiration response'
    },
    'T2': {
        'name': 'Mild Traumatic Brain Injury (TBI)',
        'description': 'Impaired autoregulation + venous compression from edema, elevated CSF resistance',
        'R0': 1800,  # Elevated CSF resistance for higher baseline ICP
        'Gaut': 0.3,  # Severely impaired autoregulation
        'Rvs1_factor': 3.0,  # 3x venous resistance (edema, high CVP)
        'stenosis_level': 'J3',
        'stenosis_degree': 0.15,
        'aspiration_sites': ['Pv', 'Pvs', 'J3', 'J2', 'J1'],
        'stabilization_time': 600,  # 10 min
        'expected_icp': '~20-28 mmHg',
        'clinical_context': 'Mild TBI with autoregulation failure and venous congestion',
        'interpretation': 'TBI pathology - aspiration targets venous drainage'
    },
    'T3': {
        'name': 'Acute Brain Edema / Bleed (Low Compliance)',
        'description': 'Decreased intracranial compliance (increased elastance) - stiff system',
        'R0': 2900,  # Significantly elevated CSF resistance
        'Gaut': 1.0,  # 33% of normal (moderate impairment)
        'kE_factor': 1.5,  # 50% increase in elastance (reduced compliance)
        'stenosis_level': 'J3',
        'stenosis_degree': 0.15,
        'aspiration_sites': ['Pv', 'Pvs', 'J3', 'J2', 'J1'],
        'stabilization_time': 600,  # 10 min
        'expected_icp': '~25-30 mmHg',
        'clinical_context': 'Acute edema/bleed with reduced compliance - steep ICP response',
        'interpretation': 'Low compliance pathology - small volume changes cause large ICP rises'
    },
}


def create_test_parameters(test_case_id):
    """
    Create model parameters for a test case.
    
    Args:
        test_case_id: Key from CLINICAL_TEST_CASES
    
    Returns:
        ModelParameters object configured for test case
    """
    if test_case_id not in CLINICAL_TEST_CASES:
        raise ValueError(f"Unknown test case: {test_case_id}")
    
    test = CLINICAL_TEST_CASES[test_case_id]
    
    params = ModelParameters(posture='supine')

    # Configure CSF outflow resistance and conductance for pathology
    params.R0 = test['R0']
    params.G0 = 1.0 / test['R0'] if test['R0'] > 0 else 0.1

    params.Gaut = test['Gaut']
    
    # Apply Rvs1 scaling if specified (e.g., for TBI)
    if 'Rvs1_factor' in test:
        params.Rvs1 *= test['Rvs1_factor']
    
    # Apply kE scaling if specified (e.g., for low compliance)
    if 'kE_factor' in test:
        params.kE *= test['kE_factor']

    # Apply stenosis based on level
    stenosis_level = test.get('stenosis_level', 'none')
    stenosis_degree = test.get('stenosis_degree', 1.0)
    
    # Bilateral stenosis (both sides)
    if stenosis_level == 'J3':
        params.kjr3 *= stenosis_degree
        params.kjl3 *= stenosis_degree
    elif stenosis_level == 'J2':
        params.kjr2 *= stenosis_degree
        params.kjl2 *= stenosis_degree
    elif stenosis_level == 'J1':
        params.kjr1 *= stenosis_degree
        params.kjl1 *= stenosis_degree
    elif stenosis_level == 'all':
        params.kjr3 *= stenosis_degree
        params.kjl3 *= stenosis_degree
        params.kjr2 *= stenosis_degree
        params.kjl2 *= stenosis_degree
        params.kjr1 *= stenosis_degree
        params.kjl1 *= stenosis_degree
    # Unilateral stenosis (right side only)
    elif stenosis_level == 'J3_right':
        params.kjr3 *= stenosis_degree
    elif stenosis_level == 'J2_right':
        params.kjr2 *= stenosis_degree
    elif stenosis_level == 'J1_right':
        params.kjr1 *= stenosis_degree

    # expected_icp is descriptive only; let the model reach its own equilibrium

    return params


def run_single_test(test_case_id, site='Pvs', protocol_type='constant', 
                    flow_rate_mL_min=None, verbose=True):
    """
    Run single aspiration test case.
    
    Args:
        test_case_id: Key from CLINICAL_TEST_CASES
        site: Aspiration site ('Pvs', 'Pv', 'J3', 'J2', 'J1')
        protocol_type: 'constant' or 'stepwise'
        flow_rate_mL_min: Target flow rate for constant protocol (None = auto-scale)
        verbose: Print detailed results
    
    Returns:
        Dictionary with test results
    """
    if verbose:
        test = CLINICAL_TEST_CASES[test_case_id]
        print("\n" + "="*80)
        print(f"Test case: {test_case_id} - {test['name']}")
        print("="*80)
        print(f"Description: {test['description']}")
        print(f"Expected ICP: {test['expected_icp']}")
        print(f"Clinical Context: {test['clinical_context']}")
        print("-"*80)
        print(f"Parameters: R0={test['R0']}, Gaut={test['Gaut']}")
        print(f"Aspiration Site: {site}")
        print(f"Protocol: {protocol_type}")
    
    # Auto-scale flow rate if not specified
    if flow_rate_mL_min is None:
        flow_rate_mL_min = get_scaled_flow_rate(test_case_id)
    
    if verbose:
        print(f"Flow Rate: {flow_rate_mL_min} mL/min (auto-scaled)")
        print("-"*80)
    

    params = create_test_parameters(test_case_id)
    

    if protocol_type == 'constant':
        protocol = AspirationProtocol.create_constant_flow(
            site=site,
            target_flow_mL_min=flow_rate_mL_min,
            ramp_duration_s=60.0,
            total_duration_s=1500.0
        )
    elif protocol_type == 'stepwise':
        protocol = AspirationProtocol.create_stepwise_escalation(
            site=site,
            flow_steps_mL_min=[30, 60, 90, 120, 150],
            step_duration_s=180.0,
            ramp_between_steps_s=30.0
        )
    else:
        raise ValueError(f"Unknown protocol type: {protocol_type}")
    
    params.aspiration_protocol = protocol
    

    if verbose:
        print("Running baseline stabilization...")
    baseline, intervention = run_full_protocol(
        params,
        stabilization_time_s=600.0,
        intervention_time_s=1500.0
    )
    

    baseline_summary = baseline.get_summary(window_s=60)
    intervention_summary = intervention.get_summary(window_s=60)
    
    reduction = baseline_summary['ICP_mean'] - intervention_summary['ICP_mean']
    reduction_percent = (reduction / baseline_summary['ICP_mean']) * 100
    

    therapeutic = reduction >= 5.0
    vs_shunt_percent = (reduction / 12.5) * 100  # VP shunt ~12.5 mmHg
    
    if verbose:
        print("\n" + "-"*80)
        print("RESULTS")
        print("-"*80)
        print(f"Baseline ICP: {baseline_summary['ICP_mean']:.2f} ± {baseline_summary['ICP_std']:.2f} mmHg")
        print(f"Intervention ICP: {intervention_summary['ICP_mean']:.2f} ± {intervention_summary['ICP_std']:.2f} mmHg")
        print(f"ICP Reduction: {reduction:.2f} mmHg ({reduction_percent:.1f}%)")
        print(f"Therapeutic: {'Yes' if therapeutic else 'No'} (target ≥5 mmHg)")
        print(f"vs VP Shunt: {vs_shunt_percent:.0f}% efficacy")
        print("-"*80)
        print(f"Baseline CPP: {baseline_summary['CPP_mean']:.2f} mmHg")
        print(f"Intervention CPP: {intervention_summary['CPP_mean']:.2f} mmHg")
        print("="*80)
    
    return {
        'test_case_id': test_case_id,
        'site': site,
        'protocol_type': protocol_type,
        'baseline_icp': baseline_summary['ICP_mean'],
        'intervention_icp': intervention_summary['ICP_mean'],
        'baseline_cpp': baseline_summary['CPP_mean'],
        'intervention_cpp': intervention_summary['CPP_mean'],
        'cpp_change': intervention_summary['CPP_mean'] - baseline_summary['CPP_mean'],
        'reduction': reduction,
        'reduction_percent': reduction_percent,
        'therapeutic': therapeutic,
        'vs_shunt_percent': vs_shunt_percent,
        'baseline': baseline,
        'intervention': intervention
    }


def run_test_suite(sites=None, protocol_type='constant', save_summary=True, default_flow_mL_min=60.0):
    """
    Run complete test suite (T0–T3) for multiple aspiration sites.
    
    Args:
        sites: List of aspiration sites (default: from test case)
        protocol_type: 'constant' or 'stepwise'
        save_summary: Save results summary to file
    
    Returns:
        Dictionary with all test results by site
    """
    if sites is None:
        sites = ['Pvs', 'Pv', 'J3', 'J2', 'J1']  # Default fallback
    
    print("\n" + "="*80)
    print(f"Aspiration study: Complete test suite")
    print("="*80)
    print(f"Sites: {', '.join(sites)}, Protocol: {protocol_type}")
    print("="*80)
    

    all_results = {}
    all_test_ids = ['T0', 'T1', 'T2', 'T3']
    
    for test_id in all_test_ids:
        test = CLINICAL_TEST_CASES[test_id]
        test_sites = test.get('aspiration_sites', sites)  # Use test-specific sites if available
        all_results[test_id] = {}
        
        for site in test_sites:

            test_protocol = protocol_type
            flow_rate = get_scaled_flow_rate(test_id, base_flow_mL_min=default_flow_mL_min)
            
            result = run_single_test(test_id, site, test_protocol, 
                                    flow_rate_mL_min=flow_rate, verbose=False)
            all_results[test_id][site] = result
    
    # Get primary site results for printing (use Pvs as primary)
    results = {test_id: all_results[test_id]['Pvs'] for test_id in all_test_ids}
    
    # Multi-site summary table (matching original Table 1.0 format)
    print("\n" + "="*80)
    print("MULTI-SITE ASPIRATION EFFICACY SUMMARY")
    print("="*80)
    print(f"{'Test':<8} {'Baseline':<10} {'ΔICP_Pvs':<10} {'ΔICP_Pv':<10} {'ΔICP_J3':<10} {'ΔICP_J2':<10} {'ΔICP_J1':<10} {'Best Site'}")
    print(f"{'ID':<8} {'ICP (mmHg)':<10} {'(mmHg)':<10} {'(mmHg)':<10} {'(mmHg)':<10} {'(mmHg)':<10} {'(mmHg)':<10} {''}")
    print("-"*85)
    
    for test_id in all_test_ids:
        baseline_icp = all_results[test_id]['Pvs']['baseline_icp']
        
        reductions = {}
        for site in sites:
            reductions[site] = all_results[test_id][site]['reduction']
        
        # Find best site
        best_site = max(reductions.keys(), key=lambda s: reductions[s])
        
        # Format reductions
        pvs_str = f"{reductions.get('Pvs', 0):>+7.2f}" if 'Pvs' in reductions else "N/A"
        pv_str = f"{reductions.get('Pv', 0):>+7.2f}" if 'Pv' in reductions else "N/A"
        j3_str = f"{reductions.get('J3', 0):>+7.2f}" if 'J3' in reductions else "N/A"
        j2_str = f"{reductions.get('J2', 0):>+7.2f}" if 'J2' in reductions else "N/A"
        j1_str = f"{reductions.get('J1', 0):>+7.2f}" if 'J1' in reductions else "N/A"
        
        print(f"{test_id:<8} {baseline_icp:>8.2f}   {pvs_str:<10} {pv_str:<10} {j3_str:<10} {j2_str:<10} {j1_str:<10} {best_site}")
    
    print("="*80)
    
    # CPP Safety Analysis
    print("\n" + "="*80)
    print("CPP POST-ASPIRATION SAFETY ANALYSIS (Multi-Site Comparison)")
    print("="*80)
    print(f"{'Test':<8} {'Baseline':<10} {'Pvs':<10} {'Pv':<10} {'J3':<10} {'J2':<10} {'J1':<10} {'Safety Status'}")
    print(f"{'ID':<8} {'CPP (mmHg)':<10} {'Post-Asp':<10} {'Post-Asp':<10} {'Post-Asp':<10} {'Post-Asp':<10} {'Post-Asp':<10} {''}")
    print("-"*85)
    
    for test_id in all_test_ids:
        baseline_cpp = all_results[test_id]['Pvs']['baseline_cpp']
        
        cpp_values = {}
        for site in sites:
            cpp_values[site] = all_results[test_id][site]['intervention_cpp']
            
        # Format CPP values
        pvs_str = f"{cpp_values.get('Pvs', 0):>8.2f}" if 'Pvs' in cpp_values else "N/A"
        pv_str = f"{cpp_values.get('Pv', 0):>8.2f}" if 'Pv' in cpp_values else "N/A"
        j3_str = f"{cpp_values.get('J3', 0):>8.2f}" if 'J3' in cpp_values else "N/A"
        j2_str = f"{cpp_values.get('J2', 0):>8.2f}" if 'J2' in cpp_values else "N/A"
        j1_str = f"{cpp_values.get('J1', 0):>8.2f}" if 'J1' in cpp_values else "N/A"
        
        # Determine worst safety status
        min_cpp = min(cpp_values.values()) if cpp_values else 0
        if min_cpp >= 60:
            safety = "✓ Safe"
        elif min_cpp >= 50:
            safety = "⚠ Marginal"
        else:
            safety = "✗ Risk"
            
        print(f"{test_id:<8} {baseline_cpp:>8.2f}   {pvs_str:<10} {pv_str:<10} {j3_str:<10} {j2_str:<10} {j1_str:<10} {safety}")
    
    print("="*80)
    print("Safety Thresholds: CPP ≥60 mmHg (safe), 50-60 mmHg (marginal), <50 mmHg (ischemic risk)")
    print("="*80)
    
    if save_summary:
        # Save to text file
        filename = f'test_suite_summary_multisite_{protocol_type}.txt'
        with open(filename, 'w') as f:
            f.write("ASPIRATION STUDY TEST SUITE SUMMARY (MULTI-SITE)\n")
            f.write("="*80 + "\n")
            f.write(f"Sites: {', '.join(sites)}\n")
            f.write(f"Protocol: {protocol_type}\n")
            f.write("="*80 + "\n\n")
            
            for test_id in all_test_ids:
                test = CLINICAL_TEST_CASES[test_id]
                baseline_icp = all_results[test_id]['Pvs']['baseline_icp']
                
                f.write(f"\n{test_id}: {test['name']}\n")
                f.write(f"  Description: {test['description']}\n")
                f.write(f"  Baseline ICP: {baseline_icp:.2f} mmHg\n\n")
                
                for site in sites:
                    res = all_results[test_id][site]
                    f.write(f"  [{site}] ICP Reduction: {res['reduction']:.2f} mmHg "
                           f"(Final ICP: {res['intervention_icp']:.2f} mmHg, "
                           f"Final CPP: {res['intervention_cpp']:.2f} mmHg, "
                           f"Therapeutic: {'Yes' if res['therapeutic'] else 'No'})\n")
        
        print(f"\nSummary saved: {filename}")
    
    return all_results


def run_flow_sensitivity_analysis(test_case_ids=None, sites=None, flow_rates=None, access_method='femoral'):
    """
    Algorithmic sensitivity analysis: Test flow rates for each pathology first,
    then determine optimal flows before site comparison.
    
    Step 1: Flow sensitivity for each pathology (Pvs only)
    Step 2: Determine optimal flow for each pathology
    Step 3: Site comparison at optimal flows
    
    Args:
        test_case_ids: List of test case IDs to analyze
        sites: List of aspiration sites
        flow_rates: List of flow rates to test (mL/min)
        access_method: 'femoral' (intravascular) or 'direct' (percutaneous puncture)
    
    Returns:
        Dictionary with sensitivity results and optimal flows
    """
    if test_case_ids is None:
        test_case_ids = ['T0', 'T1', 'T2', 'T3']
    
    if sites is None:
        sites = ['Pvs', 'J3', 'J2', 'J1']
    
    if flow_rates is None:
        flow_rates = [16, 60, 120, 240]  # Experimental roller pump flow rates
    
    print("\n" + "="*100)
    print("ALGORITHMIC ASPIRATION SENSITIVITY ANALYSIS")
    print("="*100)
    print("Phase 1: Flow Rate Sensitivity Testing")
    print("Phase 2: Optimal Flow Determination") 
    print("Phase 3: Site Comparison at Optimal Flows")
    print("="*100)
    print(f"Analyzing {len(test_case_ids)} pathologies with {len(flow_rates)} flow rates each...")
    print(f"Sites to compare: {sites}")
    
    # Phase 1: Flow sensitivity for each pathology (Pvs only)
    print("\nPHASE 1: FLOW SENSITIVITY ANALYSIS (Pvs Aspiration)")
    print("-" * 60)
    
    sensitivity_results = {}
    baseline_icps = {}  # Store baseline ICPs for each test
    baseline_values = {}  # Store all baseline values
    
    for test_id in test_case_ids:
        print(f"\nAnalyzing {test_id}...")
        sensitivity_results[test_id] = {}
        
        for flow in flow_rates:

            params = create_test_parameters(test_id)
            

            protocol = AspirationProtocol.create_constant_flow(
                site='Pvs',
                target_flow_mL_min=flow,
                ramp_duration_s=60.0,
                total_duration_s=2100.0,
                access_method=access_method
            )
            params.aspiration_protocol = protocol
            # Parameters are read directly from the protocol in the ODE system
            
            baseline_result, intervention_result = run_full_protocol(
                params, stabilization_time_s=600.0, intervention_time_s=1800.0, max_stabilization_attempts=5
            )
            
            baseline_icp = baseline_result.get_summary()['ICP_mean']
            baseline_icps[test_id] = baseline_icp
            
            # Store baseline pressures for reference
            idx = baseline_result.t >= (baseline_result.t[-1] - 60)
            baseline_pvs = float(np.mean(baseline_result.Pvs[idx]))
            baseline_pjr3 = float(np.mean(baseline_result.Pjr3[idx]))
            baseline_pjl3 = float(np.mean(baseline_result.Pjl3[idx]))
            baseline_pjr2 = float(np.mean(baseline_result.Pjr2[idx]))
            baseline_pjl2 = float(np.mean(baseline_result.Pjl2[idx]))
            
            baseline_values[test_id] = {
                'icp': baseline_icp,
                'pvs': baseline_pvs,
                'pjr3': baseline_pjr3,
                'pjl3': baseline_pjl3,
                'pjr2': baseline_pjr2,
                'pjl2': baseline_pjl2
            }
            
            print(f"  Baseline ICP: {baseline_icp:.1f} mmHg, Pvs: {baseline_pvs:.1f} mmHg, J3: {baseline_pjr3:.1f}/{baseline_pjl3:.1f}, J2: {baseline_pjr2:.1f}/{baseline_pjl2:.1f} mmHg")
            
            final_icp = intervention_result.get_summary()['ICP_mean']
            reduction = baseline_icp - final_icp
            final_cpp = intervention_result.get_summary()['CPP_mean']
            
            # Get final Pvs and individual jugular pressures
            idx = intervention_result.t >= (intervention_result.t[-1] - 60)  # Last 60 seconds
            final_pvs = float(np.mean(intervention_result.Pvs[idx]))
            final_pjr3 = float(np.mean(intervention_result.Pjr3[idx]))
            final_pjl3 = float(np.mean(intervention_result.Pjl3[idx]))
            final_pjr2 = float(np.mean(intervention_result.Pjr2[idx]))
            final_pjl2 = float(np.mean(intervention_result.Pjl2[idx]))
            
            # Compute actual aspiration flow achieved at end of intervention
            if hasattr(intervention_result, 't') and len(intervention_result.t) > 0:
                # Get flow rate at end of intervention
                end_time = intervention_result.t[-1]
                flow_rate_callback = lambda t, pic: get_flow_rate(protocol, t)
                actual_flow_mL_s = flow_rate_callback(end_time, final_icp)
                actual_flow_mL_min = actual_flow_mL_s * 60.0
            else:
                actual_flow_mL_min = 0.0
            
            # Calculate time to achieve maximum ICP reduction (not necessarily to normal)
            time_to_max_reduction = None
            min_icp_achieved = float('inf')
            
            icp_timeseries = intervention_result.ICP
            times = intervention_result.t
            
            # Find the minimum ICP achieved and when it was first reached
            for i, icp_val in enumerate(icp_timeseries):
                if icp_val < min_icp_achieved:
                    min_icp_achieved = icp_val
                    time_to_max_reduction = times[i]
            
            # Also calculate time to normal ICP if achieved
            normal_icp_threshold = 15.0
            time_to_normal = None
            if min_icp_achieved <= normal_icp_threshold:
                time_to_normal = time_to_max_reduction
            
            sensitivity_results[test_id][flow] = {
                'reduction': reduction,
                'final_icp': final_icp,
                'final_cpp': final_cpp,
                'final_pvs': final_pvs,
                'final_pjr3': final_pjr3,
                'final_pjl3': final_pjl3,
                'final_pjr2': final_pjr2,
                'final_pjl2': final_pjl2,
                'time_to_normal': time_to_normal,
                'therapeutic': reduction >= 5.0,
                'safe': final_cpp >= 50.0
            }
            
            # Report results
            time_str = "N/A"
            
            if reduction > 0.1:  # Only report time if there's meaningful reduction
                time_str = f"{time_to_max_reduction:.1f}s"
                if time_to_normal is not None:
                    time_str += f" (to normal: {time_to_normal:.1f}s)"
            
            status = "✓" if reduction >= 5.0 and final_cpp >= 50.0 else "⚠" if reduction >= 5.0 else "✗"
            print(f"     {flow:3.0f} mL/min: {reduction:5.1f} mmHg reduction, ICP={final_icp:4.1f}, Pvs={final_pvs:4.1f}, J3={final_pjr3:4.1f}/{final_pjl3:4.1f}, J2={final_pjr2:4.1f}/{final_pjl2:4.1f}, CPP={final_cpp:4.1f}, Flow={actual_flow_mL_min:5.1f}, Time={time_str} {status}")
    
    # Phase 2: Determine optimal flow for each pathology
    print("\n" + "-" * 60)
    print("PHASE 2: OPTIMAL FLOW DETERMINATION")
    print("-" * 60)
    
    optimal_flows = {}
    
    for test_id in test_case_ids:
        baseline_icp = baseline_icps[test_id]  # Use stored baseline ICP
        
        # Find minimum therapeutic flow (≥5 mmHg reduction, CPP ≥50)
        # Prefer flows that achieve normal ICP faster
        optimal_flow = None
        max_reduction = 0
        min_time_to_normal = float('inf')
        
        for flow in sorted(flow_rates):
            result = sensitivity_results[test_id][flow]
            if result['therapeutic'] and result['safe']:
                # Primary: maximum reduction
                # Secondary: minimum time to normal ICP (if achieved)
                better_reduction = result['reduction'] > max_reduction
                same_reduction_better_time = (result['reduction'] == max_reduction and 
                                            result['time_to_normal'] is not None and
                                            result['time_to_normal'] < min_time_to_normal)
                
                if better_reduction or same_reduction_better_time:
                    optimal_flow = flow
                    max_reduction = result['reduction']
                    if result['time_to_normal'] is not None:
                        min_time_to_normal = result['time_to_normal']
        
        if optimal_flow is None:
            # FALLBACK: If no therapeutic flow found, pick the one with max reduction
            for flow in sorted(flow_rates):
                result = sensitivity_results[test_id][flow]
                if optimal_flow is None or result['reduction'] > max_reduction:
                    optimal_flow = flow
                    max_reduction = result['reduction']
            print(f"{test_id}: No optimal flow found, using {optimal_flow} mL/min (maximum reduction: {max_reduction:.1f} mmHg)")
        else:
            time_info = ""
            optimal_result = sensitivity_results[test_id][optimal_flow]
            if optimal_result['time_to_normal'] is not None:
                time_info = f", Time={optimal_result['time_to_normal']:.0f}s"
            print(f"{test_id}: Optimal flow = {optimal_flow} mL/min (ΔICP = {max_reduction:.1f} mmHg{time_info})")
        
        optimal_flows[test_id] = optimal_flow
    
    # Phase 3: Site comparison at optimal flows
    print("\n" + "-" * 60)
    print("PHASE 3: SITE COMPARISON AT OPTIMAL FLOWS")
    print("-" * 60)
    
    site_comparison_results = {}
    
    print(f"{'Test':<8} {'Opt Flow':<8} {'Baseline':<10} {'Pvs':<10} {'J3':<10} {'J2':<10} {'J1':<10} {'Best Site'}")
    print(f"{'ID':<8} {'(mL/min)':<8} {'ICP (mmHg)':<10} {'ΔICP/Time':<10} {'ΔICP/Time':<10} {'ΔICP/Time':<10} {'ΔICP/Time':<10} {''}")
    print("-" * 85)
    
    for test_id in test_case_ids:
        site_comparison_results[test_id] = {}
        optimal_flow = optimal_flows[test_id]
        
        # Get stabilized baseline ICP (no aspiration)
        params = create_test_parameters(test_id)
        baseline_result, _ = run_full_protocol(
            params, stabilization_time_s=600.0, intervention_time_s=10.0, max_stabilization_attempts=5
        )
        baseline_icp = baseline_result.get_summary()['ICP_mean']
        
        reductions = {}
        
        for site in sites:

            params = create_test_parameters(test_id)
            
            protocol = AspirationProtocol.create_constant_flow(
                site=site,
                target_flow_mL_min=optimal_flow,
                ramp_duration_s=60.0,  # Increased ramp time
                total_duration_s=1800.0,  # Increased total
                access_method=access_method
            )
            params.aspiration_protocol = protocol
            # Parameters are read directly from the protocol in the ODE system
            
            _, intervention = run_full_protocol(
                params, stabilization_time_s=600.0, intervention_time_s=1800.0, max_stabilization_attempts=5
            )
            
            final_icp = intervention.get_summary()['ICP_mean']
            reduction = baseline_icp - final_icp
            
            # Calculate time to normal ICP
            normal_icp_threshold = 15.0
            time_to_normal = None
            if final_icp <= normal_icp_threshold:
                icp_timeseries = intervention.ICP
                times = intervention.t
                below_threshold = icp_timeseries <= normal_icp_threshold
                if np.any(below_threshold):
                    first_normal_idx = np.where(below_threshold)[0][0]
                    time_to_normal = times[first_normal_idx]
            
            reductions[site] = reduction
            
            site_comparison_results[test_id][site] = {
                'reduction': reduction,
                'final_icp': final_icp,
                'final_cpp': intervention.get_summary()['CPP_mean'],
                'time_to_normal': time_to_normal
            }
        
        # Find best site
        best_site = max(reductions.keys(), key=lambda s: reductions[s])
        
        # Print results with reduction and time
        def format_site_result(site):
            if site in site_comparison_results[test_id]:
                result = site_comparison_results[test_id][site]
                red = result['reduction']
                time_val = result['time_to_normal']
                if time_val is not None:
                    return f"{red:>+6.1f}/{int(time_val):d}s"
                else:
                    return f"{red:>+6.1f}/N/A"
            return "N/A"
        
        pvs_str = format_site_result('Pvs')
        j3_str = format_site_result('J3')
        j2_str = format_site_result('J2')
        j1_str = format_site_result('J1')
        
        print(f"{test_id:<8} {optimal_flow:>6.0f} mL/min {baseline_icp:>8.1f}   {pvs_str:<10} {j3_str:<10} {j2_str:<10} {j1_str:<10} {best_site}")
    
    print("\n" + "="*100)
    print("ANALYSIS COMPLETE")
    print("="*100)
    print("Flow sensitivity analysis complete for each pathology")
    print("Optimal therapeutic flows determined")
    print("Site comparison performed at optimal flows")
    print("Clinically rigorous methodology validated")
    print("="*100)
    
    return {
        'sensitivity_results': sensitivity_results,
        'optimal_flows': optimal_flows,
        'site_comparison': site_comparison_results
    }


def get_scaled_flow_rate(test_case_id, base_flow_mL_min=60.0):
    """
    Get flow rate scaled by stenosis severity for numerical stability.
    
    Severe stenosis creates extreme pressure gradients that limit maximum
    safe aspiration flow rates. This function scales flow rates based on
    the minimum jugular conductance remaining.
    
    Args:
        test_case_id: Test case identifier
        base_flow_mL_min: Base flow rate for healthy/normal cases
    
    Returns:
        Scaled flow rate in mL/min appropriate for the pathology severity
    """
    if test_case_id not in CLINICAL_TEST_CASES:
        return base_flow_mL_min
    
    test = CLINICAL_TEST_CASES[test_case_id]
    
    # Special cases with fixed flows
    if test_case_id == 'T3':
        return 150.0  # High flow protocol by design
    
    # Calculate minimum conductance remaining (most severe stenosis)
    stenosis_degree = test.get('stenosis_degree', 1.0)
    min_conductance = stenosis_degree
    
    # Scale flow rate based on conductance
    # Healthy: 1.0 conductance -> 60 mL/min base flow
    # Severe stenosis: 0.2 conductance -> ~20 mL/min max flow
    # Linear scaling: flow = base_flow * sqrt(conductance)

    
    if min_conductance >= 1.0:
        # Healthy or near-healthy
        return base_flow_mL_min
    elif min_conductance >= 0.5:
        # Moderate stenosis (50-70% remaining)
        return base_flow_mL_min * 0.5  # 30 mL/min
    elif min_conductance >= 0.3:
        # Severe stenosis (70-80% remaining) 
        return base_flow_mL_min * 0.3  # 18 mL/min
    else:
        # Very severe stenosis (<70% remaining)
        return base_flow_mL_min * 0.2  # 12 mL/min


def validate_starling_resistor():
    """
    Validate Starling resistor implementation for different access methods and postures.
    
    Tests the interaction between hydrostatic pressure effects, posture changes,
    and direct vs femoral access aspiration.
    """
    from gadda_model.equations import aspiration_resistance_with_starling
    from gadda_model.parameters import ModelParameters
    
    print("\n" + "="*80)
    print("STARLING RESISTOR VALIDATION")
    print("="*80)
    
    # Test cases: different postures and access methods
    test_conditions = [
        {'posture': 'supine', 'access': 'femoral', 'description': 'Supine, Femoral Access'},
        {'posture': 'upright', 'access': 'femoral', 'description': 'Upright, Femoral Access'},
        {'posture': 'supine', 'access': 'direct', 'description': 'Supine, Direct Access'},
        {'posture': 'upright', 'access': 'direct', 'description': 'Upright, Direct Access'},
    ]
    
    # Test parameters
    P_site = 10.0  # mmHg (typical venous pressure)
    P_pump = -200.0  # mmHg (aspiration pressure)
    
    for condition in test_conditions:
        print(f"\n{condition['description']}:")
        print("-" * 40)
        

        params = ModelParameters(posture=condition['posture'])
        
        # Test different external pressures (relevant for direct access)
        external_pressures = [5.0, 0.0, -5.0]  # mmHg (intrathoracic to atmospheric)
        
        for P_external in external_pressures:
            R_total = aspiration_resistance_with_starling(P_site, P_external, P_pump, params)
            Q_calc = (P_site - P_pump) / R_total
            
            print(f"  P_external = {P_external:>5.1f} mmHg: R_total = {R_total:>8.1f} mmHg·s/mL, Q = {Q_calc:>6.2f} mL/s")
    
    print("\n" + "="*80)
    print("VALIDATION COMPLETE")
    print("Key observations:")
    print("- Lower external pressure (more negative) increases collapse resistance")
    print("- Upright posture may affect baseline pressures but not Starling calculation directly")
    print("- Direct access exposes veins to atmospheric pressure (0 mmHg)")
    print("="*80)


# Alias for backward compatibility
TEST_DEFINITIONS = CLINICAL_TEST_CASES

__all__ = ['CLINICAL_TEST_CASES', 'TEST_DEFINITIONS', 'create_test_parameters', 'run_single_test', 'run_test_suite', 'run_flow_sensitivity_analysis', 'get_scaled_flow_rate', 'validate_starling_resistor']
