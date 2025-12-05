"""
Plotting Script for Aspiration Study Results
Generates publication-quality figures for T1 and T2 test cases
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from aspiration_study.test_cases import create_test_parameters
from gadda_model.aspiration import AspirationProtocol
from adaptive_stabilization import run_adaptive_stabilization
from scipy.integrate import solve_ivp
from gadda_model.solver import SimulationResult
from gadda_model.equations import gadda_ode_system_FIXED

def smooth_signal(data, window_length=51, polyorder=3):
    """
    Smooth signal using Savitzky-Golay filter
    
    Args:
        data: Input signal array
        window_length: Length of filter window (must be odd)
        polyorder: Order of polynomial fit
    
    Returns:
        Smoothed signal
    """
    if len(data) < window_length:
        window_length = len(data) if len(data) % 2 == 1 else len(data) - 1
        if window_length < polyorder + 2:
            return data
    
    return savgol_filter(data, window_length, polyorder)


def run_full_simulation_t1():
    """
    Run complete T1 simulation (baseline + aspiration)
    
    Returns:
        baseline: SimulationResult for baseline phase
        intervention: SimulationResult for aspiration phase
        baseline_icp: Baseline ICP value (mmHg)
        final_icp: Final ICP value (mmHg)
    """
    print("\n" + "="*80)
    print("Running T1 Simulation (Mild IIH)")
    print("="*80)
    
    params = create_test_parameters('T1')
    
    # Phase 1: Baseline stabilization
    print("\nPhase 1: Baseline Stabilization...")
    max_time = 14400 if params.R0 > 1000 else 7200
    baseline = run_adaptive_stabilization(
        params,
        max_time_s=max_time,
        check_interval_s=300,
        verbose=True
    )
    
    baseline_icp = baseline.get_summary()['ICP_mean']
    print(f"✓ Baseline ICP: {baseline_icp:.2f} mmHg")
    
    # Phase 2: Aspiration
    print("\nPhase 2: Aspiration (240 mL/min)...")
    flow_rate_mL_min = 240.0
    aspiration_duration_s = 300.0
    ramp_duration_s = 60.0
    target_flow_mL_s = flow_rate_mL_min / 60.0
    
    protocol = AspirationProtocol.create_constant_flow(
        site='Pv',
        target_flow_mL_min=flow_rate_mL_min,
        ramp_duration_s=ramp_duration_s,
        total_duration_s=aspiration_duration_s
    )
    params.aspiration_protocol = protocol
    
    y0 = baseline.get_final_state()
    num_eval_points = max(100, int(aspiration_duration_s / 5))
    t_eval = np.linspace(0, aspiration_duration_s, num_eval_points)
    
    def ode_with_callback(t, y):
        return gadda_ode_system_FIXED(t, y, params, flow_rate_callback=lambda t, pic: protocol.get_flow_rate(t))
    
    sol = solve_ivp(
        ode_with_callback,
        t_span=(0, aspiration_duration_s),
        y0=y0,
        method='BDF',
        t_eval=t_eval,
        max_step=5.0,
        rtol=1e-3,
        atol=1e-5,
        dense_output=False
    )
    
    class MinimalSol:
        def __init__(self, t, y, success, message):
            self.t = t
            self.y = y
            self.success = success
            self.message = message
    
    minimal_sol = MinimalSol(sol.t, sol.y, sol.success, sol.message)
    intervention = SimulationResult(minimal_sol, params)
    
    final_icp = intervention.get_summary()['ICP_mean']
    print(f"✓ Final ICP: {final_icp:.2f} mmHg")
    print(f"✓ Reduction: {baseline_icp - final_icp:.2f} mmHg")
    
    return baseline, intervention, baseline_icp, final_icp


def run_full_simulation_t2(flow_rate_mL_min):
    """
    Run complete T2 simulation (baseline + aspiration) for given flow rate
    
    Args:
        flow_rate_mL_min: Aspiration flow rate (mL/min)
    
    Returns:
        baseline: SimulationResult for baseline phase
        intervention: SimulationResult for aspiration phase
        baseline_icp: Baseline ICP value (mmHg)
        final_icp: Final ICP value (mmHg)
    """
    print(f"\n{'='*80}")
    print(f"Running T2 Simulation (Severe TBI) - {flow_rate_mL_min} mL/min")
    print("="*80)
    
    # Create parameters
    params = create_test_parameters('T2')
    
    # Phase 1: Baseline stabilization
    print("\nPhase 1: Baseline Stabilization...")
    max_time = 14400 if params.R0 > 1000 else 7200
    baseline = run_adaptive_stabilization(
        params,
        max_time_s=max_time,
        check_interval_s=300,
        verbose=True
    )
    
    baseline_icp = baseline.get_summary()['ICP_mean']
    print(f"✓ Baseline ICP: {baseline_icp:.2f} mmHg")
    
    # Phase 2: Aspiration
    print(f"\nPhase 2: Aspiration ({flow_rate_mL_min} mL/min)...")
    aspiration_duration_s = 300.0
    ramp_duration_s = 60.0
    target_flow_mL_s = flow_rate_mL_min / 60.0
    
    protocol = AspirationProtocol.create_constant_flow(
        site='Pv',
        target_flow_mL_min=flow_rate_mL_min,
        ramp_duration_s=ramp_duration_s,
        total_duration_s=aspiration_duration_s
    )
    params.aspiration_protocol = protocol
    
    y0 = baseline.get_final_state()
    num_eval_points = max(100, int(aspiration_duration_s / 5))
    t_eval = np.linspace(0, aspiration_duration_s, num_eval_points)
    
    def ode_with_callback(t, y):
        return gadda_ode_system_FIXED(t, y, params, flow_rate_callback=lambda t, pic: protocol.get_flow_rate(t))
    
    sol = solve_ivp(
        ode_with_callback,
        t_span=(0, aspiration_duration_s),
        y0=y0,
        method='BDF',
        t_eval=t_eval,
        max_step=5.0,
        rtol=1e-3,
        atol=1e-5,
        dense_output=False
    )
    
    class MinimalSol:
        def __init__(self, t, y, success, message):
            self.t = t
            self.y = y
            self.success = success
            self.message = message
    
    minimal_sol = MinimalSol(sol.t, sol.y, sol.success, sol.message)
    intervention = SimulationResult(minimal_sol, params)
    
    final_icp = intervention.get_summary()['ICP_mean']
    print(f"✓ Final ICP: {final_icp:.2f} mmHg")
    print(f"✓ Reduction: {baseline_icp - final_icp:.2f} mmHg")
    
    return baseline, intervention, baseline_icp, final_icp


def plot_t1_timecourse(baseline, intervention, save_path='output/T1_ICP_timecourse.png'):
    """
    Plot T1 ICP time-course (baseline + aspiration)
    
    Args:
        baseline: Baseline SimulationResult
        intervention: Intervention SimulationResult
        save_path: Path to save figure
    """
    # Extract data
    t_baseline = baseline.t
    icp_baseline = baseline.ICP
    
    # Shift intervention time to absolute time
    baseline_end = t_baseline[-1]
    t_intervention = intervention.t + baseline_end
    icp_intervention = intervention.ICP
    
    # Combine for continuous plot
    t_combined = np.concatenate([t_baseline, t_intervention])
    icp_combined = np.concatenate([icp_baseline, icp_intervention])
    
    # Smooth signal
    icp_smooth = smooth_signal(icp_combined, window_length=51, polyorder=3)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot ICP
    ax.plot(t_combined / 60, icp_smooth, 'b-', linewidth=2, label='ICP')
    
    # Mark aspiration start
    ax.axvline(baseline_end / 60, color='r', linestyle='--', linewidth=2, 
               label='Aspiration Start', alpha=0.7)
    
    # Add shaded region for aspiration period
    ax.axvspan(baseline_end / 60, t_combined[-1] / 60, alpha=0.2, color='red', 
               label='Aspiration Period (240 mL/min)')
    
    # Formatting
    ax.set_xlabel('Time (minutes)', fontsize=14, fontweight='bold')
    ax.set_ylabel('ICP (mmHg)', fontsize=14, fontweight='bold')
    ax.set_title('T1 (Mild IIH): ICP Time-Course During Venous Aspiration', 
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)
    
    # Add annotations
    baseline_icp = np.mean(icp_baseline[-100:])
    final_icp = np.mean(icp_intervention[-20:])
    reduction = baseline_icp - final_icp
    
    ax.text(0.02, 0.98, f'Baseline ICP: {baseline_icp:.2f} mmHg\nFinal ICP: {final_icp:.2f} mmHg\nReduction: {reduction:.2f} mmHg ({reduction/baseline_icp*100:.1f}%)',
            transform=ax.transAxes, fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {save_path}")
    plt.show()


def plot_t2_timecourse_comparison(baseline_300, intervention_300, 
                                   baseline_360, intervention_360,
                                   baseline_480, intervention_480,
                                   save_path='output/T2_ICP_timecourse_comparison.png'):
    """
    Plot T2 ICP time-course comparison for multiple flow rates
    
    Args:
        baseline_120, baseline_240, baseline_360: Baseline SimulationResult for each flow rate
        intervention_120, intervention_240, intervention_360: Intervention SimulationResult for each flow rate
        save_path: Path to save figure
    """
    # Use first baseline as reference (they should be identical)
    t_baseline = baseline_300.t
    icp_baseline = baseline_300.ICP
    baseline_end = t_baseline[-1]
    
    # Prepare data for each flow rate
    flow_rates = [300, 360, 480]
    interventions = [intervention_300, intervention_360, intervention_480]
    line_styles = ['-', '--', '-.']
    line_markers = ['', '', '*']
    colors = ['b', 'g', 'r']
    labels = ['300 mL/min', '360 mL/min', '480 mL/min']
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 7))
    
    # Plot shared baseline (only once)
    icp_baseline_smooth = smooth_signal(icp_baseline, window_length=51, polyorder=3)
    ax.plot(t_baseline / 60, icp_baseline_smooth, 'k-', linewidth=2.5, 
            label='Shared Baseline', alpha=0.7)
    
    # Plot each flow rate's intervention phase
    for i, (flow, interv, style, marker, color, label) in enumerate(
        zip(flow_rates, interventions, line_styles, line_markers, colors, labels)):
        
        # Shift intervention time to absolute time
        t_intervention = interv.t + baseline_end
        icp_intervention = interv.ICP
        
        # Smooth
        icp_smooth = smooth_signal(icp_intervention, window_length=21, polyorder=3)
        
        # Plot
        if marker:
            ax.plot(t_intervention / 60, icp_smooth, linestyle=style, color=color,
                   linewidth=2.5, marker=marker, markersize=8, markevery=5,
                   label=f'Aspiration: {label}', alpha=0.9)
        else:
            ax.plot(t_intervention / 60, icp_smooth, linestyle=style, color=color,
                   linewidth=2.5, label=f'Aspiration: {label}', alpha=0.9)
    
    # Mark aspiration start
    ax.axvline(baseline_end / 60, color='orange', linestyle=':', linewidth=3, 
               label='Aspiration Start', alpha=0.8)
    
    # Add shaded region for aspiration period
    ax.axvspan(baseline_end / 60, (baseline_end + 300) / 60, alpha=0.15, color='orange')
    
    # Formatting
    ax.set_xlabel('Time (minutes)', fontsize=14, fontweight='bold')
    ax.set_ylabel('ICP (mmHg)', fontsize=14, fontweight='bold')
    ax.set_title('T2 (Severe TBI): ICP Time-Course Comparison - Multi-Flow Aspiration', 
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=11, loc='upper right', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)
    
    # Add annotations for each flow rate
    baseline_icp = np.mean(icp_baseline[-100:])
    final_300 = np.mean(intervention_300.ICP[-20:])
    final_360 = np.mean(intervention_360.ICP[-20:])
    final_480 = np.mean(intervention_480.ICP[-20:])
    
    annotation_text = f'Baseline ICP: {baseline_icp:.2f} mmHg\n\n'
    annotation_text += f'300 mL/min: {final_300:.2f} mmHg (Δ{baseline_icp-final_300:.2f} mmHg)\n'
    annotation_text += f'360 mL/min: {final_360:.2f} mmHg (Δ{baseline_icp-final_360:.2f} mmHg)\n'
    annotation_text += f'480 mL/min: {final_480:.2f} mmHg (Δ{baseline_icp-final_480:.2f} mmHg)'
    
    ax.text(0.02, 0.98, annotation_text,
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {save_path}")
    plt.show()


def plot_combined_histogram(t1_baseline, t1_final, 
                            t2_baseline, t2_final_300, t2_final_360, t2_final_480,
                            save_path='output/Combined_ICP_histogram.png'):
    """
    Plot combined histogram comparing baseline and final ICP values
    
    Args:
        t1_baseline: T1 baseline ICP (mmHg)
        t1_final: T1 final ICP (mmHg)
        t2_baseline: T2 baseline ICP (mmHg)
        t2_final_120, t2_final_240, t2_final_360: T2 final ICP for each flow rate (mmHg)
        save_path: Path to save figure
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # --- T1 Histogram ---
    categories_t1 = ['Baseline', 'Final\n(240 mL/min)']
    values_t1 = [t1_baseline, t1_final]
    colors_t1 = ['lightcoral', 'lightgreen']
    
    bars_t1 = ax1.bar(categories_t1, values_t1, color=colors_t1, edgecolor='black', 
                      linewidth=2, width=0.6, alpha=0.8)
    
    # Add value labels on bars
    for bar, val in zip(bars_t1, values_t1):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.2f} mmHg',
                ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    # Add reduction annotation
    reduction_t1 = t1_baseline - t1_final
    ax1.annotate('', xy=(0, t1_final), xytext=(0, t1_baseline),
                arrowprops=dict(arrowstyle='<->', color='blue', lw=2))
    ax1.text(0.15, (t1_baseline + t1_final)/2, 
            f'Δ{reduction_t1:.2f} mmHg\n({reduction_t1/t1_baseline*100:.1f}%)',
            fontsize=11, color='blue', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    
    ax1.set_ylabel('ICP (mmHg)', fontsize=14, fontweight='bold')
    ax1.set_title('T1 (Mild IIH): Baseline vs. Post-Aspiration', 
                  fontsize=14, fontweight='bold')
    ax1.set_ylim(0, max(values_t1) * 1.3)
    ax1.grid(axis='y', alpha=0.3)
    ax1.tick_params(labelsize=11)
    
    # --- T2 Histogram ---
    categories_t2 = ['Baseline', '300\nmL/min', '360\nmL/min', '480\nmL/min']
    values_t2 = [t2_baseline, t2_final_300, t2_final_360, t2_final_480]
    colors_t2 = ['lightcoral', 'lightblue', 'lightgreen', 'lightyellow']
    
    bars_t2 = ax2.bar(categories_t2, values_t2, color=colors_t2, edgecolor='black', 
                      linewidth=2, width=0.6, alpha=0.8)
    
    # Add value labels on bars
    for bar, val in zip(bars_t2, values_t2):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.2f} mmHg',
                ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    # Add reduction annotations for each flow rate
    reductions_t2 = [t2_baseline - t2_final_300, 
                     t2_baseline - t2_final_360, 
                     t2_baseline - t2_final_480]
    
    annotation_text = 'Reductions:\n'
    for i, (flow, red) in enumerate(zip([300, 360, 480], reductions_t2)):
        annotation_text += f'{flow} mL/min: Δ{red:.2f} mmHg ({red/t2_baseline*100:.1f}%)\n'
    
    ax2.text(0.98, 0.98, annotation_text.strip(),
            transform=ax2.transAxes, fontsize=10, verticalalignment='top',
            horizontalalignment='right', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    ax2.set_ylabel('ICP (mmHg)', fontsize=14, fontweight='bold')
    ax2.set_title('T2 (Severe TBI): Baseline vs. Post-Aspiration', 
                  fontsize=14, fontweight='bold')
    ax2.set_ylim(0, max(values_t2) * 1.3)
    ax2.grid(axis='y', alpha=0.3)
    ax2.tick_params(labelsize=11)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {save_path}")
    plt.show()


def main():
    """
    Main execution function - generates all plots
    """
    print("\n" + "="*80)
    print("ASPIRATION STUDY: PLOTTING RESULTS")
    print("="*80)
    
    # Create output directory if needed
    import os
    os.makedirs('output', exist_ok=True)
    

    print("\n[1/4] Running T1 simulation...")
    baseline_t1, intervention_t1, t1_baseline_icp, t1_final_icp = run_full_simulation_t1()
    

    print("\n[2/4] Running T2 simulations...")
    baseline_t2_300, intervention_t2_300, t2_baseline_300, t2_final_300 = run_full_simulation_t2(300)
    baseline_t2_360, intervention_t2_360, t2_baseline_360, t2_final_360 = run_full_simulation_t2(360)
    baseline_t2_480, intervention_t2_480, t2_baseline_480, t2_final_480 = run_full_simulation_t2(480)
    
    # Use first baseline as reference (they should be identical)
    t2_baseline_icp = t2_baseline_300
    

    print("\n[3/4] Generating T1 time-course plot...")
    plot_t1_timecourse(baseline_t1, intervention_t1)
    
    print("\n[4/4] Generating T2 time-course comparison plot...")
    plot_t2_timecourse_comparison(baseline_t2_300, intervention_t2_300,
                                   baseline_t2_360, intervention_t2_360,
                                   baseline_t2_480, intervention_t2_480)
    
    print("\n[5/5] Generating combined histogram...")
    plot_combined_histogram(t1_baseline_icp, t1_final_icp,
                           t2_baseline_icp, t2_final_300, t2_final_360, t2_final_480)
    
    print("\n" + "="*80)
    print("Publication figures generation complete")
    print("="*80)
    print("Output files:")
    print("  - output/T1_ICP_timecourse.png")
    print("  - output/T2_ICP_timecourse_comparison.png")
    print("  - output/Combined_ICP_histogram.png")
    print("="*80)


if __name__ == '__main__':
    main()
