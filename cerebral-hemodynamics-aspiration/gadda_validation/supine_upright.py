"""
Supine vs upright posture validation

Reproduces published postural effects on cerebral hemodynamics for validation.
"""

import numpy as np
import matplotlib.pyplot as plt
from gadda_model import ModelParameters, run_simulation


def compare_postures(save_plot=True, show_plot=False):
    """
    Compare supine vs upright posture flow distributions.
    
    Reproduces published supine/upright posture reference (2015).
    
    Args:
        save_plot: Save figure to file
        show_plot: Display figure interactively
    
    Returns:
        Dictionary with supine and upright results
    """
    print("\n" + "="*70)
    print("VALIDATION: SUPINE VS UPRIGHT POSTURE")
    print("="*70)
    print("Reproducing published supine/upright posture reference (2015)")
    print("-"*70)
    

    print("\nSimulating supine posture...")
    params_supine = ModelParameters(posture='supine')
    result_supine = run_simulation(params_supine, duration_s=300)
    

    print("Simulating upright posture...")
    params_upright = ModelParameters(posture='upright')
    result_upright = run_simulation(params_upright, duration_s=300)
    
    # Extract final values
    def calculate_flows(result, params, posture):
        """Calculate venous outflow distribution."""
        # Get final state
        Pic = result.Pic[-1]
        Ppa = result.Ppa[-1]
        Pv = result.Pv[-1]
        Pvs = result.Pvs[-1]
        Pjr3 = result.Pjr3[-1]
        Pjl3 = result.Pjl3[-1]
        Pc3 = result.Pc3[-1]
        Pvv = result.Pvv[-1]
        Cpa = result.Cpa[-1]
        
        # Apply hydrostatic corrections for upright
        if posture == 'upright':
            h_hydro_j3 = params.rho * params.g * params.h_j3 / 1333.0
            Pjr3_int = Pjr3 - h_hydro_j3
            Pjl3_int = Pjl3 - h_hydro_j3
        else:
            Pjr3_int = Pjr3
            Pjl3_int = Pjl3
        
        # Jugular conductances (collapsibility)
        from gadda_model.equations import jugular_conductance
        Gjr3 = jugular_conductance(Pjr3_int, params.Pj3ext, params.kjr3, params.A)
        Gjl3 = jugular_conductance(Pjl3_int, params.Pj3ext, params.kjl3, params.A)
        
        # Flows
        Q_jugular = (Pvs - Pjr3) * Gjr3 + (Pvs - Pjl3) * Gjl3
        Q_vertebral = (Pvs - Pvv) * params.Gvvr + (Pvs - Pvv) * params.Gvvl
        Q_collateral = (Pvs - Pc3) * params.Gc3
        
        # Cerebral blood flow
        from gadda_model.equations import cerebral_arterioles_resistance
        Rpa = cerebral_arterioles_resistance(Cpa, Ppa, Pic, params)
        Q_cerebral = (params.Pa - Ppa) / (params.Rla + Rpa / 2.0)
        
        return {
            'ICP': Pic,
            'Pvs': Pvs,
            'Q_jugular': Q_jugular,
            'Q_vertebral': Q_vertebral,
            'Q_collateral': Q_collateral,
            'Q_cerebral': Q_cerebral
        }
    
    supine_flows = calculate_flows(result_supine, params_supine, 'supine')
    upright_flows = calculate_flows(result_upright, params_upright, 'upright')
    
    # Print results
    print("\n" + "-"*70)
    print("RESULTS")
    print("-"*70)
    print(f"{'Variable':<20} {'Supine':<15} {'Upright':<15}")
    print("-"*70)
    print(f"{'ICP (mmHg)':<20} {supine_flows['ICP']:>10.2f}     {upright_flows['ICP']:>10.2f}")
    print(f"{'Pvs (mmHg)':<20} {supine_flows['Pvs']:>10.2f}     {upright_flows['Pvs']:>10.2f}")
    print(f"{'Q_jugular (mL/s)':<20} {supine_flows['Q_jugular']:>10.2f}     {upright_flows['Q_jugular']:>10.2f}")
    print(f"{'Q_vertebral (mL/s)':<20} {supine_flows['Q_vertebral']:>10.2f}     {upright_flows['Q_vertebral']:>10.2f}")
    print(f"{'Q_collateral (mL/s)':<20} {supine_flows['Q_collateral']:>10.2f}     {upright_flows['Q_collateral']:>10.2f}")
    print(f"{'Q_cerebral (mL/s)':<20} {supine_flows['Q_cerebral']:>10.2f}     {upright_flows['Q_cerebral']:>10.2f}")
    print("-"*70)
    
    # Create comparison plot (Figure 4 style)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x = np.arange(4)
    width = 0.35
    
    # Reference data from published supine/upright study (2015)
    reference_supine = [11.74, 0.79, 0, 12.5]
    reference_upright = [1.38, 3.39, 7.7, 12.5]
    
    # Model data
    model_supine = [
        supine_flows['Q_jugular'],
        supine_flows['Q_vertebral'],
        supine_flows['Q_collateral'],
        supine_flows['Q_cerebral']
    ]
    model_upright = [
        upright_flows['Q_jugular'],
        upright_flows['Q_vertebral'],
        upright_flows['Q_collateral'],
        upright_flows['Q_cerebral']
    ]
    
    labels = ['Jugular\n(Qj3)', 'Vertebral\n(Qvv)', 'Collateral\n(Qc3)', 'Cerebral\n(Q)']
    
    # Reference data (bars)
    bars1 = ax.bar(x - width/2, reference_supine, width, 
                   label='Reference (2015) - Supine', 
                   color='lightblue', edgecolor='black')
    bars2 = ax.bar(x + width/2, reference_upright, width, 
                   label='Reference (2015) - Upright', 
                   color='lightcoral', edgecolor='black')
    
    # Model data (scatter overlay)
    supine_x = x - width/2 - 0.06
    upright_x = x + width/2 + 0.06
    ax.scatter(supine_x, model_supine, s=140, c='navy', marker='o',
               edgecolors='white', linewidths=1.5, label='Model - Supine', zorder=5)
    ax.scatter(upright_x, model_upright, s=140, c='darkred', marker='s',
               edgecolors='white', linewidths=1.5, label='Model - Upright', zorder=5)
    
    ax.set_ylabel('Flow (mL/s)', fontsize=12, fontweight='bold')
    ax.set_title('Validation: Flow Distribution vs. published reference (2015)', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim([0, 15])
    
    plt.tight_layout()
    
    if save_plot:
        filename = 'validation_supine_upright.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved: {filename}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    print("="*70)
    
    return {
        'supine': supine_flows,
        'upright': upright_flows,
        'reference_supine': reference_supine,
        'reference_upright': reference_upright
    }


def simulate_posture_change(t_change=300, t_total=600, save_plot=True, show_plot=False):
    """
    Simulate dynamic posture change from supine to upright.
    
    Args:
        t_change: Time of posture change (seconds)
        t_total: Total simulation time (seconds)
        save_plot: Save figure to file
        show_plot: Display figure interactively
    
    Returns:
        Dictionary with combined time series
    """
    print("\n" + "="*70)
    print("VALIDATION: POSTURE CHANGE (SUPINE → UPRIGHT)")
    print("="*70)
    print(f"Posture change at t = {t_change} seconds")
    print(f"Total simulation time = {t_total} seconds")
    print("-"*70)
    
    # Phase 1: Supine
    print("\nPhase 1: Supine posture...")
    params_supine = ModelParameters(posture='supine')
    result_supine = run_simulation(params_supine, duration_s=t_change)
    
    # Phase 2: Upright (use final state from supine as initial condition)
    print("Phase 2: Upright posture...")
    params_upright = ModelParameters(posture='upright')
    y0_upright = result_supine.get_final_state()
    result_upright = run_simulation(params_upright, duration_s=t_total-t_change, y0=y0_upright)
    
    # Combine results
    t_combined = np.concatenate([result_supine.t, result_upright.t + t_change])
    ICP_combined = np.concatenate([result_supine.ICP, result_upright.ICP])
    Pvs_combined = np.concatenate([result_supine.Pvs, result_upright.Pvs])
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    ax1.plot(t_combined, ICP_combined, 'b-', linewidth=2)
    ax1.axvline(t_change, color='red', linestyle='--', linewidth=2, label='Posture Change')
    ax1.set_ylabel('ICP (mmHg)', fontsize=12, fontweight='bold')
    ax1.set_title('Posture Change: Supine → Upright', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.plot(t_combined, Pvs_combined, 'g-', linewidth=2)
    ax2.axvline(t_change, color='red', linestyle='--', linewidth=2, label='Posture Change')
    ax2.set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Pvs (mmHg)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    
    if save_plot:
        import os
        os.makedirs('validation_plots', exist_ok=True)
        filename = f'validation_plots/validation_posture_change_{t_change}s.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved: {filename}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    print("="*70)
    
    return {
        't': t_combined,
        'ICP': ICP_combined,
        'Pvs': Pvs_combined,
        't_change': t_change
    }


__all__ = ['compare_postures', 'simulate_posture_change']
