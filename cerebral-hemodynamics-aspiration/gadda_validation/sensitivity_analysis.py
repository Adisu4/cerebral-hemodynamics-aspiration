"""
Sensitivity Analysis Module
============================

Comprehensive sensitivity analysis of cerebral venous drainage pathways.
Tests model robustness under bilateral and unilateral venous blockages.

Analyzes:
- Venous sinus pressure (Pvs) response to drainage blockages
- Jugular outflow (Qsvc) response to pathway occlusions  
- Vertebral outflow (Qvv) response to stenosis
- Mass balance verification across all conditions
"""

import numpy as np
import matplotlib.pyplot as plt
from gadda_model import ModelParameters, run_simulation
from gadda_model.equations import (
    jugular_conductance
)


def sensitivity_pvs_bilateral(save_plot=True, show_plot=False):
    """
    Bilateral blockage sensitivity: Pvs response to drainage pathway occlusions.
    
    Tests complete bilateral blockage of:
    - J3 (upper jugular)
    - J2 (middle jugular)
    - J1 (lower jugular)
    - Vertebral veins
    - Collateral pathways
    
    Returns:
        tuple: (supine_pvs_values, upright_pvs_values)
    """
    print("\n" + "="*80)
    print("SENSITIVITY ANALYSIS: Pvs (BILATERAL BLOCKAGE)")
    print("="*80)
    print("Analyzing complete occlusion of bilateral drainage pathways")
    print("-"*80)
    
    conditions = ['Basal', 'J3=0', 'J2=0', 'J1=0', 'VV=0', 'C3=0']
    pvs_supine = []
    pvs_upright = []
    
    for condition in conditions:
        print(f"\n{condition}:")
        
        # Supine
        params_sup = ModelParameters(posture='supine')
        if condition == 'J3=0':
            params_sup.kjr3 = params_sup.kjl3 = 0.0
        elif condition == 'J2=0':
            params_sup.kjr2 = params_sup.kjl2 = 0.0
        elif condition == 'J1=0':
            params_sup.kjr1 = params_sup.kjl1 = 0.0
        elif condition == 'VV=0':
            params_sup.Gvvr = params_sup.Gvvl = 0.0
        elif condition == 'C3=0':
            params_sup.Gc3 = 0.0
        
        result_sup = run_simulation(params_sup, duration_s=500)
        summary_sup = result_sup.get_summary(window_s=60)
        pvs_supine.append(summary_sup['Pvs_mean'])
        print(f"  Supine Pvs:  {summary_sup['Pvs_mean']:.2f} mmHg")
        
        # Upright
        params_upr = ModelParameters(posture='upright')
        if condition == 'J3=0':
            params_upr.kjr3 = params_upr.kjl3 = 0.0
        elif condition == 'J2=0':
            params_upr.kjr2 = params_upr.kjl2 = 0.0
        elif condition == 'J1=0':
            params_upr.kjr1 = params_upr.kjl1 = 0.0
        elif condition == 'VV=0':
            params_upr.Gvvr = params_upr.Gvvl = 0.0
        elif condition == 'C3=0':
            params_upr.Gc3 = 0.0
        
        result_upr = run_simulation(params_upr, duration_s=500)
        summary_upr = result_upr.get_summary(window_s=60)
        pvs_upright.append(summary_upr['Pvs_mean'])
        print(f"  Upright Pvs: {summary_upr['Pvs_mean']:.2f} mmHg")
    

    print("\n" + "="*80)
    print(f"{'Condition':<12} {'Supine Pvs':<15} {'Upright Pvs':<15} {'ΔPvs':<12}")
    print("-"*80)
    for i, cond in enumerate(conditions):
        delta = pvs_upright[i] - pvs_supine[i]
        print(f"{cond:<12} {pvs_supine[i]:>10.2f}     {pvs_upright[i]:>10.2f}     {delta:>+8.2f}")
    print("="*80)
    
    # Plot
    if save_plot or show_plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(conditions))
        width = 0.35
        
        ax.bar(x - width/2, pvs_supine, width, label='Supine',
               color='black', edgecolor='black', linewidth=1.5)
        ax.bar(x + width/2, pvs_upright, width, label='Upright',
               color='white', edgecolor='black', linewidth=1.5, hatch='//')
        
        ax.set_xlabel('Drainage Pathway', fontsize=12, fontweight='bold')
        ax.set_ylabel('Pvs (mmHg)', fontsize=12, fontweight='bold')
        ax.set_title('Bilateral Blockage: Pvs Sensitivity', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(conditions, fontsize=10)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        if save_plot:
            import os
            os.makedirs('validation_plots', exist_ok=True)
            plt.savefig('validation_plots/validation_pvs_sensitivity_bilateral.png', dpi=300, bbox_inches='tight')
            print("\nPlot saved: validation_plots/validation_pvs_sensitivity_bilateral.png")
        if show_plot:
            plt.show()
        else:
            plt.close()
    
    return pvs_supine, pvs_upright


def sensitivity_pvs_unilateral(save_plot=True, show_plot=False):
    """
    Unilateral blockage sensitivity: Pvs response to right-side occlusions.
    
    Tests right-side only blockage for asymmetric stenosis modeling.
    
    Returns:
        tuple: (supine_pvs_values, upright_pvs_values)
    """
    print("\n" + "="*80)
    print("SENSITIVITY ANALYSIS: Pvs (UNILATERAL BLOCKAGE - Right)")
    print("="*80)
    print("Analyzing right-side only occlusion (asymmetric stenosis)")
    print("-"*80)
    
    conditions = ['Basal', 'Jr3=0', 'Jr2=0', 'Jr1=0', 'VVr=0', 'C3=0']
    pvs_supine = []
    pvs_upright = []
    
    for condition in conditions:
        print(f"\n{condition}:")
        
        # Supine
        params_sup = ModelParameters(posture='supine')
        if condition == 'Jr3=0':
            params_sup.kjr3 = 0.0
        elif condition == 'Jr2=0':
            params_sup.kjr2 = 0.0
        elif condition == 'Jr1=0':
            params_sup.kjr1 = 0.0
        elif condition == 'VVr=0':
            params_sup.Gvvr = 0.0
        elif condition == 'C3=0':
            params_sup.Gc3 = 0.0
        
        result_sup = run_simulation(params_sup, duration_s=500)
        summary_sup = result_sup.get_summary(window_s=60)
        pvs_supine.append(summary_sup['Pvs_mean'])
        print(f"  Supine Pvs:  {summary_sup['Pvs_mean']:.2f} mmHg")
        
        # Upright
        params_upr = ModelParameters(posture='upright')
        if condition == 'Jr3=0':
            params_upr.kjr3 = 0.0
        elif condition == 'Jr2=0':
            params_upr.kjr2 = 0.0
        elif condition == 'Jr1=0':
            params_upr.kjr1 = 0.0
        elif condition == 'VVr=0':
            params_upr.Gvvr = 0.0
        elif condition == 'C3=0':
            params_upr.Gc3 = 0.0
        
        result_upr = run_simulation(params_upr, duration_s=500)
        summary_upr = result_upr.get_summary(window_s=60)
        pvs_upright.append(summary_upr['Pvs_mean'])
        print(f"  Upright Pvs: {summary_upr['Pvs_mean']:.2f} mmHg")
    
    # Summary
    print("\n" + "="*80)
    print(f"{'Condition':<12} {'Supine Pvs':<15} {'Upright Pvs':<15} {'ΔPvs':<12}")
    print("-"*80)
    for i, cond in enumerate(conditions):
        delta = pvs_upright[i] - pvs_supine[i]
        print(f"{cond:<12} {pvs_supine[i]:>10.2f}     {pvs_upright[i]:>10.2f}     {delta:>+8.2f}")
    print("="*80)
    
    # Plot
    if save_plot or show_plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(conditions))
        width = 0.35
        
        ax.bar(x - width/2, pvs_supine, width, label='Supine',
               color='black', edgecolor='black', linewidth=1.5)
        ax.bar(x + width/2, pvs_upright, width, label='Upright',
               color='white', edgecolor='black', linewidth=1.5, hatch='//')
        
        ax.set_xlabel('Drainage Pathway (Right)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Pvs (mmHg)', fontsize=12, fontweight='bold')
        ax.set_title('Unilateral Blockage: Pvs Sensitivity', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(conditions, fontsize=10)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        if save_plot:
            import os
            os.makedirs('validation_plots', exist_ok=True)
            plt.savefig('validation_plots/validation_pvs_sensitivity_unilateral.png', dpi=300, bbox_inches='tight')
            print("\nPlot saved: validation_plots/validation_pvs_sensitivity_unilateral.png")
        if show_plot:
            plt.show()
        else:
            plt.close()
    
    return pvs_supine, pvs_upright


def calculate_jugular_outflow(result, params):
    """
    Calculate total jugular outflow (Qsvc) from J2 to SVC.
    
    Args:
        result: SimulationResult object
        params: ModelParameters object
    
    Returns:
        float: Jugular outflow in mL/s
    """

    Pjr2 = result.Pjr2[-1]
    Pjl2 = result.Pjl2[-1]
    Psvc = result.Psvc[-1]
    

    Gjr1 = jugular_conductance(Pjr2, params.Pj1ext, params.kjr1, params.A)
    Gjl1 = jugular_conductance(Pjl2, params.Pj1ext, params.kjl1, params.A)
    
    # Total jugular outflow
    Qsvc = (Pjr2 - Psvc) * Gjr1 + (Pjl2 - Psvc) * Gjl1
    
    return Qsvc


def sensitivity_qsvc_bilateral(save_plot=True, show_plot=False):
    """
    Bilateral blockage sensitivity: jugular outflow response.
    
    Returns:
        tuple: (supine_qsvc_values, upright_qsvc_values)
    """
    print("\n" + "="*80)
    print("SENSITIVITY ANALYSIS: Jugular Outflow (BILATERAL BLOCKAGE)")
    print("="*80)
    
    conditions = ['Basal', 'J3=0', 'J2=0', 'J1=0', 'VV=0', 'C3=0']
    qsvc_supine = []
    qsvc_upright = []
    
    for condition in conditions:
        print(f"\n{condition}:")
        
        # Supine
        params_sup = ModelParameters(posture='supine')
        if condition == 'J3=0':
            params_sup.kjr3 = params_sup.kjl3 = 0.0
        elif condition == 'J2=0':
            params_sup.kjr2 = params_sup.kjl2 = 0.0
        elif condition == 'J1=0':
            params_sup.kjr1 = params_sup.kjl1 = 0.0
        elif condition == 'VV=0':
            params_sup.Gvvr = params_sup.Gvvl = 0.0
        elif condition == 'C3=0':
            params_sup.Gc3 = 0.0
        
        result_sup = run_simulation(params_sup, duration_s=500)
        qsvc_sup = calculate_jugular_outflow(result_sup, params_sup)
        qsvc_supine.append(qsvc_sup)
        print(f"  Supine Qsvc:  {qsvc_sup:.2f} mL/s")
        
        # Upright
        params_upr = ModelParameters(posture='upright')
        if condition == 'J3=0':
            params_upr.kjr3 = params_upr.kjl3 = 0.0
        elif condition == 'J2=0':
            params_upr.kjr2 = params_upr.kjl2 = 0.0
        elif condition == 'J1=0':
            params_upr.kjr1 = params_upr.kjl1 = 0.0
        elif condition == 'VV=0':
            params_upr.Gvvr = params_upr.Gvvl = 0.0
        elif condition == 'C3=0':
            params_upr.Gc3 = 0.0
        
        result_upr = run_simulation(params_upr, duration_s=500)
        qsvc_upr = calculate_jugular_outflow(result_upr, params_upr)
        qsvc_upright.append(qsvc_upr)
        print(f"  Upright Qsvc: {qsvc_upr:.2f} mL/s")
    
    # Summary
    print("\n" + "="*80)
    print(f"{'Condition':<12} {'Supine Qsvc':<15} {'Upright Qsvc':<15} {'ΔQsvc':<12}")
    print("-"*80)
    for i, cond in enumerate(conditions):
        delta = qsvc_upright[i] - qsvc_supine[i]
        print(f"{cond:<12} {qsvc_supine[i]:>10.2f}     {qsvc_upright[i]:>10.2f}     {delta:>+8.2f}")
    print("="*80)
    
    # Plot
    if save_plot or show_plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(conditions))
        width = 0.35
        
        ax.bar(x - width/2, qsvc_supine, width, label='Supine',
               color='black', edgecolor='black', linewidth=1.5)
        ax.bar(x + width/2, qsvc_upright, width, label='Upright',
               color='white', edgecolor='black', linewidth=1.5, hatch='//')
        
        ax.set_xlabel('Drainage Pathway', fontsize=12, fontweight='bold')
        ax.set_ylabel('Qsvc (mL/s)', fontsize=12, fontweight='bold')
        ax.set_title('Bilateral Blockage: Jugular Outflow Sensitivity', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(conditions, fontsize=10)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        if save_plot:
            import os
            os.makedirs('validation_plots', exist_ok=True)
            plt.savefig('validation_plots/validation_qsvc_sensitivity_bilateral.png', dpi=300, bbox_inches='tight')
            print("\nPlot saved: validation_plots/validation_qsvc_sensitivity_bilateral.png")
        if show_plot:
            plt.show()
        else:
            plt.close()
    
    return qsvc_supine, qsvc_upright


def sensitivity_qsvc_unilateral(save_plot=True, show_plot=False):
    """
    Unilateral blockage sensitivity: jugular outflow response (right side).
    
    Returns:
        tuple: (supine_qsvc_values, upright_qsvc_values)
    """
    print("\n" + "="*80)
    print("SENSITIVITY ANALYSIS: Jugular Outflow (UNILATERAL - Right)")
    print("="*80)
    
    conditions = ['Basal', 'Jr3=0', 'Jr2=0', 'Jr1=0', 'VVr=0', 'C3=0']
    qsvc_supine = []
    qsvc_upright = []
    
    for condition in conditions:
        print(f"\n{condition}:")
        
        # Supine
        params_sup = ModelParameters(posture='supine')
        if condition == 'Jr3=0':
            params_sup.kjr3 = 0.0
        elif condition == 'Jr2=0':
            params_sup.kjr2 = 0.0
        elif condition == 'Jr1=0':
            params_sup.kjr1 = 0.0
        elif condition == 'VVr=0':
            params_sup.Gvvr = 0.0
        elif condition == 'C3=0':
            params_sup.Gc3 = 0.0
        
        result_sup = run_simulation(params_sup, duration_s=500)
        qsvc_sup = calculate_jugular_outflow(result_sup, params_sup)
        qsvc_supine.append(qsvc_sup)
        print(f"  Supine Qsvc:  {qsvc_sup:.2f} mL/s")
        
        # Upright
        params_upr = ModelParameters(posture='upright')
        if condition == 'Jr3=0':
            params_upr.kjr3 = 0.0
        elif condition == 'Jr2=0':
            params_upr.kjr2 = 0.0
        elif condition == 'Jr1=0':
            params_upr.kjr1 = 0.0
        elif condition == 'VVr=0':
            params_upr.Gvvr = 0.0
        elif condition == 'C3=0':
            params_upr.Gc3 = 0.0
        
        result_upr = run_simulation(params_upr, duration_s=500)
        qsvc_upr = calculate_jugular_outflow(result_upr, params_upr)
        qsvc_upright.append(qsvc_upr)
        print(f"  Upright Qsvc: {qsvc_upr:.2f} mL/s")
    
    # Summary
    print("\n" + "="*80)
    print(f"{'Condition':<12} {'Supine Qsvc':<15} {'Upright Qsvc':<15} {'ΔQsvc':<12}")
    print("-"*80)
    for i, cond in enumerate(conditions):
        delta = qsvc_upright[i] - qsvc_supine[i]
        print(f"{cond:<12} {qsvc_supine[i]:>10.2f}     {qsvc_upright[i]:>10.2f}     {delta:>+8.2f}")
    print("="*80)
    
    # Plot
    if save_plot or show_plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(conditions))
        width = 0.35
        
        ax.bar(x - width/2, qsvc_supine, width, label='Supine',
               color='black', edgecolor='black', linewidth=1.5)
        ax.bar(x + width/2, qsvc_upright, width, label='Upright',
               color='white', edgecolor='black', linewidth=1.5, hatch='//')
        
        ax.set_xlabel('Drainage Pathway (Right)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Qsvc (mL/s)', fontsize=12, fontweight='bold')
        ax.set_title('Unilateral Blockage: Jugular Outflow Sensitivity', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(conditions, fontsize=10)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        if save_plot:
            import os
            os.makedirs('validation_plots', exist_ok=True)
            plt.savefig('validation_plots/validation_qsvc_sensitivity_unilateral.png', dpi=300, bbox_inches='tight')
            print("\nPlot saved: validation_plots/validation_qsvc_sensitivity_unilateral.png")
        if show_plot:
            plt.show()
        else:
            plt.close()
    
    return qsvc_supine, qsvc_upright


__all__ = [
    'sensitivity_pvs_bilateral',
    'sensitivity_pvs_unilateral',
    'sensitivity_qsvc_bilateral',
    'sensitivity_qsvc_unilateral',
    'sensitivity_qvv_bilateral',
    'sensitivity_qvv_unilateral'
]


def calculate_vertebral_outflow(result, params):
    """
    Calculate total vertebral outflow from Pvs to vertebral veins.
    
    Args:
        result: SimulationResult object
        params: ModelParameters object
    
    Returns:
        float: Vertebral outflow in mL/s
    """
    Pvs = result.Pvs[-1]
    Pvv = result.Pvv[-1]
    
    # Total vertebral flow (bilateral)
    Qvv = (Pvs - Pvv) * params.Gvvr + (Pvs - Pvv) * params.Gvvl
    
    return Qvv


def sensitivity_qvv_bilateral(save_plot=True, show_plot=False):
    """
    Bilateral blockage sensitivity: vertebral outflow response.
    
    Returns:
        tuple: (supine_qvv_values, upright_qvv_values)
    """
    print("\n" + "="*80)
    print("SENSITIVITY ANALYSIS: Vertebral Outflow (BILATERAL BLOCKAGE)")
    print("="*80)
    
    conditions = ['Basal', 'J3=0', 'J2=0', 'J1=0', 'VV=0', 'C3=0']
    qvv_supine = []
    qvv_upright = []
    
    for condition in conditions:
        print(f"\n{condition}:")
        
        # Supine
        params_sup = ModelParameters(posture='supine')
        if condition == 'J3=0':
            params_sup.kjr3 = params_sup.kjl3 = 0.0
        elif condition == 'J2=0':
            params_sup.kjr2 = params_sup.kjl2 = 0.0
        elif condition == 'J1=0':
            params_sup.kjr1 = params_sup.kjl1 = 0.0
        elif condition == 'VV=0':
            params_sup.Gvvr = params_sup.Gvvl = 0.0
        elif condition == 'C3=0':
            params_sup.Gc3 = 0.0
        
        result_sup = run_simulation(params_sup, duration_s=500)
        qvv_sup = calculate_vertebral_outflow(result_sup, params_sup)
        qvv_supine.append(qvv_sup)
        print(f"  Supine Qvv:  {qvv_sup:.2f} mL/s")
        
        # Upright
        params_upr = ModelParameters(posture='upright')
        if condition == 'J3=0':
            params_upr.kjr3 = params_upr.kjl3 = 0.0
        elif condition == 'J2=0':
            params_upr.kjr2 = params_upr.kjl2 = 0.0
        elif condition == 'J1=0':
            params_upr.kjr1 = params_upr.kjl1 = 0.0
        elif condition == 'VV=0':
            params_upr.Gvvr = params_upr.Gvvl = 0.0
        elif condition == 'C3=0':
            params_upr.Gc3 = 0.0
        
        result_upr = run_simulation(params_upr, duration_s=500)
        qvv_upr = calculate_vertebral_outflow(result_upr, params_upr)
        qvv_upright.append(qvv_upr)
        print(f"  Upright Qvv: {qvv_upr:.2f} mL/s")
    
    # Summary
    print("\n" + "="*80)
    print(f"{'Condition':<12} {'Supine Qvv':<15} {'Upright Qvv':<15} {'ΔQvv':<12}")
    print("-"*80)
    for i, cond in enumerate(conditions):
        delta = qvv_upright[i] - qvv_supine[i]
        print(f"{cond:<12} {qvv_supine[i]:>10.2f}     {qvv_upright[i]:>10.2f}     {delta:>+8.2f}")
    print("="*80)
    
    # Plot
    if save_plot or show_plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(conditions))
        width = 0.35
        
        ax.bar(x - width/2, qvv_supine, width, label='Supine',
               color='black', edgecolor='black', linewidth=1.5)
        ax.bar(x + width/2, qvv_upright, width, label='Upright',
               color='white', edgecolor='black', linewidth=1.5, hatch='//')
        
        ax.set_xlabel('Drainage Pathway', fontsize=12, fontweight='bold')
        ax.set_ylabel('Qvv (mL/s)', fontsize=12, fontweight='bold')
        ax.set_title('Bilateral Blockage: Vertebral Outflow Sensitivity', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(conditions, fontsize=10)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        if save_plot:
            import os
            os.makedirs('validation_plots', exist_ok=True)
            plt.savefig('validation_plots/validation_qvv_sensitivity_bilateral.png', dpi=300, bbox_inches='tight')
            print("\nPlot saved: validation_plots/validation_qvv_sensitivity_bilateral.png")
        if show_plot:
            plt.show()
        else:
            plt.close()
    
    return qvv_supine, qvv_upright


def sensitivity_qvv_unilateral(save_plot=True, show_plot=False):
    """
    Unilateral blockage sensitivity: vertebral outflow response (right side).
    
    Returns:
        tuple: (supine_qvv_values, upright_qvv_values)
    """
    print("\n" + "="*80)
    print("SENSITIVITY ANALYSIS: Vertebral Outflow (UNILATERAL - Right)")
    print("="*80)
    
    conditions = ['Basal', 'Jr3=0', 'Jr2=0', 'Jr1=0', 'VVr=0', 'C3=0']
    qvv_supine = []
    qvv_upright = []
    
    for condition in conditions:
        print(f"\n{condition}:")
        
        # Supine
        params_sup = ModelParameters(posture='supine')
        if condition == 'Jr3=0':
            params_sup.kjr3 = 0.0
        elif condition == 'Jr2=0':
            params_sup.kjr2 = 0.0
        elif condition == 'Jr1=0':
            params_sup.kjr1 = 0.0
        elif condition == 'VVr=0':
            params_sup.Gvvr = 0.0
        elif condition == 'C3=0':
            params_sup.Gc3 = 0.0
        
        result_sup = run_simulation(params_sup, duration_s=500)
        qvv_sup = calculate_vertebral_outflow(result_sup, params_sup)
        qvv_supine.append(qvv_sup)
        print(f"  Supine Qvv:  {qvv_sup:.2f} mL/s")
        
        # Upright
        params_upr = ModelParameters(posture='upright')
        if condition == 'Jr3=0':
            params_upr.kjr3 = 0.0
        elif condition == 'Jr2=0':
            params_upr.kjr2 = 0.0
        elif condition == 'Jr1=0':
            params_upr.kjr1 = 0.0
        elif condition == 'VVr=0':
            params_upr.Gvvr = 0.0
        elif condition == 'C3=0':
            params_upr.Gc3 = 0.0
        
        result_upr = run_simulation(params_upr, duration_s=500)
        qvv_upr = calculate_vertebral_outflow(result_upr, params_upr)
        qvv_upright.append(qvv_upr)
        print(f"  Upright Qvv: {qvv_upr:.2f} mL/s")
    
    # Summary
    print("\n" + "="*80)
    print(f"{'Condition':<12} {'Supine Qvv':<15} {'Upright Qvv':<15} {'ΔQvv':<12}")
    print("-"*80)
    for i, cond in enumerate(conditions):
        delta = qvv_upright[i] - qvv_supine[i]
        print(f"{cond:<12} {qvv_supine[i]:>10.2f}     {qvv_upright[i]:>10.2f}     {delta:>+8.2f}")
    print("="*80)
    
    # Plot
    if save_plot or show_plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(conditions))
        width = 0.35
        
        ax.bar(x - width/2, qvv_supine, width, label='Supine',
               color='black', edgecolor='black', linewidth=1.5)
        ax.bar(x + width/2, qvv_upright, width, label='Upright',
               color='white', edgecolor='black', linewidth=1.5, hatch='//')
        
        ax.set_xlabel('Drainage Pathway (Right)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Qvv (mL/s)', fontsize=12, fontweight='bold')
        ax.set_title('Unilateral Blockage: Vertebral Outflow Sensitivity', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(conditions, fontsize=10)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        if save_plot:
            import os
            os.makedirs('validation_plots', exist_ok=True)
            plt.savefig('validation_plots/validation_qvv_sensitivity_unilateral.png', dpi=300, bbox_inches='tight')
            print("\nPlot saved: validation_plots/validation_qvv_sensitivity_unilateral.png")
        if show_plot:
            plt.show()
        else:
            plt.close()
    
    return qvv_supine, qvv_upright
