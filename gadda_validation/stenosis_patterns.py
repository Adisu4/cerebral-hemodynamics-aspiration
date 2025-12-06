"""
Stenosis Pattern Validation
============================

Tests various venous stenosis patterns to validate model behavior
with obstructed outflow pathways.
"""

import matplotlib.pyplot as plt
from gadda_model import ModelParameters, run_simulation


# Stenosis pattern definitions
STENOSIS_PATTERNS = {
    'Baseline': {
        'description': 'Normal venous anatomy (no stenosis)',
        'modifications': {}
    },
    'PatternA': {
        'description': 'Proximal azygos blocked + Left jugular blocked',
        'modifications': {
            'Gazy2': 0.0,
            'kjl3': 0.0,
            'kjl2': 0.0,
            'kjl1': 0.0
        }
    },
    'PatternA': {
        'description': 'Both jugular veins + Proximal azygos blocked (complete blockage)',
        'modifications': {
            'kjr3': 0.0,
            'kjl3': 0.0,
            'kjr2': 0.0,
            'kjl2': 0.0,
            'Gazy2': 0.0
        }
    },
    'PatternC': {
        'description': 'Both jugular veins blocked (azygos patent)',
        'modifications': {
            'kjr3': 0.0,
            'kjl3': 0.0,
            'kjr2': 0.0,
            'kjl2': 0.0
        }
    },
    'PatternD': {
        'description': 'Halved conductance of lower jugular segments',
        'modifications': {
            'kjr2': 0.5,
            'kjl2': 0.5
        }
    }
}


def apply_stenosis_pattern(params, pattern_name):
    """
    Apply predefined stenosis pattern to model parameters.
    
    Args:
        params: ModelParameters object to modify
        pattern_name: Name of stenosis pattern from STENOSIS_PATTERNS
    
    Returns:
        Modified params object
    """
    if pattern_name not in STENOSIS_PATTERNS:
        raise ValueError(f"Unknown stenosis pattern: {pattern_name}")
    
    pattern = STENOSIS_PATTERNS[pattern_name]
    mods = pattern['modifications']
    
    print(f"  Applying {pattern_name}: {pattern['description']}")
    
    for param, value in mods.items():
        if param in ['kjr2', 'kjl2'] and pattern_name == 'PatternD':
            # For PatternD, multiply (halve) existing values
            current = getattr(params, param)
            setattr(params, param, current * value)
            print(f"    {param} *= {value}")
        else:
            # For other patterns, set directly
            setattr(params, param, value)
            print(f"    {param} = {value}")
    
    return params


def test_stenosis_patterns(posture='supine', save_plot=True, show_plot=False):
    """
    Test all stenosis patterns and compare ICP elevation.
    
    Args:
        posture: 'supine' or 'upright'
        save_plot: Save figure to file
        show_plot: Display figure interactively
    
    Returns:
        Dictionary with results for each pattern
    """
    print("\n" + "="*70)
    print(f"VALIDATION: STENOSIS PATTERNS ({posture.upper()})")
    print("="*70)
    print("Analyzing venous outflow obstruction patterns")
    print("-"*70)
    
    results = {}
    
    for pattern_name in STENOSIS_PATTERNS.keys():
        print(f"\n{'='*70}")
        print(f"Pattern: {pattern_name}")
        print(f"{'='*70}")
        
        # Create parameters and apply stenosis
        params = ModelParameters(posture=posture)
        apply_stenosis_pattern(params, pattern_name)
        

        print(f"  Running simulation...")
        result = run_simulation(params, duration_s=900)
        
        # Extract final values
        summary = result.get_summary(window_s=60)
        
        results[pattern_name] = {
            'ICP': summary['ICP_mean'],
            'ICP_std': summary['ICP_std'],
            'Pvs': summary['Pvs_mean'],
            'Pvs_std': summary['Pvs_std'],
            'CPP': summary['CPP_mean'],
            'description': STENOSIS_PATTERNS[pattern_name]['description']
        }
        
        print(f"  Final ICP: {summary['ICP_mean']:.2f} ± {summary['ICP_std']:.2f} mmHg")
        print(f"  Final Pvs: {summary['Pvs_mean']:.2f} ± {summary['Pvs_std']:.2f} mmHg")
        print(f"  Final CPP: {summary['CPP_mean']:.2f} ± {summary['CPP_std']:.2f} mmHg")
    
    # Summary comparison
    print("\n" + "="*70)
    print("STENOSIS PATTERN COMPARISON")
    print("="*70)
    print(f"{'Pattern':<15} {'ICP (mmHg)':<15} {'Pvs (mmHg)':<15} {'ΔICP':<10}")
    print("-"*70)
    
    baseline_icp = results['Baseline']['ICP']
    for pattern_name, data in results.items():
        delta_icp = data['ICP'] - baseline_icp
        print(f"{pattern_name:<15} {data['ICP']:>10.2f}     {data['Pvs']:>10.2f}     {delta_icp:>+7.2f}")
    print("="*70)
    
    # Create bar chart
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    patterns = list(results.keys())
    icp_values = [results[p]['ICP'] for p in patterns]
    pvs_values = [results[p]['Pvs'] for p in patterns]
    
    colors = ['green', 'orange', 'red', 'darkred', 'purple']
    
    # ICP comparison
    bars1 = ax1.bar(patterns, icp_values, color=colors, edgecolor='black', linewidth=1.5)
    ax1.axhline(baseline_icp, color='black', linestyle='--', linewidth=2, label='Baseline')
    ax1.set_ylabel('ICP (mmHg)', fontsize=12, fontweight='bold')
    ax1.set_title('ICP Elevation by Stenosis Pattern', fontsize=14, fontweight='bold')
    ax1.set_ylim([0, max(icp_values) * 1.2])
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.legend()
    
    # Add value labels on bars
    for bar, val in zip(bars1, icp_values):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f}', ha='center', va='bottom', fontweight='bold')
    
    # Pvs comparison
    bars2 = ax2.bar(patterns, pvs_values, color=colors, edgecolor='black', linewidth=1.5)
    baseline_pvs = results['Baseline']['Pvs']
    ax2.axhline(baseline_pvs, color='black', linestyle='--', linewidth=2, label='Baseline')
    ax2.set_ylabel('Pvs (mmHg)', fontsize=12, fontweight='bold')
    ax2.set_title('Venous Sinus Pressure by Stenosis Pattern', fontsize=14, fontweight='bold')
    ax2.set_ylim([0, max(pvs_values) * 1.2])
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.legend()
    
    # Add value labels on bars
    for bar, val in zip(bars2, pvs_values):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.1f}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    
    if save_plot:
        import os
        os.makedirs('validation_plots', exist_ok=True)
        filename = f'validation_plots/validation_stenosis_patterns_{posture}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved: {filename}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    return results


__all__ = ['STENOSIS_PATTERNS', 'apply_stenosis_pattern', 'test_stenosis_patterns']
