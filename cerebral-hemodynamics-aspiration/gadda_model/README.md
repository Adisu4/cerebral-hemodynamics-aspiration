# Cerebral Hemodynamics Model Package

Computational model of intracranial pressure (ICP) dynamics and cerebral venous outflow based on Gadda et al. (2015).

## Overview

This package implements a 15-state ordinary differential equation (ODE) system modeling:
- Intracranial pressure dynamics
- Cerebral venous outflow pathways
- CSF formation and absorption
- Cerebral autoregulation
- Venous aspiration interventions

## Installation

```python
# Package is located in the project directory
# Import directly from gadda_model/
from gadda_model import ModelParameters, run_simulation
```

## Basic Usage

### 1. Healthy Baseline Simulation

```python
from gadda_model import ModelParameters, run_simulation

# Create default parameters (healthy baseline)
params = ModelParameters(posture='supine')

# Run 15-minute simulation
result = run_simulation(params, duration_s=900)

# Access results
print(f"Final ICP: {result.ICP[-1]:.1f} mmHg")
print(f"Mean ICP: {result.ICP.mean():.1f} mmHg")
```

### 2. IIH Phenotype with Pathologies

```python
from gadda_model import (
    ModelParameters,
    run_simulation,
    apply_csf_pathology,
    apply_venous_stenosis,
    apply_autoregulation_impairment
)

# Start with healthy parameters
params = ModelParameters(posture='supine')

# Apply pathologies
apply_csf_pathology(params, R0=1000)  # Elevated CSF resistance
apply_venous_stenosis(params, severity=0.3, location='bilateral_J3')  # 70% stenosis
apply_autoregulation_impairment(params, Gaut=0.0)  # Disabled autoregulation

# Run simulation
result = run_simulation(params, duration_s=900)
print(f"IIH ICP: {result.ICP[-1]:.1f} mmHg")
```

### 3. Aspiration Protocol

```python
from gadda_model import (
    ModelParameters,
    AspirationProtocol,
    run_full_protocol,
    apply_venous_stenosis
)

# Create IIH phenotype
params = ModelParameters(posture='supine')
params.R0 = 1000
params.Gaut = 0.0
apply_venous_stenosis(params, severity=0.3, location='bilateral_J3')

# Define aspiration protocol
protocol = AspirationProtocol.create_constant_flow(
    site='Pvs',  # Venous sinus
    target_flow_mL_min=120.0,
    ramp_duration_s=120.0,
    total_duration_s=900.0
)
params.aspiration_protocol = protocol

# Run two-stage protocol (baseline + intervention)
baseline, intervention = run_full_protocol(
    params,
    stabilization_time_s=900.0,
    intervention_time_s=900.0
)

# Compare results
baseline_icp = baseline.ICP[-60:].mean()
intervention_icp = intervention.ICP[-60:].mean()
print(f"Baseline ICP: {baseline_icp:.1f} mmHg")
print(f"Intervention ICP: {intervention_icp:.1f} mmHg")
print(f"ICP Reduction: {baseline_icp - intervention_icp:.1f} mmHg")
```

## Core Components

### ModelParameters
Stores all model parameters including:
- Hemodynamic constants (pressures, resistances, capacitances)
- Physiological parameters (autoregulation, CSF dynamics)
- Aspiration protocol configuration
- Body posture (supine/upright with hydrostatic corrections)

### Simulation Functions
- `run_simulation()` - Single simulation run
- `run_full_protocol()` - Two-stage protocol (baseline + intervention)

### Aspiration Protocols
- `AspirationProtocol.create_constant_flow()` - Steady flow rate
- `AspirationProtocol.create_stepwise_escalation()` - Gradual increase

### Pathology Functions
- `apply_csf_pathology()` - Modify CSF outflow resistance
- `apply_venous_stenosis()` - Apply jugular stenosis
- `apply_autoregulation_impairment()` - Modify autoregulation gain
- `create_iih_phenotype()` - Combined pathologies

## Aspiration Sites

- **Pvs**: Venous sinuses (transverse/sigmoid)
- **Pv**: Cerebral veins
- **J3**: Upper jugular vein (C1-C3 level)
- **J2**: Middle jugular vein (C3-C5 level)
- **J1**: Lower jugular vein (C5-C7 level)

## Validated Parameters

Based on Gadda et al. (2015) and clinical validation:

- **R0 = 526.3 mmHg·s/mL** - CSF outflow resistance (normal)
- **Rf = 2380 mmHg·s/mL** - CSF formation resistance
- **Rla = 0.6 mmHg·s/mL** - Large artery resistance
- **Gaut = 3.0** - Autoregulation gain (normal)
- **Pa = 100 mmHg** - Arterial pressure
- **Pvp = 5.0 mmHg** - Venous peripheral pressure

## Clinical Benchmarks

### IIH Severity Classification
- **Normal**: ICP < 15 mmHg
- **Borderline**: ICP 15-20 mmHg
- **Mild IIH**: ICP 20-25 mmHg
- **Moderate IIH**: ICP 25-35 mmHg
- **Severe IIH**: ICP 35-45 mmHg
- **Critical IIH**: ICP > 45 mmHg

### Therapeutic Targets
- **Minimum reduction**: ≥5 mmHg (prevent vision loss)
- **VP shunt efficacy**: 10-15 mmHg reduction
- **Diamox efficacy**: 4-8 mmHg reduction

### Safety Limits
- **Minimum ICP**: 7 mmHg (cerebral hypoperfusion risk)
- **Minimum CPP**: 50 mmHg (CPP = MAP - ICP)
- **Minimum jugular pressure**: 3-4 mmHg (venous collapse risk)

## Example: Complete Clinical Test

```python
from gadda_model import *

# Setup IIH patient (moderate severity)
params = ModelParameters(posture='supine')
params.R0 = 800  # Mild CSF pathology
params.Gaut = 0.0  # Disabled autoregulation
params.kjr3 *= 0.4  # 60% bilateral stenosis
params.kjl3 *= 0.4
params.kjr2 *= 0.4
params.kjl2 *= 0.4

# Aspiration protocol
protocol = AspirationProtocol.create_constant_flow(
    site='J3',
    target_flow_mL_min=120.0,
    ramp_duration_s=120.0,
    total_duration_s=1500.0
)
params.aspiration_protocol = protocol

# Run protocol
baseline, intervention = run_full_protocol(params)

# Analyze efficacy
summary_baseline = baseline.get_summary()
summary_intervention = intervention.get_summary()

reduction = summary_baseline['ICP_mean'] - summary_intervention['ICP_mean']
print(f"\nResults:")
print(f"  Baseline ICP: {summary_baseline['ICP_mean']:.1f} ± {summary_baseline['ICP_std']:.2f} mmHg")
print(f"  Intervention ICP: {summary_intervention['ICP_mean']:.1f} ± {summary_intervention['ICP_std']:.2f} mmHg")
print(f"  Reduction: {reduction:.1f} mmHg")
print(f"  Therapeutic: {'Yes' if reduction >= 5.0 else 'No'}")
```

## Reference

Gadda G, et al. (2015). "A new hemodynamic model for the study of intracranial hypertension." Fluids and Barriers of the CNS, 12:15.
Ursino M, Lodi CA. (1997). A simple mathematical model of the interaction between intracranial pressure and cerebral hemodynamics. Journal of Applied Physiology, 82(4):1256–1269.
Marmarou A, et al. (1975). A nonlinear analysis of the CSF system and intracranial pressure dynamics. Journal of Neurosurgery.
L. O. Müller and E. F. Toro, “Enhanced global mathematical model for studying cerebral venous blood flow,” J. Biomech., vol. 47, no. 13, pp. 3361–3372, Oct. 2014, doi: 10.1016/j.jbiomech.2014.08.005.
N. Westerhof, J.-W. Lankhaar, and B. E. Westerhof, “The arterial Windkessel,” Med. Biol. Eng. Comput., vol. 47, no. 2, pp. 131–141, Feb. 2009, doi: 10.1007/s11517-008-0359-2.
M. Ursino, “A mathematical study of human intracranial hydrodynamics part 1—The cerebrospinal fluid pulse pressure,” Ann. Biomed. Eng., vol. 16, no. 4, pp. 379–401, July 1988, doi: 10.1007/BF02364625. 
Ghitti, R., Toro, E.F., Müller, L.O.Nonlinear lumped-parameter models for blood flow simulations in networks of vessels. Appl Math Model 108, 270–303 (2022).https://doi.org/10.1016/j.apm.2022.03.013



## Notes

- All pressures in mmHg
- All flows in mL/s (protocols specified in mL/min for convenience)
- All resistances in mmHg·s/mL
- All capacitances in mL/mmHg
- Time in seconds
- Use `max_step=1.0` for stability in `run_simulation()`
