## Model Validation & Reproduction

To reproduce all validation analyses (baseline, posture, sensitivity, stenosis, mass balance) as in the original Gadda et al. (2015) paper, run:

```bash
python -m gadda_validation.model_validation_reproduction
```

This must be run from the root of the repository (the folder containing `gadda_model/` and `gadda_validation/`).

This will print all results and generate plots matching the original validation script. If you see `ModuleNotFoundError: No module named 'gadda_model'`, make sure you are running from the root and using the `-m` flag as shown above.

# Intracranial Pressure Model with Venous Aspiration

A computational model of intracranial pressure (ICP) dynamics and cerebral venous outflow with venous aspiration capabilities, adapted from published venous hemodynamics literature.

## Overview

This project implements a 15-state ODE system modeling cerebral hemodynamics and evaluates venous aspiration as a therapeutic intervention for elevated ICP in conditions like Idiopathic Intracranial Hypertension (IIH) and traumatic brain injury (TBI).

**Key Features:**
- Multi-compartment cerebral circulation model (arterial, capillary, venous, CSF)
- Cerebral autoregulation with impairment modeling
- Starling resistor mechanics for collapsible veins
- CSF dynamics (formation, absorption, compliance)
- Venous stenosis modeling (transverse sinus, jugular veins)
- Aspiration protocol implementation with multiple access sites
- Adaptive convergence detection for baseline stabilization

## Installation

### Prerequisites
- Python 3.7+
- NumPy, SciPy, Matplotlib

### Setup

```bash
# Clone the repository
git clone https://github.com/adisumengesha/cerebral-hemodynamics-aspiration.git
cd cerebral-hemodynamics-aspiration

# Create virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install package in development mode
pip install -e .
```

## Quick Start

### Run Aspiration Simulations
```bash
python adaptive_stabilization.py
```

### Generate Publication Figures
```bash
python plot_aspiration_results.py
```

## Project Structure

```
├── gadda_model/              # Core model package
│   ├── equations.py          # ODE system (15 states)
│   ├── parameters.py         # Model parameters
│   ├── solver.py             # Numerical integration
│   └── aspiration.py         # Aspiration protocol management
├── aspiration_study/         # Clinical test cases
│   └── test_cases.py         # Phenotypes T0–T3 and test suite helpers
├── gadda_validation/         # Model validation
│   ├── supine_upright.py     # Posture change validation
│   ├── stenosis_patterns.py  # Venous obstruction patterns
│   └── sensitivity_analysis.py # Parameter sensitivity
├── adaptive_stabilization.py # Main simulation runner with convergence detection
└── plot_aspiration_results.py # Generate figures for aspiration studies
```

## Usage Examples

### 1. Baseline ICP Simulation

```python
from gadda_model import ModelParameters, run_simulation

# Healthy baseline
params = ModelParameters(posture='supine')
result = run_simulation(params, duration_s=900)
print(f"ICP: {result.ICP.mean():.1f} mmHg")
```

### 2. IIH Pathology with Aspiration

```python
from aspiration_study.test_cases import create_test_parameters
from gadda_model import AspirationProtocol, run_full_protocol

# Create T1 (Mild IIH) parameters
params = create_test_parameters('T1')

# Add aspiration protocol
protocol = AspirationProtocol.create_constant_flow(
    site='J3',                   # Right jugular J3
    target_flow_mL_min=60.0,     # 60 mL/min
    ramp_duration_s=600.0,       # 10-min ramp
    total_duration_s=1500.0,     # 25-min total
    catheter_resistance=10.0
)
params.aspiration_protocol = protocol

# Run baseline + intervention
baseline, intervention = run_full_protocol(params)
print(f"Baseline ICP: {baseline.ICP.mean():.1f} mmHg")
print(f"Final ICP: {intervention.ICP[-1]:.1f} mmHg")
print(f"Reduction: {baseline.ICP.mean() - intervention.ICP[-1]:.1f} mmHg")
```

### 3. Multi-Site Comparison (built into test suite)

```python
from aspiration_study.test_cases import run_test_suite

results = run_test_suite(
    sites=['Pvs', 'Pv', 'J3', 'J2', 'J1'],
    protocol_type='constant'
)
```

## Running Simulations

### Adaptive Stabilization (Main Simulation Runner)

```python
from adaptive_stabilization import run_two_phase_aspiration_test

# Run T1 (Mild IIH) with aspiration
result = run_two_phase_aspiration_test(
    test_case_id='T1',
    site='J3',
    flow_rate_mL_min=60.0,
    aspiration_duration_s=1500.0
)
```

### Generate Plots

```bash
# Generate all publication figures (T1 and T2 results)
python plot_aspiration_results.py
```

## Clinical Test Cases

| Test ID | Condition | R0 (mmHg·s/mL) | Gaut | Description |
|---------|-----------|----------------|------|-------------|
| T0 | Healthy control | 526.3 | 3.0 | Normal anatomy, reference |
| T1 | Mild IIH | 2000 | 0.5 | Bilateral J3 stenosis, impaired autoregulation |
| T2 | Severe TBI | 1800 | 0.3 | Venous congestion with impaired autoregulation |
| T3 | Acute edema/bleed (low compliance) | 2900 | 1.0 | Elevated elastance, J3 stenosis |

**Aspiration Sites:**
- **Pvs** - Superior sagittal sinus
- **Pv** - Confluence of sinuses
- **J3** - Right internal jugular (J3 level)
- **J2** - Right internal jugular (J2 level)
- **J1** - Right internal jugular (J1 level)

## Model Validation

The model has been validated against:
1. **Gadda et al. (2015)** - Baseline ICP and venous pressures in supine/upright postures
2. **Stenosis patterns** - ICP elevation with varying degrees of venous obstruction
3. **Autoregulation** - CPP maintenance under MAP changes
4. **CSF dynamics** - Formation/absorption balance

## Output Files

All results are saved to:
- `output/` - CSV files with simulation results
- `validation_plots/` - Publication-quality figures (PNG)

## Key Findings

- **Therapeutic threshold**: ≥5 mmHg ICP reduction
- **Optimal flow rates**: 60-90 mL/min (pathology-dependent)
- **Best sites**: J3/J2 show superior efficacy vs. upstream sites
- **Safety**: CPP ≥60 mmHg maintained across protocols
- **Dose-response**: Saturation at high flow rates due to collapsible vein mechanics

## Model Details

### State Variables (15)
1. `Pic` - Intracranial pressure
2. `Ppa` - Pial arterial pressure
3. `Pv` - Venous outflow pressure
4. `Pvs` - Superior sagittal sinus pressure
5. `Pjr3` - Right jugular J3
6. `Pjl3` - Left jugular J3
7. `Pjr2` - Right jugular J2
8. `Pjl2` - Left jugular J2
9. `Pc3` - Collateral c3 pressure
10. `Pc2` - Collateral c2 pressure
11. `Pvv` - Vertebral venous pressure
12. `Pazy` - Azygos pressure
13. `Psvc` - Superior vena cava pressure
14. `xaut` - Autoregulation state
15. `Cpa` - Pial arterial compliance

### Key Equations

**ICP Balance:**
```
dPic/dt = (1/Cic) * [Cpa*dPpa/dt + Cv*dPv_trans/dt + Ccsf*dPcsf/dt]
```

**Venous Aspiration:**
```
Q_asp = (Pv - P_external) / (R_vein + R_catheter)
dPv/dt = dPv_trans/dt + dPic/dt
dPv_trans/dt = (1/Cv) * (Q_in - Q_out - Q_asp)
```

**Starling Resistor (Collapsible Veins):**
```
R_collapse = R0 * exp(-kE * (P_internal - P_external))
```

## Citation

If you use this model in your research, please cite:

**Original Model:**
Gadda G, et al. (2015) "A new hemodynamic model shows that temporal venous stenosis can cause idiopathic intracranial hypertension." *Acta Neurochirurgica*

**Aspiration Extension:**
If you use the aspiration extension, please cite:
Mengesha Assefa, A., "Venous Aspiration Protocols for Intracranial Pressure Modulation: Computational Model Extension," University of Nebraska at Omaha, 2025. [In preparation]

## License

MIT License - See LICENSE file for details

## Contact

**Author:** Adisu Mengesha Assefa  
**Role:** PhD Student & Graduate Assistant  
**Institution:** Department of Biomechanics, University of Nebraska at Omaha  
**Email:** aassefa@unomaha.edu  



## Acknowledgments
This work was carried out under the supervision of Prof. Majid Jadidi, Department of Biomechanics, University of Nebraska at Omaha. I sincerely appreciate his guidance and support throughout the development of this project.

- Original model: Gadda et al. (2015)
- Cerebral autoregulation: Ursino & Lodi (1997)
- CSF dynamics: Marmarou et al. (1975)
