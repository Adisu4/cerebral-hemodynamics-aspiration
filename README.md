# Cerebral hemodynamics during controlled venous aspiration

A lumped-parameter model of intracranial pressure (ICP), cerebral circulation,
extracranial venous drainage, and controlled venous aspiration. The model is
adapted from published cerebral and venous hemodynamics literature and examines
how aspiration location and flow rate affect ICP in a constructed TBI-like
parameter state.

> This is a computational research model. It has not been validated for patient
> treatment, device operation, or clinical decision-making.

## Overview

The model contains 14 dynamic states representing intracranial and
extracranial pressures together with cerebral autoregulation. Pial arterial
compliance is calculated from the autoregulatory state rather than integrated
as an additional state.

Key features include:

- intracranial arterial, capillary, venous, and cerebrospinal-fluid dynamics;
- bilateral jugular, vertebral, collateral, azygos, and caval pathways;
- pressure-dependent jugular conductance;
- cerebral autoregulation and pressure-dependent intracranial compliance;
- aspiration at the cerebral-vein, venous-sinus, J3, and J2 nodes;
- matched baseline and intervention simulations;
- explicit convergence tests and terminal-window averaging; and
- dose-response, site-comparison, reinfusion, solver, and sensitivity analyses.

Only the supine configuration is implemented. Primary aspiration is modeled as
external withdrawal. Reinfusion is evaluated as a separate comparator.

## Reference reproduction

The reference parameter set reproduces selected **supine baseline** pressures
and flows reported by Gadda et al. (2015). This is an implementation-
reproduction check, not independent biological validation of the aspiration
model. Upright posture is not included in the maintained implementation.

The repository retains the exact result-producing source and its SHA-256
manifest under `reference_implementation/`. Regression tests compare the
maintained modules with that archived implementation and verify every saved run
report.

## Installation

Python 3.10–3.13 is supported.

```bash
git clone https://github.com/Adisu4/cerebral-hemodynamics-aspiration.git
cd cerebral-hemodynamics-aspiration

python -m venv .venv

# Linux or macOS
source .venv/bin/activate

# Windows PowerShell
.venv\Scripts\Activate.ps1

python -m pip install -e ".[test,visualization]"
python -m pytest
```

## Quick start

Run the predefined baseline, dose-response, site-comparison, and sensitivity
experiments:

```bash
cerebral-hemodynamics-run --output-dir results/simulations
```

Generate the scientific figures from the saved reference results:

```bash
cerebral-hemodynamics-plot \
  --results-dir reference_results/figure_data \
  --output-dir figures
```

The plotting command does not rerun the simulations. It produces editable
vector PDF files and 600-dpi PNG files using Arial and a color-vision-safe
palette.

## Python example

```python
from pathlib import Path
import numpy as np

from cerebral_hemodynamics_aspiration import integrate, make_tbi_parameters

parameters = make_tbi_parameters()
output_dir = Path("results/example")

baseline = integrate(
    parameters,
    parameters.initial_state(),
    label="tbi_baseline",
    output_dir=output_dir,
)

baseline_state = np.asarray(baseline["terminal_window_mean_state"], dtype=float)
intervention = integrate(
    parameters,
    baseline_state,
    label="pv_240",
    site="Pv",
    flow_ml_min=240.0,
    output_dir=output_dir,
)

delta_icp = (
    baseline["summary"]["icp_mean_mmhg"]
    - intervention["summary"]["icp_mean_mmhg"]
)
print(f"Predicted ICP reduction: {delta_icp:.2f} mmHg")
```

## Repository structure

```text
src/cerebral_hemodynamics_aspiration/
    model.py            Governing equations and flow calculations
    parameters.py       Physiological and numerical parameter sets
    simulation.py       Integration, convergence, and result serialization
    experiments.py      Dose-response, site-comparison, and sensitivity runs
    visualization.py    Scientific visualization of saved results

reference_implementation/  Archived result-producing source and source manifest
reference_results/         Saved simulation reports and aggregate data
figures/                   Vector PDF and 600-dpi PNG figures
tests/                     Equation, provenance, and regression tests
docs/                      Model scope, provenance, and release notes
```

Development and reuse should use the modules under `src/`. The historical
filenames in `reference_implementation/` are retained unchanged because their
hashes are embedded in the saved results.

## Model states

| State | Description |
|---|---|
| `Pic` | Intracranial pressure |
| `Ppa` | Pial arterial pressure |
| `Pv` | Cerebral-vein pressure |
| `Pvs` | Venous-sinus pressure |
| `Pjr3`, `Pjl3` | Right and left jugular J3 pressure |
| `Pjr2`, `Pjl2` | Right and left jugular J2 pressure |
| `Pc3`, `Pc2` | Collateral pathway pressures |
| `Pvv` | Vertebral venous pressure |
| `Pazy` | Azygos pressure |
| `Psvc` | Superior vena cava pressure |
| `xaut` | Cerebral autoregulatory state |

The J1 confluence pressure is solved algebraically. Pial arterial compliance
`Cpa` is derived from `xaut` at each model evaluation.

## Governing relationships

The intracranial pressure balance is evaluated as

```text
dPic/dt = [Cpa·d(Ppa-Pic)/dt + Cvi·d(Pv-Pic)/dt
           + dCpa/dt·(Ppa-Pic) + Qf - Q0] / Cic
```

Cerebral-vein aspiration enters the venous mass balance as an external sink:

```text
d(Pv-Pic)/dt = (Qin,v - Qout,v - Qasp) / Cvi
```

Jugular conductance follows the pressure-dependent sigmoid used by Gadda et
al.:

```text
Gj = kj [1 + (2/π) arctan((Pupstream - Pexternal)/A)]²
```

See `model.py` and `parameters.py` for the complete implemented equations,
units, branch conditions, and parameter provenance.

## Reference results

For the implemented TBI-like parameter state, the converged simulations give:

| Aspiration condition | ICP reduction (mmHg) |
|---|---:|
| Cerebral-vein node, 240 mL/min | 3.0682 |
| Cerebral-vein node, 480 mL/min | 6.2045 |
| Venous-sinus node, 240 mL/min | 0.2689 |
| Bilateral J3 nodes, 240 mL/min | 0.1450 |
| Bilateral J2 nodes, 240 mL/min | 0.1212 |

These are deterministic model predictions rather than measured biological
effects. The constructed TBI-like state and aspiration predictions require
independent experimental validation.

## Reproducibility

The test suite checks:

- equation-level parity with the archived result-producing implementation;
- mass balance and model-domain conditions;
- convergence and terminal pressure residuals;
- source and input fingerprints;
- consistency of 91 saved simulation reports; and
- independently solved local equilibria for the principal results.

See [model scope and limitations](docs/MODEL_SCOPE_AND_LIMITATIONS.md) and
[source provenance](docs/PROVENANCE.md) before interpreting or citing the
results.

## Scientific basis and citation

The cerebral and venous model is based principally on:

- Gadda G, et al. (2015), “A new hemodynamic model shows that temporal venous
  stenosis can cause idiopathic intracranial hypertension,” *Acta
  Neurochirurgica*.
- Ursino and Lodi (1997), cerebral autoregulation modeling.
- Marmarou et al. (1975), intracranial pressure-volume and CSF dynamics.

The venous-aspiration extension is being prepared for publication by Adisu
Mengesha Assefa and collaborators. Add a formal software citation only after
the author list and manuscript citation have been approved by all contributors.

## Contact

**Adisu Mengesha Assefa**

- PhD Student and Graduate Assistant
- Department of Biomechanics, University of Nebraska at Omaha
- Email: aassefa@unomaha.edu

## Acknowledgments

This work was conducted under the supervision of Prof. Majid Jadidi,
Department of Biomechanics, University of Nebraska at Omaha. The model builds
on the cerebral hemodynamics, autoregulation, and CSF literature cited above.

## License

Licensed under the MIT License. The license permits software reuse; it does not
certify clinical performance or fitness for medical use.
