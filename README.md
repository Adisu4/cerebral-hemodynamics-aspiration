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
- dose-response, site-comparison, solver, and sensitivity analyses.

Only the supine configuration is implemented. In every primary aspiration run,
flow extracted at the selected site is returned at the same instantaneous rate
to the dynamic lower-SVC compartment (`return_fraction=1.0`). Arterial and
central venous boundary pressures remain fixed. The model does not simulate an
extracorporeal device or a systemic physiological response.

## Reference reproduction

The reference parameter set reproduces selected **supine baseline** pressures
and flows reported by Gadda et al. (2015). This is an implementation-
reproduction check, not independent biological validation of the aspiration
model. Upright posture is not included in the maintained implementation.

The maintained implementation is in `src/cerebral_hemodynamics_aspiration/`.
Saved reports record the numerical source hashes and every run input.
Three frozen fixtures under `tests/reference/` support equation and trajectory
regression checks; earlier study scripts and figures are retained in Git history.

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

Run the complete baseline, dose-response, site-comparison, matched-sensitivity,
and solver-check study with equal lower-SVC return. The command writes 89
simulation reports plus the study manifest and summary:

```bash
cerebral-hemodynamics-run --output-dir results/simulations
```

Generate the manuscript and supplementary figures and tables from the versioned
reference results:

```bash
cerebral-hemodynamics-plot \
  --results-dir data/reports \
  --output-dir figures \
  --tables-dir tables
```

The command reads the four plotted trajectories from
`data/trajectories/`. To build from a fresh run, use
`--results-dir results/simulations`; the trajectory directory then defaults to
the same run directory. Plotting does not rerun the model. It produces vector
PDFs and 600-dpi PNGs with Arial labels; current figures are versioned in `figures/`, and numerical tables in `tables/`. Main Figure 1 (model schematic) and Table 1 (model
parameters) are manuscript artwork and are not produced by this plotting
command. See [reproduction instructions](docs/reproduction.md) for the complete display mapping.

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
    return_fraction=1.0,
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
    model.py          Governing equations and flow calculations
    parameters.py     Physiological and post-traumatic parameter sets
    simulation.py     Integration, convergence, and result serialization
    experiments.py    Baseline, aspiration, sensitivity, and solver runs
    figures.py        Figures and tables from saved results

data/
    reports/          89 converged simulation reports
    trajectories/     Four trajectories used for Figures 2 and S3
figures/
    supplementary/    Reference reproduction, sensitivity, and resistance diagnostics
tables/
    supplementary/    Complete aspiration results and model comparisons
tests/
    reference/        Frozen fixtures for regression tests
docs/                 Reproduction, model scope, and provenance
```

The numerical modules retain their source hashes. Filenames describe their
scientific role; figure and table numbers match the revised manuscript and
supplement. Only the current displays are present on the main branch.

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

Cerebral-vein extraction enters its venous mass balance as a local sink:

```text
d(Pv-Pic)/dt = (Qin,v - Qout,v - Qasp) / Cvi
```

An equal instantaneous `Qasp` enters the lower-SVC pressure balance as
`Qreturn`. Baseline and control runs have zero extraction and zero return.

Jugular conductance follows the pressure-dependent sigmoid used by Gadda et
al.:

```text
Gj = kj [1 + (2/π) arctan((Pupstream - Pexternal)/A)]²
```

See `model.py` and `parameters.py` for the complete implemented equations,
units, branch conditions, and parameter provenance.

## Reference results

For the implemented TBI-like parameter state, the converged simulations give:

| Aspiration site | ΔICP at 240 mL/min (mmHg) | ΔICP at 480 mL/min (mmHg) |
|---|---:|---:|
| Cerebral veins | 3.033 | 6.133 |
| Venous sinus | 0.234 | 0.469 |
| Bilateral J3 | 0.110 | 0.221 |
| Bilateral J2 | 0.086 | 0.174 |

Baseline ICP is 22.698 mmHg. Table 2 reports these nominal and maximum-rate
results; supplementary Table S7 contains final ICP and ΔICP at all seven rates.
Every intervention uses equal lower-SVC return and the final 120-s mean after
convergence. Supplementary Figure S3 shows pressures and terminal resistance
before and during aspiration; Table S9 reports actual convergence diagnostics
for the baseline and all 28 primary interventions. The site/rate bar and line
charts have been removed to avoid duplicating the table.

These are deterministic model predictions rather than measured biological
effects. The constructed TBI-like state and aspiration predictions require
independent experimental validation.

## Reproducibility

The test suite checks:

- equation-level parity with the frozen reference fixtures;
- mass balance and model-domain conditions;
- convergence and terminal pressure residuals;
- source and input fingerprints;
- consistency of 89 saved simulation reports; and
- independently solved local equilibria for the principal results.

See [model scope and limitations](docs/model_scope.md) and
[source provenance](docs/provenance.md) before interpreting or citing the
results.

## Scientific basis and citation

The cerebral and venous model is based principally on:

- Gadda G, et al. (2015), “A new hemodynamic model for the study of cerebral
  venous outflow,” *American Journal of Physiology–Heart and Circulatory
  Physiology* 308:H217–H231. [doi:10.1152/ajpheart.00469.2014](https://doi.org/10.1152/ajpheart.00469.2014).
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
