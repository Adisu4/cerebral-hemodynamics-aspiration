# Source and result provenance

The maintained numerical implementation is in
`src/cerebral_hemodynamics_aspiration/`. Each of the 89 saved reports records
the SHA-256 hashes of `model.py`, `parameters.py`, `simulation.py`, and
`experiments.py`, along with the parameter set, solver settings, initial state,
return fraction, convergence diagnostics, and runtime environment.

The repository reorganization changes paths and presentation only. These four
numerical modules and all saved report/trajectory bytes are unchanged. The
source-hash tests verify this directly. All aspiration runs in the released
study use `return_fraction=1.0` and return flow to `lower_SVC_state`.

## Frozen regression fixtures

`tests/reference/` holds three byte-preserved files from the audited source.
The manifest maps their concise filenames to their original names and hashes.
They are test fixtures, not an alternative study workflow. Regression tests
compare parameter sets, equation evaluations, and short trajectories.

The complete earlier source snapshot and prior figure layouts are preserved at
[commit daaa0ff](https://github.com/Adisu4/cerebral-hemodynamics-aspiration/tree/daaa0ff817ba75f58e6abe6cffabea261d5157e0).
Obsolete study drivers and superseded figures have been removed from the
current tree. Earlier no-return results remain in Git history and are not the
primary analysis.

## Numerical and display data

`data/reports/` contains the current converged reports. The study manifest and
summary in `data/` describe the same reinfusion study. `data/trajectories/`
contains the four time courses used in Figures 2 and S3 and their checksums.

`figures.py` generates the figures and tables without altering scientific
inputs. Its provenance manifest records the precise reports, trajectories,
generator, and display files used. The main results appear in Table 2; the
complete rate-by-site data appear in supplementary Table S7.

The supine reference reproduction checks implementation against published
values; it is not independent biological validation of aspiration predictions.
