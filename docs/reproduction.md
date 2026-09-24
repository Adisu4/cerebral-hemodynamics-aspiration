# Reproducing the study

Install the package with `python -m pip install -e ".[test,visualization]"` and
run `python -m pytest`. Python 3.10–3.13 is supported. To reproduce the saved
numerical environment, use Python 3.13 and the NumPy/SciPy versions in
`requirements-reproduction.txt`.

## Run the simulations

```bash
cerebral-hemodynamics-run --output-dir results/simulations
```

The study writes 89 simulation reports, trajectories, a summary CSV, and a
manifest. Existing runs are reused only after input and source checks; use
`--force` to recompute them. All primary and comparator aspiration simulations
return the extracted flow to the lower SVC. The manuscript endpoint is the
matched baseline ICP minus the final 120-s mean ICP after convergence.

## Generate figures and tables

From the versioned results:

```bash
cerebral-hemodynamics-plot --results-dir data/reports --output-dir figures --tables-dir tables
```

From a new simulation run:

```bash
cerebral-hemodynamics-plot --results-dir results/simulations --output-dir results/figures --tables-dir results/tables
```

The generator reads saved outputs without rerunning the model. It checks the
return condition, aspiration site, rate, primary solver, and convergence flags.
Figures are vector PDFs and 600-dpi PNGs. The specified font is Arial; install
Arial for an exact typographic match, otherwise Matplotlib uses a substitute.
CSV tables retain nine decimal places; manuscript tables round to three.

## Manuscript displays

| Display | File | Content |
|---|---|---|
| Figure 2 | `figures/figure_2_icp_response` | Elevated baseline, aspiration onset, and convergence |
| Figure 3 | `figures/figure_3_venous_pressure` | Venous pressure reductions at 240 mL/min |
| Table 2 | `tables/table_2_icp_reduction` | Four sites at 240 and 480 mL/min |
| Figure S1 | `figures/supplementary/figure_s1_reference_hemodynamics` | Published supine reference comparison |
| Figure S2 | `figures/supplementary/figure_s2_parameter_sensitivity` | Matched-baseline one-at-a-time sensitivity |
| Table S7 | `tables/supplementary/table_s7_aspiration_response` | All seven rates and four sites, final ICP and ΔICP |
| Table S8 | `tables/supplementary/table_s8_model_comparison` | Primary, physiological, and fixed-resistance comparisons |

Figure 1 is the model schematic retained in the manuscript. Main Table 1 and
supplementary Tables S1–S6 document parameters, initial states, and site
definitions in the Word documents. They are not outputs of the plotting
command. The supplement is cited as Additional file 1.

The rate-response graph is intentionally omitted: its results are reported in
Table 2 and Table S7. Numerical uncertainty is not inferred from repeated
deterministic output samples. The four trajectory files in `data/trajectories/`
are sufficient for Figure 2; all other generated displays use JSON reports.

`figures/provenance.json` records the generator and input/output hashes.
`SHA256SUMS.txt` in each output directory checks the distributed artifacts.
The test suite checks all saved reports against the current equations and
independently solves the principal equilibria.
