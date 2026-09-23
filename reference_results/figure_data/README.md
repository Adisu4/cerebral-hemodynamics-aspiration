# Saved inputs for manuscript figures and tables

This directory contains the four trajectories used by Figure 3. The 89
converged run reports are in `reference_results/reports/`; the study summary and
manifest are one directory above. The visualization module combines these
files to regenerate Figures 2–6 and Tables 2–3 without rerunning the model.
Every nonzero-flow report specifies `return_fraction=1.0` and
`return_site=lower_SVC_state`. These are deterministic model predictions, not
clinical observations or population uncertainty estimates.
