# Manuscript and Extended Data displays

`cerebral-hemodynamics-plot` generates the displays below from the full-return
reports in `reference_results/reports/` and the plotted trajectories in
`reference_results/figure_data/`. Main Figure 1 (model schematic) and Table 1
(model parameters) are manuscript artwork and are not generated here.

| Manuscript display | Generated source file (PDF/PNG, or CSV/Markdown) |
|---|---|
| Figure 2: ICP time course | `figure_02_icp_time_course_svc_return` |
| Figure 3: ICP versus extraction and return rate | `figure_03_flow_response_svc_return` |
| Figure 4: venous-pressure reductions, Pv and Pvs | `figure_04_venous_pressure_reductions_svc_return` |
| Table 2: model comparisons at 240 mL/min | `table_02_model_comparisons_svc_return` |
| Extended Data Figure 1: supine Gadda reproduction | `extended_data_figure_01_reference_reproduction` |
| Extended Data Figure 2: parameter sensitivity | `extended_data_figure_02_parameter_sensitivity_svc_return` |
| Extended Data Table 1: complete dose response | `extended_data_table_01_dose_response_svc_return` |

The flow-response plot displays the complete dose series; its numerical values
are also available in Extended Data Table 1. Figure 4 displays only venous
pressures because Figure 3 already presents final ICP. The PDFs are vector
graphics, and the PNGs are 600 dpi. `artifact_provenance.json` records SHA-256
hashes for the visualization code, inputs, and outputs.

The Gadda comparison is a supine implementation-reproduction benchmark. It is
not independent validation of the intervention model or evidence of clinical
efficacy.
