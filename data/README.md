# Reference simulation data

The 89 JSON files in `reports/` are the converged runs of the equal lower-SVC
return study. Each report includes inputs, source hashes, numerical settings,
state summaries, and convergence diagnostics. `study_manifest.json` and
`study_summary.csv` summarize that study; `key_results.json` provides the
principal regression values.

`trajectories/` contains the elevated baseline and three intervention time
courses used in Figures 2 and S3. Every nonzero-flow report uses return fraction 1.0.
The files are preserved byte for byte from the verified study. Their original
run metadata is retained rather than rewritten after moving the directory.

Use `cerebral-hemodynamics-plot` to generate the current figures and tables.
Use `cerebral-hemodynamics-run --output-dir results/simulations` for a new run.
