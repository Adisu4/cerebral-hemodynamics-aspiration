# Scientific release checklist

- Verify that all five files under `reference_implementation/` match
  `SOURCE_MANIFEST.json`.
- Install the package in a new environment and run `python -m pytest`.
- Run the complete study with `--force`; do not reuse filename-only caches.
- Compare the new summary with `reference_results/study_summary.csv` and record
  the numerical tolerance used.
- Record Python, NumPy, SciPy, operating-system, Git commit, and tag identifiers.
- Obtain coauthor approval of the model scope, interpretation, figures, and
  manuscript before creating an archival release.
- Do not call the source-reference reproduction an independent validation.
- Do not describe model outputs as clinical efficacy or safety evidence.
- Add `CITATION.cff` only after all contributors approve authorship and the
  preferred software citation.
