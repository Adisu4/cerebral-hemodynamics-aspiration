# Source and result provenance

The result-producing source was recovered and audited outside Git. Every one of
the 91 run reports and both aggregate JSON records that contain source hashes
identify the same five files.

| Archived source file | SHA-256 |
|---|---|
| `publication_parameters.py` | `d008b39646a71c93d4a15c76c1eb26a84bad3f3cb9cebbeeae6f2bed56ef0b88` |
| `publication_model.py` | `ab027031e2aca97252ba5d9f3f67daf3f17c55aeed0c0ab76528d733dc181f7f` |
| `run_publication.py` | `033f1a77ad05146d0eadfc63d22807d9fa0bfba379d345da35cb717956a12140` |
| `run_publication_study.py` | `fcbcd3112644f56505b5d6f3eb0f78e49762d84166e2c364d3bb0e3ff341a7a0` |
| `test_publication_model.py` | `995add083ec484051582f2e7e6ff26ec0d1f596b9029ad459b64d43b32959b32` |

The complete immutable review archive has SHA-256
`392ebb87b6eb0c9b4df770ffabdad053b35b5ce1bf3c8231770cdae9a92b2d5d`.
It is a review artifact, not a clinical or publication certification.

`reference_implementation/` preserves the five source files exactly. Their
historical filenames are retained because changing them would invalidate the
recorded hashes. The installable package uses scientific module names in a
conventional `src/` layout. Regression tests compare parameters and equation
evaluations with the archived source and inspect every saved report.

The previous repository history is not provenance for these results. Neither
the earlier main branch nor the pre-existing `v1.0.0` tag contains this audited
source set.
