# Code verification results — 17 October 2026

Passed locally:

- `python3 tests/self_test.py`: all staged Python files parsed without execution; all shell files passed `bash -n`; C++ scanner compiled and passed synthetic forward/reverse, ambiguous-base, Q20/Q30, per-read repeated-hit, invalid-pattern, and truncated-FASTQ tests. Empty-stage preflight and missing-root guard behaved as expected.
- `python -B tests/reconstruction_test.py` using the already-installed original analysis environment: synthetic 31-base anchor generation, reverse complement, reference-mismatch rejection, top-20/panel union, and region merging passed. Only synthetic temporary files were written; original environment and project files were not modified.
- Read-only `preflight.py --stage annotation` against the original project reported all listed prerequisite paths present and required command-line tools discoverable. System Python lacked the four analysis packages; use the configured analysis environment. This was a presence check, not a rerun or integrity check.
- SHA-256 comparison confirmed all 15 original saved scripts remained unchanged.
- Allowlisted staging copied only 15 original code files and two software-version text receipts, then added code/documentation/tests/provenance. No raw data, genotype tables, original configs, phenotype inputs, clinical document, broad result files, caches, logs, external helpers or credentials were copied.
- Text scan found no personal home-directory paths, original sample filenames, HPO identifier literals, Hugging Face token patterns or AWS access-key patterns. This scan supplements the allowlist; it is not a general secret-detection guarantee.

These verification runs used existing dependencies and synthetic inputs. No network calls or large genomic analyses were run. Full end-to-end reproduction and scientific validation remain unperformed. See WORKFLOW.md for the required curated inputs and reference constraints. The code retains gene-specific scientific logic.

Final full preflight was also run with the existing analysis interpreter and no private configuration values set: every listed static prerequisite path was present, all five checked command-line tools were available, and the four direct analysis package versions matched requirements.txt. Exit 2 correctly reported seven distinct unset configuration variables across clinical extraction/curation, QC, normalization, phenotype, review-gene selection and anchor generation. The analysis environment lacks optional python-docx; the clinical synthetic test passed in the existing bundled interpreter with python-docx 1.2.0. These reviewer-controlled settings and the appropriate dependency environment must be supplied before reproduction; preflight does not inspect private content.

After incorporating the earlier records, all three test suites passed again. `tests/clinical_test.py` additionally verified synthetic DOCX table extraction and reviewer-selected family-history assignment using the already available python-docx interpreter. New network stages were syntax-checked only; reconstructed gnomAD stdlib transport was not contacted. FastQC/ClinVar/authenticated metadata historical commands are now present; original private inputs and mutable reference identity remain external.

## Independent critical review

The isolated revision additionally passes `tests/integrity_test.py`: same-size modified input, mismatched receipt revision, incomplete download verification, unwritable scanner outputs, changed anchor names, tampered count outputs and a failed rerun with a stale prior receipt. Existing scanner/syntax, reconstruction and clinical synthetic suites also passed independently. Source-script hashes were updated to describe this revised code, preserving the original hashes.

The old implementation's saved-checksum trust, unchecked output streams and stale-evidence acceptance were real defects. Fixes improve future execution but do not retroactively validate historical data or outputs. No real scan, whole-dataset hash, live API retrieval or clinical analysis was run. These results describe synthetic code checks, not a reanalysis of case data. See DEPENDENCIES_AND_RIGHTS.md for reuse and rights limits.
