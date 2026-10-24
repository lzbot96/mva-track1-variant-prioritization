# Small-variant prioritization and evidence checks for the MVA Hackathon

This repository contains analysis code for Track 1 of *Rare Disease, Real Kid: The 2026 MVA Hackathon*. The work was directed by zoul6. The stages cover input integrity and quality checks, variant normalization and annotation, phenotype-based gene ranking, public-reference review, coding-sequence reconstruction and exact read-anchor searches. The associated submission presents a provisional research hypothesis, not a clinical diagnosis.

**These are partially reproducible, configured stages, not a complete one-command pipeline.** Authorized case inputs, reviewer decisions, archived reference snapshots and compatible installed tools are required. No controlled inputs, populated private configuration, candidate submission files or analysis outputs are included here. Gene/transcript-specific logic and four-lane assumptions remain explicit.

## Run and inspect

Start with [WORKFLOW.md](WORKFLOW.md) for stage commands and their prerequisites, then [INPUT_SCHEMAS.md](INPUT_SCHEMAS.md) for configuration. Set `MVA_PROJECT_ROOT` to a separate controlled run directory. Keep that directory private. The scripts create controlled outputs there. Do not run Python with `-O`; existing checks use assertions.

The workflow does not choose a final candidate, confidence estimate, phase or diagnosis. Exact anchors provide raw-sequence support; they do not establish aligned, duplicate-aware coverage, phase or causality. The local CSV validator checks format, reference and source-call consistency; it does not score biological correctness or calibrate confidence.

## Verify the code

Run these synthetic checks from this directory with the indicated dependencies available:

```sh
python3 tests/self_test.py
python3 -B tests/integrity_test.py
python3 -B tests/reconstruction_test.py
python3 -B tests/clinical_test.py
python3 -B tests/genotype_carrier_test.py
python3 -B tests/phenotype_sensitivity_test.py
python3 -B tests/phase_recheck_test.py
python3 -B tests/population_reference_test.py
```

The first suite parses Python, checks shell syntax, builds the C++ scanner in a temporary directory and exercises synthetic reads. Integrity tests use small temporary files. Reconstruction tests require pysam and Biopython. Clinical tests require python-docx. These checks do not access case data or validate clinical findings. [VALIDATION.md](VALIDATION.md) records the dated results and their limits.

## Environment and provenance

The analysis environment inspected during preparation used Python 3.9.6. The clinical synthetic test used Python 3.12.14 with python-docx 1.2.0. These observations are not a complete historical lockfile or proof of whole-pipeline compatibility. `requirements.txt` records four direct analysis dependencies; `requirements-clinical.txt` records the optional clinical dependency. External tools include bcftools/samtools, SeqKit, FastQC, a C++17 compiler and zlib. The optional selected-read regional check uses BWA (0.7.19-r1273 observed); it accepts a PATH executable or --bwa override. Historical compiler/zlib versions remain unavailable. Exomiser downloads were incomplete and Exomiser was not run; VEP was not run.

[DEPENDENCIES_AND_RIGHTS.md](DEPENDENCIES_AND_RIGHTS.md) records dependency versions, license observations and reuse limits. The `provenance/` directory records software versions, original script hashes and documented transformations. Public references are not bundled. The historical acquisition command/date for the public transcript input is unavailable; its schema, inspected identity and fingerprint are documented.

The code includes fixes for three verified execution defects: trusting stale download checksums, ignoring scanner output failures, and accepting stale scan evidence. Current bytes are hashed; output failures reject a scan; a version-1 completion receipt binds anchors, input metadata and output hashes. These fixes do not retroactively validate historical outputs. Old scan receipts require regeneration by a future authorized run. FASTQ size/mtime metadata is not a cryptographic content checksum. Run only one scanner invocation per controlled run directory.

## Attribution and license

Copyright 2026 zoul6. Original contributions are licensed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/); see [LICENSE](LICENSE). Third-party software, references and controlled source data retain their own terms. No dependencies or binaries are bundled.

[AI assistance](AI_DISCLOSURE.md) · [Dataset reference and acknowledgment](DATASET_REFERENCE.md) · [Rights notice](NOTICE.md)
