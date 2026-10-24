# Dependencies and rights

Copyright 2026 zoul6. The original contributions are licensed under Creative Commons Attribution 4.0 International; see [LICENSE](LICENSE). This license does not grant rights to controlled source data or third-party materials. Source hashes establish provenance, not ownership of external materials. Dependency terms remain separate; see [NOTICE.md](NOTICE.md).

## Direct dependencies inspected locally

| Dependency | Installed version | License label/evidence in the installed distribution |
|---|---|---|
| huggingface_hub | 1.8.0 | Apache-2.0 metadata; LICENSE listed |
| pysam | 0.24.1 | MIT License metadata; COPYING and AUTHORS listed |
| Biopython | 1.85 | Installed LICENSE says Biopython License Agreement by default; some explicitly marked files also offer BSD 3-Clause. Do not label the entire distribution BSD-only. |
| py2bit | 1.0.1 | MIT metadata; LICENSE.txt listed |
| python-docx, optional clinical extraction | 1.2.0 | MIT metadata; LICENSE listed |

These are observations of installed package metadata, not a license for this repository. Dependency implementations and binary distributions are not included in this code archive. The scripts import installed modules or call separately installed tools. A future redistribution of dependencies, containers or binaries requires review of their own notices and bundled components; transitive packages and bundled libraries are not exhaustively audited here.

The analysis environment observed during this review uses Python 3.9.6. The separate clinical synthetic test uses the existing bundled Python 3.12.14 environment. These observations are not a historical lockfile or a guarantee that every stage works in either environment. requirements.txt captures four direct analysis versions; requirements-clinical.txt records the current optional clinical dependency. The older download-environment receipt is a historical package list, not a tested universal installation specification.

External requirements include bcftools/samtools, SeqKit, FastQC, a C++17 compiler and zlib. Reported analysis versions are retained in provenance/analysis_software_versions.txt. Their binaries and sources are not bundled; their licenses and transitive components have not been comprehensively audited. Public reference snapshots and their terms are also not bundled. Use authorized snapshots and source attribution, not an assumption that publicly downloadable means unrestricted redistribution.

## Limits of reuse

This remains a collection of configured analysis stages, not a one-command generic disease pipeline. It preserves gene/transcript-specific logic, four-lane assumptions and reviewer choices. A new authorized run needs external private inputs, pinned manifests/references, the documented public transcript input and compatible installed tools. No controlled case data or populated private configuration is distributed here. No full pipeline replay or live API validation was performed during packaging or critical review.
