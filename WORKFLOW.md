# Analysis stages and required inputs

Use a separate controlled run directory and inputs that you are authorized to access. Set CODE to the absolute directory containing this README, MVA_PROJECT_ROOT to an existing writable controlled run directory, and PYTHON to an interpreter with requirements.txt installed. These commands describe the configured stages. A complete replay of the workflow has not been performed.

```sh
export CODE='/absolute/path/to/code'
export MVA_PROJECT_ROOT='/absolute/path/to/controlled-run'
export PYTHON=python3
mkdir -p "$MVA_PROJECT_ROOT"/data/raw "$MVA_PROJECT_ROOT"/data/references \
  "$MVA_PROJECT_ROOT"/work "$MVA_PROJECT_ROOT"/logs \
  "$MVA_PROJECT_ROOT"/results/qc "$MVA_PROJECT_ROOT"/results/phenotype \
  "$MVA_PROJECT_ROOT"/results/candidates "$MVA_PROJECT_ROOT"/results/read_validation
```

Read-only preflight (no content reads, writes, downloads or analysis):

```sh
python3 "$CODE/scripts/preflight.py" --root "$MVA_PROJECT_ROOT"
# Optional stage selection: --stage annotation (see --help)
```

Exit 2 means stage inputs/configuration are absent. Tool/package availability is reported separately, not included in stage readiness. Presence is not integrity, schema, or biological validation. Use authorized archived reference receipts to check hashes before reproduction.

These stage commands follow saved scripts and successful recovered historical command records. Do not execute a stage until its listed inputs exist. No automatic orchestrator masks missing inputs.

| Stage | Required inputs not bundled | Invocation | Main output |
|---|---|---|---|
| Public references (network) | No patient input; compare downloads to archived hashes | `bash "$CODE/scripts/acquire_public_references.sh"` | UCSC, HPO, Ensembl GFF3 |
| Download (network; gated access required) | `sources/dataset_manifest.json`, official account access | `"$PYTHON" "$CODE/scripts/download_data.py"` | `data/raw/`, download receipt |
| Authenticated metadata (network) | original pinned public manifest and locally configured gated access | `"$PYTHON" "$CODE/scripts/authenticated_manifest.py"` | authenticated file metadata |
| Verify | `data/authenticated_manifest.json`, downloaded files/receipt | `"$PYTHON" "$CODE/scripts/verify_downloads.py"` | integrity report |
| FASTQ statistics + FastQC | raw FASTQs; FastQC/SeqKit | `bash "$CODE/scripts/fastq_qc.sh"` | fastq_stats.tsv, FastQC outputs |
| VCF QC | Set `MVA_INPUT_VCF` to relative VCF path | `"$PYTHON" "$CODE/scripts/input_assessment.py"` | QC JSON |
| Pair checks | Set `MVA_FASTQ_TEMPLATE` to relative pattern with `{lane:03d}` and `{mate}`; four lanes | `"$PYTHON" "$CODE/scripts/check_fastq_pairs.py"` | lane/summary JSON |
| Reference conversion | `data/references/ucsc/hg38.2bit` of expected size | `"$PYTHON" "$CODE/scripts/prepare_ucsc_reference.py"` | primary FASTA, provenance |
| Index FASTA | converted reference | `samtools faidx "$MVA_PROJECT_ROOT/data/references/ucsc/GRCh38.primary.no_chr.fa"` | FASTA index |
| Normalize source | indexed source VCF, FASTA/index; `MVA_INPUT_VCF` | `bash "$CODE/scripts/normalize_primary.sh"` | normalized source VCF/index |
| Dated ClinVar (network) | remote availability | `bash "$CODE/scripts/acquire_clinvar.sh"` then `"$PYTHON" "$CODE/scripts/verify_clinvar.py"` | dated source, index, checksum receipt |
| Normalize ClinVar | dated ClinVar source VCF, FASTA/index | `bash "$CODE/scripts/normalize_clinvar.sh"` | normalized ClinVar/index |
| Clinical extraction | private DOCX; set `MVA_CLINICAL_DOCX`; python-docx | `"$PYTHON" "$CODE/scripts/extract_clinical.py"` | private extraction JSON |
| Clinical table curation | extraction; set reviewer-chosen `MVA_FAMILY_HISTORY_HPO` | `"$PYTHON" "$CODE/scripts/curate_hpo.py"` | curated HPO JSON/TSV |
| Phenotype baseline | `data/references/hpo/hp.obo`, `genes_to_phenotype.txt`, controlled `results/phenotype/curated_hpo.json`; set `MVA_EXCLUDED_HPO` | `"$PYTHON" "$CODE/scripts/phenotype_rank.py"` | semantic rankings |
| Select review genes | three ranking TSVs, refGene GTF; set `MVA_PANEL_JSON` to private JSON with `genes` | `"$PYTHON" "$CODE/scripts/build_review_gene_manifest.py"` | review manifest/regions |
| Consequence + ClinVar | normalized/indexed source VCF, FASTA/index, Ensembl release 116 GFF3, normalized/indexed ClinVar 2026-09-28 | `bash "$CODE/scripts/annotate_primary.sh"` | annotated primary VCF |
| Public frequencies (network) | controlled review gene manifest; optional `GNOMAD_HELPER` | `"$PYTHON" "$CODE/scripts/fetch_gnomad_review_genes.py"` | gene-level reference JSON |
| Candidate summary | annotated VCF, gene manifest, rankings, separately acquired public frequency JSON | `"$PYTHON" "$CODE/scripts/summarize_candidates.py"` | candidate tables and summary |
| Gene-specific CDS check | saved public Ensembl transcript JSON, FASTA/index, normalized VCF/index | `"$PYTHON" "$CODE/scripts/reconstruct_bub1b.py"` | reconstruction JSON |
| Compile anchor scanner | compiler and zlib | `c++ -std=c++17 -O3 "$CODE/scripts/scan_allele_kmers.cpp" -lz -o "$MVA_PROJECT_ROOT/work/scan_allele_kmers"` | local executable |
| Build anchors | FASTA/index, reviewer SNV TSV; set `MVA_ANCHOR_CANDIDATES` | `"$PYTHON" "$CODE/scripts/build_allele_kmers.py"` | private 31-base anchor table |
| Anchor scan | private `work/allele_kmers.tsv`, raw FASTQs, compiled scanner | `"$PYTHON" "$CODE/scripts/run_kmer_validation.py"` | matched reads, counts, completion marker |
| Anchor aggregation | scan outputs and `results/qc/fastq_stats.tsv` | `"$PYTHON" "$CODE/scripts/aggregate_read_validation.py"` | read-support summary |
| Reference specificity | private anchor table, hg38.2bit | `"$PYTHON" "$CODE/scripts/kmer_reference_specificity.py"` | reference-match counts |
| Draft validation | private `results/track1_submission_draft.csv`, normalized VCF/index, FASTA/index | `"$PYTHON" "$CODE/scripts/validate_submission.py"` | local schema/reference/genotype receipt |

## Missing steps and constrained reproducibility

Reconstructed clinical extraction and first-table transcription generate private HPO inputs after explicit reviewer configuration. The code does not choose the reviewer panel, subject assignment, sensitivity exclusion, final CSV, confidence or narrative conclusions. These reviewed decisions/results are deliberately excluded. Reconstructed gene-selection and anchor-generation scripts now produce the intermediate manifest and anchor table from private reviewer inputs, using recovered historical logic. Phenotype input must provide `all_challenge_hpo` and `proband_only_hpo`; gene manifest must provide `genes`. Anchor-generator TSV input requires headers `name`, `chrom`, `pos`, `ref`, `alt`; positions are 1-based and contigs match the no-chr FASTA. It accepts SNVs only and verifies reference alleles. The scanner anchor input is two whitespace-delimited columns: unique name and a 31-base A/C/G/T sequence. Both orientations must be supplied explicitly; scanner does not reverse-complement automatically. Names ending in an orientation suffix are grouped by the aggregator's last underscore. No synthetic example represents the actual case.

Recovered successful records supply exact primary normalization (`-c e -m -any --old-rec-tag ORIG`, primary nuclear plus mitochondrial contigs), ClinVar normalization (`-c x -m -any`, primary nuclear only), index commands, public UCSC/HPO/GFF3 acquisition, top-20 union review-gene selection and 31-base SNV anchors. Newly reconstructed scripts preserve these flags. ClinVar reference mismatches are excluded by `-c x`; inspect its log. FASTA indexing is recovered as `samtools faidx`.

Earlier recovered records additionally supply exact SeqKit statistics/FastQC commands, dated ClinVar download/upstream MD5 verification, authenticated metadata retrieval, and clinical extraction. Curation was the first Python block of a shell command whose later unrelated bcftools step failed; pre-existing curated outputs corroborate that first block, but the whole command is not labeled successful. The reconstructed transcription explicitly validates the expected table structure and requires the family's HPO designation as reviewer input.

Remaining reproducibility constraints are private reviewed input values, the original pinned dataset manifest, archived reference snapshots/fingerprints, and public API availability. Mutable HPO/UCSC URLs cannot guarantee historical byte identity. The aggregator expects tab-separated `file` and `num_seqs` columns produced by the recovered SeqKit command. The original gnomAD helper is optional. A reconstructed stdlib transport in gnomad_public.py preserves the same public gene-level query and expected raw JSON shape; it has not been network-tested. The initial gene-frequency file is separately acquired with that CLI using the reviewer-selected public gene and output path required by summarize_candidates.py. The loop preserves historical exclusions; it does not fabricate missing records.

The converted UCSC primary reference is complete according to the original project report. The separate redundant Ensembl reference and Exomiser reference downloads were paused incomplete. The GFF3 and normalized ClinVar paths above are inputs to the completed local annotation, not evidence that all proposed references were completed. Reference manifests/hashes need checking in the controlled workspace before rerun; code-package hashes cover only code.

The saved validator checks at most ten rows, primary findings, positive EPCR/order, reference matching, source heterozygosity and duplicate alleles. It does **not** invoke the official evaluator, assess scientific correctness, calibrate confidence or establish phase. Its assertion-based schema checks are narrower than a full official form validation. Any independent official-evaluator acceptance must be recorded separately.

Exact anchors are sequence evidence, not alignment, deduplicated molecular coverage, phasing or clinical confirmation. The scanner counts each pattern once per read and uses the first hit's central-base quality if a pattern repeats. It neither performs mapping nor pairing. The aggregator's read-name grouping is not molecular deduplication. Missing frequency remains unknown, not evidence of rarity. Local consequence annotation does not resolve phase. Original full-alignment, structural/mosaic, noncoding and functional gaps remain.

Input schemas and all environment settings are documented in INPUT_SCHEMAS.md. The workflow table is stage-oriented; satisfy dependencies before execution (in particular ClinVar acquisition/verification precede its normalization).

## Integrity and scan completion checks

The verifier now reads and hashes current file bytes against authenticated upstream metadata and checks receipt repository/revision identity. An incomplete check returns nonzero and clears prior success. Running it on the real dataset is a substantial I/O operation; it was tested only with tiny synthetic files in this review. September integrity claims in the report remain historical receipts, not a new whole-dataset recheck.

The scanner now fails for output open/write/close errors. The runner invalidates prior completion, stages new artifacts and atomically publishes a version-1 receipt only after success. The aggregator verifies anchor SHA-256, exact pattern identities, FASTQ size/mtime inventory, count totals and output hashes. Metadata inventory is not a replacement for raw-read content checksums. Historical completion receipts are not compatible with this stronger evidence binding and need regeneration in a future authorized run; historical scientific outputs were not rerun. Keep one scanner invocation per run directory. The new `cds` preflight stage checks the documented public transcript input and indexed reference/VCF prerequisites.
