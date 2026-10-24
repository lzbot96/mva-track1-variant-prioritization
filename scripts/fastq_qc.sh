#!/bin/bash
# Reconstructed from successful early records 11 and 13. Potentially long-running.
set -euo pipefail
cd "${MVA_PROJECT_ROOT:?Set controlled run directory}"
mkdir -p results/qc/fastqc logs
seqkit stats --all --tabular --threads 4 data/raw/*.fastq.gz > results/qc/fastq_stats.tsv 2> logs/seqkit_stats.log
fastqc --version
seqkit version
fastqc --threads 4 --memory 1024 --outdir results/qc/fastqc data/raw/*.fastq.gz > logs/fastqc.log 2>&1
