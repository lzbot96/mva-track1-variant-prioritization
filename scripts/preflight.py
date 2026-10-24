"""Read-only stage readiness inventory. Never executes a stage or reads content."""
import argparse, importlib.metadata, json, os, shutil
from pathlib import Path
STAGES = {
 'clinical_extract': ['results/phenotype'],
 'clinical_curate': ['results/phenotype/clinical_document_extracted.json'],
 'qc': ['results/qc', 'data/raw'],
 'reference': ['data/references/ucsc/hg38.2bit'],
 'normalize': ['data/references/ucsc/GRCh38.primary.no_chr.fa','data/references/ucsc/GRCh38.primary.no_chr.fa.fai'],
 'phenotype': ['data/references/hpo/hp.obo','data/references/hpo/genes_to_phenotype.txt','results/phenotype/curated_hpo.json'],
 'review_genes': ['results/phenotype/ranking_all_eight.tsv','results/phenotype/ranking_proband_seven.tsv','results/phenotype/ranking_without_cancer.tsv','data/references/ucsc/hg38.refGene.gtf.gz'],
 'clinvar': ['data/references/clinvar/clinvar_20260928.vcf.gz','data/references/ucsc/GRCh38.primary.no_chr.fa','data/references/ucsc/GRCh38.primary.no_chr.fa.fai'],
 'cds': ['data/references/bub1b_ensembl.json','data/references/ucsc/GRCh38.primary.no_chr.fa','data/references/ucsc/GRCh38.primary.no_chr.fa.fai','work/primary.normalized.vcf.gz','work/primary.normalized.vcf.gz.csi'],
 'annotation': ['data/references/ucsc/GRCh38.primary.no_chr.fa','data/references/ucsc/GRCh38.primary.no_chr.fa.fai','work/primary.normalized.vcf.gz','work/primary.normalized.vcf.gz.csi','data/references/ensembl/Homo_sapiens.GRCh38.116.gff3.gz','data/references/clinvar/clinvar_20260928.primary.normalized.vcf.gz','data/references/clinvar/clinvar_20260928.primary.normalized.vcf.gz.csi'],
 'candidates': ['work/primary.csq.clinvar.vcf.gz','results/candidates/review_gene_manifest.json','results/phenotype/semantic_ranking.json','data/references/bub1b_gnomad.json','data/references/gnomad_review_genes'],
 'anchors': ['data/references/ucsc/GRCh38.primary.no_chr.fa','data/references/ucsc/GRCh38.primary.no_chr.fa.fai'],
 'scan': ['work/allele_kmers.tsv','work/scan_allele_kmers','data/raw'],
 'aggregate': ['results/qc/fastq_stats.tsv','results/read_validation/complete.json','work/allele_kmers.tsv'],
 'validate': ['data/references/ucsc/GRCh38.primary.no_chr.fa','work/primary.normalized.vcf.gz.csi','results/track1_submission_draft.csv','work/primary.normalized.vcf.gz','data/references/ucsc/GRCh38.primary.no_chr.fa.fai'],
}
ENVS={'clinical_extract':['MVA_CLINICAL_DOCX'],'clinical_curate':['MVA_FAMILY_HISTORY_HPO'],'qc':['MVA_INPUT_VCF','MVA_FASTQ_TEMPLATE'],'normalize':['MVA_INPUT_VCF'],'phenotype':['MVA_EXCLUDED_HPO'],'review_genes':['MVA_PANEL_JSON'],'anchors':['MVA_ANCHOR_CANDIDATES']}
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True);p.add_argument('--stage',choices=list(STAGES));a=p.parse_args()
r=Path(a.root).expanduser().resolve()
report={'root_exists':r.is_dir(),'stages':{},'tools':{t:shutil.which(t) is not None for t in ['bcftools','samtools','seqkit','fastqc','c++']},'python_packages':{}}
for name in [a.stage] if a.stage else STAGES:
    missing=[x for x in STAGES[name] if not (r/x).exists()]
    envs=[e for e in ENVS.get(name,[]) if not os.environ.get(e)]
    for e in ENVS.get(name,[]):
        value=os.environ.get(e)
        if value and e in ['MVA_INPUT_VCF','MVA_PANEL_JSON','MVA_ANCHOR_CANDIDATES','MVA_CLINICAL_DOCX'] and not (r/value).is_file(): missing.append(e+' configured file')
    if name=='scan' and not list((r/'data/raw').glob('*.fastq.gz')): missing.append('FASTQ files')
    report['stages'][name]={'required_paths_present':not missing,'missing':missing,'unset_configuration':envs}
for pkg in ['pysam','py2bit','biopython','huggingface_hub','python-docx']:
    try: report['python_packages'][pkg]=importlib.metadata.version(pkg)
    except importlib.metadata.PackageNotFoundError: report['python_packages'][pkg]=None
print(json.dumps(report,indent=2))
# Inventory is not a correctness test; missing prerequisites return nonzero.
raise SystemExit(0 if r.is_dir() and all(x['required_paths_present'] and not x['unset_configuration'] for x in report['stages'].values()) else 2)
