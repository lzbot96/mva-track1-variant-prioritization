"""Synthetic fixtures for reconstructed anchors and review-gene union. Needs pysam/Bio."""
import csv, gzip, json, os, subprocess, sys, tempfile
from pathlib import Path
import pysam
R=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='mva-reconstruction-test-') as td:
    t=Path(td)
    for d in ['data/references/ucsc','work','results/phenotype','results/candidates']: (t/d).mkdir(parents=True)
    fa=t/'data/references/ucsc/GRCh38.primary.no_chr.fa'
    seq='ACGT'*30
    fa.write_text('>synthetic\n'+seq+'\n');pysam.faidx(str(fa))
    candidates=t/'candidates.tsv';candidates.write_text('name\tchrom\tpos\tref\talt\nexample\tsynthetic\t40\tT\tA\n')
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','MVA_PROJECT_ROOT':td,'MVA_ANCHOR_CANDIDATES':'candidates.tsv','MVA_PANEL_JSON':'panel.json'}
    def run(script): return subprocess.run([sys.executable,str(R/'scripts'/script)],env=env,capture_output=True,text=True)
    result=run('build_allele_kmers.py');assert result.returncode==0,result.stderr
    lines=[l.split() for l in (t/'work/allele_kmers.tsv').read_text().splitlines()]
    assert len(lines)==4 and all(len(v)==31 for _,v in lines)
    assert lines[0][1][15]=='T' and lines[2][1][15]=='A'
    assert lines[1][1]==lines[0][1].translate(str.maketrans('ACGT','TGCA'))[::-1]
    candidates.write_text('name\tchrom\tpos\tref\talt\nexample\tsynthetic\t40\tG\tA\n')
    assert run('build_allele_kmers.py').returncode!=0
    (t/'panel.json').write_text(json.dumps({'genes':['PANEL_ONLY']}))
    for label in ['all_eight','proband_seven','without_cancer']:
        (t/f'results/phenotype/ranking_{label}.tsv').write_text('gene\n'+''.join(f'G{i}\n' for i in range(21)))
    with gzip.open(t/'data/references/ucsc/hg38.refGene.gtf.gz','wt') as f:
        f.write('# synthetic\nchr1\ttest\ttranscript\t10\t20\t.\t+\t.\tgene_id "G0";\n')
        f.write('chr1\ttest\ttranscript\t8\t25\t.\t+\t.\tgene_id "G0";\n')
    result=run('build_review_gene_manifest.py');assert result.returncode==0,result.stderr
    manifest=json.loads((t/'results/candidates/review_gene_manifest.json').read_text())
    assert len(manifest['genes'])==21 and 'G20' not in manifest['genes'] and 'PANEL_ONLY' in manifest['genes']
    assert (t/'work/review_gene_regions.tsv').read_text()=='1\t8\t25\tG0\n'
print('PASS: synthetic anchor orientation/reference/mismatch and top-20 panel union/region merging.')
