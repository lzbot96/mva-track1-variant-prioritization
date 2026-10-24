"""Reconstruct saved top-20 union selection; panel supplied privately by reviewer."""
import csv, gzip, json, re
from project_paths import project_root, required_env
R = project_root()
panel = json.loads((R/required_env('MVA_PANEL_JSON')).read_text())['genes']
genes = set(panel)
for label in ['all_eight', 'proband_seven', 'without_cancer']:
    with (R/f'results/phenotype/ranking_{label}.tsv').open() as f:
        genes.update(row['gene'] for row in list(csv.DictReader(f, delimiter='\t'))[:20])
regions = {}
with gzip.open(R/'data/references/ucsc/hg38.refGene.gtf.gz', 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        c = line.split('\t')
        if len(c) < 9 or c[2] != 'transcript': continue
        match = re.search(r'gene_id "([^"]+)"', c[8])
        if not match: continue
        g = match[1]
        if g not in genes or c[0] not in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']: continue
        key=(g,c[0][3:]); start=int(c[3]); end=int(c[4])
        previous=regions.get(key,[start,end])
        regions[key]=[min(start,previous[0]),max(end,previous[1])]
with (R/'work/review_gene_regions.tsv').open('w') as out:
    for (g,c),(s,e) in sorted(regions.items()): out.write(f'{c}\t{s}\t{e}\t{g}\n')
(R/'results/candidates/review_gene_manifest.json').write_text(json.dumps({
    'genes':sorted(genes), 'missing_genes':sorted(genes-{g for g,c in regions}),
    'selection':'All genes in top 20 of three phenotype-only ranks plus independently sourced reviewer panel; old-symbol aliases included.'
},indent=2)+'\n')
print('Review gene selection complete; controlled outputs saved locally.')
