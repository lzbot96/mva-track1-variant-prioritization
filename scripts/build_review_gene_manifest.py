"""Build a review gene manifest, retaining ties at each phenotype cutoff.

The reviewer panel and generated manifest are controlled local inputs/outputs.
A fixed alphabetical slice must not decide which equally scored genes are reviewed.
"""
import csv
import gzip
import json
import math
import re


def select_boundary_ties(rows, limit=20, tolerance=1e-12):
    """Retain the top limit plus ties on both original ranking metrics."""
    if not isinstance(limit, int) or limit < 1:
        raise ValueError('limit must be a positive integer')
    scored = []
    for row in rows:
        score = float(row['mean_resnik'])
        exact = int(row['exact_query_terms'])
        if not math.isfinite(score) or exact < 0:
            raise ValueError('invalid phenotype score or exact-query count')
        scored.append((score, exact, row['gene']))
    scored.sort(key=lambda x: (-x[0], -x[1], x[2]))
    if not scored:
        return []
    cutoff_score, cutoff_exact, _ = scored[min(limit, len(scored)) - 1]
    # Preserve the raw-sort prefix: tolerance-based tie expansion must never
    # remove a gene already admitted by the original top-N selection.
    original_top = {gene for _, _, gene in scored[:limit]}
    return [gene for score, exact, gene in scored
            if gene in original_top
            or score > cutoff_score + tolerance
            or (abs(score - cutoff_score) <= tolerance and exact >= cutoff_exact)]


def main():
    from project_paths import project_root, required_env
    root = project_root()
    panel = json.loads((root / required_env('MVA_PANEL_JSON')).read_text())['genes']
    genes = set(panel)
    selections = {}
    for label in ['all_eight', 'proband_seven', 'without_cancer']:
        with (root / f'results/phenotype/ranking_{label}.tsv').open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        selected = select_boundary_ties(rows, 20)
        genes.update(selected)
        selections[label] = {'requested_top_n': 20, 'selected_with_boundary_ties': len(selected)}
    regions = {}
    with gzip.open(root / 'data/references/ucsc/hg38.refGene.gtf.gz', 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue
            match = re.search(r'gene_id "([^"]+)"', fields[8])
            if not match:
                continue
            gene = match[1]
            if gene not in genes or fields[0] not in ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']:
                continue
            key = (gene, fields[0][3:])
            start, end = int(fields[3]), int(fields[4])
            previous = regions.get(key, [start, end])
            regions[key] = [min(start, previous[0]), max(end, previous[1])]
    with (root / 'work/review_gene_regions.tsv').open('w') as out:
        for (gene, chrom), (start, end) in sorted(regions.items()):
            out.write(f'{chrom}\t{start}\t{end}\t{gene}\n')
    (root / 'results/candidates/review_gene_manifest.json').write_text(json.dumps({
        'genes': sorted(genes),
        'missing_genes': sorted(genes - {g for g, c in regions}),
        'selection': 'Top 20 of each phenotype ranking, including boundary ties in mean Resnik score and exact-query count, plus the separately supplied reviewer panel.',
        'phenotype_selections': selections,
        'score_tolerance': 1e-12,
    }, indent=2) + '\n')
    print('Review gene selection complete; controlled outputs saved locally.')


if __name__ == '__main__':
    main()
