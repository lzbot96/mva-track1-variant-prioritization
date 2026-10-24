"""Audit cached phenotype rankings with leave-one-term-out and tie-aware ranks.

Requires the cached per-term score JSON and the matching public HPO gene annotation.
This is robustness analysis, not causal-probability calibration. Output includes
query terms and gene rankings and must remain in the controlled workspace.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

TOLERANCE = 1e-12


def rank_rows(cached_rows, terms, exact_annotations):
    if not terms or len(set(terms)) != len(terms):
        raise ValueError('query terms must be nonempty and unique')
    rows = []
    for item in cached_rows:
        scores = [float(item['per_term'][term]) for term in terms]
        if any(not math.isfinite(score) for score in scores):
            raise ValueError('cached scores must be finite')
        rows.append({'gene': item['gene'], 'mean_resnik': sum(scores) / len(scores),
                     'exact_query_terms': sum(term in exact_annotations.get(item['gene'], set()) for term in terms)})
    return sorted(rows, key=lambda item: (-item['mean_resnik'], -item['exact_query_terms'], item['gene']))


def rank_interval(rows, gene, tolerance=TOLERANCE):
    item = next(row for row in rows if row['gene'] == gene)
    higher = sum(row['mean_resnik'] > item['mean_resnik'] + tolerance for row in rows)
    tied = [row for row in rows if abs(row['mean_resnik'] - item['mean_resnik']) <= tolerance]
    exact_higher = sum(row['exact_query_terms'] > item['exact_query_terms'] for row in tied)
    exact_ties = sum(row['exact_query_terms'] == item['exact_query_terms'] for row in tied)
    return {'gene': gene, 'score': item['mean_resnik'],
            'mean_rank_min': higher + 1, 'mean_rank_max': higher + len(tied),
            'composite_rank_min': higher + exact_higher + 1,
            'composite_rank_max': higher + exact_higher + exact_ties}


def cutoff_summary(rows, limit, tolerance=TOLERANCE):
    if limit < 1 or not rows:
        raise ValueError('positive limit and nonempty ranking required')
    cutoff = rows[min(limit, len(rows)) - 1]
    # A numerical tolerance may expand the boundary but cannot shrink top-N.
    original_top = {row['gene'] for row in rows[:limit]}
    selected = [row['gene'] for row in rows
                if row['gene'] in original_top
                or row['mean_resnik'] > cutoff['mean_resnik'] + tolerance
                or (abs(row['mean_resnik'] - cutoff['mean_resnik']) <= tolerance
                    and row['exact_query_terms'] >= cutoff['exact_query_terms'])]
    return {'requested_top_n': limit, 'boundary_score': cutoff['mean_resnik'],
            'boundary_exact_count': cutoff['exact_query_terms'], 'selected_genes': selected,
            'count_including_boundary_ties': len(selected)}


def analyze(cache, annotations, baseline, tracked):
    scenarios = cache['rankings']
    source = scenarios[baseline]
    terms = list(source[0]['per_term'])
    if len(terms) < 2:
        raise ValueError('leave-one-out needs at least two query terms')
    rebuilt = rank_rows(source, terms, annotations)
    old = {row['gene']: row for row in source}
    if len(old) != len(source):
        raise ValueError('duplicate genes in cached baseline')
    for row in rebuilt:
        previous = old[row['gene']]
        if abs(row['mean_resnik'] - previous['mean_resnik']) > TOLERANCE or row['exact_query_terms'] != previous['exact_query_terms']:
            raise ValueError('cached baseline disagrees with scores or supplied annotation')
    if [row['gene'] for row in rebuilt] != [row['gene'] for row in source]:
        raise ValueError('cached baseline ordering not reproduced')
    outputs = {}
    for label, omitted in [('baseline', None)] + [(f'leave_one_out_{i}', term) for i, term in enumerate(terms, 1)]:
        selected = terms if omitted is None else [term for term in terms if term != omitted]
        rows = rank_rows(source, selected, annotations)
        outputs[label] = {'omitted_term': omitted, 'query_term_count': len(selected),
                          'tracked_gene_ranks': [rank_interval(rows, gene) for gene in tracked],
                          'top20': cutoff_summary(rows, 20), 'top50': cutoff_summary(rows, 50)}
    return {'baseline': baseline, 'gene_count': len(source), 'baseline_reproduced': True,
            'tolerance': TOLERANCE, 'scenarios': outputs, 'calibration_performed': False,
            'interpretation': 'Phenotype-score robustness only; no genotype inference, held-out validation or calibrated causal probability.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--annotations', type=Path, required=True)
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--track-gene', action='append', default=[])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cache = json.loads(args.cache.read_text())
    annotations = {}
    with args.annotations.open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            annotations.setdefault(row['gene_symbol'], set()).add(row['hpo_id'])
    result = analyze(cache, annotations, args.baseline, args.track_gene)
    result['input_sha256'] = {'cache': hashlib.sha256(args.cache.read_bytes()).hexdigest(),
                              'annotations': hashlib.sha256(args.annotations.read_bytes()).hexdigest()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('Sensitivity analysis complete; controlled results saved locally. No calibrated probability produced.')


if __name__ == '__main__':
    main()
