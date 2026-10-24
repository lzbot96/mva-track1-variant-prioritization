"""Local exact-allele joins only. No network imports or calls.

Inputs are private candidate CSVs and whole public-gene gnomAD JSON references.
Missing records remain unknown. Group AC/AN maxima are not FAF95, calibrated
causal probabilities, or ancestry matched to a person. Keep joined results private.
"""
import argparse
import csv
import json
from pathlib import Path


def frequency(record):
    if record is None:
        return {'status': 'unknown_no_exact_match', 'max_overall_af': None,
                'raw_group_max_af': None, 'cohorts': {}}
    cohorts = {}
    for name in ('exome', 'genome'):
        v = record.get(name)
        if not v:
            cohorts[name] = None
            continue
        groups = [{**g, 'af': g['ac']/g['an']} for g in v.get('populations', [])
                  if g.get('ac') is not None and g.get('an') and g['an'] > 0]
        maximum = max((g['af'] for g in groups), default=None)
        cohorts[name] = {'ac': v.get('ac'), 'an': v.get('an'), 'af': v.get('af'),
            'raw_group_max_af': maximum, 'maximum_groups': [g for g in groups if g['af'] == maximum]}
    overall = [v['af'] for v in cohorts.values() if v and v['af'] is not None]
    group = [v['raw_group_max_af'] for v in cohorts.values() if v and v['raw_group_max_af'] is not None]
    return {'status': 'matched_public_record', 'max_overall_af': max(overall, default=None),
            'raw_group_max_af': max(group, default=None), 'cohorts': cohorts}


def join(reference_dir, candidate_csv):
    index = {}
    for path in sorted(Path(reference_dir).glob('*.json')):
        data = json.loads(path.read_text())
        if not isinstance(data, dict):
            continue
        gene = (data.get('data') or {}).get('gene')
        if not isinstance(gene, dict):
            continue
        for v in gene.get('variants', []):
            if v['variant_id'] in index and frequency(index[v['variant_id']]) != frequency(v):
                raise ValueError('Conflicting reference snapshots; reconcile before joining')
            index[v['variant_id']] = v
    result = []
    with Path(candidate_csv).open(newline='') as stream:
        reader = csv.DictReader(stream)
        compound = 'chrom_1' in (reader.fieldnames or [])
        for row in reader:
            for suffix in ('_1', '_2') if compound else ('',):
                if not row.get('chrom'+suffix):
                    continue
                key = '-'.join([row['chrom'+suffix].removeprefix('chr'), row['pos'+suffix],
                                row['ref'+suffix], row['alt'+suffix]])
                result.append({'private_variant_key': key, **frequency(index.get(key))})
    return {'scope': 'Private local exact-allele join; no data transmitted.',
        'limits': 'Missing is unknown; raw group maximum is not FAF95 or ancestry matched.', 'results': result}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--references', required=True)
    p.add_argument('--candidates', required=True)
    p.add_argument('--output', required=True)
    a = p.parse_args()
    result = join(a.references, a.candidates)
    output = Path(a.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2)+'\n')
    print('Saved private local joins:', len(result['results']), 'alleles')
