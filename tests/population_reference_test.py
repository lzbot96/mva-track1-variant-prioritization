"""Synthetic checks for gene-only requests and unknown-safe local joins."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

scripts = Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0, str(scripts))
import gnomad_public as fetcher
import join_population_references as m


class Tests(unittest.TestCase):
    def test_missing_stays_unknown(self):
        self.assertIsNone(m.frequency(None)['max_overall_af'])
        self.assertEqual(m.frequency(None)['status'], 'unknown_no_exact_match')

    def test_group_denominator_and_null_genome(self):
        r=m.frequency({'exome': {'ac': 2, 'an': 100, 'af': .02,
            'populations': [{'id': 'a', 'ac': 2, 'an': 20}, {'id': 'b', 'ac': 0, 'an': 0}]}, 'genome': None})
        self.assertEqual(r['raw_group_max_af'], .1)
        self.assertIsNone(r['cohorts']['genome'])

    def test_request_gene_only_and_legacy_default(self):
        self.assertNotIn('populations', fetcher.gene_query('GENE'))
        self.assertIn('populations { id ac an }', fetcher.gene_query('GENE', populations=True))
        for v in ('chr1-100-A-G', '1-100-A-G', '../input'):
            with self.assertRaises(ValueError): fetcher.gene_query(v)

    def test_local_join_exact_allele_no_network(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'ref.json').write_text(json.dumps({'data': {'gene': {'variants': [
                {'variant_id': '1-100-A-G','exome':{'ac':1,'an':100,'af':.01,'populations':[]},'genome':None}]}}}))
            (root/'input.csv').write_text('chrom,pos,ref,alt\nchr1,100,A,G\nchr1,100,A,T\n')
            with patch.object(fetcher.urllib.request,'urlopen',side_effect=AssertionError('No network')):
                r=m.join(root,root/'input.csv')['results']
            self.assertEqual(r[0]['max_overall_af'],.01)
            self.assertIsNone(r[1]['max_overall_af'])


if __name__=='__main__': unittest.main()
