"""Synthetic regressions: no patient data or network access."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from build_review_gene_manifest import select_boundary_ties
from phenotype_sensitivity import rank_rows, rank_interval, cutoff_summary, analyze


class PhenotypeSensitivityTests(unittest.TestCase):
    def test_alphabetical_boundary_does_not_exclude_equal_gene(self):
        rows=[{'gene':g,'mean_resnik':s,'exact_query_terms':e} for g,s,e in [('A',3,2),('B',2,1),('Z',2,1),('C',2,0),('D',1,2)]]
        self.assertEqual(select_boundary_ties(rows,2),['A','B','Z'])
        self.assertEqual(cutoff_summary(rows,2)['selected_genes'],['A','B','Z'])

    def test_tolerance_never_removes_original_top_n(self):
        rows = [{'gene': 'A', 'mean_resnik': 1 + 5e-13, 'exact_query_terms': 0},
                {'gene': 'B', 'mean_resnik': 1, 'exact_query_terms': 1}]
        self.assertEqual(select_boundary_ties(rows, 2), ['A', 'B'])
        self.assertEqual(cutoff_summary(rows, 2)['selected_genes'], ['A', 'B'])

    def test_ties_preserve_rank_interval(self):
        rows=[{'gene':g,'mean_resnik':s,'exact_query_terms':e} for g,s,e in [('A',3,2),('B',2,1),('Z',2,1),('C',2,0)]]
        self.assertEqual(rank_interval(rows,'Z')['mean_rank_min'],2)
        self.assertEqual(rank_interval(rows,'Z')['mean_rank_max'],4)
        self.assertEqual(rank_interval(rows,'Z')['composite_rank_max'],3)

    def test_leave_one_out_can_reverse_ranking(self):
        cached=[{'gene':'A','per_term':{'T1':6,'T2':0}}, {'gene':'B','per_term':{'T1':2,'T2':2}}]
        annotations={'A':{'T1'},'B':{'T2'}}
        self.assertEqual(rank_rows(cached,['T1','T2'],annotations)[0]['gene'],'A')
        self.assertEqual(rank_rows(cached,['T2'],annotations)[0]['gene'],'B')

    def test_stale_annotations_reject_baseline(self):
        cached=[{'gene':'A','mean_resnik':3,'exact_query_terms':1,'per_term':{'T1':6,'T2':0}}]
        with self.assertRaises(ValueError): analyze({'rankings':{'test':cached}}, {}, 'test', [])

    def test_invalid_score_rejected(self):
        with self.assertRaises(ValueError):select_boundary_ties([{'gene':'A','mean_resnik':'nan','exact_query_terms':1}],1)
        with self.assertRaises(ValueError):select_boundary_ties([],0)


if __name__ == '__main__':
    unittest.main()
