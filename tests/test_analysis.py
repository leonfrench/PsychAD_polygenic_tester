import unittest
import numpy as np
import pandas as pd
from analysis import run_analysis, adjust_bh, parse_genes, load_data

class AnalysisTests(unittest.TestCase):
    def sample(self, scores):
        return pd.DataFrame({'coef': ['CERAD'] * len(scores), 'assay': ['Micro'] * len(scores), 'ID': list('ABCDE')[:len(scores)], 'statistic': scores, 'estimate': scores})

    def test_auc_pairwise_ties(self):
        d = self.sample([2., 3., 2., 1., 4.])
        r = run_analysis(d, ['A', 'B']).iloc[0]
        expected = np.mean([float(x > y) + .5 * (x == y) for x in [2, 3] for y in [2, 1, 4]])
        self.assertAlmostEqual(r.AUROC, expected)
        self.assertEqual(r['Comparison genes (n)'], 3)

    def test_universe_and_missing(self):
        d = self.sample([1., 2., 3., 4., 5.])
        r = run_analysis(d, ['A', 'E'], ['A', 'B', 'C']).iloc[0]
        self.assertEqual(r.AUROC, 0)
        self.assertEqual(r['Target genes (n)'], 1)
        self.assertEqual(r['Target genes not used'], 'E')

    def test_all_tied(self):
        r = run_analysis(self.sample([1., 1., 1.]), ['A']).iloc[0]
        self.assertEqual(r.AUROC, .5)
        self.assertEqual(r['p-value'], 1.)

    def test_untestable_and_nonfinite(self):
        d = self.sample([1., np.nan, 3.])
        r = run_analysis(d, ['A', 'C']).iloc[0]
        self.assertTrue(np.isnan(r.AUROC))
        self.assertEqual(r.Status, 'No comparison genes')

    def test_bh(self):
        np.testing.assert_allclose(adjust_bh([.01, .04, .03, .2]), [.04, .0533333333, .0533333333, .2])

    def test_parser(self):
        self.assertEqual(parse_genes('A,A;B\nC\tB'), ['A','B','C'])

    def test_real_data(self):
        d = load_data()
        r = run_analysis(d, ['APOE', 'TREM2', 'BIN1', 'CLU'])
        self.assertEqual(len(r), 108)
        self.assertTrue(r.AUROC.between(0, 1).all())
        self.assertGreater(r['Available genes (n)'].nunique(), 1)
        self.assertTrue((r['FDR (BH)'] >= r['p-value'] - 1e-15).all())

if __name__ == '__main__':
    unittest.main()
