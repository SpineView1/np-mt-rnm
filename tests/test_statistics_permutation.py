"""permutation_pvalues: the vectorised perm_p_two_sided used for the rescue figures."""
import numpy as np

from np_mt_rnm.statistics import permutation_pvalues


def test_permutation_pvalues_separated_and_identical_columns():
    rng = np.random.RandomState(0)
    x = np.column_stack([np.ones(50), np.full(50, 0.3), rng.random_sample(50)])
    y = np.column_stack([np.zeros(50), np.full(50, 0.3), rng.random_sample(50)])
    p = permutation_pvalues(x, y, n_perm=200, rng=np.random.RandomState(4001))
    assert p[0] == 1 / 201          # fully separated: no permutation as extreme
    assert p[1] == 1.0              # identical constants: every permutation ties
    assert 0.0 < p[2] <= 1.0
