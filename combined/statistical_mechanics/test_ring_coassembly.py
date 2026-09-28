"""Checks against explicit rings and the thesis's limiting cases."""

import itertools
import math
import unittest
from collections import Counter

import numpy as np

from combined.statistical_mechanics import (
    activity_distribution,
    composition_distribution,
    contact_classes,
    fit_hetero_energy,
)


class RingCoassemblyTests(unittest.TestCase):
    def test_contact_counts_match_explicit_rings(self):
        for size in range(3, 11):
            explicit = Counter()
            for ring in itertools.product("AB", repeat=size):
                n_a = ring.count("A")
                n_aa = sum(ring[i] == ring[(i + 1) % size] == "A"
                           for i in range(size))
                n_bb = sum(ring[i] == ring[(i + 1) % size] == "B"
                           for i in range(size))
                n_ab = size - n_aa - n_bb
                explicit[(n_a, n_aa, n_ab, n_bb)] += 1
            grouped = {
                (item.n_a, item.n_aa, item.n_ab, item.n_bb): item.multiplicity
                for item in contact_classes(size)
            }
            self.assertEqual(grouped, dict(explicit))

    def test_unbiased_is_binomial_at_equal_and_unequal_mixing(self):
        size = 12
        for fraction in (0.2, 0.5, 0.8):
            actual = composition_distribution(size, 0.0, mole_fraction_a=fraction)
            expected = np.array([
                math.comb(size, i) * fraction**i * (1 - fraction)**(size - i)
                for i in range(size + 1)
            ])
            np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-14)

    def test_energy_biases_favor_expected_arrangements(self):
        size = 12
        self_favored = composition_distribution(size, 12.5)
        hetero_favored = composition_distribution(size, -12.5)
        self.assertGreater(self_favored[0] + self_favored[-1], 0.99)
        self.assertGreater(hetero_favored[size // 2], 0.99)
        self.assertAlmostEqual(self_favored.sum(), 1)
        self.assertAlmostEqual(hetero_favored.sum(), 1)

    def test_fixed_mean_fraction_with_unequal_homotypic_energies(self):
        size = 20
        distribution = composition_distribution(
            size, 1.5, g_aa=-1.0, g_bb=0.5, mole_fraction_a=1/3
        )
        self.assertAlmostEqual(np.dot(np.arange(size + 1), distribution) / size,
                               1/3, places=12)

    def test_joint_fit_recovers_synthetic_energy_and_activities(self):
        size = 12
        energy = 1.98
        log_ratios = (math.log(0.5), 0.0, math.log(2.0))
        observed = [activity_distribution(size, energy, log_activity_ratio=ratio)
                    for ratio in log_ratios]
        fitted = fit_hetero_energy(observed)
        self.assertTrue(fitted.success)
        self.assertAlmostEqual(fitted.hetero_energy_kj_mol, energy, places=8)
        np.testing.assert_allclose(fitted.log_activity_ratios, log_ratios,
                                   rtol=0, atol=1e-8)
        self.assertLess(fitted.root_mean_square_error, 1e-10)


if __name__ == "__main__":
    unittest.main()
