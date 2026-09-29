"""Detailed balance, numerical stability, and number-fraction semantics."""

import unittest

import numpy as np

from single.mass_distribution import mass_distribution
from single.oligomer_distribution import oligomer_distribution


class OligomerDistributionTests(unittest.TestCase):
    def test_detailed_balance_and_rate_rescaling(self):
        rates = (20.0, 10.0, 1.0)
        sizes, p = oligomer_distribution(*rates, max_size=60)
        for size in range(2, 61):
            off = size if size % 2 == 0 else size - 1 + 10
            self.assertAlmostEqual(p[size - 1] / p[size - 2], 20 / off)
        _, rescaled = oligomer_distribution(60, 30, 3)
        np.testing.assert_allclose(p, rescaled, atol=1e-14)
        self.assertEqual(int(sizes[np.argmax(p)]), 16)
        self.assertAlmostEqual(float(sizes @ p), 16.0126065179, places=8)

    def test_large_on_rate_stays_normalized(self):
        sizes, p = oligomer_distribution(1000, 1, 1, max_size=1500)
        self.assertTrue(np.all(np.isfinite(p)))
        self.assertAlmostEqual(float(p.sum()), 1)
        self.assertAlmostEqual(float(sizes @ p), 1000, places=7)

    def test_monomer_relative_weights(self):
        _, weights = oligomer_distribution(20, 10, 1, normalise=False)
        _, p = oligomer_distribution(20, 10, 1)
        self.assertEqual(weights[0], 1)
        np.testing.assert_allclose(p, weights / weights.sum(), atol=1e-14)
        with self.assertRaises(OverflowError):
            oligomer_distribution(1000, 1, 1, max_size=1500, normalise=False)

    def test_invalid_inputs_are_rejected(self):
        for bad in (0, -1, float("nan"), float("inf")):
            for index in range(3):
                rates = [20, 10, 1]
                rates[index] = bad
                with self.assertRaises(ValueError):
                    oligomer_distribution(*rates)
        for size in (0, -1, 2.5, True):
            with self.assertRaises(ValueError):
                oligomer_distribution(20, 10, 1, max_size=size)

    def test_mass_axis_preserves_number_fractions(self):
        sizes, p = oligomer_distribution(20, 10, 1)
        masses, mass_axis_p = mass_distribution(20, 10, 1, 20000)
        np.testing.assert_array_equal(masses, sizes * 20000)
        np.testing.assert_array_equal(p, mass_axis_p)
        for bad in (0, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                mass_distribution(20, 10, 1, bad)


if __name__ == "__main__":
    unittest.main()
