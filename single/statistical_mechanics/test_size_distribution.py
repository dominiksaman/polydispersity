"""Closed-form rate-model check and the equal-contact scaffold identity."""

import math
import unittest

import numpy as np

from single.oligomer_distribution import oligomer_distribution
from single.statistical_mechanics import (
    energies_from_rates,
    site_count_distribution,
)
from single.statistical_mechanics.geometry import (
    cube, dimer_ring, octahedron,
)


class SizeDistributionTests(unittest.TestCase):
    def test_site_count_closed_form_matches_kinetic_recursion(self):
        for rates in ((20, 10, 1), (24, 12, 1), (24, 1, 1), (8, 2, 0.5),
                      (0.5, 2, 3), (20000, 1e8, 1), (1e100, 1e200, 1)):
            energies = energies_from_rates(*rates)
            _, kinetic = oligomer_distribution(*rates, max_size=80)
            statistical = site_count_distribution(
                energies.edge_kj_mol, energies.dimer_kj_mol, max_size=80
            )
            np.testing.assert_allclose(statistical, kinetic, rtol=1e-12,
                                       atol=1e-14)

    def test_equal_off_rates_give_positive_size_poisson(self):
        alpha = 8.0
        energies = energies_from_rates(alpha, 1.0, 1.0)
        actual = site_count_distribution(
            energies.edge_kj_mol, energies.dimer_kj_mol, max_size=40
        )
        expected = np.array([alpha ** (size - 1) / math.factorial(size)
                             for size in range(1, 41)])
        expected /= expected.sum()
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)

    def test_24mer_scaffolds_have_same_contact_counts(self):
        scaffolds = (dimer_ring(12), cube(), octahedron())
        self.assertEqual([s.vertex_degrees[0] for s in scaffolds], [2, 3, 4])
        for scaffold in scaffolds:
            self.assertEqual(scaffold.monomers, 24)
            self.assertEqual(len(scaffold.dimer_contacts), 12)
            self.assertEqual(len(scaffold.c_terminal_contacts), 24)
            contacts = scaffold.c_terminal_contacts
            self.assertEqual(set(donor for donor, _ in contacts), set(range(24)))
            self.assertEqual(set(receiver for _, receiver in contacts), set(range(24)))
            self.assertEqual(scaffold.contact_free_energy(-1, -2), -60)


if __name__ == "__main__":
    unittest.main()
