"""Thermodynamic reformulation and contact-model comparison checks."""

import math
import unittest

import numpy as np

from single.oligomer_distribution import oligomer_distribution
from single.statistical_mechanics import (
    additive_contact_distribution,
    energies_from_rates,
    fit_ideal_cluster,
    site_count_distribution,
)
from single.statistical_mechanics.geometry import (
    cube, dimer_ring, octahedron, ring_distribution,
)


class SizeDistributionTests(unittest.TestCase):
    def test_site_count_closed_form_matches_kinetic_recursion(self):
        for rates in ((24, 12, 1), (24, 1, 1), (8, 2, 0.5), (0.5, 2, 3)):
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

    def test_simple_contacts_do_not_preserve_rate_mapping(self):
        energies = energies_from_rates(24, 12, 1)
        _, kinetic = oligomer_distribution(24, 12, 1)
        no_entropy = additive_contact_distribution(
            energies.edge_kj_mol, energies.dimer_kj_mol
        )
        factorial = additive_contact_distribution(
            energies.edge_kj_mol, energies.dimer_kj_mol,
            factorial_entropy=True,
        )
        self.assertEqual(int(np.argmax(kinetic)) + 1, 18)
        self.assertEqual(int(np.argmax(no_entropy)) + 1, 60)
        self.assertEqual(int(np.argmax(factorial)) + 1, 6)

    def test_refitted_factorial_model_is_similar_but_not_identical(self):
        _, kinetic = oligomer_distribution(24, 12, 1)
        energies = energies_from_rates(24, 12, 1)
        fitted = fit_ideal_cluster(kinetic)
        self.assertTrue(fitted.success)
        self.assertLess(fitted.total_variation, 0.06)
        self.assertGreater(abs(fitted.dimer_kj_mol - energies.dimer_kj_mol), 1)
        self.assertEqual(int(np.argmax(fitted.distribution)) + 1, 18)

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

    def test_constant_ring_energy_cannot_give_an_interior_peak(self):
        for delta_g, expected_mode in ((-1, 60), (1, 6)):
            sizes, probabilities = ring_distribution(delta_g)
            self.assertEqual(int(sizes[np.argmax(probabilities)]), expected_mode)
            ratios = probabilities[1:] / probabilities[:-1]
            np.testing.assert_allclose(ratios, ratios[0], rtol=1e-14)
        _, neutral = ring_distribution(0)
        np.testing.assert_allclose(neutral, np.full(28, 1 / 28))


if __name__ == "__main__":
    unittest.main()
