"""Mass balance, chemical potentials, analytic limits, and reservoir sampling."""

import unittest

import numpy as np

from .reservoir_thermodynamics import (
    GlobularModel, closed_equilibrium, equilibrium_at_activity,
    globular_equilibrium, sample_reservoir,
)


class ReservoirTests(unittest.TestCase):
    def test_monomer_dimer_limit_matches_quadratic_mass_balance(self):
        binding, total, c0 = 4.0, 0.2, 1.0
        result = closed_equilibrium([0, -binding], total, standard_concentration=c0)
        k = np.exp(binding) / c0
        monomer = 2 * total / (1 + np.sqrt(1 + 8 * k * total))
        np.testing.assert_allclose(result.concentration, [monomer, k * monomer**2], rtol=1e-12)
        pure = closed_equilibrium([0], total, standard_concentration=c0)
        self.assertAlmostEqual(pure.concentration[0], total)

    def test_shared_chemical_potential_and_mass_conservation(self):
        energy = GlobularModel().free_energies()
        for total in (1e-5, 0.03, 3, 30):
            r = closed_equilibrium(energy, total)
            self.assertAlmostEqual(float(r.sizes @ r.concentration) / total, 1, places=11)
            self.assertAlmostEqual(float(r.number_fraction.sum()), 1)
            self.assertAlmostEqual(float(r.subunit_fraction.sum()), 1)
            present = r.concentration > 1e-200
            np.testing.assert_allclose(np.log(r.concentration[present]) + energy[present],
                                       r.sizes[present] * r.log_activity, atol=1e-12)
            reopened = equilibrium_at_activity(energy, r.log_activity)
            np.testing.assert_allclose(reopened.concentration, r.concentration)

    def test_standard_state_conversion_preserves_physical_concentrations(self):
        energy = GlobularModel().free_energies()
        n = np.arange(1, len(energy) + 1)
        before = closed_equilibrium(energy, 3, standard_concentration=1)
        after = closed_equilibrium(energy - (n - 1) * np.log(1000), 3,
                                   standard_concentration=1000)
        np.testing.assert_allclose(before.concentration, after.concentration, rtol=1e-11, atol=1e-14)

    def test_size_support_extends_without_imposing_a_24mer_cap(self):
        small = globular_equilibrium(GlobularModel(), 3, max_size=12)
        large = globular_equilibrium(GlobularModel(), 3, max_size=240)
        self.assertGreater(len(small.sizes), 24)
        self.assertAlmostEqual(small.log_activity, large.log_activity, places=10)
        self.assertGreater(small.subunit_fraction[small.sizes > 24].sum(), 0)
        self.assertEqual(int(small.sizes[np.argmax(small.conditional_number_fraction())]), 16)
        with self.assertRaises(ValueError):
            globular_equilibrium(GlobularModel(), 3, max_size=12, size_limit=12)

    def test_pool_buffers_growth_while_large_clusters_take_up_protein(self):
        low = globular_equilibrium(GlobularModel(), 0.003)
        high = globular_equilibrium(GlobularModel(), 3)
        self.assertLess(low.subunit_fraction[2:].sum(), 0.01)
        self.assertGreater(high.subunit_fraction[2:].sum(), 0.9)
        self.assertGreater(high.concentration[0], low.concentration[0])
        self.assertLess(high.concentration[0] / low.concentration[0], 100)

    def test_grand_canonical_snapshots_obey_count_statistics(self):
        r = globular_equilibrium(GlobularModel(), 3)
        scale, snapshots = 1000, 4000
        counts = sample_reservoir(r, standard_state_particles=scale, snapshots=snapshots)
        expected = scale * r.concentration / r.standard_concentration
        error = np.abs(counts.mean(axis=0) - expected)
        self.assertTrue(np.all(error < 6 * np.sqrt(expected / snapshots) + 0.005))
        total_counts = counts @ r.sizes
        self.assertGreater(total_counts.std(), 0)  # Open reservoir, not fixed N.
        self.assertLess(abs(total_counts.mean() - 3 * scale), 6 * total_counts.std() / np.sqrt(snapshots))

    def test_invalid_inputs(self):
        for total in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                closed_equilibrium([0, -4], total)
        for energy in ([1, 2], [], [0, float("nan")]):
            with self.assertRaises(ValueError):
                closed_equilibrium(energy, 1)
        with self.assertRaises(ValueError):
            GlobularModel(packing_kbt=0)


if __name__ == "__main__":
    unittest.main()
