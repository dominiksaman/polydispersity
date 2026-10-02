"""Independent analytic limits and concentration-series equilibrium checks."""

import unittest

import numpy as np
from scipy.special import lambertw

from single.kinetic_equilibrium import KineticModel, anchor_kinetic_model, kinetic_equilibrium
from single.oligomer_distribution import oligomer_distribution
from single.compare_concentrations import predict_series


class KineticEquilibriumTests(unittest.TestCase):
    def test_equal_off_rates_have_analytic_mass_balance(self):
        # For equal off-rates: c_n=c_1*alpha**(n-1)/n!, C=c_1*exp(alpha).
        # Hence alpha = W(k_plus*C/k_off), independently of the root solver.
        model = KineticModel(7, 2, 2)
        for total in (1e-8, .03, 3, 300):
            result = kinetic_equilibrium(model, total)
            alpha = float(lambertw(7 * total / 2).real)
            self.assertAlmostEqual(result.effective_on / 2, alpha, places=10)
            self.assertAlmostEqual(result.free_monomer, 2 * alpha / 7, places=10)
            self.assertAlmostEqual(result.sizes @ result.concentration / total, 1, places=10)

    def test_reference_recovery_and_detailed_balance(self):
        model = anchor_kinetic_model()
        anchor = kinetic_equilibrium(model, 3)
        _, expected = oligomer_distribution(20, 10, 1, max_size=len(anchor.sizes))
        np.testing.assert_allclose(anchor.number_fraction, expected, atol=1e-13)
        self.assertAlmostEqual(anchor.effective_on, 20, places=10)
        for total in (.003, 3, 30):
            r = kinetic_equilibrium(model, total)
            np.testing.assert_allclose(r.number_fraction.sum(), 1, atol=1e-12)
            np.testing.assert_allclose(r.subunit_fraction.sum(), 1, atol=1e-12)
            n = r.sizes[1:]
            off = np.where(n % 2 == 0, n, n - 1 + 10)
            np.testing.assert_allclose(np.diff(r.log_concentration),
                                       np.log(r.effective_on) - np.log(off), atol=1e-12)
            self.assertAlmostEqual(r.sizes @ r.concentration, total, places=10)

    def test_concentration_units_and_time_rescaling(self):
        model = anchor_kinetic_model()
        r = kinetic_equilibrium(model, 3)
        # µM -> nM multiplies all concentrations by 1000, divides k_plus by 1000.
        converted = kinetic_equilibrium(KineticModel(model.k_plus / 1000, 10, 1), 3000)
        rescaled_time = kinetic_equilibrium(KineticModel(model.k_plus * 9, 90, 9), 3)
        np.testing.assert_allclose(converted.concentration, r.concentration * 1000, atol=1e-11)
        np.testing.assert_allclose(rescaled_time.concentration, r.concentration, atol=1e-12)

    def test_adaptive_support_and_large_weights(self):
        model = anchor_kinetic_model()
        small = kinetic_equilibrium(model, 30, max_size=12)
        large = kinetic_equilibrium(model, 30, max_size=240)
        self.assertGreater(len(small.sizes), 24)
        np.testing.assert_allclose(small.concentration, large.concentration[:len(small.sizes)], atol=1e-12)
        r = kinetic_equilibrium(KineticModel(1e200, 1, 1), 1)
        self.assertTrue(np.all(np.isfinite(r.number_fraction)))
        self.assertAlmostEqual(r.effective_on, float(lambertw(1e200).real), places=6)
        self.assertAlmostEqual(r.sizes @ r.concentration, 1, places=7)
        with self.assertRaisesRegex(ValueError, "not converged"):
            kinetic_equilibrium(model, 30, max_size=12, size_limit=12)

    def test_series_uses_one_fixed_parameter_set(self):
        series = predict_series([.003, 3, 30])
        for total, kinetic, thermo in zip(series.totals, series.kinetic, series.thermodynamic):
            self.assertAlmostEqual(kinetic.sizes @ kinetic.concentration, total, places=10)
            self.assertAlmostEqual(thermo.sizes @ thermo.concentration, total, places=10)
            np.testing.assert_array_equal(thermo.free_energy_kbt,
                                         series.thermodynamic_model.free_energies(len(thermo.sizes)))
            self.assertAlmostEqual(kinetic.effective_on,
                                   series.kinetic_model.k_plus * kinetic.free_monomer)
        self.assertAlmostEqual(series.thermodynamic[1].subunit_fraction[:2].sum(), .2, places=10)
        self.assertGreater(series.thermodynamic[0].subunit_fraction[:2].sum(), .99)
        self.assertLess(series.thermodynamic[2].subunit_fraction[:2].sum(), .03)
        self.assertLess(series.kinetic[0].effective_on, 20)
        self.assertGreater(series.kinetic[2].effective_on, 20)

    def test_invalid_inputs(self):
        for bad in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                KineticModel(bad)
            with self.assertRaises(ValueError):
                kinetic_equilibrium(KineticModel(1), bad)
            with self.assertRaises(ValueError):
                anchor_kinetic_model(bad)
        for kwargs in (dict(max_size=2), dict(max_size=True), dict(max_size=2.5),
                       dict(size_limit=12), dict(tail_tolerance=0), dict(tail_tolerance=float("nan"))):
            with self.assertRaises(ValueError):
                kinetic_equilibrium(KineticModel(1), 1, **kwargs)
        for bad in ([], [0], [float("nan")], [[1, 2]]):
            with self.assertRaises(ValueError):
                predict_series(bad)
        result = kinetic_equilibrium(KineticModel(1), 1)
        for bad in (0, True, 2.5, len(result.sizes) + 1):
            with self.assertRaises(ValueError):
                result.conditional_number_fraction(bad)


if __name__ == "__main__":
    unittest.main()
