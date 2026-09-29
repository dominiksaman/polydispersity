"""Scientific checks for geometry-conditioned partition weights."""

import unittest

from single.statistical_mechanics.geometry_ensemble import (
    EnergyParameters, conditional_contributions, default_candidates,
    split_octahedron,
)


class GeometryEnsembleTests(unittest.TestCase):
    def test_every_supported_size_normalizes_and_has_unique_names(self):
        candidates = default_candidates(40)
        for size in range(5, 41):
            with self.subTest(size=size):
                ensemble = conditional_contributions(size, candidates=candidates)
                self.assertTrue(ensemble)
                self.assertAlmostEqual(sum(item.share for item in ensemble), 1)
                self.assertEqual(len({item.geometry for item in ensemble}), len(ensemble))

    def test_even_24mer_contact_energies_cancel(self):
        for epsilon_d, epsilon_c in ((0, 0), (1, 1), (20, -3)):
            ensemble = conditional_contributions(
                24, EnergyParameters(epsilon_d, epsilon_c, 0, 4))
            self.assertEqual({item.geometry for item in ensemble},
                             {"single ring (12 dimers)", "cube", "octahedron"})
            self.assertTrue(all(item.dimer_contacts == 12 and
                                item.c_terminal_contacts == 24 for item in ensemble))
            for item in ensemble:
                self.assertAlmostEqual(item.share, 1 / 3)

    def test_degree_preference_and_offsets_change_shares(self):
        params = EnergyParameters(1, 1, 0.3, 4)
        shares = {item.geometry: item.share
                  for item in conditional_contributions(24, params)}
        self.assertGreater(shares["octahedron"], shares["cube"])
        self.assertGreater(shares["cube"], shares["single ring (12 dimers)"])
        offset_shares = {item.geometry: item.share for item in
                         conditional_contributions(24, params,
                             offsets_kbt={"cube": -5})}
        self.assertGreater(offset_shares["cube"], shares["cube"])

    def test_odd_defects_count_sites_and_contacts(self):
        ensemble = conditional_contributions(23, EnergyParameters(0, 0, 0, 4))
        by_name = {item.geometry: item for item in ensemble}
        self.assertEqual(by_name["single ring (12 dimers)"].microstates, 24)
        self.assertEqual(by_name["single ring (12 dimers)"].c_terminal_contacts, 22)
        self.assertEqual(by_name["cube"].c_terminal_contacts, 23)
        self.assertEqual(by_name["octahedron"].c_terminal_contacts, 23)
        for item in ensemble:
            self.assertEqual(item.dimer_contacts, 11)
            self.assertAlmostEqual(item.share, 1 / 3)

    def test_vertex_splits_have_intended_capacity(self):
        one, two = split_octahedron(), split_octahedron(True)
        self.assertEqual((one.monomers, one.vertex_degrees),
                         (26, (3, 3, 4, 4, 4, 4, 4)))
        self.assertEqual((two.monomers, two.vertex_degrees),
                         (28, (3, 3, 3, 3, 4, 4, 4, 4)))


if __name__ == "__main__":
    unittest.main()
