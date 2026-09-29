"""Interface compatibility, explicit entropy, and conditional partition sums."""

import math
import unittest

from single.statistical_mechanics.geometry import cube, dimer_ring, octahedron, tetrahedron
from single.statistical_mechanics.geometry_catalogue import monomer_candidates
from single.statistical_mechanics.geometry_weights import (
    BondPattern, GeometryState, conditional_geometry_shares, fully_bound_pattern,
    open_chain_pattern, remove_monomer, scaffold_pattern,
)


class GeometryWeightsTests(unittest.TestCase):
    def test_joint_assignments_on_compact_graphs(self):
        for geometry in monomer_candidates(12):
            if geometry.family != "compact polyhedron":
                continue
            with self.subTest(geometry=geometry.name):
                bonds = fully_bound_pattern(geometry)
                self.assertEqual(len(bonds.dimer_pairs), geometry.monomers // 2)
                self.assertEqual(len(bonds.c_terminal_contacts), geometry.monomers)
                available = {frozenset(edge) for edge in geometry.edges}
                dimers = {frozenset(edge) for edge in bonds.dimer_pairs}
                terminal = {frozenset(edge) for edge in bonds.c_terminal_contacts}
                self.assertTrue(dimers <= available and terminal <= available)
                self.assertFalse(dimers & terminal)
                self.assertEqual({a for a, _ in bonds.c_terminal_contacts},
                                 set(range(geometry.monomers)))
                self.assertEqual({b for _, b in bonds.c_terminal_contacts},
                                 set(range(geometry.monomers)))

    def test_odd_ring_cannot_assume_full_valence_with_paired_neighbours(self):
        geometries = {item.name: item for item in monomer_candidates(12)}
        for size in (3, 5, 11):
            with self.assertRaises(ValueError):
                fully_bound_pattern(geometries[f"closed ring ({size})"])
        even = fully_bound_pattern(geometries["closed ring (10)"])
        self.assertEqual((len(even.dimer_pairs), len(even.c_terminal_contacts)), (5, 10))

    def test_open_chain_counts_depend_on_interface_type(self):
        expected = {1: (0, 0), 2: (1, 0), 3: (1, 2), 4: (2, 2),
                    5: (2, 4), 6: (3, 4)}
        for size, counts in expected.items():
            bonds = open_chain_pattern(size)
            self.assertEqual((len(bonds.dimer_pairs), len(bonds.c_terminal_contacts)), counts)

    def test_equal_contacts_cancel_and_residual_energy_controls_shares(self):
        a, b = scaffold_pattern(cube()), scaffold_pattern(octahedron())
        for eps_d, eps_c in ((0, 0), (2, 3), (-2, 1)):
            equal = conditional_geometry_shares(
                [GeometryState("cube", a, 0, 0), GeometryState("octahedron", b, 0, 0)],
                epsilon_dimer=eps_d, epsilon_c_terminal=eps_c)
            self.assertAlmostEqual(equal[0].probability, 0.5)
        shifted = conditional_geometry_shares(
            [GeometryState("cube", a, 2, 0), GeometryState("octahedron", b, 0, 0)],
            epsilon_dimer=2, epsilon_c_terminal=3)
        self.assertAlmostEqual(shifted[1].probability, 1 / (1 + math.exp(-2)))
        entropy = conditional_geometry_shares(
            [GeometryState("cube", a, 0, math.log(3)), GeometryState("octahedron", b, 0, 0)],
            epsilon_dimer=2, epsilon_c_terminal=3)
        self.assertAlmostEqual(entropy[0].probability, 0.75)

    def test_states_are_summed_and_duplicates_and_mixed_sizes_rejected(self):
        a, b = scaffold_pattern(cube()), scaffold_pattern(octahedron())
        first = GeometryState("cube", a, 0, 0)
        probabilities = conditional_geometry_shares(
            [first, GeometryState("cube", a, math.log(2), 0),
             GeometryState("octahedron", b, 0, 0)],
            epsilon_dimer=0, epsilon_c_terminal=0)
        self.assertAlmostEqual(probabilities[0].probability, 0.6)
        for invalid in ([first, first],
                        [first, GeometryState("tetrahedron", scaffold_pattern(tetrahedron()), 0, 0)]):
            with self.assertRaises(ValueError):
                conditional_geometry_shares(invalid, epsilon_dimer=1, epsilon_c_terminal=1)

    def test_fixed_defect_breaks_both_incident_terminal_contacts(self):
        for parent in (dimer_ring(12), octahedron()):
            defect = remove_monomer(scaffold_pattern(parent), 0)
            self.assertEqual((defect.monomers, len(defect.dimer_pairs),
                              len(defect.c_terminal_contacts)), (23, 11, 22))

    def test_invalid_valences_and_fragments_are_rejected(self):
        cases = (
            (3, ((0, 1), (1, 2)), ()),
            (2, ((0, 1),), ((0, 1),)),
            (3, (), ((0, 1), (0, 2))),
            (4, ((0, 1), (2, 3)), ()),
        )
        for size, dimers, terminal in cases:
            with self.assertRaises(ValueError):
                BondPattern("invalid", size, dimers, terminal)


if __name__ == "__main__":
    unittest.main()
