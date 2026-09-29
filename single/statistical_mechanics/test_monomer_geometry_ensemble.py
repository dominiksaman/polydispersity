"""Check the small-size geometry hypothesis and its contact bookkeeping."""

from functools import lru_cache
import unittest

from single.statistical_mechanics.monomer_geometry_ensemble import (
    MonomerEnergy, monomer_candidates, monomer_contributions,
)


def matching_size(size, edges):
    neighbours = [0] * size
    for a, b in edges:
        neighbours[a] |= 1 << b
        neighbours[b] |= 1 << a

    @lru_cache(None)
    def solve(mask):
        if not mask:
            return 0
        bit = mask & -mask
        vertex = bit.bit_length() - 1
        rest = mask ^ bit
        best = solve(rest)
        options = neighbours[vertex] & rest
        while options:
            partner_bit = options & -options
            best = max(best, 1 + solve(rest ^ partner_bit))
            options ^= partner_bit
        return best

    return solve((1 << size) - 1)


def has_hamiltonian_cycle(size, edges):
    neighbours = [set() for _ in range(size)]
    for a, b in edges:
        neighbours[a].add(b)
        neighbours[b].add(a)

    @lru_cache(None)
    def visit(current, mask):
        if mask == (1 << size) - 1:
            return 0 in neighbours[current]
        return any(visit(next_vertex, mask | (1 << next_vertex))
                   for next_vertex in neighbours[current]
                   if not mask & (1 << next_vertex))

    return visit(0, 1)


class MonomerGeometryTests(unittest.TestCase):
    def test_small_default_catalogue_and_optional_compact_shapes(self):
        default = monomer_candidates(12)
        for size in (3, 4, 5):
            self.assertEqual({item.family for item in default if item.monomers == size},
                             {"open chain", "closed ring"})
        six = {item.name for item in default if item.monomers == 6}
        self.assertIn("octahedron (6 monomers)", six)
        self.assertIn("triangular prism (6 monomers)", six)
        eleven = [item for item in default if item.monomers == 11]
        self.assertEqual({item.family for item in eleven},
                         {"open chain", "closed ring", "compact polyhedron"})
        self.assertEqual(len([item for item in eleven if item.family ==
                              "compact polyhedron"]), 2)
        optional = monomer_candidates(5, include_small_polyhedra=True)
        self.assertIn("tetrahedron (4 monomers)", {item.name for item in optional})
        self.assertIn("triangular dipyramid (5 monomers)",
                      {item.name for item in optional})

    def test_graphs_admit_assumed_dimer_and_terminal_contacts(self):
        for item in monomer_candidates(12, include_small_polyhedra=True):
            with self.subTest(geometry=item.name):
                self.assertEqual(matching_size(item.monomers, item.edges),
                                 item.dimer_contacts)
                if item.closed:
                    self.assertTrue(has_hamiltonian_cycle(item.monomers, item.edges))
                self.assertEqual(item.c_terminal_contacts,
                                 item.monomers if item.closed
                                 else max(0, item.monomers - 1))

    def test_even_closed_geometries_share_contact_counts(self):
        values = monomer_contributions(6, MonomerEnergy(0, 0, 0))
        closed = [item for item in values if item.family != "open chain"]
        self.assertGreater(len(closed), 2)
        self.assertTrue(all(item.dimer_contacts == 3 and
                            item.c_terminal_contacts == 6 for item in closed))
        self.assertAlmostEqual(sum(item.share for item in values), 1)

    def test_eleven_vertex_polyhedra_are_distinct_and_normalized(self):
        compact = [item for item in monomer_candidates(12)
                   if item.monomers == 11 and item.family == "compact polyhedron"]
        self.assertEqual({len(item.edges) for item in compact}, {20, 25})
        self.assertTrue(all(item.dimer_contacts == 5 and
                            item.c_terminal_contacts == 11 for item in compact))
        values = monomer_contributions(11, MonomerEnergy(1, 1, 0.1))
        self.assertAlmostEqual(sum(item.share for item in values), 1)
        self.assertEqual(len(values), 4)

    def test_contact_and_shape_energies_change_shares(self):
        no_terminal = monomer_contributions(3, MonomerEnergy(1, 0, 0))
        self.assertAlmostEqual(no_terminal[0].share, no_terminal[1].share)
        terminal = monomer_contributions(3, MonomerEnergy(1, 1, 0))
        self.assertGreater(terminal[1].share, terminal[0].share)
        no_penalty = {item.geometry: item.share for item in
                      monomer_contributions(6, MonomerEnergy(1, 1, 0))}
        penalty = {item.geometry: item.share for item in
                   monomer_contributions(6, MonomerEnergy(1, 1, 0.2))}
        self.assertLess(penalty["octahedron (6 monomers)"],
                        no_penalty["octahedron (6 monomers)"])


if __name__ == "__main__":
    unittest.main()
