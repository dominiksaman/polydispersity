"""Checks for the rate-free equilibrium scaffold partition function."""

from collections import Counter
import unittest

import numpy as np

from single.statistical_mechanics.equilibrium_scaffold import (
    EquilibriumParameters, enumerate_states,
)
from single.statistical_mechanics.geometry import cube, dimer_ring, octahedron


class EquilibriumScaffoldTests(unittest.TestCase):
    def test_ring_has_one_distinct_connected_state_per_size(self):
        catalogue = enumerate_states(dimer_ring(12))
        by_size = Counter()
        for (n, _d, _c, _v), count in catalogue.counts.items():
            by_size[n] += count
        self.assertEqual([by_size[n] for n in range(1, 25)], [1] * 24)

    def test_zero_energies_give_topology_count_distribution(self):
        catalogue = enumerate_states(octahedron())
        by_size = Counter()
        for (n, _d, _c, _v), count in catalogue.counts.items():
            by_size[n] += count
        counts = np.array([by_size[n] for n in range(1, 25)])
        probabilities = catalogue.distribution(EquilibriumParameters(0, 0, 0))
        np.testing.assert_allclose(probabilities, counts / counts.sum())
        self.assertEqual(by_size[24], 1)

    def test_example_parameters_select_a_polyhedral_interior_size(self):
        parameters = EquilibriumParameters(-1, 4, 0.05, 1.1)
        modes = [int(np.argmax(enumerate_states(scaffold).distribution(parameters))) + 1
                 for scaffold in (dimer_ring(12), cube(), octahedron())]
        self.assertEqual(modes, [24, 18, 18])


if __name__ == "__main__":
    unittest.main()
