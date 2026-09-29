"""Test tetrahedral parent scaffolds against the requested kinetic example.

Run: python -m single.statistical_mechanics.compare_tetrahedra
"""

from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from single.oligomer_distribution import oligomer_distribution
from single.statistical_mechanics.equilibrium_scaffold import enumerate_states
from single.statistical_mechanics.fit_scaffold_to_kinetics import (
    BASE_RATES, MAX_SIZE, fit_catalogue, total_variation,
)
from single.statistical_mechanics.geometry import (
    octahedron, subdivided_tetrahedron, tetrahedron,
)


def skewness(probabilities: np.ndarray) -> float:
    sizes = np.arange(1, len(probabilities) + 1)
    mean = float(sizes @ probabilities)
    centred = sizes - mean
    variance = float((centred ** 2) @ probabilities)
    return float((centred ** 3) @ probabilities) / variance ** 1.5


def main() -> None:
    sizes, kinetic = oligomer_distribution(*BASE_RATES, max_size=MAX_SIZE)
    scaffolds = (tetrahedron(), subdivided_tetrahedron(), octahedron())
    colours = {"tetrahedron": "#b57931", "subdivided tetrahedron": "#8172b2",
               "octahedron": "#55a868"}
    fig, (ax_prob, ax_count) = plt.subplots(
        2, 1, figsize=(9, 7.5), constrained_layout=True)
    ax_prob.bar(sizes, kinetic, color="#b8d9ed", width=0.9,
                label="Kinetic target (20, 10, 1)")
    print(f"Kinetic mean: {sizes @ kinetic:.3f}; skewness: {skewness(kinetic):+.3f}")
    print("Parent                 capacity  P(target > capacity)  TV full  mode   skew")
    for scaffold in scaffolds:
        catalogue = enumerate_states(scaffold)
        capacity = catalogue.max_size
        parameters, probabilities = fit_catalogue(catalogue, kinetic[:capacity])
        padded = np.zeros(MAX_SIZE)
        padded[:capacity] = probabilities
        abundance_above = float(kinetic[capacity:].sum())
        colour = colours[scaffold.name]
        local_sizes = np.arange(1, capacity + 1)
        ax_prob.plot(local_sizes, probabilities, marker="o", markersize=2.8,
                     linewidth=1.6, color=colour, label=scaffold.name)
        by_size = Counter()
        for (n, _d, _c, _v), count in catalogue.counts.items():
            by_size[n] += count
        ax_count.plot(local_sizes, [by_size[n] for n in local_sizes],
                      marker="o", markersize=2.8, linewidth=1.6, color=colour,
                      label=scaffold.name)
        print(f"{scaffold.name:23} {capacity:2d}           "
              f"{abundance_above:.4f}            "
              f"{total_variation(padded, kinetic):.4f}     "
              f"{int(np.argmax(probabilities)) + 1:2d}    "
              f"{skewness(probabilities):+.3f}")
    ax_prob.set(xlim=(0, 31), ylabel="Fraction of oligomer molecules",
                title="Tetrahedral comparisons: six-edge tetrahedron is a 12-mer")
    ax_prob.legend(frameon=False, fontsize=8)
    ax_count.set(xlim=(0.5, 24.5), yscale="log",
                 xlabel="Oligomer size (monomers)",
                 ylabel="Patterns within each parent template")
    ax_count.set_xticks(np.arange(2, 25, 2))
    output = Path(__file__).with_name("tetrahedral_comparison.png")
    fig.savefig(output, dpi=180)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
