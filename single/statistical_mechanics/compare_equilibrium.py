"""Plot a rate-free scaffold partition function beside the kinetic reference.

Run from the repository root:

    python -m single.statistical_mechanics.compare_equilibrium
"""

from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from single.oligomer_distribution import oligomer_distribution
from single.statistical_mechanics.equilibrium_scaffold import (
    EquilibriumParameters, enumerate_states,
)
from single.statistical_mechanics.geometry import cube, dimer_ring, octahedron


def main() -> None:
    # These dimensionless values illustrate a finite-size polyhedral peak.
    # They are not measured interface energies or a fit to native MS data.
    parameters = EquilibriumParameters(
        log_monomer_activity=-1.0,
        dimer_stabilization=4.0,
        c_terminal_stabilization=0.05,
        saturated_vertex_penalty=1.1,
    )
    colours = {"12-dimer ring": "#c44e52", "cube": "#8172b2",
               "octahedron": "#55a868"}
    catalogues = [enumerate_states(scaffold)
                  for scaffold in (dimer_ring(12), cube(), octahedron())]
    sizes, kinetic = oligomer_distribution(24, 12, 1, max_size=24)

    fig, (ax_count, ax_poly, ax_ring) = plt.subplots(
        3, 1, figsize=(9, 9), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [1.2, 1.2, 0.65]})
    for catalogue in catalogues:
        multiplicity = Counter()
        for (n, _d, _c, _v), count in catalogue.counts.items():
            multiplicity[n] += count
        count_by_size = np.array([multiplicity[n] for n in sizes])
        probability = catalogue.distribution(parameters)
        colour = colours[catalogue.scaffold_name]
        ax_count.plot(sizes, count_by_size, marker="o", markersize=3,
                      linewidth=1.5, label=catalogue.scaffold_name,
                      color=colour)
        ax = ax_ring if catalogue.scaffold_name.endswith("ring") else ax_poly
        ax.plot(sizes, probability, marker="o", markersize=3,
                linewidth=1.5, label=catalogue.scaffold_name, color=colour)
        print(f"{catalogue.scaffold_name}: {sum(multiplicity.values())} "
              f"occupied parent-template patterns; illustrative mode "
              f"{int(sizes[np.argmax(probability)])}-mer")
    ax_count.set(yscale="log", ylabel="Patterns within each parent template",
                 title="Restricted scaffold model: occupied-template counts and size weights")
    ax_count.legend(frameon=False, ncol=3, fontsize=9)
    ax_poly.bar(sizes, kinetic, width=0.8, color="#b8d9ed", alpha=0.65,
                label="Kinetic reference, conditional on sizes ≤24")
    ax_poly.set(ylabel="Polyhedral probability", ylim=(0, 0.25))
    ax_poly.legend(frameon=False, fontsize=8)
    ax_ring.set(xlabel="Oligomer size (monomers)",
                ylabel="Ring probability", xlim=(0.5, 24.5), ylim=(0, 0.92))
    ax_ring.set_xticks(np.arange(1, 25, 2))
    ax_ring.legend(frameon=False, fontsize=8, loc="upper left")
    output = Path(__file__).with_name("equilibrium_scaffold_comparison.png")
    fig.savefig(output, dpi=180)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
