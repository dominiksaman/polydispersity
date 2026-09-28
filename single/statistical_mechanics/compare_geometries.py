"""Compare contact topology and variable-size ring weights with kinetics.

Run: python -m single.statistical_mechanics.compare_geometries
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from single.oligomer_distribution import oligomer_distribution
from single.statistical_mechanics.geometry import (
    cube, dimer_ring, octahedron, ring_distribution,
)


def main() -> None:
    scaffolds = (dimer_ring(12), cube(), octahedron())
    print("24-mer scaffold         vertices   degrees   dimer bonds   C-terminal contacts")
    for scaffold in scaffolds:
        print(f"{scaffold.name:23} {len(scaffold.vertex_degrees):2d}"
              f"       {set(scaffold.vertex_degrees)}"
              f"          {len(scaffold.dimer_contacts):2d}"
              f"              {len(scaffold.c_terminal_contacts):2d}")

    sizes, kinetic = oligomer_distribution(24, 12, 1, max_size=60)
    fig, ax = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    ax.bar(sizes, kinetic, width=0.9, color="#b8d9ed",
           label="Kinetic model (includes odd sizes)")
    for delta_g, color, description in (
        (-1.0, "#c44e52", "favoured growth"),
        (0.0, "#8172b2", "neutral growth"),
        (+1.0, "#55a868", "costly growth"),
    ):
        ring_sizes, probabilities = ring_distribution(delta_g, max_dimers=30)
        ax.plot(ring_sizes, probabilities, marker="o", markersize=2.5,
                linewidth=1.4, color=color,
                label=f"Closed ring: ΔG/dimer = {delta_g:+g} kJ/mol ({description})")
    ax.set(xlabel="Oligomer size (monomers)",
           ylabel="Fraction of oligomer molecules within each model",
           xlim=(0, 61), ylim=(0, 0.36))
    ax.set_xticks(np.arange(0, 61, 5))
    ax.legend(frameon=False, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.5, 1.0))
    ax.set_title("Constant contact energies in closed dimer rings do not select a finite size")
    output = Path(__file__).with_name("geometry_comparison.png")
    fig.savefig(output, dpi=180)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
