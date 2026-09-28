"""Plot statistical weights against the existing single-protein kinetics.

Run from the repository root:

    python -m single.statistical_mechanics.compare_models
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from single.oligomer_distribution import oligomer_distribution
from single.statistical_mechanics import (
    additive_contact_distribution,
    energies_from_rates,
    fit_ideal_cluster,
    site_count_distribution,
)


def main() -> None:
    rates = (24.0, 12.0, 1.0)
    sizes, kinetic = oligomer_distribution(*rates, max_size=60)
    energies = energies_from_rates(*rates)
    site_count = site_count_distribution(
        energies.edge_kj_mol, energies.dimer_kj_mol, max_size=60
    )
    contacts_only = additive_contact_distribution(
        energies.edge_kj_mol, energies.dimer_kj_mol, max_size=60
    )
    factorial = additive_contact_distribution(
        energies.edge_kj_mol, energies.dimer_kj_mol,
        max_size=60, factorial_entropy=True,
    )
    fitted = fit_ideal_cluster(kinetic)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    comparisons = (
        (site_count, "Site-count partition weight", "Exact match"),
        (contacts_only, "Additive contacts only", "Peaks at the size cutoff"),
        (factorial, "Contacts + 1/n! (rate-derived energies)",
         "Wrong size and parity bias"),
        (fitted.distribution, "Contacts + 1/n! (refitted energies)",
         f"Similar shape; TV distance = {fitted.total_variation:.3f}"),
    )
    for ax, (other, title, subtitle) in zip(axes.flat, comparisons):
        ax.bar(sizes, kinetic, color="#b8d9ed", width=0.9,
               label="Kinetic solution")
        ax.plot(sizes, other, color="#d95f0e", linewidth=1.8,
                marker="o", markersize=2.5, label="Statistical weight")
        ax.set_title(title, fontsize=11)
        ax.text(0.98, 0.95, subtitle, ha="right", va="top",
                transform=ax.transAxes, fontsize=8)
        ax.set(xlim=(0, 61), xlabel="Oligomer size (monomers)",
               ylabel="Fraction of oligomers")
        ax.set_xticks(np.arange(0, 61, 10))
    axes[0, 0].legend(frameon=False, loc="upper left", fontsize=8)
    fig.suptitle("Single-protein assembly: 24 effective on, 12 unpaired off, 1 paired off")
    output = Path(__file__).with_name("size_distribution_comparison.png")
    fig.savefig(output, dpi=180)
    print(f"Wrote {output}")
    print(f"Site-count max difference: {np.max(np.abs(site_count - kinetic)):.2e}")
    print(f"Refitted factorial TV distance: {fitted.total_variation:.4f}")
    print(f"Rate-derived energies (edge, dimer): {energies.edge_kj_mol:.3f}, "
          f"{energies.dimer_kj_mol:.3f} kJ/mol")
    print(f"Descriptive refit (edge, dimer): {fitted.edge_kj_mol:.3f}, "
          f"{fitted.dimer_kj_mol:.3f} kJ/mol")


if __name__ == "__main__":
    main()
