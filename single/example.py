#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""One reference figure. Run from the repository root: python -m single.example."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .mass_distribution import mass_distribution
from .oligomer_distribution import oligomer_distribution
from .statistical_mechanics import energies_from_rates, site_count_distribution


def main() -> None:
    rates = (20.0, 10.0, 1.0)  # effective on, unpaired off, paired off
    sizes, p = oligomer_distribution(*rates, max_size=60)
    energies = energies_from_rates(*rates)
    closed_form = site_count_distribution(
        energies.edge_kj_mol, energies.dimer_kj_mol, max_size=60)
    masses, mass_p = mass_distribution(*rates, monomer_mass=20000, max_size=60)

    fig, (ax_size, ax_mass) = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    colors = np.where(sizes % 2 == 0, "#2c7fb8", "#d95f0e")
    ax_size.bar(sizes, p, color=colors, width=0.88, label="Rate-model equilibrium")
    ax_size.plot(sizes, closed_form, color="#263238", linewidth=1.1,
                 label="Same model, closed-form weights")
    ax_size.set(xlim=(0.5, 32.5), xlabel="Monomers per oligomer",
                ylabel="Oligomer number fraction")
    ax_size.legend(frameon=False, fontsize=8)
    ax_mass.bar(masses / 1000, mass_p, color=colors, width=17.6)
    ax_mass.set(xlim=(10, 650), xlabel="Neutral mass (kDa; 20 kDa per monomer)",
                ylabel="Oligomer number fraction")
    fig.suptitle("Effective on = 20; unpaired off = 10; paired off = 1")
    path = Path(__file__).with_name("size_distribution.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"Wrote {path}")
    print(f"Mode {sizes[np.argmax(p)]}; mean {sizes @ p:.6f}; "
          f"P(n > 24) {p[sizes > 24].sum():.6f}; "
          f"closed-form maximum error {np.max(np.abs(p - closed_form)):.2e}")


if __name__ == "__main__":
    main()
