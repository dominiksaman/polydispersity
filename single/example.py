#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runnable example: produce the oligomer-size and mass distributions and save the
two figures shown in the README.

Run with:  python example.py
"""

import matplotlib.pyplot as plt

from oligomer_distribution import oligomer_distribution, plot_oligomer_distribution
from mass_distribution import mass_distribution, plot_mass_distribution

# Three rate constants for a polydisperse, even-biased protein:
#   k_on / k_off_dimer = 24  -> peaks near a ~20-mer
#   k_off_monomer (12) > k_off_dimer (1) -> even sizes preferred
K_ON, K_OFF_MONOMER, K_OFF_DIMER = 24.0, 12.0, 1.0
MONOMER_MASS = 20_000.0  # Da

# ---- oligomer-size distribution -----------------------------------------
sizes, abundance = oligomer_distribution(K_ON, K_OFF_MONOMER, K_OFF_DIMER,
                                         max_size=60)
fig, ax = plt.subplots(figsize=(7, 4))
plot_oligomer_distribution(sizes, abundance, K_ON, K_OFF_MONOMER, K_OFF_DIMER,
                           ax=ax, show=False)
fig.tight_layout()
fig.savefig("example_oligomer.png", dpi=120)

# ---- mass distribution ---------------------------------------------------
masses, abundance = mass_distribution(K_ON, K_OFF_MONOMER, K_OFF_DIMER,
                                      MONOMER_MASS, max_size=60)
fig, ax = plt.subplots(figsize=(7, 4))
plot_mass_distribution(masses, abundance, MONOMER_MASS, ax=ax, show=False)
fig.tight_layout()
fig.savefig("example_mass.png", dpi=120)

print("wrote example_oligomer.png and example_mass.png")
