#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mass distribution of a single self-assembling protein.

This builds directly on :mod:`oligomer_distribution`: it takes the same three
rate constants, computes the equilibrium oligomer-size distribution, and then
converts the size axis (number of monomers) into a mass axis using the mass of a
single monomer.  The result is what a native mass-spectrometry or mass-photometry
experiment effectively measures: the abundance of each oligomer plotted against
its molecular mass.

Because every oligomer of size ``i`` has mass ``i * monomer_mass``, the
conversion is just a relabelling of the x-axis -- the abundances are unchanged.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from oligomer_distribution import oligomer_distribution, _save_csv


def mass_distribution(
    k_on: float,
    k_off_monomer: float,
    k_off_dimer: float,
    monomer_mass: float,
    max_size: int = 60,
    normalise: bool = True,
    csv_path: Optional[str] = None,
    plot: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the equilibrium mass distribution of a self-assembling protein.

    Parameters
    ----------
    k_on, k_off_monomer, k_off_dimer : float
        The three rate constants, exactly as in
        :func:`oligomer_distribution.oligomer_distribution` (see that function
        for what each one controls).
    monomer_mass : float
        Mass of a single monomer, in whatever units you like (Da, kDa, ...).
        The output masses inherit these units.
    max_size : int, optional
        Largest oligomer size to include.  Default 60.
    normalise : bool, optional
        If True (default) the abundances sum to 1.
    csv_path : str, optional
        If given, write a two-column CSV ``mass,abundance`` to this path.
    plot : bool, optional
        If True, draw and show a mass-vs-abundance plot.

    Returns
    -------
    masses : ndarray of float, shape (max_size,)
        Oligomer masses ``monomer_mass, 2*monomer_mass, ...``.
    abundance : ndarray of float, shape (max_size,)
        Relative abundance of each mass (identical to the size-distribution
        abundances; only the x-axis has changed).

    Raises
    ------
    ValueError
        If ``monomer_mass`` is non-positive (other inputs are checked by
        :func:`oligomer_distribution`).
    """
    if monomer_mass <= 0:
        raise ValueError("monomer_mass must be positive.")

    # Reuse the size distribution, then map size -> mass.
    sizes, abundance = oligomer_distribution(
        k_on, k_off_monomer, k_off_dimer,
        max_size=max_size, normalise=normalise,
    )
    masses = sizes * float(monomer_mass)

    if csv_path is not None:
        _save_csv(masses, abundance, csv_path, header="mass,abundance")

    if plot:
        plot_mass_distribution(masses, abundance, monomer_mass)

    return masses, abundance


def plot_mass_distribution(
    masses: np.ndarray,
    abundance: np.ndarray,
    monomer_mass: Optional[float] = None,
    mass_unit: str = "Da",
    ax=None,
    show: bool = True,
):
    """Plot a mass distribution (mass vs. abundance) as a stem/peak spectrum.

    Returns the matplotlib ``Axes`` so the caller can customise or save it.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))

    # A stem plot reads like a mass spectrum: one peak per oligomer.
    markerline, stemlines, baseline = ax.stem(masses, abundance, basefmt=" ")
    plt.setp(stemlines, color="#2c7fb8", linewidth=1.2)
    plt.setp(markerline, color="#2c7fb8", markersize=3)

    ax.set_xlabel(f"mass ({mass_unit})")
    ax.set_ylabel("relative abundance")
    title = "Equilibrium mass distribution"
    if monomer_mass is not None:
        title += f"\nmonomer mass = {monomer_mass:g} {mass_unit}"
    ax.set_title(title)
    ax.margins(x=0.01)

    if show:
        plt.tight_layout()
        plt.show()
    return ax


if __name__ == "__main__":
    # Example: same even-biased ~24-mer as in oligomer_distribution, with a
    # 20 kDa monomer (so the most abundant species sits near 480 kDa).
    masses, abundance = mass_distribution(
        k_on=24.0, k_off_monomer=12.0, k_off_dimer=1.0,
        monomer_mass=20_000.0, max_size=60,
    )
    peak_mass = masses[np.argmax(abundance)]
    print(f"most abundant mass: {peak_mass:,.0f} Da "
          f"({peak_mass / 20_000:.0f}-mer)")
    # Uncomment to save the spectrum to CSV and/or pop up a plot:
    # mass_distribution(24.0, 12.0, 1.0, 20_000.0, csv_path="mass.csv", plot=True)
