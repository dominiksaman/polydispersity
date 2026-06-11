#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Equilibrium oligomer-size distribution for a single self-assembling protein.

Many proteins -- small heat-shock proteins (sHsps) being a classic example --
are *polydisperse*: rather than forming a single well-defined complex, they
populate a whole range of oligomer sizes at once.  This module computes the
equilibrium size distribution of such a protein from just three rate constants,
using the "modified helical polymerisation" model of Baldwin et al.

The model
---------
The protein assembles and disassembles one monomer at a time,

    P_1 + P_1  <->  P_2
    P_2 + P_1  <->  P_3
    ...
    P_(i-1) + P_1  <->  P_i

where ``P_i`` is an oligomer of ``i`` monomers and ``P_1`` is a free monomer.
At equilibrium the abundance of each successive size is fixed by the ratio of
the on- and off-rates for that step,

    [P_i] / [P_(i-1)]  =  k_on / (effective off-rate of P_i).

The one subtlety is that monomers prefer to pair up into dimers, so an oligomer
with an *even* number of monomers (all paired) is more stable than one with an
*odd* number (one unpaired "spare" monomer).  This is captured by using two
different off-rates, giving:

    even i :  [P_i] / [P_(i-1)]  =  k_on / ( i * k_off_dimer )
    odd  i :  [P_i] / [P_(i-1)]  =  k_on / ( (i-1) * k_off_dimer + k_off_monomer )

Starting from the monomer (``[P_1] = 1``) and applying this recursion gives the
full relative distribution, which we then normalise.

What each rate does (intuitively)
---------------------------------
* ``k_on``          -- how readily a monomer adds on.  Sets the OVERALL SIZE:
                       the distribution peaks near ``k_on / k_off_dimer``.
                       Bigger ``k_on`` -> bigger oligomers.
* ``k_off_dimer``   -- how readily a monomer leaves a fully-paired (even) part
                       of the oligomer.  Breaking a dimer is hard, so this is
                       usually the SMALLEST rate.  Bigger -> smaller oligomers.
* ``k_off_monomer`` -- how readily the single unpaired monomer leaves an
                       odd-sized oligomer.  Sets the EVEN:ODD BIAS: when it is
                       much larger than ``k_off_dimer`` the spare monomer falls
                       off quickly, so even sizes dominate.

Note on ``k_on``
----------------
Physically the on-step is bimolecular (oligomer + free monomer), so the true
rate is ``k_on_true * [P_1]``.  Here we fold the free-monomer concentration into
a single effective ``k_on`` (a "pseudo-first-order" rate), so the three rate
constants fully determine the *shape* of the distribution.  Only their ratios
matter, so e.g. (k_on, k_off_monomer, k_off_dimer) = (24, 5, 1) and
(48, 10, 2) give the same distribution.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


def oligomer_distribution(
    k_on: float,
    k_off_monomer: float,
    k_off_dimer: float,
    max_size: int = 60,
    normalise: bool = True,
    csv_path: Optional[str] = None,
    plot: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the equilibrium oligomer-size distribution.

    Parameters
    ----------
    k_on : float
        Effective monomer association rate ``k+`` (controls the overall size;
        the distribution peaks near ``k_on / k_off_dimer``).
    k_off_monomer : float
        Off-rate of the single unpaired monomer from an odd-sized oligomer
        (``k-_m``); larger values favour even sizes.
    k_off_dimer : float
        Off-rate of a monomer from a dimer in an even-sized oligomer (``k-_d``);
        usually the smallest of the three rates.
    max_size : int, optional
        Largest oligomer size (number of monomers) to compute.  Default 60.
    normalise : bool, optional
        If True (default) the abundances sum to 1 (a probability distribution
        over sizes).  If False, abundances are relative to the monomer
        (``[P_1] = 1``).
    csv_path : str, optional
        If given, write a two-column CSV ``oligomer_size,abundance`` to this path.
    plot : bool, optional
        If True, draw and show a bar plot of the distribution (even and odd
        sizes coloured separately to make the even:odd bias visible).

    Returns
    -------
    sizes : ndarray of int, shape (max_size,)
        Oligomer sizes ``1, 2, ..., max_size`` (number of monomers).
    abundance : ndarray of float, shape (max_size,)
        Relative abundance of each size.

    Raises
    ------
    ValueError
        If any rate is non-positive or ``max_size < 1``.
    """
    # ---- validate inputs -------------------------------------------------
    if min(k_on, k_off_monomer, k_off_dimer) <= 0:
        raise ValueError("All rate constants must be positive.")
    if max_size < 1:
        raise ValueError("max_size must be at least 1.")

    # ---- build the distribution by recursion ----------------------------
    sizes = np.arange(1, max_size + 1)
    abundance = np.empty(max_size, dtype=float)
    abundance[0] = 1.0  # [P_1]: the monomer is the reference point

    for i in range(2, max_size + 1):
        if i % 2 == 0:
            # even oligomer: every monomer is paired -> dimer off-rate, with a
            # statistical factor i for the number of equivalent positions.
            effective_off_rate = i * k_off_dimer
        else:
            # odd oligomer: one spare monomer (k_off_monomer) plus the paired
            # part (the (i-1) already-paired monomers leave at k_off_dimer).
            effective_off_rate = (i - 1) * k_off_dimer + k_off_monomer
        # abundance[i-1] is [P_i]; abundance[i-2] is [P_(i-1)]
        abundance[i - 1] = abundance[i - 2] * k_on / effective_off_rate

    # ---- normalise -------------------------------------------------------
    if normalise:
        abundance = abundance / abundance.sum()

    # ---- optional outputs ------------------------------------------------
    if csv_path is not None:
        _save_csv(sizes, abundance, csv_path,
                  header="oligomer_size,abundance")

    if plot:
        plot_oligomer_distribution(sizes, abundance,
                                   k_on, k_off_monomer, k_off_dimer)

    return sizes, abundance


def _save_csv(x: np.ndarray, y: np.ndarray, path: str, header: str) -> None:
    """Write two columns ``x, y`` to ``path`` as CSV with a header line."""
    data = np.column_stack([x, y])
    np.savetxt(path, data, delimiter=",", header=header, comments="",
               fmt=["%d", "%.8e"] if np.issubdtype(x.dtype, np.integer)
               else "%.8e")


def plot_oligomer_distribution(
    sizes: np.ndarray,
    abundance: np.ndarray,
    k_on: Optional[float] = None,
    k_off_monomer: Optional[float] = None,
    k_off_dimer: Optional[float] = None,
    ax=None,
    show: bool = True,
):
    """Bar plot of an oligomer-size distribution.

    Even and odd sizes are coloured differently so the even:odd preference is
    immediately visible.  Returns the matplotlib ``Axes`` so the caller can
    customise or save the figure.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4))

    even = sizes % 2 == 0
    ax.bar(sizes[even], abundance[even], width=0.9,
           color="#2c7fb8", label="even (all monomers paired)")
    ax.bar(sizes[~even], abundance[~even], width=0.9,
           color="#d95f0e", label="odd (one unpaired monomer)")

    ax.set_xlabel("oligomer size (number of monomers)")
    ax.set_ylabel("relative abundance")
    title = "Equilibrium oligomer-size distribution"
    if None not in (k_on, k_off_monomer, k_off_dimer):
        title += (f"\n$k_{{on}}$={k_on:g}, $k^-_m$={k_off_monomer:g}, "
                  f"$k^-_d$={k_off_dimer:g}")
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.margins(x=0.01)

    if show:
        plt.tight_layout()
        plt.show()
    return ax


if __name__ == "__main__":
    # Example: a polydisperse, even-biased oligomer that peaks near a 24-mer.
    # k_on/k_off_dimer = 24 sets the size; k_off_monomer > k_off_dimer makes
    # even sizes dominate.
    sizes, abundance = oligomer_distribution(
        k_on=24.0, k_off_monomer=12.0, k_off_dimer=1.0, max_size=60,
    )
    peak = sizes[np.argmax(abundance)]
    even_odd = abundance[sizes % 2 == 0].sum() / abundance[sizes % 2 == 1].sum()
    print(f"most abundant size: {peak}-mer")
    print(f"even:odd abundance ratio: {even_odd:.2f}:1")
    # Uncomment to save the distribution to CSV and/or pop up a plot:
    # oligomer_distribution(24.0, 12.0, 1.0, csv_path="distribution.csv", plot=True)
