"""Recreate the thesis's three-panel 12-mer counting-model figure.

Reference: Ver3-2/text/ch1-intro.tex, lines 206-230, and the figure at
Ver3-2/figures/modelling/counting_model.eps. Run from the repository root:

    python -m combined.statistical_mechanics.replicate_counting_model

The original figure has no numerical source table. We reproduce the model
curves from its caption and the thesis equations. Its horizontal energy
convention is negative for self-assembly and positive for AB-favoured
co-assembly, so x = (g_homo - g_AB) / (RT) = -g_AB/(RT) when g_AA=g_BB=0.
"""

from math import comb
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .ring_coassembly import R_KJ_MOL_K, composition_distribution


def probabilities(energy_difference_kt: float, size: int = 12) -> np.ndarray:
    """Composition fractions in the thesis's energy-sign convention."""
    temperature = 298.15
    g_ab = -energy_difference_kt * R_KJ_MOL_K * temperature
    return composition_distribution(
        size, g_ab, temperature=temperature, mole_fraction_a=0.5
    )


def main() -> None:
    size = 12
    compositions = np.arange(size + 1)
    integer_energies = np.arange(-5, 6)
    integer_distributions = np.array([probabilities(e) for e in integer_energies])
    dense_energies = np.linspace(-5, 5, 201)
    heatmap = np.array([probabilities(e) for e in dense_energies]).T

    fig = plt.figure(figsize=(12.2, 4.0), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=(1, 1.15, 1.32))
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])

    # Panel a in the thesis is scaled to its modal height of one.
    unbiased = probabilities(0)
    binomial = np.array([comb(size, i) / 2**size for i in compositions])
    ax_a.bar(compositions, unbiased / unbiased.max(), color="#1976b8",
             width=0.86, label="Ring model")
    ax_a.plot(compositions, binomial / binomial.max(), color="#d9252b",
              lw=2.0, label="Binomial")
    ax_a.set(xlim=(-0.55, 12.55), ylim=(0, 1.07),
             xlabel="Number of A subunits", ylabel="Relative abundance")
    ax_a.set_xticks(np.arange(0, 13, 2))
    ax_a.legend(frameon=False, loc="upper right", fontsize=8)
    ax_a.set_title("a  Unbiased 12-mer", loc="left", fontweight="bold")

    # The two homomer endpoints are pooled, as required for their combined
    # abundance to approach one in the strongly self-assembling limit.
    homomer = integer_distributions[:, 0] + integer_distributions[:, -1]
    obligate = integer_distributions[:, size // 2]
    ax_b.bar(integer_energies, homomer, width=0.84,
             color="#1976b8", label="Pure 12:0 + 0:12")
    ax_b.bar(integer_energies, obligate, width=0.84,
             color="#f49b21", label="6:6")
    ax_b.set(xlim=(-5.6, 5.6), ylim=(0, 1.07),
             xlabel="Energy difference per contact (kT)",
             ylabel="Fraction of 12-mers")
    ax_b.set_xticks(np.arange(-5, 6, 2))
    ax_b.legend(frameon=False, loc="upper center", fontsize=8)
    ax_b.set_title("b  Assembly extremes", loc="left", fontweight="bold")

    image = ax_c.imshow(
        heatmap, origin="lower", aspect="auto", cmap="PuRd", vmin=0, vmax=1,
        extent=(-5, 5, -0.5, 12.5), interpolation="nearest",
    )
    ax_c.set(xlabel="Energy difference per contact (kT)",
             ylabel="Number of A subunits")
    ax_c.set_xticks(np.arange(-5, 6, 2))
    ax_c.set_yticks(np.arange(0, 13, 2))
    ax_c.set_title("c  All compositions", loc="left", fontweight="bold")
    fig.colorbar(image, ax=ax_c, fraction=0.046, pad=0.03,
                 label="Fraction of 12-mers")

    output = Path(__file__).with_name("counting_model_replication.png")
    fig.savefig(output, dpi=200)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
