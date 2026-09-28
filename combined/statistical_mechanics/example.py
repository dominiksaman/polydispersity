"""Plot fixed-size circular co-assembly for three contact-energy biases.

Run from the repository root: ``python -m combined.statistical_mechanics.example``.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .ring_coassembly import composition_distribution


def main() -> None:
    size = 12
    temperature = 298.15
    figure, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, fraction in zip(axes, (0.5, 1/3)):
        for energy, label, color in (
            (4.0, "self contacts favoured", "#6a51a3"),
            (0.0, "unbiased", "#3182bd"),
            (-4.0, "AB contacts favoured", "#e6550d"),
        ):
            distribution = composition_distribution(
                size, energy, temperature=temperature, mole_fraction_a=fraction
            )
            ax.plot(np.arange(size + 1), distribution, marker="o", markersize=3,
                    label=label, color=color)
        ax.set_title(f"Mean A fraction = {fraction:.2f}")
        ax.set_xlabel("Number of A subunits in a 12-mer")
        ax.set_xticks(range(size + 1, 2))
    axes[0].set_ylabel("Fraction of 12-mers")
    axes[0].legend(frameon=False, fontsize=8)
    figure.suptitle("Circular A/B co-assembly at 298.15 K")
    figure.tight_layout()
    output = Path(__file__).with_name("example_ring.png")
    figure.savefig(output, dpi=160)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
