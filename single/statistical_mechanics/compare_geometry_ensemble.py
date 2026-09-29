"""Write conditional geometry shares for every modelled size and a figure.

Run: python -m single.statistical_mechanics.compare_geometry_ensemble
"""

from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .geometry_ensemble import (EnergyParameters, conditional_contributions,
                                default_candidates)


def main() -> None:
    output_dir = Path(__file__).parent
    candidates = default_candidates(40)
    scenarios = (0.0, 0.1, 0.3)
    csv_path = output_dir / "geometry_contributions.csv"
    rows = []
    for kappa in scenarios:
        params = EnergyParameters(epsilon_dimer=1.0, epsilon_c_terminal=1.0,
                                  degree_penalty=kappa, preferred_degree=4)
        for size in range(5, 41):
            for contribution in conditional_contributions(
                    size, params, candidates=candidates):
                rows.append({"degree_penalty_kbt": kappa, **asdict(contribution)})
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    fig, (ax_all, ax_24) = plt.subplots(2, 1, figsize=(11, 8.2),
                                      layout="constrained")
    family_order = ("ring", "degree-3 polyhedron", "mixed-degree polyhedron",
                    "degree-4 polyhedron")
    colors = ("#6e8297", "#e7a36c", "#b49acb", "#58aa91")
    sizes = np.arange(5, 41)
    params = EnergyParameters(1.0, 1.0, 0.1, 4)
    family_shares = np.zeros((len(family_order), len(sizes)))
    for j, size in enumerate(sizes):
        for item in conditional_contributions(int(size), params, candidates=candidates):
            family_shares[family_order.index(item.family), j] += item.share
    bottoms = np.zeros(len(sizes))
    for family, color, shares in zip(family_order, colors, family_shares):
        ax_all.bar(sizes, shares, bottom=bottoms, width=0.88,
                   label=family, color=color)
        bottoms += shares
    ax_all.set(xlim=(5, 40), ylim=(0, 1), ylabel="Conditional share within size",
               title="Illustrative geometry ensemble, degree penalty = 0.1 kBT")
    ax_all.set_xticks(np.arange(6, 41, 2))
    ax_all.legend(loc="lower left", ncol=2, fontsize=8, frameon=False)
    ax_all.grid(axis="y", alpha=0.18)

    penalty_grid = np.linspace(0, 0.5, 101)
    geometry_names = ("single ring (12 dimers)", "cube", "octahedron")
    geometry_colors = (colors[0], colors[1], colors[3])
    for name, color in zip(geometry_names, geometry_colors):
        shares = []
        for kappa in penalty_grid:
            ensemble = conditional_contributions(
                24, EnergyParameters(1.0, 1.0, float(kappa), 4),
                candidates=candidates)
            shares.append(next(item.share for item in ensemble if item.geometry == name))
        ax_24.plot(penalty_grid, shares, label=name, color=color, linewidth=2)
    ax_24.set(xlim=(0, 0.5), ylim=(0, 1), xlabel="Degree mismatch penalty (kBT per squared degree)",
              ylabel="P(geometry | 24-mer)",
              title="The 24-mer shares depend on an unknown shape free energy")
    ax_24.legend(frameon=False, ncol=3, fontsize=9)
    ax_24.grid(alpha=0.18)
    image_path = output_dir / "geometry_contributions.png"
    fig.savefig(image_path, dpi=180)
    plt.close(fig)

    print(f"Wrote {csv_path}")
    print(f"Wrote {image_path}")
    print("\n24-mer: candidate shares at each degree penalty")
    for kappa in scenarios:
        ensemble = conditional_contributions(
            24, EnergyParameters(1.0, 1.0, kappa, 4), candidates=candidates)
        shares = ", ".join(f"{item.geometry}: {item.share:.3f}" for item in ensemble)
        print(f"  kappa={kappa:.1f}: {shares}")


if __name__ == "__main__":
    main()
