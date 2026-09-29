"""Export the distinct monomer-at-vertex geometry hypothesis.

Run: python -m single.statistical_mechanics.compare_monomer_geometries
"""

from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .monomer_geometry_ensemble import (MonomerEnergy, monomer_candidates,
                                        monomer_contributions)


def main() -> None:
    directory = Path(__file__).parent
    candidates = monomer_candidates(12)
    penalties = (0.0, 0.1, 0.3)
    rows = []
    for penalty in penalties:
        for size in range(1, 13):
            for item in monomer_contributions(
                    size, MonomerEnergy(1.0, 1.0, penalty), candidates=candidates):
                rows.append({"shape_penalty_kbt": penalty, **asdict(item)})
    csv_path = directory / "monomer_geometry_contributions.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    fig, (ax_all, ax_eleven) = plt.subplots(2, 1, figsize=(10.5, 7.8),
                                       layout="constrained")
    families = ("open chain", "closed ring", "compact polyhedron")
    colors = ("#718da6", "#dc9971", "#72ad96")
    sizes = np.arange(1, 13)
    family_shares = np.zeros((3, len(sizes)))
    for index, size in enumerate(sizes):
        for item in monomer_contributions(
                int(size), MonomerEnergy(1, 1, 0.1), candidates=candidates):
            family_shares[families.index(item.family), index] += item.share
    bottoms = np.zeros(len(sizes))
    for family, color, shares in zip(families, colors, family_shares):
        ax_all.bar(sizes, shares, bottom=bottoms, width=0.82,
                   color=color, label=family)
        bottoms += shares
    ax_all.set(xlim=(0.5, 12.5), ylim=(0, 1), xticks=sizes,
               ylabel="Conditional share within size",
               title="Monomer-at-vertex examples (shape penalty 0.1 kBT)")
    ax_all.legend(frameon=False, ncol=3, loc="lower left")
    ax_all.grid(axis="y", alpha=0.17)

    penalty_grid = np.linspace(0, 0.15, 121)
    eleven_names = [item.name for item in candidates if item.monomers == 11]
    line_colors = ("#718da6", "#dc9971", "#ae88bf", "#72ad96")
    for name, color in zip(eleven_names, line_colors):
        shares = []
        for penalty in penalty_grid:
            values = monomer_contributions(
                11, MonomerEnergy(1, 1, float(penalty)), candidates=candidates)
            shares.append(next(item.share for item in values if item.geometry == name))
        ax_eleven.plot(penalty_grid, shares, linewidth=2, label=name, color=color)
    ax_eleven.set(xlim=(0, 0.15), ylim=(0, 1),
               xlabel="Assumed compact-shape penalty (kBT per squared excess degree)",
               ylabel="P(geometry | 11-mer)",
               title="Two distinct eleven-monomer polyhedral candidates")
    ax_eleven.legend(frameon=False, ncol=2, fontsize=8)
    ax_eleven.grid(alpha=0.17)
    image_path = directory / "monomer_geometry_contributions.png"
    fig.savefig(image_path, dpi=180)
    plt.close(fig)
    print(f"Wrote {csv_path}")
    print(f"Wrote {image_path}")
    for size in (3, 5, 6, 8, 11, 12):
        values = monomer_contributions(
            size, MonomerEnergy(1, 1, 0.1), candidates=candidates)
        print(f"{size}-mer: " + ", ".join(
            f"{item.geometry} {item.share:.3f}" for item in values))


if __name__ == "__main__":
    main()
