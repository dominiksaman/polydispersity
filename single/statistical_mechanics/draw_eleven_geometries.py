"""Draw the two illustrative 11-monomer polyhedral contact graphs.

Run: python -m single.statistical_mechanics.draw_eleven_geometries
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from .geometry_catalogue import monomer_candidates
from .geometry_weights import BondPattern, fully_bound_pattern


def pentagon(radius: float, z: float, rotation: float = 0) -> np.ndarray:
    angles = rotation + 2 * np.pi * np.arange(5) / 5
    return np.column_stack((radius * np.cos(angles),
                            radius * np.sin(angles), np.full(5, z)))


def capped_prism_coordinates() -> np.ndarray:
    return np.vstack((pentagon(1, 0.35), pentagon(1, -0.9),
                      np.array([[0, 0, 1.5]])))


def pentagonal_face_coordinates() -> np.ndarray:
    # Remaining vertices of a regular icosahedron after removing its apex.
    radius = 2 / np.sqrt(5)
    height = 1 / np.sqrt(5)
    return np.vstack((pentagon(radius, height),
                      pentagon(radius, -height, np.pi / 5),
                      np.array([[0, 0, -1]])))


def triangle_faces(edges: tuple[tuple[int, int], ...], size: int) -> list[tuple[int, ...]]:
    edge_set = {frozenset(edge) for edge in edges}
    return [vertices for vertices in combinations(range(size), 3)
            if all(frozenset(pair) in edge_set for pair in combinations(vertices, 2))]


def draw_candidate(ax, coordinates: np.ndarray, edges: tuple[tuple[int, int], ...],
                   faces: list[tuple[int, ...]], highlighted: list[tuple[int, ...]],
                   title: str, bonds: BondPattern) -> None:
    shell = Poly3DCollection([[coordinates[i] for i in face] for face in faces],
                             facecolors="#83abc2", edgecolors="none", alpha=0.06)
    ax.add_collection3d(shell)
    highlight = Poly3DCollection(
        [[coordinates[i] for i in face] for face in highlighted],
        facecolors="#f2b664", edgecolors="none", alpha=0.18)
    ax.add_collection3d(highlight)
    for a, b in edges:
        a_xyz, b_xyz = coordinates[[a, b]]
        ax.plot((a_xyz[0], b_xyz[0]), (a_xyz[1], b_xyz[1]),
                (a_xyz[2], b_xyz[2]), color="#8c99a3", linewidth=0.65,
                linestyle="--", alpha=0.6)
    for a, b in bonds.dimer_pairs:
        segment = coordinates[[a, b]]
        ax.plot(*segment.T, color="#2166ac", linewidth=3.0)
    for a, b in bonds.c_terminal_contacts:
        segment = coordinates[[a, b]]
        ax.plot(*segment.T, color="#d96818", linewidth=1.8)
        delta = coordinates[b] - coordinates[a]
        start = coordinates[a] + 0.37 * delta
        ax.quiver(*start, *(0.25 * delta), color="#d96818", arrow_length_ratio=0.5,
                  linewidth=1.0)
    ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2],
               s=57, color="#143d52", edgecolor="white", linewidth=0.65,
               depthshade=False)
    ax.set(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), zlim=(-1.55, 1.85),
           title=f"{title}\n11 monomers; 5 dimers; 11 C-terminal contacts")
    ax.set_box_aspect((1, 1, 1.1))
    ax.set_axis_off()
    ax.view_init(elev=19, azim=-61)


def main() -> None:
    candidates = {item.name: item for item in monomer_candidates(12)}
    prism = candidates["capped pentagonal prism (11 monomers)"]
    ico = candidates["pentagonal-face icosahedral polyhedron (11 monomers)"]
    prism_faces = [tuple(range(5, 10))]
    prism_faces += [(i, (i + 1) % 5, (i + 1) % 5 + 5, i + 5)
                    for i in range(5)]
    prism_cap = [(i, (i + 1) % 5, 10) for i in range(5)]
    ico_cap = [tuple(range(5))]
    ico_faces = triangle_faces(ico.edges, 11)
    if len(prism_faces) + len(prism_cap) != 11 or len(ico_faces) != 15:
        raise ValueError("unexpected face count in 11-vertex graph")

    fig = plt.figure(figsize=(12, 6.6), facecolor="white")
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.11, top=0.82, wspace=0.02)
    first = fig.add_subplot(1, 2, 1, projection="3d")
    second = fig.add_subplot(1, 2, 2, projection="3d")
    draw_candidate(first, capped_prism_coordinates(), prism.edges, prism_faces,
                   prism_cap, "Capped pentagonal prism", fully_bound_pattern(prism))
    draw_candidate(second, pentagonal_face_coordinates(), ico.edges, ico_faces,
                   ico_cap, "Pentagonal-face icosahedral graph",
                   fully_bound_pattern(ico))
    fig.suptitle("Two 11-monomer geometry hypotheses with explicit interface assignments", fontsize=15, y=0.97)
    fig.text(0.5, 0.035,
             "Blue: dimer interface. Orange arrows: C-terminal contacts. Grey: other possible adjacencies.\n"
             "One feasible graph assignment per shape; protein conformations and thermodynamic weights remain unknown.",
             ha="center", fontsize=10, color="#344d60")
    path = Path(__file__).with_name("eleven_monomer_geometries.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
