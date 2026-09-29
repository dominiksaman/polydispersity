"""Illustrative equilibrium cluster partition functions on finite scaffolds.

This is a static statistical-mechanical model within one specified parent
scaffold: each connected occupancy pattern is given a Boltzmann weight, and
patterns related by that scaffold's symmetry are counted once. These are
*not* all molecular arrangements at a given size. No association or
dissociation rate appears in the partition function. The scaffold is a
template of possible contacts, not a pre-existing empty protein shell.

One may occupy whole dimeric edges, with at most one monomer left unpaired.
At a vertex with q occupied monomer endpoints, q >= 2, the q monomers are
assumed able to form a cyclic set of q directed C-terminal contacts. This
local rearrangement rule is a coarse-grained hypothesis, not a structure
determination. Models include only one parent scaffold at a time and do not
yet include bending/closure strain or translational/rotational entropy.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import permutations, product
from math import log

import numpy as np
from scipy.special import logsumexp

from .geometry import Scaffold


def _popcount(value: int) -> int:
    # The repository's default macOS Python is 3.9 (before int.bit_count).
    return bin(value).count("1")


@dataclass(frozen=True)
class EquilibriumParameters:
    """Dimensionless chemical activity and favourable bond energies (in kBT)."""

    log_monomer_activity: float
    dimer_stabilization: float
    c_terminal_stabilization: float
    saturated_vertex_penalty: float = 0.0


@dataclass(frozen=True)
class StateCatalogue:
    """Symmetry-inequivalent occupied subgraphs of one parent scaffold."""

    scaffold_name: str
    max_size: int
    # (monomers, dimer contacts, C-terminal contacts, saturated vertices) -> count
    counts: dict[tuple[int, int, int, int], int]

    def distribution(self, parameters: EquilibriumParameters) -> np.ndarray:
        """P(n) from Z_n = sum exp[n log(z) + D eps_d + C eps_c - V kappa]."""
        values = np.array((parameters.log_monomer_activity,
                           parameters.dimer_stabilization,
                           parameters.c_terminal_stabilization,
                           parameters.saturated_vertex_penalty), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("parameters must be finite")
        log_terms: list[list[float]] = [[] for _ in range(self.max_size)]
        for (n, d, c, saturated), multiplicity in self.counts.items():
            log_terms[n - 1].append(log(multiplicity) + n * values[0]
                                    + d * values[1] + c * values[2]
                                    - saturated * values[3])
        log_weights = np.array([logsumexp(terms) if terms else -np.inf
                                for terms in log_terms])
        return np.exp(log_weights - logsumexp(log_weights))


def _vertex_automorphisms(scaffold: Scaffold) -> list[tuple[int, ...]]:
    vertices = sorted({v for edge in scaffold.edges for v in edge})
    if vertices != list(range(len(vertices))):
        raise ValueError("scaffold vertices must be numbered consecutively")
    if scaffold.name.endswith("-dimer ring"):
        m = len(vertices)
        return [tuple((direction * v + shift) % m for v in vertices)
                for direction in (1, -1) for shift in range(m)]
    if scaffold.name == "cube":
        return [tuple(sum((((v >> bit) & 1) ^ flip[bit]) << perm[bit]
                          for bit in range(3)) for v in vertices)
                for perm in permutations(range(3)) for flip in product((0, 1), repeat=3)]
    if scaffold.name == "octahedron":
        return [tuple(2 * perm[v // 2] + ((v % 2) ^ flip[v // 2])
                      for v in vertices)
                for perm in permutations(range(3)) for flip in product((0, 1), repeat=3)]
    raise ValueError("symmetry enumeration supports rings, cube, and octahedron")


def _monomer_automorphisms(scaffold: Scaffold) -> list[tuple[int, ...]]:
    edge_lookup = {frozenset(edge): i for i, edge in enumerate(scaffold.edges)}
    maps = []
    for vertex_map in _vertex_automorphisms(scaffold):
        monomer_map = []
        for a, b in scaffold.edges:
            mapped_a, mapped_b = vertex_map[a], vertex_map[b]
            new_edge_index = edge_lookup[frozenset((mapped_a, mapped_b))]
            new_a, new_b = scaffold.edges[new_edge_index]
            monomer_map.extend((2 * new_edge_index + (0 if mapped_a == new_a else 1),
                                2 * new_edge_index + (0 if mapped_b == new_a else 1)))
        maps.append(tuple(monomer_map))
    return list(dict.fromkeys(maps))


def enumerate_states(scaffold: Scaffold) -> StateCatalogue:
    """Enumerate connected full-dimer or single-defect template states.

    This is practical for a 24-mer parent (12 dimer edges). The orbit of a
    state under parent-scaffold symmetry is counted once, avoiding artificial
    multiplication by equivalent positions on an unanchored assembly.
    """
    if len(scaffold.edges) > 14:
        raise ValueError("exact enumeration is limited to 14 dimeric edges")
    n_sites = scaffold.monomers
    vertex_masks: dict[int, int] = {}
    neighbours = [1 << (i ^ 1) for i in range(n_sites)]
    for edge_index, (a, b) in enumerate(scaffold.edges):
        vertex_masks[a] = vertex_masks.get(a, 0) | (1 << (2 * edge_index))
        vertex_masks[b] = vertex_masks.get(b, 0) | (1 << (2 * edge_index + 1))
    for vertex_mask in vertex_masks.values():
        for i in range(n_sites):
            if vertex_mask & (1 << i):
                neighbours[i] |= vertex_mask & ~(1 << i)
    automorphisms = _monomer_automorphisms(scaffold)

    def connected(mask: int) -> bool:
        reached = mask & -mask
        while True:
            candidates = reached
            bits = reached
            while bits:
                bit = bits & -bits
                candidates |= neighbours[bit.bit_length() - 1] & mask
                bits -= bit
            if candidates == reached:
                return reached == mask
            reached = candidates

    def canonical(mask: int) -> int:
        minimum = mask
        occupied = []
        bits = mask
        while bits:
            bit = bits & -bits
            occupied.append(bit.bit_length() - 1)
            bits -= bit
        for automorphism in automorphisms:
            transformed = 0
            for monomer in occupied:
                transformed |= 1 << automorphism[monomer]
            minimum = min(minimum, transformed)
        return minimum

    seen: set[int] = set()
    counts: Counter[tuple[int, int, int, int]] = Counter()
    for edge_subset in range(1, 1 << len(scaffold.edges)):
        full_mask = 0
        bits = edge_subset
        while bits:
            bit = bits & -bits
            index = bit.bit_length() - 1
            full_mask |= 3 << (2 * index)
            bits -= bit
        candidates = [full_mask]
        bits = full_mask
        while bits:
            bit = bits & -bits
            candidates.append(full_mask ^ bit)
            bits -= bit
        for mask in candidates:
            if not connected(mask):
                continue
            representative = canonical(mask)
            if representative in seen:
                continue
            seen.add(representative)
            dimer_contacts = sum((mask >> (2 * i)) & 3 == 3
                                 for i in range(len(scaffold.edges)))
            c_terminal_contacts = sum(q for vertex_mask in vertex_masks.values()
                                      if (q := _popcount(mask & vertex_mask)) >= 2)
            saturated_vertices = sum(_popcount(mask & vertex_mask) == _popcount(vertex_mask)
                                     and _popcount(vertex_mask) >= 3
                                     for vertex_mask in vertex_masks.values())
            counts[(_popcount(mask), dimer_contacts,
                    c_terminal_contacts, saturated_vertices)] += 1
    return StateCatalogue(scaffold.name, n_sites, dict(counts))
