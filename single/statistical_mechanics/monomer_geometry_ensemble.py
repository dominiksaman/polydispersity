"""Exploratory oligomer geometries with one monomer at each graph vertex.

This is a different coarse graining from ``geometry_ensemble``: that module
puts a *dimer* on each scaffold edge, as in the published alphaB-crystallin
models. Here a six-vertex octahedron is a six-monomer hypothesis. Graph edges
are possible spatial adjacencies, not all simultaneously occupied bonds.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite
from typing import Mapping

from .geometry import cube, octahedron, tetrahedron
from .geometry_ensemble import antiprism, dipyramid, prism, pyramid


@dataclass(frozen=True)
class MonomerGeometry:
    name: str
    family: str
    monomers: int
    edges: tuple[tuple[int, int], ...]
    closed: bool
    evidence: str = "illustrative graph"

    def __post_init__(self) -> None:
        if self.monomers < 1:
            raise ValueError("at least one monomer is required")
        if any(a < 0 or b < 0 or a >= self.monomers or b >= self.monomers or a == b
               for a, b in self.edges):
            raise ValueError("invalid graph edge")
        if len({tuple(sorted(edge)) for edge in self.edges}) != len(self.edges):
            raise ValueError("duplicate graph edge")

    @property
    def degrees(self) -> tuple[int, ...]:
        result = [0] * self.monomers
        for a, b in self.edges:
            result[a] += 1
            result[b] += 1
        return tuple(result)

    @property
    def dimer_contacts(self) -> int:
        # Assumes the graph permits a maximum matching of this size.
        return self.monomers // 2

    @property
    def c_terminal_contacts(self) -> int:
        # One directed contact per step along the chosen path/cycle.
        return self.monomers if self.closed else max(0, self.monomers - 1)

    @property
    def excess_coordination(self) -> int:
        # Available adjacencies above a path/ring's degree two are a shape
        # descriptor. They are not counted as extra C-terminal bonds.
        return sum(max(0, degree - 2) ** 2 for degree in self.degrees)


@dataclass(frozen=True)
class MonomerEnergy:
    """Dimensionless energies in kBT; positive contact terms are favourable."""

    epsilon_dimer: float = 1.0
    epsilon_c_terminal: float = 1.0
    shape_penalty: float = 0.0

    def __post_init__(self) -> None:
        if not all(isfinite(x) for x in (self.epsilon_dimer,
                                         self.epsilon_c_terminal,
                                         self.shape_penalty)):
            raise ValueError("energies must be finite")
        if self.shape_penalty < 0:
            raise ValueError("shape penalty must be nonnegative")


@dataclass(frozen=True)
class MonomerContribution:
    size: int
    geometry: str
    family: str
    evidence: str
    dimer_contacts: int
    c_terminal_contacts: int
    excess_coordination: int
    log_weight: float
    share: float


def _path(size: int) -> MonomerGeometry:
    return MonomerGeometry(f"open chain ({size})", "open chain", size,
                           tuple((i, i + 1) for i in range(size - 1)), False,
                           "generic topology")


def _cycle(size: int) -> MonomerGeometry:
    return MonomerGeometry(f"closed ring ({size})", "closed ring", size,
                           tuple((i, (i + 1) % size) for i in range(size)), True,
                           "generic topology")


def _compact(name: str, size: int,
             edges: tuple[tuple[int, int], ...],
             evidence: str = "illustrative graph") -> MonomerGeometry:
    return MonomerGeometry(name, "compact polyhedron", size, edges, True, evidence)


def _icosahedron() -> MonomerGeometry:
    """Twelve vertices and thirty edges, with two staggered pentagonal rings."""
    edges = []
    for i in range(5):
        edges.extend(((0, 1 + i), (1 + i, 1 + (i + 1) % 5),
                      (6 + i, 6 + (i + 1) % 5), (6 + i, 11),
                      (1 + i, 6 + i), (1 + i, 6 + (i - 1) % 5)))
    return _compact("icosahedron (12 monomers)", 12, tuple(edges))


def _capped_pentagonal_prism() -> MonomerGeometry:
    """Eleven-vertex prism with a pyramid erected on one pentagonal face."""
    edges = prism(5).edges + tuple((i, 10) for i in range(5))
    return _compact("capped pentagonal prism (11 monomers)", 11, edges)


def _pentagonal_face_icosahedron() -> MonomerGeometry:
    """Delete one icosahedral apex, leaving an eleven-vertex polyhedral graph.

    The five neighbours of the deleted apex bound a pentagonal face in the
    resulting graph; no missing monomer is included in its contact counts.
    """
    edges = tuple((a - 1, b - 1) for a, b in _icosahedron().edges
                  if a != 0 and b != 0)
    return _compact("pentagonal-face icosahedral polyhedron (11 monomers)",
                    11, edges)


def monomer_candidates(max_size: int = 12, *,
                       include_small_polyhedra: bool = False) -> tuple[MonomerGeometry, ...]:
    """One open path and one closed cycle per size, plus selected compact graphs.

    Small compact graphs are mathematically possible at four/five monomers,
    but omitted by default to make the user's proposed six-monomer onset an
    explicit hypothesis. No claim is made that the library is exhaustive.
    """
    if not isinstance(max_size, int) or isinstance(max_size, bool) or max_size < 1:
        raise ValueError("max_size must be a positive integer")
    candidates = [_path(size) for size in range(1, max_size + 1)]
    candidates += [_cycle(size) for size in range(3, max_size + 1)]
    compact = [
        _compact("tetrahedron (4 monomers)", 4, tetrahedron().edges),
        _compact("square pyramid (5 monomers)", 5, pyramid(4).edges),
        _compact("triangular dipyramid (5 monomers)", 5, dipyramid(3).edges),
        _compact("triangular prism (6 monomers)", 6, prism(3).edges),
        _compact("octahedron (6 monomers)", 6, octahedron().edges),
        _compact("cube (8 monomers)", 8, cube().edges),
        _capped_pentagonal_prism(),
        _pentagonal_face_icosahedron(),
        _icosahedron(),
    ]
    for sides in range(5, 9):
        compact.append(_compact(f"{sides}-gonal pyramid ({sides + 1} monomers)",
                                sides + 1, pyramid(sides).edges))
    for sides in range(5, 9):
        compact.append(_compact(f"{sides}-gonal dipyramid ({sides + 2} monomers)",
                                sides + 2, dipyramid(sides).edges))
    for sides in range(4, 7):
        if sides != 4:  # four-sided prism is the cube
            compact.append(_compact(f"{sides}-gonal prism ({2 * sides} monomers)",
                                    2 * sides, prism(sides).edges))
        compact.append(_compact(f"{sides}-gonal antiprism ({2 * sides} monomers)",
                                2 * sides, antiprism(sides).edges))
    minimum = 4 if include_small_polyhedra else 6
    candidates += [item for item in compact
                   if minimum <= item.monomers <= max_size]
    return tuple(candidates)


def monomer_contributions(
    size: int, energy: MonomerEnergy = MonomerEnergy(), *,
    candidates: tuple[MonomerGeometry, ...] | None = None,
    offsets_kbt: Mapping[str, float] | None = None,
    log_multiplicities: Mapping[str, float] | None = None,
) -> tuple[MonomerContribution, ...]:
    """Boltzmann shares conditional on size and the stated graph catalogue.

    One representative maximal dimer matching and Hamiltonian path/cycle is
    assumed per graph. Their unknown conformational and symmetry factors may
    be supplied through ``log_multiplicities``; defaults are all equal.
    """
    if not isinstance(size, int) or isinstance(size, bool) or size < 1:
        raise ValueError("size must be a positive integer")
    candidates = monomer_candidates(max(12, size)) if candidates is None else candidates
    selected = [item for item in candidates if item.monomers == size]
    if not selected:
        return ()
    names = {item.name for item in candidates}
    if len(names) != len(candidates):
        raise ValueError("candidate names must be unique")
    offsets_kbt = {} if offsets_kbt is None else offsets_kbt
    log_multiplicities = {} if log_multiplicities is None else log_multiplicities
    if (set(offsets_kbt) | set(log_multiplicities)) - names:
        raise ValueError("unknown candidate in offset or multiplicity mapping")
    weighted = []
    for item in selected:
        offset = offsets_kbt.get(item.name, 0.0)
        log_mult = log_multiplicities.get(item.name, 0.0)
        if not isfinite(offset) or not isfinite(log_mult):
            raise ValueError("offsets and log multiplicities must be finite")
        logw = (energy.epsilon_dimer * item.dimer_contacts
                + energy.epsilon_c_terminal * item.c_terminal_contacts
                - energy.shape_penalty * item.excess_coordination
                - offset + log_mult)
        weighted.append((item, logw))
    top = max(logw for _, logw in weighted)
    denominator = sum(exp(logw - top) for _, logw in weighted)
    return tuple(MonomerContribution(
        size, item.name, item.family, item.evidence, item.dimer_contacts,
        item.c_terminal_contacts, item.excess_coordination, logw,
        exp(logw - top) / denominator)
        for item, logw in weighted)
