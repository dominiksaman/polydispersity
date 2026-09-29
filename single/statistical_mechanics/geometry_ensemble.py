"""Conditional equilibrium weights for alternative oligomer architectures.

An edge represents an intact dimer and its two endpoints are monomers. Closed
scaffolds have one directed C-terminal contact per monomer. An odd oligomer is
modelled as one missing monomer on an otherwise closed even parent, with the
remaining vertex contacts allowed to reconnect. These are coarse contact
graphs, not atomic structures or measured solution populations.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite, log
from typing import Mapping

from .geometry import Scaffold, cube, dimer_ring, octahedron, tetrahedron


@dataclass(frozen=True)
class GeometryCandidate:
    name: str
    family: str
    scaffold: Scaffold
    evidence: str = "illustrative graph"


@dataclass(frozen=True)
class EnergyParameters:
    """All energies are dimensionless, in units of k_B*T.

    Positive epsilon values stabilize contacts. ``degree_penalty`` charges
    each vertex by kappa*(occupied_degree - preferred_degree)^2. This is an
    illustrative shape/conformational free energy, not a measured bond energy.
    """

    epsilon_dimer: float = 1.0
    epsilon_c_terminal: float = 1.0
    degree_penalty: float = 0.0
    preferred_degree: int = 4

    def __post_init__(self) -> None:
        if not all(isfinite(x) for x in (self.epsilon_dimer,
                                         self.epsilon_c_terminal,
                                         self.degree_penalty)):
            raise ValueError("energies must be finite")
        if self.degree_penalty < 0 or self.preferred_degree < 2:
            raise ValueError("invalid degree penalty or preferred degree")


@dataclass(frozen=True)
class Contribution:
    size: int
    geometry: str
    family: str
    evidence: str
    parent_size: int
    microstates: int
    dimer_contacts: int | None
    c_terminal_contacts: int | None
    log_weight: float
    share: float


def prism(sides: int) -> Scaffold:
    if sides < 3:
        raise ValueError("prism needs at least three sides")
    edges = [(i, (i + 1) % sides) for i in range(sides)]
    edges += [(sides + i, sides + (i + 1) % sides) for i in range(sides)]
    edges += [(i, sides + i) for i in range(sides)]
    return Scaffold(f"{sides}-gonal prism", tuple(edges))


def antiprism(sides: int) -> Scaffold:
    if sides < 3:
        raise ValueError("antiprism needs at least three sides")
    edges = [(i, (i + 1) % sides) for i in range(sides)]
    edges += [(sides + i, sides + (i + 1) % sides) for i in range(sides)]
    for i in range(sides):
        edges.extend(((i, sides + i), (i, sides + (i - 1) % sides)))
    return Scaffold(f"{sides}-gonal antiprism", tuple(edges))


def pyramid(sides: int) -> Scaffold:
    if sides < 3:
        raise ValueError("pyramid needs at least three sides")
    edges = [(i, (i + 1) % sides) for i in range(sides)]
    edges += [(i, sides) for i in range(sides)]
    return Scaffold(f"{sides}-gonal pyramid", tuple(edges))


def dipyramid(sides: int) -> Scaffold:
    if sides < 3:
        raise ValueError("dipyramid needs at least three sides")
    edges = [(i, (i + 1) % sides) for i in range(sides)]
    edges += [(i, sides) for i in range(sides)]
    edges += [(i, sides + 1) for i in range(sides)]
    return Scaffold(f"{sides}-gonal dipyramid", tuple(edges))


def split_octahedron(split_opposites: bool = False) -> Scaffold:
    """A vertex-split octahedral graph, without a claim of atomic geometry."""
    edges = list(octahedron().edges)
    for old, new in ((0, 6), (1, 7))[:2 if split_opposites else 1]:
        neighbours = sorted(b if a == old else a for a, b in edges
                            if a == old or b == old)
        # Split four original contacts into two pairs and add one new edge.
        for neighbour in neighbours[2:]:
            edge = (min(old, neighbour), max(old, neighbour))
            edges.remove(edge)
            edges.append((min(new, neighbour), max(new, neighbour)))
        edges.append((old, new))
    return Scaffold("two vertex-split octahedron" if split_opposites
                    else "vertex-split octahedron", tuple(edges))


def default_candidates(max_size: int = 40) -> tuple[GeometryCandidate, ...]:
    """A small, explicit library; one graph topology per named architecture.

    Rings above 16 dimers are marked as extrapolations. Aliases such as the
    four-sided prism/cube and three-sided antiprism/octahedron are omitted.
    Multiple-ring cages need a separate connected contact model and are not
    silently counted here.
    """
    candidates = [GeometryCandidate(f"single ring ({m} dimers)", "ring",
                                    dimer_ring(m),
                                    "extrapolated" if m > 16 else "catalogued")
                  for m in range(3, max_size // 2 + 1)]
    candidates += [
        GeometryCandidate("tetrahedron", "degree-3 polyhedron", tetrahedron(), "catalogued"),
        GeometryCandidate("square pyramid", "mixed-degree polyhedron", pyramid(4)),
        GeometryCandidate("pentagonal pyramid", "mixed-degree polyhedron", pyramid(5)),
        GeometryCandidate("triangular prism", "degree-3 polyhedron", prism(3), "catalogued"),
        GeometryCandidate("triangular dipyramid", "mixed-degree polyhedron", dipyramid(3)),
        GeometryCandidate("cube", "degree-3 polyhedron", cube(), "catalogued"),
        GeometryCandidate("octahedron", "degree-4 polyhedron", octahedron(), "catalogued"),
        GeometryCandidate("vertex-split octahedron", "mixed-degree polyhedron",
                          split_octahedron(), "illustrative vertex split"),
        GeometryCandidate("two vertex-split octahedron", "mixed-degree polyhedron",
                          split_octahedron(True), "illustrative vertex split"),
        GeometryCandidate("pentagonal prism", "degree-3 polyhedron", prism(5), "catalogued"),
        GeometryCandidate("pentagonal dipyramid", "mixed-degree polyhedron", dipyramid(5)),
        GeometryCandidate("square antiprism", "degree-4 polyhedron", antiprism(4), "catalogued"),
        GeometryCandidate("hexagonal prism", "degree-3 polyhedron", prism(6), "catalogued"),
        GeometryCandidate("pentagonal antiprism", "degree-4 polyhedron", antiprism(5), "catalogued"),
    ]
    return tuple(c for c in candidates if c.scaffold.monomers <= max_size)


def _state_counts(scaffold: Scaffold, removed: int | None) -> tuple[int, int, tuple[int, ...]]:
    incident = scaffold._incidence()
    dimer_count = len(scaffold.edges) - int(removed is not None)
    terminal_count = 0
    occupied_degrees = []
    for monomers in incident.values():
        degree = len(monomers) - int(removed in monomers)
        occupied_degrees.append(degree)
        if degree >= 2:
            terminal_count += degree
    return dimer_count, terminal_count, tuple(occupied_degrees)


def _log_weight(scaffold: Scaffold, size: int, params: EnergyParameters,
                offset: float, log_multiplicity: float) -> tuple[float, int, int | None, int | None]:
    removed_sites = (None,) if size == scaffold.monomers else range(scaffold.monomers)
    states = []
    counts = []
    for removed in removed_sites:
        dimers, terminals, degrees = _state_counts(scaffold, removed)
        logw = (params.epsilon_dimer * dimers
                + params.epsilon_c_terminal * terminals
                - params.degree_penalty * sum((d - params.preferred_degree) ** 2
                                              for d in degrees)
                - offset + log_multiplicity)
        states.append(logw)
        counts.append((dimers, terminals))
    top = max(states)
    total = top + log(sum(exp(value - top) for value in states))
    unique_counts = set(counts)
    d, c = next(iter(unique_counts)) if len(unique_counts) == 1 else (None, None)
    return total, len(states), d, c


def conditional_contributions(
    size: int, params: EnergyParameters = EnergyParameters(),
    *, candidates: tuple[GeometryCandidate, ...] | None = None,
    offsets_kbt: Mapping[str, float] | None = None,
    log_multiplicities: Mapping[str, float] | None = None,
) -> tuple[Contribution, ...]:
    """Return P(geometry | size) for the stated candidate catalogue.

    Odd sizes sum over possible missing-monomer sites of the next even parent.
    Equal default architecture factors are a convention, not measured priors.
    Monomer activity cancels because all candidate structures have this size.
    """
    if not isinstance(size, int) or isinstance(size, bool) or size < 5:
        raise ValueError("size must be an integer at least five")
    candidates = default_candidates(max(40, size + (size % 2))) if candidates is None else candidates
    offsets_kbt = {} if offsets_kbt is None else offsets_kbt
    log_multiplicities = {} if log_multiplicities is None else log_multiplicities
    known = {candidate.name for candidate in candidates}
    if len(known) != len(candidates):
        raise ValueError("candidate names must be unique")
    if (set(offsets_kbt) | set(log_multiplicities)) - known:
        raise ValueError("unknown candidate in offset or multiplicity mapping")
    selected = [candidate for candidate in candidates
                if candidate.scaffold.monomers == size + size % 2]
    if not selected:
        return ()
    intermediate = []
    for candidate in selected:
        offset = offsets_kbt.get(candidate.name, 0.0)
        log_mult = log_multiplicities.get(candidate.name, 0.0)
        if not isfinite(offset) or not isfinite(log_mult):
            raise ValueError("offsets and log multiplicities must be finite")
        logw, states, dimers, terminals = _log_weight(
            candidate.scaffold, size, params, offset, log_mult)
        intermediate.append((candidate, logw, states, dimers, terminals))
    top = max(item[1] for item in intermediate)
    normalizer = sum(exp(item[1] - top) for item in intermediate)
    return tuple(Contribution(
        size, candidate.name, candidate.family, candidate.evidence,
        candidate.scaffold.monomers, states, dimers, terminals, logw,
        exp(logw - top) / normalizer)
        for candidate, logw, states, dimers, terminals in intermediate)
