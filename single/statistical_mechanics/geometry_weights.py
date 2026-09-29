"""Conditional Boltzmann weights for explicit, connected contact patterns.

No shape free energy or physical degeneracy is inferred from a graph. Those
quantities must be supplied for each distinct state before fractions can be
calculated. A graph feasibility check supplies one representative wiring;
it does not enumerate conformations or establish a protein structure.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Iterator, Sequence

import numpy as np

from .geometry import Scaffold
from .geometry_catalogue import MonomerGeometry


def _connected(size: int, edges: Sequence[tuple[int, int]]) -> bool:
    neighbours = [set() for _ in range(size)]
    for a, b in edges:
        neighbours[a].add(b)
        neighbours[b].add(a)
    seen, pending = {0}, [0]
    while pending:
        vertex = pending.pop()
        for neighbour in neighbours[vertex] - seen:
            seen.add(neighbour)
            pending.append(neighbour)
    return len(seen) == size


@dataclass(frozen=True)
class BondPattern:
    """Actual interface contacts between monomers numbered 0 through n-1.

    Dimer contacts form a matching. Each C-terminal contact is directed from
    one donor to one receiver, with at most one of each per monomer. Here
    C-terminal contacts join different dimer partners, as an explicit sHSP
    interface assumption. Reciprocal contacts are two distinct donor sites.
    """

    name: str
    monomers: int
    dimer_pairs: tuple[tuple[int, int], ...]
    c_terminal_contacts: tuple[tuple[int, int], ...]

    def __post_init__(self) -> None:
        if (not isinstance(self.monomers, int) or isinstance(self.monomers, bool)
                or self.monomers < 1):
            raise ValueError("monomers must be a positive integer")
        edges = self.dimer_pairs + self.c_terminal_contacts
        if any(not all(isinstance(v, int) and not isinstance(v, bool) for v in (a, b))
               or a == b or min(a, b) < 0 or max(a, b) >= self.monomers
               for a, b in edges):
            raise ValueError("invalid monomer contact")
        paired = [v for pair in self.dimer_pairs for v in pair]
        if len(paired) != len(set(paired)):
            raise ValueError("a monomer can have at most one dimer partner")
        donors = [a for a, _ in self.c_terminal_contacts]
        receivers = [b for _, b in self.c_terminal_contacts]
        if len(donors) != len(set(donors)) or len(receivers) != len(set(receivers)):
            raise ValueError("a monomer can donate and receive at most one C-terminal contact")
        dimers = {frozenset(pair) for pair in self.dimer_pairs}
        if any(frozenset(pair) in dimers for pair in self.c_terminal_contacts):
            raise ValueError("this model places C-terminal contacts outside the dimer interface")
        if not _connected(self.monomers, edges):
            raise ValueError("occupied contacts must form one connected oligomer")


@dataclass(frozen=True)
class GeometryState:
    """One distinct contact/conformer state with explicitly supplied thermodynamics.

    ``shape_free_energy_kbt`` includes the residual shape and conformational
    free energy relative to a common reference. ``log_degeneracy`` is an
    explicitly specified additional statistical factor, including symmetry
    where appropriate. Do not count the same entropy in both terms.
    """

    architecture: str
    bonds: BondPattern
    shape_free_energy_kbt: float
    log_degeneracy: float

    def __post_init__(self) -> None:
        if not self.architecture or not all(isfinite(x) for x in
                (self.shape_free_energy_kbt, self.log_degeneracy)):
            raise ValueError("an architecture name and finite state weights are required")


@dataclass(frozen=True)
class GeometryShare:
    architecture: str
    log_weight: float
    probability: float


def conditional_geometry_shares(
    states: Sequence[GeometryState], *, epsilon_dimer: float,
    epsilon_c_terminal: float,
) -> tuple[GeometryShare, ...]:
    """P(architecture | n), summing explicitly provided distinct states.

    log(q_s) = epsilon_dimer*D + epsilon_c_terminal*C - F_shape/kBT + log(g_s).
    Positive epsilon values favour contacts. A monomer activity factor
    cancels at fixed size; states of different sizes cannot be mixed here.
    """
    if not states or not all(isfinite(x) for x in (epsilon_dimer, epsilon_c_terminal)):
        raise ValueError("provide states and finite contact energies")
    if len({state.bonds.monomers for state in states}) != 1:
        raise ValueError("conditional shares require one oligomer size")
    keys = [(state.architecture, tuple(sorted(tuple(sorted(pair))
             for pair in state.bonds.dimer_pairs)),
             tuple(sorted(state.bonds.c_terminal_contacts)), state.shape_free_energy_kbt)
            for state in states]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate states must be combined in their statistical factor")
    grouped: dict[str, float] = {}
    for state in states:
        log_weight = (epsilon_dimer * len(state.bonds.dimer_pairs)
                      + epsilon_c_terminal * len(state.bonds.c_terminal_contacts)
                      - state.shape_free_energy_kbt + state.log_degeneracy)
        if not isfinite(log_weight):
            raise ValueError("state weights exceed the numerical range")
        grouped[state.architecture] = float(np.logaddexp(
            grouped.get(state.architecture, -np.inf), log_weight))
    names = tuple(grouped)
    log_weights = np.array([grouped[name] for name in names])
    weights = np.exp(log_weights - log_weights.max())
    probabilities = weights / weights.sum()
    return tuple(GeometryShare(name, grouped[name], float(probability))
                 for name, probability in zip(names, probabilities))


def scaffold_pattern(scaffold: Scaffold) -> BondPattern:
    """One explicit cyclic C-terminal wiring on a full dimer-on-edge scaffold."""
    return BondPattern(scaffold.name, scaffold.monomers,
                       scaffold.dimer_contacts, scaffold.c_terminal_contacts)


def open_chain_pattern(size: int) -> BondPattern:
    """Alternating dimers and reciprocal interdimer contacts along a path.

    This concrete interface convention gives C=n-2 for even n and C=n-1
    for odd n, rather than assigning one generic bond to every chain edge.
    """
    if not isinstance(size, int) or isinstance(size, bool) or size < 1:
        raise ValueError("size must be a positive integer")
    dimers = tuple((i, i + 1) for i in range(0, size - 1, 2))
    terminal = tuple(pair for i in range(1, size - 1, 2)
                     for pair in ((i, i + 1), (i + 1, i)))
    return BondPattern(f"open chain ({size})", size, dimers, terminal)


def _pairings(remaining: tuple[int, ...], neighbours: list[set[int]],
              target: int, pairs: tuple[tuple[int, int], ...] = ()) -> Iterator[tuple[tuple[int, int], ...]]:
    if len(pairs) == target:
        yield pairs
        return
    if len(pairs) + len(remaining) // 2 < target:
        return
    a, rest = remaining[0], remaining[1:]
    for b in rest:
        if b in neighbours[a]:
            yield from _pairings(tuple(v for v in rest if v != b), neighbours,
                                 target, pairs + ((a, b),))
    if len(remaining) % 2:
        yield from _pairings(rest, neighbours, target, pairs)


def fully_bound_pattern(geometry: MonomerGeometry) -> BondPattern:
    """Find one JOINT dimer/C-terminal assignment on a small adjacency graph.

    Each monomer has one C-terminal donor and receiver; floor(n/2) dimers
    are paired. The C-terminal cycles may be separate, provided their union
    with dimers is connected. A pairing and cycle checked separately do not
    establish this compatibility. A failure excludes only this saturation
    assumption, not every possible partly bound state on the geometry.
    """
    size = geometry.monomers
    if size > 12:
        raise ValueError("the feasibility search is limited to twelve monomers")
    neighbours = [set() for _ in range(size)]
    for a, b in geometry.edges:
        neighbours[a].add(b)
        neighbours[b].add(a)
    for dimers in _pairings(tuple(range(size)), neighbours, size // 2):
        excluded = {frozenset(pair) for pair in dimers}
        permitted = [tuple(sorted(b for b in neighbours[a]
                                 if frozenset((a, b)) not in excluded))
                     for a in range(size)]
        if any(not choices for choices in permitted):
            continue
        order = sorted(range(size), key=lambda a: len(permitted[a]))
        terminal: list[tuple[int, int]] = []

        def assign(position: int, used_receivers: set[int]):
            if position == size:
                contacts = tuple(terminal)
                return contacts if _connected(size, dimers + contacts) else None
            donor = order[position]
            for receiver in permitted[donor]:
                if receiver not in used_receivers:
                    terminal.append((donor, receiver))
                    result = assign(position + 1, used_receivers | {receiver})
                    if result is not None:
                        return result
                    terminal.pop()
            return None

        contacts = assign(0, set())
        if contacts is not None:
            return BondPattern(geometry.name, size, dimers, contacts)
    raise ValueError(f"{geometry.name} has no fully saturated joint contact assignment under these rules")


def remove_monomer(pattern: BondPattern, site: int) -> BondPattern:
    """Remove actual incident bonds, without silently rewiring a defect.

    Raises if removal produces separate fragments; that product would need a
    multiple-cluster model. A relaxed defect requires another explicit state.
    """
    if (pattern.monomers == 1 or not isinstance(site, int) or isinstance(site, bool)
            or not 0 <= site < pattern.monomers):
        raise ValueError("site must identify a removable monomer")

    def retained(edges):
        return tuple((a - int(a > site), b - int(b > site))
                     for a, b in edges if site not in (a, b))
    return BondPattern(f"{pattern.name}, site {site} removed", pattern.monomers - 1,
                       retained(pattern.dimer_pairs), retained(pattern.c_terminal_contacts))
