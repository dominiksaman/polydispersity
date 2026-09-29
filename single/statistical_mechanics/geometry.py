"""Contact topology for dimer-ring and polyhedral single-protein assemblies.

Each scaffold edge is one *dimer*, containing a monomer at each endpoint.
At a scaffold vertex, the incident monomers are linked in a directed cycle:
each donates one C-terminal contact and receives one. This connectivity is a
minimal graph idealisation of the structural models, not an atomic structure.
It does not supply geometry-dependent strain or conformational entropy.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite


@dataclass(frozen=True)
class Scaffold:
    name: str
    edges: tuple[tuple[int, int], ...]

    def __post_init__(self) -> None:
        if not self.edges or any(a == b or a < 0 or b < 0 for a, b in self.edges):
            raise ValueError("scaffold requires distinct nonnegative vertices")
        if len({tuple(sorted(edge)) for edge in self.edges}) != len(self.edges):
            raise ValueError("scaffold edges must be unique")
        if any(len(incident) < 2 for incident in self._incidence().values()):
            raise ValueError("each vertex needs at least two incident dimers")

    def _incidence(self) -> dict[int, list[int]]:
        incident: dict[int, list[int]] = {}
        for dimer_index, (a, b) in enumerate(self.edges):
            incident.setdefault(a, []).append(2 * dimer_index)
            incident.setdefault(b, []).append(2 * dimer_index + 1)
        return incident

    @property
    def monomers(self) -> int:
        return 2 * len(self.edges)

    @property
    def dimer_contacts(self) -> tuple[tuple[int, int], ...]:
        return tuple((2 * i, 2 * i + 1) for i in range(len(self.edges)))

    @property
    def c_terminal_contacts(self) -> tuple[tuple[int, int], ...]:
        """One directed donor-to-receiver link per monomer.

        The arbitrary cyclic ordering at each vertex affects detailed
        connectivity but not the contact count. It is not a conformation fit.
        """
        contacts = []
        for incident in self._incidence().values():
            contacts.extend((monomer, incident[(i + 1) % len(incident)])
                            for i, monomer in enumerate(incident))
        return tuple(contacts)

    @property
    def vertex_degrees(self) -> tuple[int, ...]:
        return tuple(sorted(len(incident) for incident in self._incidence().values()))

    def contact_free_energy(self, dimer_kj_mol: float,
                            c_terminal_kj_mol: float) -> float:
        """Additive energy only; no chemical potential, degeneracy, or strain."""
        if not all(isfinite(x) for x in (dimer_kj_mol, c_terminal_kj_mol)):
            raise ValueError("contact energies must be finite")
        return (len(self.dimer_contacts) * dimer_kj_mol
                + len(self.c_terminal_contacts) * c_terminal_kj_mol)


def dimer_ring(dimers: int) -> Scaffold:
    """A closed single ring of at least three dimeric edges."""
    if not isinstance(dimers, int) or isinstance(dimers, bool) or dimers < 3:
        raise ValueError("dimers must be an integer of at least three")
    return Scaffold(f"{dimers}-dimer ring",
                    tuple((i, (i + 1) % dimers) for i in range(dimers)))


def cube() -> Scaffold:
    """Eight degree-three vertices, twelve dimeric edges: a 24-mer."""
    return Scaffold("cube", tuple((v, v ^ (1 << bit))
                                  for v in range(8) for bit in range(3)
                                  if v < (v ^ (1 << bit))))


def octahedron() -> Scaffold:
    """Six degree-four vertices, twelve dimeric edges: a 24-mer."""
    return Scaffold("octahedron", tuple((a, b) for a in range(6)
                                         for b in range(a + 1, 6)
                                         if a // 2 != b // 2))


def tetrahedron() -> Scaffold:
    """Four degree-three vertices, six dimeric edges: a 12-mer."""
    return Scaffold("tetrahedron", tuple((a, b) for a in range(4)
                                         for b in range(a + 1, 4)))


def _validate_sides(sides: int) -> None:
    if not isinstance(sides, int) or isinstance(sides, bool) or sides < 3:
        raise ValueError("sides must be an integer of at least three")


def prism(sides: int) -> Scaffold:
    _validate_sides(sides)
    edges = [(i, (i + 1) % sides) for i in range(sides)]
    edges += [(sides + i, sides + (i + 1) % sides) for i in range(sides)]
    edges += [(i, sides + i) for i in range(sides)]
    return Scaffold(f"{sides}-gonal prism", tuple(edges))


def antiprism(sides: int) -> Scaffold:
    _validate_sides(sides)
    edges = [(i, (i + 1) % sides) for i in range(sides)]
    edges += [(sides + i, sides + (i + 1) % sides) for i in range(sides)]
    for i in range(sides):
        edges.extend(((i, sides + i), (i, sides + (i - 1) % sides)))
    return Scaffold(f"{sides}-gonal antiprism", tuple(edges))


def pyramid(sides: int) -> Scaffold:
    _validate_sides(sides)
    edges = [(i, (i + 1) % sides) for i in range(sides)]
    edges += [(i, sides) for i in range(sides)]
    return Scaffold(f"{sides}-gonal pyramid", tuple(edges))


def dipyramid(sides: int) -> Scaffold:
    _validate_sides(sides)
    edges = [(i, (i + 1) % sides) for i in range(sides)]
    edges += [(i, sides) for i in range(sides)]
    edges += [(i, sides + 1) for i in range(sides)]
    return Scaffold(f"{sides}-gonal dipyramid", tuple(edges))
