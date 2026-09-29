"""Mathematical geometry hypotheses with one monomer at each graph vertex.

Edges are possible spatial adjacencies. No interface energies, contact counts
or equilibrium populations are inferred from vertex degree. This convention
differs from the dimer-on-edge scaffolds in ``geometry.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from .geometry import antiprism, cube, dipyramid, octahedron, prism, pyramid, tetrahedron


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
                       include_small_polyhedra: bool = True) -> tuple[MonomerGeometry, ...]:
    """One open path and one closed cycle per size, plus selected compact graphs.

    Compact graphs are included from four monomers. Excluding them at four/
    five with ``include_small_polyhedra=False`` is an optional physical
    hypothesis, not a geometric threshold. The library is not exhaustive.
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
