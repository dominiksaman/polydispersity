"""Statistical-mechanical co-assembly of A and B subunits on a ring.

Implements the circular-arrangement model in the thesis's Chapter 2,
``text/ch1-intro.tex`` lines 206-245. A fixed-size ring has one contact
between each neighbouring pair, including the last and first subunits.
Rotations of an arrangement contribute with their full multiplicity; mirror
images are distinct if they differ as labelled arrangements. Arrangements
with the same composition and contact counts are grouped exactly, so no
2**size enumeration is needed for normal calculations.

Energies are in kJ/mol per contact. With a fitted activity ratio, the
identifiable contact contrast is g_AB - (g_AA + g_BB)/2; g_AA-g_BB is
confounded with that activity ratio. Absolute interface energies are not
identified by a fixed-size composition distribution.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import comb, isfinite, log
from typing import Optional, Sequence, Tuple

import numpy as np

R_KJ_MOL_K = 0.00831446261815324


@dataclass(frozen=True)
class ContactClass:
    """All labelled ring arrangements with these composition/contact counts."""

    n_a: int
    n_b: int
    n_aa: int
    n_ab: int
    n_bb: int
    multiplicity: int


@dataclass(frozen=True)
class FitResult:
    """A shared AB contact energy and one fitted activity ratio per mixture."""

    hetero_energy_kj_mol: float
    log_activity_ratios: Tuple[float, ...]
    fitted_mole_fractions_a: Tuple[float, ...]
    predicted: Tuple[np.ndarray, ...]
    root_mean_square_error: float
    success: bool


@lru_cache(maxsize=None)
def contact_classes(size: int) -> Tuple[ContactClass, ...]:
    """Count every labelled A/B ring exactly, grouped by composition and contacts.

For a mixed ring with i A subunits and r runs of A (hence r runs of B),
there are 2r AB contacts, i-r AA contacts, and size-i-r BB contacts. Its
labelled multiplicity is size/r * C(i-1,r-1) * C(size-i-1,r-1). This is
equivalent to summing the rotational multiplicities of the thesis's ring
arrangements, including periodic arrangements such as ABABAB.
"""
    if not isinstance(size, int) or isinstance(size, bool) or size < 3:
        raise ValueError("size must be an integer of at least 3")

    classes = [ContactClass(0, size, 0, 0, size, 1)]
    for n_a in range(1, size):
        n_b = size - n_a
        for runs in range(1, min(n_a, n_b) + 1):
            numerator = size * comb(n_a - 1, runs - 1) * comb(n_b - 1, runs - 1)
            multiplicity, remainder = divmod(numerator, runs)
            assert remainder == 0
            classes.append(ContactClass(
                n_a, n_b, n_a - runs, 2 * runs, n_b - runs, multiplicity
            ))
    classes.append(ContactClass(size, 0, size, 0, 0, 1))
    return tuple(classes)


def _validate_energies(g_aa: float, g_ab: float, g_bb: float, temperature: float) -> None:
    if not all(isfinite(value) for value in (g_aa, g_ab, g_bb, temperature)):
        raise ValueError("energies and temperature must be finite")
    if temperature <= 0:
        raise ValueError("temperature must be positive in kelvin")


def _class_weights(
    size: int,
    g_aa: float,
    g_ab: float,
    g_bb: float,
    temperature: float,
    log_activity_ratio: float,
) -> Tuple[Tuple[ContactClass, ...], np.ndarray]:
    _validate_energies(g_aa, g_ab, g_bb, temperature)
    if not isfinite(log_activity_ratio):
        raise ValueError("log_activity_ratio must be finite")
    classes = contact_classes(size)
    beta = 1.0 / (R_KJ_MOL_K * temperature)
    log_weights = np.array([
        log(item.multiplicity) + item.n_a * log_activity_ratio
        - beta * (item.n_aa * g_aa + item.n_ab * g_ab + item.n_bb * g_bb)
        for item in classes
    ])
    weights = np.exp(log_weights - log_weights.max())
    weights /= weights.sum()
    return classes, weights


def activity_distribution(
    size: int,
    g_ab: float,
    *,
    g_aa: float = 0.0,
    g_bb: float = 0.0,
    temperature: float = 298.15,
    log_activity_ratio: float = 0.0,
) -> np.ndarray:
    """Return probabilities for 0, 1, ..., ``size`` A subunits.

    ``log_activity_ratio`` is ln(a/b). The default of zero gives equal
    activities, which gives a 1:1 mixture when AA and BB energies are equal.
    """
    classes, weights = _class_weights(
        size, g_aa, g_ab, g_bb, temperature, log_activity_ratio
    )
    result = np.zeros(size + 1)
    for item, weight in zip(classes, weights):
        result[item.n_a] += weight
    return result


def composition_distribution(
    size: int,
    g_ab: float,
    *,
    g_aa: float = 0.0,
    g_bb: float = 0.0,
    temperature: float = 298.15,
    mole_fraction_a: float = 0.5,
) -> np.ndarray:
    """Return the distribution with mean A fraction ``mole_fraction_a``.

    Solves the relative activity a/b so that the mean ring composition is
    the requested mixing fraction. This is a fixed-size equilibrium model;
    it does not model assembly kinetics or oligomer-size distributions.
    """
    if not isfinite(mole_fraction_a) or not 0 <= mole_fraction_a <= 1:
        raise ValueError("mole_fraction_a must lie between 0 and 1")
    _validate_energies(g_aa, g_ab, g_bb, temperature)
    contact_classes(size)
    if mole_fraction_a == 0:
        return np.eye(1, size + 1, 0, dtype=float)[0]
    if mole_fraction_a == 1:
        return np.eye(1, size + 1, size, dtype=float)[0]

    def mean_fraction(log_ratio: float) -> float:
        probabilities = activity_distribution(
            size, g_ab, g_aa=g_aa, g_bb=g_bb, temperature=temperature,
            log_activity_ratio=log_ratio,
        )
        return float(np.dot(np.arange(size + 1), probabilities) / size)

    lower, upper = -1.0, 1.0
    while mean_fraction(lower) > mole_fraction_a:
        lower *= 2
    while mean_fraction(upper) < mole_fraction_a:
        upper *= 2
    for _ in range(80):
        midpoint = (lower + upper) / 2
        if mean_fraction(midpoint) < mole_fraction_a:
            lower = midpoint
        else:
            upper = midpoint
    return activity_distribution(
        size, g_ab, g_aa=g_aa, g_bb=g_bb, temperature=temperature,
        log_activity_ratio=(lower + upper) / 2,
    )


def fit_hetero_energy(
    observed: Sequence[Sequence[float]],
    *,
    g_aa: float = 0.0,
    g_bb: float = 0.0,
    temperature: float = 298.15,
    initial_mole_fractions_a: Optional[Sequence[float]] = None,
    energy_bounds: Tuple[float, float] = (-30.0, 30.0),
) -> FitResult:
    """Fit one AB contact energy and a separate a/b ratio for each mixture.

    Each observed series has ``size + 1`` nonnegative abundances ordered by
    A count (0 through size). It is normalized before fitting. The activity
    ratios are nuisance parameters, as in the thesis's simultaneous fit of
    three mixing ratios. AA and BB are held fixed: the identifiable energy
    is the AB contrast against their mean, while their difference can be
    absorbed into the activity ratio.

    Requires SciPy. The fit minimizes unweighted squared differences in
    *normalized* composition probabilities; experimental uncertainties are
    not inferred or used.
    """
    try:
        from scipy.optimize import least_squares
    except ImportError as exc:
        raise ImportError("fit_hetero_energy requires scipy") from exc

    series = [np.asarray(row, dtype=float) for row in observed]
    if not series or any(row.ndim != 1 or len(row) < 4 for row in series):
        raise ValueError("provide at least one 1D series for a ring of size >= 3")
    size = len(series[0]) - 1
    if any(len(row) != size + 1 for row in series):
        raise ValueError("all observed series must have the same size")
    if any(not np.all(np.isfinite(row)) or np.any(row < 0) or row.sum() <= 0
           for row in series):
        raise ValueError("observed abundances must be finite, nonnegative, and nonzero")
    series = [row / row.sum() for row in series]
    _validate_energies(g_aa, 0.0, g_bb, temperature)
    if (len(energy_bounds) != 2 or not all(isfinite(x) for x in energy_bounds)
            or energy_bounds[0] >= energy_bounds[1]):
        raise ValueError("energy_bounds must be finite and increasing")

    if initial_mole_fractions_a is None:
        fractions = [float(np.dot(np.arange(size + 1), row) / size) for row in series]
    else:
        fractions = list(initial_mole_fractions_a)
        if len(fractions) != len(series):
            raise ValueError("give one initial mole fraction per series")
    if any(not isfinite(p) or not 0 < p < 1 for p in fractions):
        raise ValueError("initial mole fractions must lie strictly between 0 and 1")

    initial_energy = (0.0 if energy_bounds[0] < 0 < energy_bounds[1]
                      else (energy_bounds[0] + energy_bounds[1]) / 2)
    initial = [initial_energy] + [log(p / (1 - p)) for p in fractions]
    lower = [energy_bounds[0]] + [-50.0] * len(series)
    upper = [energy_bounds[1]] + [50.0] * len(series)

    def residual(parameters: np.ndarray) -> np.ndarray:
        return np.concatenate([
            activity_distribution(
                size, parameters[0], g_aa=g_aa, g_bb=g_bb,
                temperature=temperature, log_activity_ratio=parameters[j + 1],
            ) - row
            for j, row in enumerate(series)
        ])

    fit = least_squares(residual, initial, bounds=(lower, upper),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12)
    predictions = tuple(activity_distribution(
        size, fit.x[0], g_aa=g_aa, g_bb=g_bb, temperature=temperature,
        log_activity_ratio=value,
    ) for value in fit.x[1:])
    fitted_fractions = tuple(
        float(np.dot(np.arange(size + 1), row) / size) for row in predictions
    )
    return FitResult(
        float(fit.x[0]), tuple(float(x) for x in fit.x[1:]),
        fitted_fractions, predictions,
        float(np.sqrt(np.mean(residual(fit.x) ** 2))), bool(fit.success),
    )
