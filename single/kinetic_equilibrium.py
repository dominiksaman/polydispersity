"""Mass-conserving equilibrium of the monomer-addition rate reference.

This solves equilibrium concentrations, not a time trajectory. k_plus is a
fixed bimolecular association constant; the effective on-rate k_plus*c_1
changes with total concentration through the free-monomer mass balance.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, log

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from .oligomer_distribution import _log_relative_weights


@dataclass(frozen=True)
class KineticModel:
    """k_plus in concentration^-1 time^-1, both off-rates in time^-1."""

    k_plus: float
    k_off_monomer: float = 10.0
    k_off_dimer: float = 1.0

    def __post_init__(self) -> None:
        if not all(isfinite(x) and x > 0 for x in
                   (self.k_plus, self.k_off_monomer, self.k_off_dimer)):
            raise ValueError("rate constants must be finite and positive")


@dataclass(frozen=True)
class KineticEquilibrium:
    sizes: np.ndarray
    concentration: np.ndarray
    log_concentration: np.ndarray
    number_fraction: np.ndarray
    subunit_fraction: np.ndarray
    free_monomer: float
    effective_on: float
    total_concentration: float

    def conditional_number_fraction(self, minimum_size: int = 3) -> np.ndarray:
        """Number fractions conditioned on n >= minimum_size, on the full axis."""
        if (not isinstance(minimum_size, int) or isinstance(minimum_size, bool)
                or not 1 <= minimum_size <= len(self.sizes)):
            raise ValueError("minimum_size must lie in the computed size range")
        weights = self.log_concentration.copy()
        weights[self.sizes < minimum_size] = -np.inf
        return np.exp(weights - logsumexp(weights))


def _validate_support(max_size: int, size_limit: int, tail_tolerance: float) -> None:
    if (not isinstance(max_size, int) or isinstance(max_size, bool) or max_size < 3
            or not isinstance(size_limit, int) or isinstance(size_limit, bool)
            or size_limit < max_size):
        raise ValueError("support requires integer size_limit >= max_size >= 3")
    if not isfinite(tail_tolerance) or not 0 < tail_tolerance < 1:
        raise ValueError("tail_tolerance must lie between zero and one")


def _boundary_converged(log_weights: np.ndarray, tolerance: float) -> bool:
    sizes = np.arange(1, len(log_weights) + 1)
    log_mass = np.log(sizes) + log_weights
    window = min(len(sizes), max(5, len(sizes) // 10))
    fraction = np.exp(logsumexp(log_mass[-window:]) - logsumexp(log_mass))
    return bool(fraction < tolerance and np.all(np.diff(log_weights[-window:]) < 0))


def kinetic_equilibrium(model: KineticModel, total_concentration: float, *,
                        max_size: int = 120, tail_tolerance: float = 1e-10,
                        size_limit: int = 3840) -> KineticEquilibrium:
    """Solve C_total = sum_n n*c_1*w_n(k_plus*c_1) in log concentration.

    All concentrations use the same unit as k_plus. The support expands until
    the declining boundary has negligible subunit fraction, with no 24-mer cap.
    """
    if not isfinite(total_concentration) or total_concentration <= 0:
        raise ValueError("total concentration must be finite and positive")
    _validate_support(max_size, size_limit, tail_tolerance)
    log_total = log(total_concentration)
    while True:
        sizes = np.arange(1, max_size + 1)
        # c_n = exp(n*log(c_1) + log_coefficient_n). Products of off-rates
        # supply the coefficient; its numerical value uses the chosen unit.
        coefficients = _log_relative_weights(log(model.k_plus), model.k_off_monomer,
                                              model.k_off_dimer, max_size)

        def residual(log_monomer):
            return logsumexp(np.log(sizes) + sizes * log_monomer + coefficients) - log_total

        # c_1 <= C_total. At the lower bound each mass term is below C_total/M.
        upper = log_total
        lower = float(np.min((log_total - log(max_size) - np.log(sizes)
                              - coefficients) / sizes)) - 1.0
        log_monomer = brentq(residual, lower, upper, xtol=1e-13)
        log_concentration = sizes * log_monomer + coefficients
        if _boundary_converged(log_concentration, tail_tolerance):
            mass_total = logsumexp(np.log(sizes) + log_concentration)
            return KineticEquilibrium(
                sizes, np.exp(log_concentration), log_concentration,
                np.exp(log_concentration - logsumexp(log_concentration)),
                np.exp(np.log(sizes) + log_concentration - mass_total),
                float(np.exp(log_monomer)),
                float(np.exp(log(model.k_plus) + log_monomer)),
                float(np.exp(mass_total)),
            )
        if max_size == size_limit:
            raise ValueError("size range has not converged; increase size_limit")
        max_size = min(2 * max_size, size_limit)


def anchor_kinetic_model(total_concentration: float = 3.0, effective_on: float = 20.0,
                         k_off_monomer: float = 10.0, k_off_dimer: float = 1.0, *,
                         max_size: int = 120, tail_tolerance: float = 1e-10,
                         size_limit: int = 3840) -> KineticModel:
    """Choose one k_plus reproducing a reference curve at one chosen total.

    c_1,anchor = C_anchor / sum(n*w_n), k_plus = effective_on/c_1,anchor.
    This sets a theoretical concentration scale, not a measured rate or a fit
    to the thermodynamic model. The returned parameters stay fixed thereafter.
    """
    if not all(isfinite(x) and x > 0 for x in
               (total_concentration, effective_on, k_off_monomer, k_off_dimer)):
        raise ValueError("concentration and rates must be finite and positive")
    _validate_support(max_size, size_limit, tail_tolerance)
    while True:
        weights = _log_relative_weights(log(effective_on), k_off_monomer,
                                        k_off_dimer, max_size)
        if _boundary_converged(weights, tail_tolerance):
            log_monomer = log(total_concentration) - logsumexp(
                np.log(np.arange(1, max_size + 1)) + weights)
            log_plus = log(effective_on) - log_monomer
            if log_plus > log(np.finfo(float).max):
                raise OverflowError("anchored k_plus exceeds floating-point range")
            return KineticModel(float(np.exp(log_plus)), k_off_monomer, k_off_dimer)
        if max_size == size_limit:
            raise ValueError("reference size range has not converged; increase size_limit")
        max_size = min(2 * max_size, size_limit)
