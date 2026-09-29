"""Closed-form equilibrium weights implied by the specified rate model.

The log-gamma formula checks the kinetic recurrence independently as a
calculation. It uses the same physical assumptions; it supplies neither a
new geometric mechanism nor independently measured contact free energies.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite, lgamma, log

import numpy as np

R_KJ_MOL_K = 0.00831446261815324


@dataclass(frozen=True)
class ContactEnergies:
    """Rate-derived addition parameters in kJ/mol, including monomer activity.

    ``edge_kj_mol`` is not a standard-state C-terminal interface energy.
    """

    edge_kj_mol: float
    dimer_kj_mol: float


def _validate(max_size: int, edge_kj_mol: float, dimer_kj_mol: float,
              temperature: float) -> None:
    if not isinstance(max_size, int) or isinstance(max_size, bool) or max_size < 1:
        raise ValueError("max_size must be a positive integer")
    if not all(isfinite(value) for value in
               (edge_kj_mol, dimer_kj_mol, temperature)):
        raise ValueError("energies and temperature must be finite")
    if temperature <= 0:
        raise ValueError("temperature must be positive in kelvin")


def _normalise(log_weights: np.ndarray) -> np.ndarray:
    weights = np.exp(log_weights - np.max(log_weights))
    return weights / weights.sum()


def energies_from_rates(
    effective_on_rate: float,
    unpaired_off_rate: float,
    paired_off_rate: float,
    temperature: float = 298.15,
) -> ContactEnergies:
    """Map the kinetic rates onto the thesis's specific-site free energies.

    ``effective_on_rate`` is k+ times free-monomer concentration. It is not
    the bimolecular k+ alone. The paired off-rate is k_d and the unpaired
    off-rate is k_m in ``single/oligomer_distribution.py``.
    """
    if not all(isfinite(value) and value > 0 for value in
               (effective_on_rate, unpaired_off_rate, paired_off_rate,
                temperature)):
        raise ValueError("rates and temperature must be finite and positive")
    rt = R_KJ_MOL_K * temperature
    return ContactEnergies(
        edge_kj_mol=-rt * (log(effective_on_rate) - log(unpaired_off_rate)),
        dimer_kj_mol=-rt * (log(unpaired_off_rate) - log(paired_off_rate)),
    )


def site_count_distribution(
    edge_kj_mol: float,
    dimer_kj_mol: float,
    *,
    temperature: float = 298.15,
    max_size: int = 60,
) -> np.ndarray:
    """Exact equilibrium size weights with the oligomer dissociation counts.

    Write alpha = exp[-(G_edge + G_dimer)/RT] = k_on_eff/k_d and
    r = exp[-G_dimer/RT] = k_m/k_d. If w_1=1, detailed balance gives

        w_i/w_(i-1) = alpha/i              (i even)
        w_i/w_(i-1) = alpha/(i-1+r)        (i odd).

    The implementation evaluates the product in closed form using log-gamma
    identities rather than repeating the kinetic recursion. This is an
    *equilibrium reformulation* of the same rate model, not independent
    thermodynamic evidence for it.
    """
    _validate(max_size, edge_kj_mol, dimer_kj_mol, temperature)
    rt = R_KJ_MOL_K * temperature
    log_alpha = -(edge_kj_mol + dimer_kj_mol) / rt
    log_r = -dimer_kj_mol / rt
    if not all(isfinite(x) for x in (log_alpha, log_r)):
        raise ValueError("energy/temperature ratios exceed the numerical range")
    # Large r causes catastrophic cancellation between log-gamma values.
    # Evaluate the same rising-factorial ratio as a finite product in logs.
    if log_r > log(1e5):
        factors = np.logaddexp(log_r - log(2), np.log(np.arange(1, max_size // 2 + 1)))
        rising = np.concatenate(([0.0], np.cumsum(factors)))
        r = None
    else:
        r = exp(log_r)
        gamma_reference = lgamma(1.0 + r / 2.0)
    log_weights = np.empty(max_size)
    for index in range(max_size):
        size = index + 1
        pairs = size // 2
        if r is None:
            gamma_difference = rising[pairs - int(size % 2 == 0)]
        elif size % 2 == 0:
            gamma_difference = lgamma(pairs + r / 2.0) - gamma_reference
        else:
            gamma_difference = lgamma(pairs + 1.0 + r / 2.0) - gamma_reference
        log_weights[index] = (
            (size - 1) * (log_alpha - log(2.0))
            - lgamma(pairs + 1.0) - gamma_difference
        )
    if not np.all(np.isfinite(log_weights)):
        raise ValueError("partition weights exceed the numerical range")
    return _normalise(log_weights)
