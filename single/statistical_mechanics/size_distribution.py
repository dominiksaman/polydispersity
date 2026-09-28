"""Statistical weights for a single protein assembling into variable sizes.

There are two distinct calculations here:

* ``site_count_distribution`` is the exact equilibrium partition weight
  implied by the thesis/Baldwin model. It includes the number of possible
  dissociation pathways at each size. Its closed form uses gamma functions,
  so it is an independent numerical check of the kinetic recursion, but it
  is not an independent physical model.
* ``additive_contact_distribution`` is a simple grand-canonical contact
  model. It is an exploratory alternative, not the Baldwin partition
  function. Without size-dependent entropy it grows to the chosen cutoff;
  an optional 1/n! factor makes it peaked but does not preserve the
  microscopic energy-to-rate mapping when paired and unpaired off-rates
  differ.

Each output is normalized over oligomer sizes 1 through ``max_size``. The
effective monomer concentration is already included in the edge free energy.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite, lgamma, log
from typing import Sequence

import numpy as np

R_KJ_MOL_K = 0.00831446261815324


@dataclass(frozen=True)
class ContactEnergies:
    """Microscopic specific-site addition energies in kJ/mol."""

    edge_kj_mol: float
    dimer_kj_mol: float


@dataclass(frozen=True)
class IdealClusterFit:
    """A descriptive fit of the factorial contact model to a target curve."""

    edge_kj_mol: float
    dimer_kj_mol: float
    distribution: np.ndarray
    total_variation: float
    success: bool


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
        edge_kj_mol=-rt * log(effective_on_rate / unpaired_off_rate),
        dimer_kj_mol=-rt * log(unpaired_off_rate / paired_off_rate),
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
    if log_r > 700:
        raise ValueError("dimer energy implies an unrepresentably large rate ratio")
    r = exp(log_r)
    gamma_reference = lgamma(1.0 + r / 2.0)
    log_weights = np.empty(max_size)
    for index in range(max_size):
        size = index + 1
        pairs = size // 2
        if size % 2 == 0:
            gamma_odd = lgamma(pairs + r / 2.0)
        else:
            gamma_odd = lgamma(pairs + 1.0 + r / 2.0)
        log_weights[index] = (
            (size - 1) * (log_alpha - log(2.0))
            - lgamma(pairs + 1.0) - gamma_odd + gamma_reference
        )
    return _normalise(log_weights)


def additive_contact_distribution(
    edge_kj_mol: float,
    dimer_kj_mol: float,
    *,
    temperature: float = 298.15,
    max_size: int = 60,
    factorial_entropy: bool = False,
) -> np.ndarray:
    """Exploratory contact model with an optional ideal-cluster 1/i! factor.

    Its weights are proportional to

        exp[-((i-1)G_edge + floor(i/2)G_dimer)/RT] / i!

    when ``factorial_entropy`` is true, and omit ``i!`` otherwise. The
    factorial variant is Poisson-like at G_dimer=0, but the independent
    contact model does not encode the exact paired/unpaired exit pathways.
    """
    _validate(max_size, edge_kj_mol, dimer_kj_mol, temperature)
    rt = R_KJ_MOL_K * temperature
    sizes = np.arange(1, max_size + 1)
    log_weights = -((sizes - 1) * edge_kj_mol
                    + (sizes // 2) * dimer_kj_mol) / rt
    if factorial_entropy:
        log_weights -= np.array([lgamma(int(size) + 1) for size in sizes])
    return _normalise(log_weights)


def fit_ideal_cluster(
    target: Sequence[float], *, temperature: float = 298.15,
) -> IdealClusterFit:
    """Fit both effective contact energies of the factorial toy model.

    The fit minimizes KL divergence to the normalized target. Its energies
    are descriptive and should not be interpreted as the microscopic rates'
    specific-site free energies.
    """
    try:
        from scipy.optimize import minimize
    except ImportError as exc:
        raise ImportError("fit_ideal_cluster requires scipy") from exc

    target_array = np.asarray(target, dtype=float)
    if (target_array.ndim != 1 or len(target_array) < 3
            or not np.all(np.isfinite(target_array))
            or np.any(target_array < 0) or target_array.sum() <= 0):
        raise ValueError("target must be a nonzero, finite 1D distribution")
    if not isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature must be finite and positive")
    target_array = target_array / target_array.sum()
    sizes = np.arange(1, len(target_array) + 1)
    features = np.column_stack((sizes - 1, sizes // 2))
    log_factorial = np.array([lgamma(int(size) + 1) for size in sizes])
    target_features = target_array @ features

    def objective(parameters: np.ndarray):
        log_weights = features @ parameters - log_factorial
        maximum = float(np.max(log_weights))
        weights = np.exp(log_weights - maximum)
        log_partition = maximum + log(float(weights.sum()))
        probabilities = weights / weights.sum()
        loss = log_partition - float(target_array @ log_weights)
        gradient = probabilities @ features - target_features
        return loss, gradient

    result = minimize(objective, np.array([log(2.0), 0.0]), jac=True,
                      method="L-BFGS-B", bounds=((-20, 20), (-20, 20)),
                      options={"ftol": 1e-14, "gtol": 1e-10})
    rt = R_KJ_MOL_K * temperature
    edge_energy, dimer_energy = -rt * result.x
    distribution = additive_contact_distribution(
        float(edge_energy), float(dimer_energy), temperature=temperature,
        max_size=len(target_array), factorial_entropy=True,
    )
    total_variation = float(np.abs(distribution - target_array).sum() / 2)
    return IdealClusterFit(float(edge_energy), float(dimer_energy),
                           distribution, total_variation, bool(result.success))
