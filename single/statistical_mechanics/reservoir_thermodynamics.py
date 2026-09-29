"""Ideal cluster thermodynamics with a self-consistent small-species pool.

Inspired by Dear et al., J. Phys. Chem. B (2018), DOI 10.1021/acs.jpcb.8b07805,
especially Supporting Information Eqs. 2--4 and 6. This implements equilibrium
cluster concentrations, not the paper's spatial amphiphile Monte Carlo code.
All formation free energies share the free-monomer reference and standard state.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, log

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp


@dataclass(frozen=True)
class GlobularModel:
    """An illustrative formation-free-energy family, in kBT at fixed T.

    Bulk attraction, exposed surface, and a superlinear packing/connectivity
    cost give -b*(n-1) + s*(n**(2/3)-1) + h*(n**(5/3)-1). The exponents follow
    the globular model; the chosen coefficients are not measured for sHSPs.
    A separate dimer free energy and an odd-size cost are native-protein
    extensions. The odd cost describes an unpaired-site penalty on the smooth
    envelope; it is not a dissociation rate or a count of inferred interfaces.
    """

    bulk_kbt: float = 8.0
    surface_kbt: float = 8.0
    packing_kbt: float = 0.18
    dimer_binding_kbt: float = 4.0
    unpaired_penalty_kbt: float = 0.3

    def __post_init__(self) -> None:
        values = (self.bulk_kbt, self.surface_kbt, self.packing_kbt,
                  self.dimer_binding_kbt, self.unpaired_penalty_kbt)
        if not all(isfinite(x) and x >= 0 for x in values) or self.packing_kbt == 0:
            raise ValueError("coefficients must be finite and nonnegative; packing must be positive")

    def free_energies(self, max_size: int = 120) -> np.ndarray:
        """Dimensionless standard formation energies Delta G_n^0/kBT; F_1=0."""
        if not isinstance(max_size, int) or isinstance(max_size, bool) or max_size < 2:
            raise ValueError("max_size must be an integer of at least two")
        sizes = np.arange(1, max_size + 1)
        energy = (-self.bulk_kbt * (sizes - 1)
                  + self.surface_kbt * (sizes ** (2 / 3) - 1)
                  + self.packing_kbt * (sizes ** (5 / 3) - 1))
        energy[(sizes >= 3) & (sizes % 2 == 1)] += self.unpaired_penalty_kbt
        energy[0] = 0.0
        energy[1] = -self.dimer_binding_kbt
        return energy


@dataclass(frozen=True)
class Equilibrium:
    sizes: np.ndarray
    free_energy_kbt: np.ndarray
    concentration: np.ndarray
    number_fraction: np.ndarray
    subunit_fraction: np.ndarray
    standard_concentration: float
    log_activity: float
    total_concentration: float

    def conditional_number_fraction(self, minimum_size: int = 3) -> np.ndarray:
        """Number fractions conditioned on n >= minimum_size, on the full axis."""
        if (not isinstance(minimum_size, int) or isinstance(minimum_size, bool)
                or not 1 <= minimum_size <= len(self.sizes)):
            raise ValueError("minimum_size must lie in the computed size range")
        log_weights = self.sizes * self.log_activity - self.free_energy_kbt
        log_weights[self.sizes < minimum_size] = -np.inf
        return np.exp(log_weights - logsumexp(log_weights))


def _energies(free_energy_kbt) -> np.ndarray:
    energy = np.asarray(free_energy_kbt, dtype=float)
    if (energy.ndim != 1 or len(energy) == 0 or not np.all(np.isfinite(energy))
            or energy[0] != 0):
        raise ValueError("provide finite formation energies for sizes 1..M, with F_1=0")
    return energy


def equilibrium_at_activity(free_energy_kbt, log_activity: float, *,
                            standard_concentration: float = 1.0) -> Equilibrium:
    """Open, clamped reservoir: c_n = c0 exp(n*log(z) - F_n), with z=c_1/c0.

    c0 and returned concentrations use the same user-chosen concentration unit.
    This is an ideal mixture of clusters; interactions between clusters are absent.
    Monomer and dimer chemical potentials are linked, not independently set.
    """
    energy = _energies(free_energy_kbt)
    if (not isfinite(log_activity) or not isfinite(standard_concentration)
            or standard_concentration <= 0):
        raise ValueError("provide finite log activity and positive standard concentration")
    sizes = np.arange(1, len(energy) + 1)
    log_relative = sizes * log_activity - energy
    log_concentration = log(standard_concentration) + log_relative
    log_total = logsumexp(np.log(sizes) + log_concentration)
    if not isfinite(log_total) or log_total > log(np.finfo(float).max):
        raise OverflowError("predicted concentration exceeds floating-point range")
    concentration = np.exp(log_concentration)
    return Equilibrium(
        sizes, energy.copy(), concentration,
        np.exp(log_relative - logsumexp(log_relative)),
        np.exp(np.log(sizes) + log_concentration - log_total),
        standard_concentration, log_activity, float(np.exp(log_total)),
    )


def closed_equilibrium(free_energy_kbt, total_concentration: float, *,
                       standard_concentration: float = 1.0) -> Equilibrium:
    """Solve sum(n*c_n)=C_total for the shared monomer chemical potential.

    Minimizes the ideal-mixture free energy with a mass constraint in the
    macroscopic limit. The supplied support is a numerical truncation.
    """
    energy = _energies(free_energy_kbt)
    if not all(isfinite(x) and x > 0 for x in
               (total_concentration, standard_concentration)):
        raise ValueError("total and standard concentrations must be finite and positive")
    sizes = np.arange(1, len(energy) + 1)
    log_total = log(total_concentration) - log(standard_concentration)

    def residual(log_activity):
        return logsumexp(np.log(sizes) + sizes * log_activity - energy) - log_total

    # Upper bound follows from c_1 <= C_total. The lower bound makes every
    # subunit-concentration term smaller than C_total/M.
    upper = log_total
    lower = float(np.min((log_total - log(len(energy)) - np.log(sizes) + energy)
                         / sizes)) - 1.0
    root = brentq(residual, lower, upper, xtol=1e-13)
    return equilibrium_at_activity(energy, root,
                                   standard_concentration=standard_concentration)


def globular_equilibrium(model: GlobularModel, total_concentration: float, *,
                         standard_concentration: float = 1.0, max_size: int = 120,
                         tail_tolerance: float = 1e-10,
                         size_limit: int = 3840) -> Equilibrium:
    """Extend the size range until its declining boundary has negligible mass.

    Positive packing penalizes unbounded growth. This boundary check is a
    numerical convergence criterion, not a physical maximum oligomer size.
    """
    if not isfinite(tail_tolerance) or not 0 < tail_tolerance < 1:
        raise ValueError("tail_tolerance must lie between zero and one")
    if (not isinstance(size_limit, int) or isinstance(size_limit, bool)
            or size_limit < max_size):
        raise ValueError("size_limit must be an integer at least max_size")
    while True:
        result = closed_equilibrium(model.free_energies(max_size), total_concentration,
                                    standard_concentration=standard_concentration)
        window = min(max_size, max(5, max_size // 10))
        log_weights = result.sizes * result.log_activity - result.free_energy_kbt
        if (result.subunit_fraction[-window:].sum() < tail_tolerance
                and np.all(np.diff(log_weights[-window:]) < 0)):
            return result
        if max_size == size_limit:
            raise ValueError("size range has not converged; increase size_limit")
        max_size = min(2 * max_size, size_limit)


def sample_reservoir(result: Equilibrium, *, standard_state_particles: float = 10000,
                     snapshots: int = 2000, seed: int = 7) -> np.ndarray:
    """Exact ideal grand-canonical cluster-count snapshots, shape (snapshots, M).

    standard_state_particles = N_A*V*c0 in consistent molar/volume units.
    N_n is Poisson with mean (N_A*V*c0)*exp(n*log(z)-F_n). This samples
    equilibrium fluctuations at the solved/clamped chemical potential. A finite
    snapshot's protein count fluctuates; it is not a fixed-N particle trajectory.
    """
    if not isfinite(standard_state_particles) or standard_state_particles <= 0:
        raise ValueError("standard_state_particles must be finite and positive")
    if not isinstance(snapshots, int) or isinstance(snapshots, bool) or snapshots < 1:
        raise ValueError("snapshots must be a positive integer")
    mean = standard_state_particles * result.concentration / result.standard_concentration
    return np.random.default_rng(seed).poisson(mean, size=(snapshots, len(mean)))
