"""Quantify how well the current equilibrium scaffolds mimic kinetic curves.

This is calibration to a synthetic target, not a derivation of the kinetic
law or a fit to experimental native MS measurements.

Run from the repository root:

    python -m single.statistical_mechanics.fit_scaffold_to_kinetics
"""

from math import log
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from single.oligomer_distribution import oligomer_distribution
from single.statistical_mechanics.equilibrium_scaffold import (
    EquilibriumParameters, StateCatalogue, enumerate_states,
)
from single.statistical_mechanics.geometry import cube, dimer_ring, octahedron


# Argument order: effective on-rate, unpaired off-rate, paired off-rate.
BASE_RATES = (20.0, 10.0, 1.0)
MAX_SIZE = 60
PARENT_SIZE = 24


def fit_catalogue(catalogue: StateCatalogue,
                  target: np.ndarray) -> tuple[EquilibriumParameters, np.ndarray]:
    """Fit four shared dimensionless terms, with stabilizations nonnegative."""
    if target.shape != (catalogue.max_size,) or np.any(target < 0):
        raise ValueError("target must be a nonnegative distribution on this scaffold")
    target = target / target.sum()

    def loss(values: np.ndarray) -> float:
        predicted = catalogue.distribution(EquilibriumParameters(*values))
        return -float(target @ np.log(np.maximum(predicted, 1e-300)))

    bounds = ((-10, 10), (0, 10), (0, 10), (0, 10))
    starts = ((-1, 3, 0.1, 0.1), (-2, 4, 1, 1), (0, 1, 0, 0))
    results = [minimize(loss, start, method="L-BFGS-B", bounds=bounds)
               for start in starts]
    best = min(results, key=lambda result: result.fun)
    if not best.success:
        raise RuntimeError(f"fit failed for {catalogue.scaffold_name}: {best.message}")
    parameters = EquilibriumParameters(*map(float, best.x))
    return parameters, catalogue.distribution(parameters)


def total_variation(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.abs(left - right).sum() / 2)


def main() -> None:
    sizes, full_kinetic = oligomer_distribution(*BASE_RATES, max_size=MAX_SIZE)
    conditional_target = full_kinetic[:PARENT_SIZE]
    conditional_target = conditional_target / conditional_target.sum()
    missing_tail = float(full_kinetic[PARENT_SIZE:].sum())
    catalogues = [enumerate_states(scaffold)
                  for scaffold in (dimer_ring(12), cube(), octahedron())]
    fits = [(catalogue, *fit_catalogue(catalogue, conditional_target))
            for catalogue in catalogues]

    fig, (ax_full, ax_residual) = plt.subplots(
        2, 1, figsize=(10, 7.2), constrained_layout=True)
    ax_full.bar(sizes, full_kinetic, width=0.9, color="#b8d9ed",
                label="Kinetic target")
    colours = {"12-dimer ring": "#c44e52", "cube": "#8172b2",
               "octahedron": "#55a868"}
    print(f"Kinetic target mode: {int(np.argmax(full_kinetic)) + 1}-mer")
    print(f"Kinetic probability above 24 monomers: {missing_tail:.6f}")
    print("Parent      parameters (log z, eps_d, eps_C, kappa)"
          "        TV n<=24  TV full  mode")
    for catalogue, parameters, predicted in fits:
        padded = np.zeros(MAX_SIZE)
        padded[:PARENT_SIZE] = predicted
        colour = colours[catalogue.scaffold_name]
        ax_full.plot(sizes, padded, marker="o", markersize=2.4,
                     linewidth=1.5, color=colour, label=catalogue.scaffold_name)
        ax_residual.plot(sizes[:PARENT_SIZE], predicted - conditional_target,
                         marker="o", markersize=2.4, linewidth=1.5,
                         color=colour, label=catalogue.scaffold_name)
        values = tuple(vars(parameters).values())
        print(f"{catalogue.scaffold_name:13} "
              f"{str(tuple(round(x, 3) for x in values)):38} "
              f"{total_variation(predicted, conditional_target):.4f}    "
              f"{total_variation(padded, full_kinetic):.4f}    "
              f"{int(np.argmax(predicted)) + 1}")

    ax_full.axvspan(24.5, 30.5, color="#f1e5dc", alpha=0.5)
    ax_full.set(xlim=(0, 31), ylabel="Fraction of oligomer molecules",
                title=f"Kinetic example ({BASE_RATES[0]:g}, "
                      f"{BASE_RATES[1]:g}, {BASE_RATES[2]:g}): "
                      f"{100 * (1 - missing_tail):.1f}% at sizes ≤24")
    ax_full.legend(frameon=False, fontsize=8)
    ax_residual.axhline(0, color="#666666", linewidth=0.8)
    ax_residual.set(xlim=(0.5, 24.5), xlabel="Oligomer size (monomers)",
                    ylabel="Equilibrium minus kinetic fraction\n(conditional on n ≤ 24)")
    ax_residual.set_xticks(np.arange(2, 25, 2))
    output = Path(__file__).with_name("scaffold_kinetic_fit.png")
    fig.savefig(output, dpi=180)
    print(f"Wrote {output}")

    print(f"\nHoldout test: fit at effective on-rate {BASE_RATES[0]:g}; "
          "change only log(z) by log(new on-rate / baseline)")
    print("Parent         on-rate  TV on n<=24  predicted mode  kinetic mode")
    for catalogue, parameters, _ in fits:
        for new_on_rate in (16.0, 18.0):
            _, target = oligomer_distribution(
                new_on_rate, BASE_RATES[1], BASE_RATES[2], max_size=PARENT_SIZE)
            shifted = EquilibriumParameters(
                parameters.log_monomer_activity + log(new_on_rate / BASE_RATES[0]),
                parameters.dimer_stabilization,
                parameters.c_terminal_stabilization,
                parameters.saturated_vertex_penalty,
            )
            predicted = catalogue.distribution(shifted)
            print(f"{catalogue.scaffold_name:15} {new_on_rate:5.0f}"
                  f"       {total_variation(predicted, target):.4f}"
                  f"              {int(np.argmax(predicted)) + 1:2d}"
                  f"              {int(np.argmax(target)) + 1:2d}")


if __name__ == "__main__":
    main()
