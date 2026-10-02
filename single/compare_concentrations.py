"""Theoretical concentration series with independently fixed model parameters.

Run: python -m single.compare_concentrations
No measured abundances, fitting, sampling, or time integration are used.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .kinetic_equilibrium import KineticModel, KineticEquilibrium, anchor_kinetic_model, kinetic_equilibrium
from .thermodynamics.reservoir_thermodynamics import (
    Equilibrium, GlobularModel, enrich_reservoir, globular_equilibrium,
)


@dataclass(frozen=True)
class ConcentrationSeries:
    totals: np.ndarray
    anchor_total: float
    kinetic_model: KineticModel
    thermodynamic_model: GlobularModel
    kinetic: tuple[KineticEquilibrium, ...]
    thermodynamic: tuple[Equilibrium, ...]


def predict_series(totals, *, anchor_total: float = 3.0,
                   pool_subunit_fraction: float = 0.2) -> ConcentrationSeries:
    """Set each model at one reference total, then freeze it for the whole series.

    The kinetic anchor reproduces effective_on=20, unpaired_off=10, paired_off=1.
    Thermodynamics uses the existing illustrative globular energies and one
    reservoir enrichment at the anchor. Neither model is fitted to the other.
    Concentrations are interpreted in µM; time units remain unspecified.
    """
    totals = np.asarray(totals, dtype=float)
    if (totals.ndim != 1 or len(totals) == 0 or not np.all(np.isfinite(totals))
            or np.any(totals <= 0)):
        raise ValueError("totals must be a nonempty vector of finite positive concentrations")
    kinetic_model = anchor_kinetic_model(anchor_total)
    thermodynamic_model, _ = enrich_reservoir(
        GlobularModel(), anchor_total, pool_subunit_fraction,
    )
    return ConcentrationSeries(
        totals.copy(), anchor_total, kinetic_model, thermodynamic_model,
        tuple(kinetic_equilibrium(kinetic_model, float(c)) for c in totals),
        tuple(globular_equilibrium(thermodynamic_model, float(c)) for c in totals),
    )


def _summary(result):
    p = result.conditional_number_fraction()
    mean = float(result.sizes @ p)
    return (int(result.sizes[np.argmax(p)]), mean,
            float(np.sqrt((result.sizes - mean) ** 2 @ p)),
            float(p[result.sizes > 24].sum()),
            float(result.subunit_fraction[:2].sum()))


def draw_series(series: ConcentrationSeries, output: Path) -> None:
    """One summary figure; heatmaps show conditional larger-species populations."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import FuncFormatter, PercentFormatter

    colors = ("#247b89", "#bd5335")
    families = (series.thermodynamic, series.kinetic)
    labels = ("Thermodynamic", "Rate-model equilibrium")
    fig, axes = plt.subplots(2, 3, figsize=(15.3, 9.2), layout="constrained")
    totals = series.totals
    display_max = 40  # Display only; all predictions use converged larger supports.
    shown_sizes = np.arange(3, display_max + 1)
    for ax, results, label in zip(axes[0, :2], families, labels):
        probabilities = np.zeros((len(totals), len(shown_sizes)))
        for row, result in enumerate(results):
            p = result.conditional_number_fraction()
            mask = shown_sizes <= len(p)
            probabilities[row, mask] = p[shown_sizes[mask] - 1]
        mesh = ax.pcolormesh(
            shown_sizes, totals, np.ma.masked_less(probabilities, 1e-4),
            shading="nearest", norm=LogNorm(1e-4, 1), cmap="magma",
            rasterized=True,
        )
        ax.axhline(series.anchor_total, color="white", linestyle=":", linewidth=1.3)
        ax.set_yscale("log")
        ax.set_xlabel("Monomers per oligomer (displayed sizes)")
        ax.set_ylabel("Total protein (µM of subunits)")
        ax.set_title(label + ": larger-species distribution")
        ax.set_xlim(2.5, display_max + .5)
        ax.set_ylim(totals.min(), totals.max())
    fig.colorbar(mesh, ax=list(axes[0, :2]), shrink=.88,
                 label="Number fraction conditioned on n ≥ 3 (below 10⁻⁴ hidden)")

    for results, label, color in zip(families, labels, colors):
        summaries = np.array([_summary(r) for r in results])
        _, means, widths, tails, pools = summaries.T
        axes[0, 2].loglog(totals, pools, color=color, label=label, linewidth=2)
        axes[1, 0].semilogx(totals, means, color=color, label=label, linewidth=2)
        axes[1, 0].fill_between(totals, means - widths, means + widths,
                                color=color, alpha=.16)
        axes[1, 1].semilogx(totals, tails, color=color, label=label, linewidth=2)
        axes[1, 2].loglog(totals, [r.concentration[0] for r in results],
                          color=color, linewidth=2, label=label + ": monomer")
        axes[1, 2].loglog(totals, [r.concentration[1] for r in results],
                          color=color, linestyle="--", linewidth=2, label=label + ": dimer")
    axes[0, 2].set_title("Small-species pool (all subunits)")
    axes[0, 2].set_ylabel("Subunit fraction in 1- and 2-mers")
    axes[0, 2].yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{100*y:g}%"))
    axes[0, 2].set_ylim(1e-6, 1.15)
    axes[0, 2].legend(frameon=False, fontsize=9)
    axes[1, 0].set_title("Larger-species size: mean ± standard deviation")
    axes[1, 0].set_ylabel("Monomers per oligomer, conditioned on n ≥ 3")
    axes[1, 0].legend(frameon=False, fontsize=9)
    axes[1, 1].set_title("Larger-species tail above 24")
    axes[1, 1].set_ylabel("Number fraction, conditioned on n ≥ 3")
    axes[1, 1].yaxis.set_major_formatter(PercentFormatter(1))
    axes[1, 1].set_ylim(bottom=0)
    axes[1, 2].set_title("Free monomer and dimer concentrations")
    axes[1, 2].set_ylabel("Species concentration (µM)")
    axes[1, 2].legend(frameon=False, fontsize=8.5, loc="center right")
    for ax in (axes[0, 2], *axes[1]):
        ax.axvline(series.anchor_total, color="0.5", linestyle=":", linewidth=1)
        ax.set_xlabel("Total protein (µM of subunits)")
        ax.set_xlim(totals.min(), totals.max())
        ax.grid(alpha=.15)
    fig.suptitle(
        "Theoretical equilibrium across concentration — fixed parameters\n"
        f"Reference at {series.anchor_total:g} µM: rate-model effective on = 20; "
        "thermodynamic 1/2-mer pool chosen once",
        fontsize=15,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor-total", type=float, default=3)
    parser.add_argument("--minimum-total", type=float, default=.003)
    parser.add_argument("--maximum-total", type=float, default=30)
    parser.add_argument("--pool-subunit-fraction", type=float, default=.2)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("concentration_comparison.png"))
    args = parser.parse_args()
    if (not np.isfinite(args.minimum_total) or not np.isfinite(args.maximum_total)
            or not 0 < args.minimum_total < args.maximum_total
            or not args.minimum_total <= args.anchor_total <= args.maximum_total):
        parser.error("require 0 < minimum_total < maximum_total with anchor_total in the range")
    # Include the reference exactly; the remaining points are log spaced.
    totals = np.unique(np.append(np.geomspace(args.minimum_total, args.maximum_total, 81),
                                  args.anchor_total))
    try:
        series = predict_series(totals, anchor_total=args.anchor_total,
                                pool_subunit_fraction=args.pool_subunit_fraction)
    except (ValueError, OverflowError) as error:
        parser.error(str(error))
    draw_series(series, args.output)
    print("THEORETICAL SCENARIO: no experimental data or parameter fitting")
    print("Fixed rate parameters:", series.kinetic_model)
    print("Fixed thermodynamic parameters:", series.thermodynamic_model)
    print("C_total  model     c1 (µM)     c2 (µM)    pool mass    mode>=3  mean>=3  SD>=3  tail>24>=3")
    for i in np.unique(np.linspace(0, len(totals) - 1, 5).round().astype(int)):
        for label, result in (("thermo", series.thermodynamic[i]), ("rate", series.kinetic[i])):
            mode, mean, width, tail, pool = _summary(result)
            print(f"{totals[i]:7.4g}  {label:6s}  {result.concentration[0]:11.5g}  "
                  f"{result.concentration[1]:11.5g}  {pool:10.4%}  "
                  f"{mode:7d}  {mean:7.3f}  {width:5.3f}  {tail:9.4%}")
    error = max(abs(r.sizes @ r.concentration / c - 1)
                for results in (series.kinetic, series.thermodynamic)
                for c, r in zip(totals, results))
    print(f"Maximum relative mass-balance error: {error:.3g}")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
