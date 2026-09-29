"""One reservoir figure; run python -m single.statistical_mechanics.simulate_reservoir."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from single.oligomer_distribution import oligomer_distribution
from .reservoir_thermodynamics import (
    GlobularModel, enrich_reservoir, globular_equilibrium, sample_reservoir,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--total", type=float, default=3.0, help="Total protein in subunit-equivalent micromolar")
    parser.add_argument("--standard", type=float, default=1.0, help="Standard concentration in micromolar; energies refer to this state")
    parser.add_argument("--bulk", type=float, default=8.0, help="Bulk attraction coefficient in kBT")
    parser.add_argument("--surface", type=float, default=8.0, help="Surface coefficient in kBT")
    parser.add_argument("--packing", type=float, default=0.18, help="Superlinear packing coefficient in kBT")
    parser.add_argument("--dimer", type=float, default=4.0, help="Dimer stabilization in kBT")
    parser.add_argument("--unpaired", type=float, default=0.3, help="Odd-size free-energy penalty in kBT")
    parser.add_argument("--assembly-offset", type=float, default=0.0, help="Formation offset for sizes >=3 in kBT")
    parser.add_argument("--pool-subunit-fraction", type=float, default=0.2,
                        help="Chosen monomer/dimer subunit fraction at --total; 0 uses unadjusted energies")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("reservoir_thermodynamics.png"))
    args = parser.parse_args()

    model = GlobularModel(args.bulk, args.surface, args.packing, args.dimer, args.unpaired,
                          args.assembly_offset)
    reference = globular_equilibrium(model, args.total, standard_concentration=args.standard)
    if args.pool_subunit_fraction != 0:
        model, result = enrich_reservoir(model, args.total, args.pool_subunit_fraction,
                                         standard_concentration=args.standard)
    else:
        result = reference
    n = result.sizes
    larger = result.conditional_number_fraction(3)
    snapshots = sample_reservoir(result)
    sampled = snapshots.mean(axis=0)
    sampled_total = sampled[n >= 3].sum()
    if sampled_total > 0:
        sampled /= sampled_total
    sizes, kinetic = oligomer_distribution(20, 10, 1, max_size=len(n))
    kinetic[sizes < 3] = 0
    kinetic /= kinetic.sum()
    displayed_size = max(40, int(n[np.searchsorted(np.cumsum(result.subunit_fraction), 0.99999)]) + 2)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), layout="constrained")
    first, second, third = axes
    colors = np.where(n <= 2, "#b45f22", "#277da8")
    first.bar(n, result.concentration, color=colors, width=0.82)
    first.set(yscale="log", xlim=(0.4, displayed_size + 0.5),
              ylim=(max(args.total * 1e-7, 1e-300), args.total),
              xlabel="Monomers per species", ylabel="Species concentration (µM)",
              title="Small pool and larger oligomers")
    first.plot(reference.sizes, reference.concentration, color="#777777", linewidth=1,
                linestyle="--", label="Original pool scenario")
    first.plot([], [], color="#b45f22", linewidth=6, label="Monomer / dimer pool")
    first.plot([], [], color="#277da8", linewidth=6, label="Sizes ≥3")
    first.legend(frameon=False, fontsize=8)

    second.bar(n[n >= 3], larger[n >= 3], color="#277da8", alpha=0.75,
               label="Thermodynamic model")
    second.plot(n[n >= 3], kinetic[n >= 3], color="#393939", linewidth=1.4,
                label="Rate-model equilibrium (20,10,1)")
    if sampled_total > 0:
        second.scatter(n[(n >= 3) & (larger > 1e-5)], sampled[(n >= 3) & (larger > 1e-5)],
                       s=11, color="#ba4b32", label="Exact reservoir snapshots", zorder=3)
    second.set(xlim=(2.5, displayed_size + 0.5), xlabel="Monomers per oligomer",
               ylabel="Number fraction conditioned on size ≥3",
               title="Larger-species shape preserved")
    second.legend(frameon=False, fontsize=7)

    totals = np.geomspace(args.total / 1000, args.total * 10, 80)
    equilibria = [globular_equilibrium(model, total, standard_concentration=args.standard)
                  for total in totals]
    third.semilogx(totals, [r.concentration[0] for r in equilibria],
                   color="#b45f22", label="Free monomer")
    third.semilogx(totals, [r.concentration[1] for r in equilibria],
                   color="#7d5c95", label="Free dimer")
    third.set(xlabel="Total protein, subunit equivalents (µM)",
              ylabel="Small-species concentration (µM)",
              title="Pool grows slowly after assembly")
    third.axvline(args.total, color="#666666", linestyle=":", linewidth=1)
    third.legend(frameon=False, fontsize=8, loc="upper left")
    pool_mass = float(result.subunit_fraction[:2].sum())
    pool_number = float(result.number_fraction[:2].sum())
    fig.suptitle(f"Equilibrium cluster gas: total {args.total:g} µM; {100 * pool_mass:.1f}% of subunits in monomers/dimers")
    fig.supxlabel(f"Chosen equilibrium scenario: the pool holds {100 * pool_number:.1f}% of species by number. Larger-size panel excludes sizes 1 and 2.",
                  fontsize=9, color="#4a5560")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)

    print(f"Wrote {args.output}")
    print(f"Mass balance: {np.dot(n, result.concentration):.10g} / {args.total:g} µM")
    print(f"Free monomer {result.concentration[0]:.6g} µM; dimer {result.concentration[1]:.6g} µM")
    print(f"Pool subunit fraction {pool_mass:.6f}; pool number fraction {pool_number:.6f}")
    print(f"Adjusted bulk {model.bulk_kbt:.6f} kBT; assembly offset {model.assembly_offset_kbt:.6f} kBT")
    print(f"Higher-size mode {n[np.argmax(larger)]}; mean {n @ larger:.4f}; "
          f"conditional P(n>24) {larger[n > 24].sum():.4f}")
    print(f"Subunit fraction in sizes >=3 {result.subunit_fraction[n >= 3].sum():.4f}")
    print(f"Computed up to {len(n)}; high-size boundary subunit fraction {result.subunit_fraction[-12:].sum():.2e}")


if __name__ == "__main__":
    main()
