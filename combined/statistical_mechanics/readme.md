# Circular co-assembly of two proteins

This implements the **fixed-size ring** statistical-mechanics model described
in Chapter 2 of the thesis, `text/ch1-intro.tex`, under “Description of
co-assembly between two monodisperse proteins using statistical mechanics”.
It is a reconstruction from the thesis equations and counting description;
the original analysis script was not present in the attached folder.
It is separate from the repository's `single/` model, which predicts a range
of oligomer sizes from chemical kinetics.

For a ring of `N` subunits, each arrangement has `N` neighbour contacts,
including the contact between the last and first subunits. Its weight is

`multiplicity × (a/b)^number_of_A × exp[-(nAA*gAA + nAB*gAB + nBB*gBB)/(RT)]`.

Here `gAA`, `gAB`, and `gBB` are contact energies in **kJ/mol**, and the
temperature is in **kelvin**. The common energy offset cancels when the
distribution is normalised. Positive `gAB` relative to `gAA=gBB=0` favours
self-assembly; negative `gAB` favours alternating AB arrangements.

`contact_classes(N)` counts rings by A composition and AA/AB/BB contacts.
For a mixed ring with `i` A subunits and `r` runs of each type, its exact
labelled multiplicity is

`N/r × C(i−1,r−1) × C(N−i−1,r−1)`.

This groups arrangements with identical energies without enumerating all
`2**N` strings. It retains rotational multiplicities: an alternating even
ring has multiplicity 2, while either pure ring has multiplicity 1.
Reflections contribute separately when they correspond to distinct labelled
arrangements. At zero energy bias, the composition distribution is binomial.

## Predict a composition distribution

From the repository root:

```python
from combined.statistical_mechanics import composition_distribution

# p[i] is the fraction of 12-mers containing i A subunits, i=0,...,12.
p = composition_distribution(
    12, g_ab=1.98, g_aa=0.0, g_bb=0.0,
    temperature=298.15, mole_fraction_a=0.5,
)
```

`mole_fraction_a` specifies the **mean subunit fraction** of A. The function
solves the relative activity `a/b` to achieve it. For a specified activity
ratio instead, use `activity_distribution(..., log_activity_ratio=log(a/b))`.

## Fit several measured mixing ratios together

```python
import numpy as np
from combined.statistical_mechanics import fit_hetero_energy

# The CSV has three rows (1:2, 1:1, 2:1 mixtures), each with 13
# measured abundances ordered by A count from 0 through 12.
observed = np.loadtxt("my_12mer_abundances.csv", delimiter=",")

fit = fit_hetero_energy(
    observed,
    initial_mole_fractions_a=[1/3, 1/2, 2/3],
)
print(fit.hetero_energy_kj_mol, fit.fitted_mole_fractions_a)
```

The fit shares one AB contact energy across datasets and fits an activity
ratio for each mixture, as described in the thesis. It minimises unweighted
squared differences in normalised composition fractions. Provide uncertainty
weights or a measurement model before interpreting fit confidence intervals.
The thesis's reported 1.98 kJ/mol cannot be reproduced without its original
experimental abundance table, which was not in the attached folder.

## Verify and plot

```bash
python -m unittest combined.statistical_mechanics.test_ring_coassembly
python -m combined.statistical_mechanics.example
python -m combined.statistical_mechanics.replicate_counting_model
```

The example saves `combined/statistical_mechanics/example_ring.png`. The replication script
recreates the three panels of the thesis's `counting_model.eps` as
`combined/statistical_mechanics/counting_model_replication.png`: the unbiased 12-mer binomial
distribution, the pure-homomer and 6:6 fractions across −5 to +5 kT, and the
full composition heatmap. Its plotted energy difference is
`(g_homo − g_AB)/RT`, the opposite sign of the API's `g_ab` when the homotypic
energies are zero. Panel a is divided by its mode, as in the source figure;
panels b and c use fractions summing to one. The original numerical plotting
table was not supplied, so this is a reconstruction from the model and caption.

The tests compare contact
counts to explicit enumeration for rings of size 3–10, check the unbiased
binomial limit, test both assembly extremes, and recover a known synthetic
energy from three mixtures.

## Scope

The ring has a fixed size and two subunit types. This predicts equilibrium
**composition fractions**, not size distributions or time courses. It does
not account for native MS response differences, charge-state overlap, or
mass-deconvolution uncertainty. Those belong in a measurement/fitting layer
for real spectra.
