# Scientific and software audit — 29 September 2026

## Subsequent thermodynamic extension

The [reservoir model](../single/thermodynamics/readme.md) adds
an independent ideal-mixture thermodynamic calculation with supplied formation
free energies, monomer mass balance, a small-species pool, and exact open-reservoir
count sampling. Its bulk/surface/packing family is inspired by the Knowles group's
2018 globular model; its coefficients are illustrative. It does not establish
geometry-dependent free energies or replicate the paper's spatial Monte Carlo.
Seven additional conservation, standard-state, convergence, and sampling checks
bring the suite to 28 tests. A further pool-enrichment check brings the current
suite to 29 tests: the new chosen 20% monomer/dimer subunit pool preserves the
conditional larger-species shape and mass balance at the reference total.
The original audit below remains a record
of the preceding cleanup.

## Conclusions

The single-protein equilibrium recurrence agrees with the thesis equations and
the modified-helical model's published equilibrium relations. The A/B ring
multiplicities and its thesis figure reconstruction also check out. The earlier
geometry percentages were controlled by unmeasured shape penalties and assumed
state factors. They do not establish solution populations or an independent
physical reproduction of the kinetic size distribution.

## Corrections and checks

| Item | Result |
| --- | --- |
| Requested rates: effective on 20, unpaired off 10, paired off 1 | Mode 16, mean 16.012607, number fraction above 24 = 0.032656 |
| Rate-model recurrence versus closed-form weights | Agreement at floating-point precision for the reference; tested across distinct rate ratios |
| Large-rate numerical overflow | Replaced raw recursion with log weights; a 1000-rate Poisson case remains normalized |
| Package import of mass conversion | Repaired relative import; conversion preserves oligomer number fractions |
| Ring composition counting | Exact small-ring enumeration and binomial limits pass |
| Ring energy fitting | Synthetic recovery passes; only heterotypic preference is identifiable with a free activity ratio |
| Monomer geometry compatibility | Dimer and C-terminal contacts are now assigned jointly, with connected occupied networks |
| Odd rings | Full C-terminal saturation is incompatible with maximal nearest-neighbour dimers under the chosen separate-interface rule; partial binding is still possible |
| Both compact 11-mer drawings | Joint five-dimer/eleven-C-terminal assignments exist; atomic feasibility remains untested |
| Frozen monomer defects | Two incident C-terminal contacts are lost; rewiring requires a separate specified state |
| Test suite | 21 tests pass, covering numerical, combinatorial, and contact-network constraints |

All three retained figures were regenerated and visually inspected. The neutral
mass panel is an axis conversion, not a simulated mass spectrum. Normalization
is over the computed range; a 24-mer numerical cutoff would hide the reference
tail rather than remove it physically.

## Why the presentation was simplified

Removed from the current tree: duplicate examples and dependency lists, arbitrary
factorial fits, parent-scaffold occupancy fits, degree-penalty geometry percentages,
and their comparison plots/CSVs. Counting one representative symmetry orbit of
an occupied parent does not supply its complete physical statistical weight.
Separate parent templates can also represent the same small oligomer, causing
double counting if their weights are added.

These experiments remain recoverable in the
[last exploratory commit](https://github.com/dominiksaman/polydispersity/tree/11b307c).
Their numerical fits are historical exploratory results, not validated structural
inference. The retained API requires explicit residual shape free energies and
statistical factors, and sums only supplied distinct states at a common size.
Equal-contact closed geometries remain indistinguishable by uniform contact
energies alone.

## Thesis comparison and source limits

The source checked was `/Users/dominiksaman/Downloads/Ver3-2/text/ch1-intro.tex`.
Its free-monomer time equation has a positive extra dimer-association term where
the [original kinetic paper, Eq. 2](https://baldwinlab.chem.ox.ac.uk/publications/2011%20alphaB%201.pdf)
has a negative term. This thesis typo does not enter the implemented equilibrium
recurrence. The thesis source was not edited. The original co-assembly abundance
table is absent, so the reported 1.98 kJ/mol cannot be experimentally refitted.

The [structural paper](https://baldwinlab.chem.ox.ac.uk/publications/2011%20alphaB%203.pdf)
and [supplement](https://baldwinlab.chem.ox.ac.uk/publications/2011%20alphaB%203%20supp.pdf)
motivate dimer-on-edge candidates and an ion-mobility comparison. They do not
provide Boltzmann weights for the alternative solution architectures. The
monomer-at-vertex polyhedra are additional graph hypotheses.

## Productive next steps

1. Specify and deduplicate contact/conformer states, keeping actual interfaces
   distinct from possible spatial adjacencies. Check steric feasibility.
2. Constrain architecture candidates using oligomer-resolved native IM-MS and
   orthogonal solution measurements; schematic drawing coordinates are not CCS models.
3. Supply or infer shared interface and residual shape free energies across
   concentrations, with monomer mass balance and an appropriate MS response model.
4. Test predictions on held-out concentrations and structural measurements.
   Matching one synthetic abundance curve or its mode is insufficient evidence.

A full geometry-derived equilibrium size distribution remains open. The exact
rate-derived equilibrium curve is retained with its physical assumptions explicit.
