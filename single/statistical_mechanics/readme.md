# Statistical mechanics of a single polydisperse protein

This compares equilibrium statistical weights for **variable-size
homo-oligomers** with the existing [kinetic calculation](../oligomer_distribution.py).
It follows the homo-oligomer section of Chapter 2 in the attached thesis and
the specific-site thermodynamic interpretation in
[Baldwin et al., *J. Mol. Biol.* 413 (2011)](https://baldwinlab.chem.ox.ac.uk/publications/2011%20alphaB%201.pdf),
especially its Eqs. 5–7. The code here does not fit new experimental data.

## What carries over from the circular model

For a fixed-size two-protein ring, arrangements have the same number of
subunits, so contact energies and arrangement multiplicities determine the
**composition** distribution. For a single polydisperse protein, the number
of subunits changes. We must also account for monomer activity and the
size-dependent number of ways an oligomer can lose a subunit.

Write `alpha = k_on_eff / k_d` and `r = k_m / k_d`, where `k_on_eff = k+ [P1]`,
`k_d` is the off-rate from a paired monomer, and `k_m` is that from the lone
unpaired monomer. The thesis's microscopic free energies are

```
G_edge  = -RT log(k_on_eff / k_m)
G_dimer = -RT log(k_m / k_d)
```

The **whole-oligomer** off-rate includes all exit pathways. It is `i*k_d`
for even size `i`, and `(i-1)*k_d + k_m` for odd size `i`. These counts give
an equilibrium partition weight, relative to the monomer,

```
w_1 = 1
w_i/w_(i-1) = alpha / i           (even i)
w_i/w_(i-1) = alpha / (i-1 + r)   (odd i)
P_i = w_i / sum_j w_j
```

`site_count_distribution` evaluates the product using gamma functions,
independently of the kinetic code's numerical recursion. It agrees with that
recursion to floating-point precision. When `r=1`, it reduces to a Poisson
weight proportional to `alpha^(i-1)/i!` for positive oligomer sizes.
**This exact agreement is a reformulation of the same detailed-balance
assumptions, not independent confirmation of the physical model.**

## What simple contact energies predict

An independent additive-contact toy model assigns an edge energy for every
addition and a dimer bonus for every pair:

```
w_i ∝ exp(-[(i-1)*G_edge + floor(i/2)*G_dimer]/RT)
```

With the rate-derived energies above, this has no finite interior peak: its
distribution piles up at the arbitrary `max_size` cutoff. Adding a
phenomenological `1/i!` ideal-cluster factor produces a peak, but still does
not reproduce the kinetic shape with the **same** microscopic energies. For
the repository's `(k_on_eff, k_m, k_d)=(24,12,1)` example, the kinetic and
site-count curves peak at an **18-mer**, whereas the contact-only model peaks
at the 60-mer cutoff and the factorial version peaks at a **6-mer**.

If both energies of the factorial model are fitted freely to the kinetic
curve, it gives a similar distribution with total-variation distance `0.049`
and the same 18-mer mode. The fitted dimer energy is about `-1.15 kJ/mol`,
versus `-6.16 kJ/mol` obtained from the microscopic rate ratio at 298.15 K.
These fitted toy-model energies should **not** be read as the microscopic
interface energies.

## Ring versus polyhedral geometry

The [structural study of αB-crystallin](https://baldwinlab.chem.ox.ac.uk/publications/2011%20alphaB%203.pdf)
places dimeric cores along scaffold edges and C-terminal interactions at
vertices. For a 24-mer it considers a single ring of 12 dimers, a cube, and
an octahedron, among other architectures. The study's ion-mobility comparison
favours degree-four polyhedra over single rings or degree-three polyhedra;
its proposed 24-, 26-, and 28-mers are an octahedron, augmented triangular
prism, and gyrobifastigium. Their odd neighbours can be formed by removing
one monomer. These are structural candidates, not measured Boltzmann weights.

`geometry.py` represents each scaffold edge as one dimer and links incident
monomers at each vertex in a directed C-terminal cycle. Thus each monomer
donates and receives one C-terminal contact. This is a **contact graph**;
the cyclic ordering is arbitrary and does not model atomic coordinates.

| 24-mer scaffold | Vertices | Vertex degree | Dimer contacts | Directed C-terminal contacts |
| --- | ---: | ---: | ---: | ---: |
| Single 12-dimer ring | 12 | 2 | 12 | 24 |
| Cube | 8 | 3 | 12 | 24 |
| Octahedron | 6 | 4 | 12 | 24 |

Therefore `12*g_dimer + 24*g_C-terminal` is **identical** for all three
24-mer candidates. A model with only two constant contact energies cannot
choose among their shapes. In particular, the `G_edge` inferred from an
effective monomer-addition rate must not automatically be equated with the
energy of one C-terminal contact; the reference states and contact counts
differ.

For a closed ring of `m` dimers, the simple grand-canonical model gives

```
ΔG = g_dimer + 2*g_C-terminal - 2*mu_monomer
w_m ∝ exp(-m*ΔG/RT)                  (m >= 3)
w_(m+1)/w_m = exp(-ΔG/RT)
```

`mu_monomer` is the monomer chemical potential. The constant ratio means
that the distribution is increasing, decreasing, or flat over the allowed
ring sizes; it has no selected interior mode. Closed rings also exclude odd
oligomers. `compare_geometries.py` plots these three regimes against the
kinetic solution. A size-dependent **shape entropy, closure/curvature cost,
or scaffold degeneracy**, together with explicitly modeled odd structures,
is required to predict a broad finite-size distribution from geometry. The
paper does not provide those free energies, so the current geometry code
does not fit or claim to reproduce the experimental distribution. The exact
site-count curve above remains the meaningful numerical match to the
kinetic model, because it carries over its size-dependent pathway counts.

## Run the comparison

From the repository root:

```bash
python -m unittest single.statistical_mechanics.test_size_distribution
python -m single.statistical_mechanics.compare_models
python -m single.statistical_mechanics.compare_geometries
```

The plotting commands save `size_distribution_comparison.png` and
`geometry_comparison.png` in this folder.

```python
from single.statistical_mechanics import energies_from_rates, site_count_distribution

energies = energies_from_rates(24, 12, 1, temperature=298.15)
probabilities = site_count_distribution(
    energies.edge_kj_mol, energies.dimer_kj_mol,
    temperature=298.15, max_size=60,
)
# probabilities[0] is the monomer fraction, probabilities[17] the 18-mer.
```

All curves here are fractions of **oligomer molecules**, normalized over the
chosen size range. They are not monomer mass fractions or calibrated native
MS intensities. The effective on-rate already contains free-monomer
concentration; predicting a concentration series requires solving monomer
mass balance separately.
