# Statistical mechanics of a single polydisperse protein

This compares equilibrium statistical weights for **variable-size
homo-oligomers** with the existing [kinetic calculation](../oligomer_distribution.py).
It follows the homo-oligomer section of Chapter 2 in the attached thesis and
the specific-site thermodynamic interpretation in
[Baldwin et al., *J. Mol. Biol.* 413 (2011)](https://baldwinlab.chem.ox.ac.uk/publications/2011%20alphaB%201.pdf),
especially its Eqs. 5–7. The code here does not fit new experimental data.

## A small-size, monomer-at-vertex geometry hypothesis

The first geometry ensemble below followed the structural paper's
**dimer-on-edge** convention. That leaves out elementary open chains and
makes an octahedron a 24-mer. `monomer_geometry_ensemble.py` adds a separate
exploratory convention in which **one graph vertex is one monomer**. It does
not mix its partition weights with the dimer-on-edge calculations.

The default candidate library has an open chain at every size 1–12, a closed
ring from three monomers onward, and selected compact polyhedra from six
onward. The six-monomer candidates are an open chain, closed ring, triangular
prism, octahedron, and pentagonal pyramid. At eight monomers, the catalogue
includes a cube and square antiprism; at twelve it includes an icosahedron.
These are graph hypotheses rather than experimentally established αB-crystallin
structures. A tetrahedron is mathematically possible with four monomers and
a square pyramid or triangular dipyramid with five. Set
`include_small_polyhedra=True` to include those; omitting them by default is
the proposed **physical onset assumption**, not a geometric theorem.

Graph edges here are **possible spatial adjacencies**. The model assumes
`floor(n/2)` dimer contacts and a C-terminal path with `n-1` contacts for an
open chain, or a C-terminal cycle with `n` contacts for a closed graph. The
tests verify that each listed graph permits the required dimer matching and
path/cycle counts separately. They do not establish a simultaneous atomic
conformation or a physically correct C-terminal wiring. Extra polyhedral
edges are not counted as extra chemical bonds.

The illustrative conditional weight is

```
q(g,n) = A_g exp[epsilon_d*floor(n/2) + epsilon_C*C_g
                 - kappa*Σ_v max(0, degree_g(v)-2)² - Δf_g]
P(g | n) = q(g,n) / Σ_h q(h,n)
```

The `kappa` term penalizes high *potential* coordination as a shape proxy.
This is an unmeasured assumption; the sign and magnitude of the actual
compact-shape free energy could differ. In particular, it is a **different
shape hypothesis** from the degree-four preference used in the dimer-on-edge
example below; their numerical fractions should not be compared. `A_g` and
`Δf_g` allow explicit conformer/symmetry factors and geometry offsets,
defaulting to one and zero.
At a fixed size, the dimer energy cancels. All closed graphs also have the
same C-terminal contact count, so distinguishing two closed shapes needs
shape or entropy information. The [monomer-level CSV](monomer_geometry_contributions.csv)
contains every candidate at sizes 1–12 under three shape penalties, and the
[figure](monomer_geometry_contributions.png) shows family shares and six-mer
sensitivity. At `epsilon_d=epsilon_C=1` and `kappa=0.1 kBT`, the six-mer
shares are 16.3% open chain, 44.4% ring, 24.3% triangular prism, 4.0%
octahedron, and 10.9% pentagonal pyramid. **These numbers illustrate the
assumptions; they are not inferred solution abundances.** Run
`python -m single.statistical_mechanics.compare_monomer_geometries` to
regenerate the outputs.

## Relative geometry contributions in the dimer-on-edge model

`geometry_ensemble.py` compares **alternative closed scaffold graphs at each
oligomer size**. The catalogue includes single rings; tetrahedral, prismatic,
antiprismatic, cubic and octahedral graphs; illustrative pyramids and
dipyramids; and two vertex-split octahedral graphs for 26 and 28 monomers.
An odd oligomer is a one-monomer defect of the next even closed scaffold.
The complete candidate list and an evidence tag for each graph are in
`geometry_contributions.csv`. This is a **finite hypothesis catalogue**, not
an enumeration of every possible αB-crystallin architecture. The published
multiple-ring candidates are not represented because their connecting
contacts and relative conformational weights still need a graph definition.

For geometry `g` at size `n`, the conditional partition weight is

```
q(g,n) = A_g Σ_(defect sites s) exp[epsilon_d*D_s + epsilon_C*C_s
                                  - kappa*Σ_v (d_v(s)-d_*)² - Δf_g]
P(g | n) = q(g,n) / Σ_h q(h,n)
```

Here `D` counts intact dimers, `C` counts directed C-terminal contacts,
`d_v` is the occupied degree at vertex `v`, and all energy terms are in
`k_B T`. `A_g` is an optional architecture/conformation multiplicity and
`Δf_g` an optional shape offset. The default takes `A_g=1`, `Δf_g=0` and
`d_*=4`. The penalty `kappa` is an **unmeasured illustrative shape free
energy**. For even closed structures the sum has one graph state; for odd
sizes it sums all labelled sites at which one monomer could be absent, with
the remaining contacts at that vertex allowed to rearrange. This site
counting is a modelling convention; actual conformer and symmetry weights
are unknown. A ring's degree-two defect loses two C-terminal contacts,
whereas a degree-three or degree-four vertex defect loses one.

At every even size `n`, all closed graphs here have exactly `n/2` intact
dimers and `n` C-terminal contacts. **Uniform interaction energies therefore
cancel from `P(g | n)`**. The relative shares are controlled by shape free
energies and multiplicities, which have not been measured. Monomer activity
also cancels at fixed `n`. This calculation does not predict the overall
size distribution `P(n)`; that would need a common standard state, monomer
activity, and comparable architecture partition functions across sizes.

For example, the 24-mer catalogue has a 12-dimer ring, cube and octahedron.
With `epsilon_d=epsilon_C=1` and all offsets zero, the shares are:

| Degree mismatch penalty `kappa` | Ring | Cube | Octahedron |
| ---: | ---: | ---: | ---: |
| 0 `k_B T` | 33.3% | 33.3% | 33.3% |
| 0.1 `k_B T` | 0.6% | 30.8% | 68.6% |
| 0.3 `k_B T` | <0.1% | 8.3% | 91.7% |

The [structural study](https://baldwinlab.chem.ox.ac.uk/publications/2011%20alphaB%203.pdf)
and its [candidate catalogue](https://baldwinlab.chem.ox.ac.uk/publications/2011%20alphaB%203%20supp.pdf)
motivate several of these graphs and a degree-four preference in its
ion-mobility comparison. That comparison is **not** a measurement of
`kappa` or solution-state Boltzmann fractions. In particular, 25–28-mer
shares in this small catalogue look very decisive because each size has
only a single ring and an illustrative vertex-split graph; extra plausible
architectures or a different shape offset could change them sharply.

Run `python -m single.statistical_mechanics.compare_geometry_ensemble` at
the repository root to regenerate `geometry_contributions.csv` (every size
5–40 at three assumed penalties) and `geometry_contributions.png`. Call
`conditional_contributions(size, EnergyParameters(...), offsets_kbt=...,
log_multiplicities=...)` for other interaction and shape-energy hypotheses.

## Restricted equilibrium scaffold model

`equilibrium_scaffold.py` now provides a **rate-free** statistical-mechanical
model **within one chosen parent scaffold**. A ring, cube, or octahedron
specifies possible dimeric edges and C-terminal contacts. We enumerate
connected occupied configurations with intact dimer units and at most one
unpaired monomer. Configurations related by that parent's symmetry are counted
once. For each size `n`, the restricted partition sum is

```
Z_n^(parent) = Σ_{allowed occupied states s of size n in this parent}
      exp[n*log(z) + D_s*epsilon_d + C_s*epsilon_C - V_s*kappa]
P_n^(parent) = Z_n^(parent) / Σ_m Z_m^(parent)
F_n^(parent) = -k_B*T*log(Z_n^(parent))
```

Here `z` is monomer activity, `D_s` counts intact dimer contacts, `C_s`
counts directed C-terminal contacts, and `V_s` counts fully occupied
vertices of degree at least three. The favourable contact strengths
`epsilon_d` and `epsilon_C`, and the illustrative crowded-vertex penalty
`kappa`, are dimensionless energies in units of `k_B T`. **No rate or
kinetic recurrence enters these weights.** The size dependence instead
comes from the number of allowed occupied states in the chosen parent and
their contact energies. In particular, the octahedral template has many more
partially filled patterns than nearly full ones, whereas a **single**
12-dimer parent ring has one symmetry-inequivalent pattern per size under
these restricted rules. That is not a count of all possible ring oligomers.

`compare_equilibrium.py` plots those restricted counts and a sample equilibrium
distribution against the kinetic result, conditioned on sizes `1..24` so the
supports match. The example parameters (`log(z)=-1`, `epsilon_d=4`,
`epsilon_C=0.05`, `kappa=1.1`) place the cube and octahedron modes at 18
monomers, while the ring mode is 24. They are illustrative values, **not
measured interface free energies or a fit to experiment**. This shows that
polyhedral configurational multiplicity *can* select an interior size in an
equilibrium model; it does not establish that this is the mechanism in
αB-crystallin or quantitatively reproduce its abundance curve.

The structural assumptions matter: the scaffold is a template of possible
contacts, not a pre-existing empty protein shell. At a vertex, occupied
monomers are assumed to be able to rearrange into a cyclic C-terminal
contact pattern; the model does not enumerate atomic conformations or their
entropy. The crowded-vertex penalty is a tunable shape term with no measured
value. Only one 24-monomer parent scaffold is considered at a time, so this
example cannot predict oligomers above 24. The conditional geometry library
above adds other closed sizes, but predicting a full size distribution still
needs their relative shape free energies and concentration-dependent data to
constrain `z`.

### Other arrangements at the same size

The [structural study](https://baldwinlab.chem.ox.ac.uk/publications/2011%20alphaB%203.pdf)
explicitly sketches **seven** candidate 24-mer scaffold families: a single
12-dimer ring; double 6-dimer rings; triple 4-dimer rings; quadruple
3-dimer rings; a cube; an octahedron; and a mixed-degree elongated triangular
pyramid. It proposes an augmented triangular prism for 26 monomers and a
gyrobifastigium for 28. The restricted occupied-pattern calculation below
uses only the single-ring, cube, and octahedron parents; the separate
conditional ensemble above adds other graph candidates. The paper's
ion-mobility comparison favours
degree-four polyhedra, but that does not establish zero solution population
for the other candidates.

Even within one scaffold, there may be alternative dimer pairings, C-terminal
wirings at a vertex, unpaired-monomer positions, and conformations. The code
fixes the first two, identifies symmetry-equivalent positions, and omits
conformational entropy. Its pattern counts therefore **must not be treated as
the complete physical degeneracy**. Nor can separately normalized curves for
the three parents be added: many small occupied fragments could represent
the same molecule and would be counted more than once.

The fuller equilibrium object would be

```
Z_n = z^n Σ_{architecture a at size n} q_(n,a)
q_(n,a) = Σ_{distinct conformations c of a} exp[-F_(n,a,c)/(k_B*T)]
P_n = Z_n / Σ_m Z_m
```

`q_(n,a)` must include the appropriate symmetry and conformational weights,
and architectures must be deduplicated across parent templates. The geometry
ensemble above is a small conditional step toward this object, with
illustrative shape energies; it is not a complete full-size partition sum.
IM-MS collision cross sections and other structural restraints can narrow
which `a` are plausible before their equilibrium weights are fitted.

### Calibration against the kinetic example

`fit_scaffold_to_kinetics.py` fits **four shared equilibrium parameters**
(`log(z)`, dimer stabilization, C-terminal stabilization, and saturated-vertex
penalty) to the repository's synthetic kinetic example
`(k_on_eff, k_m, k_d)=(20,10,1)`: effective on-rate 20, unpaired off-rate
10, and paired off-rate 1. Binding stabilizations and the penalty are
constrained nonnegative. The kinetic curve is a calibration target, not an
input to the partition function itself. Both curves are normalized over
sizes `1..24` for the fit, then the fitted equilibrium curve is compared with
the complete kinetic curve out to size 60.

| Parent template | Total variation, sizes ≤24 | Total variation, full kinetic curve | Fitted mode |
| --- | ---: | ---: | ---: |
| Single 12-dimer ring | 0.317 | 0.326 | 24 |
| Cube | 0.119 | 0.138 | 18 |
| Octahedron | 0.110 | 0.132 | 16 |

With these rates, the kinetic curve peaks at a 16-mer and has `0.03266`
probability above size 24. The octahedral equilibrium fit has the same mode;
its remaining mismatch is partly within the shared size range. A
concentration-proxy check fits at effective on-rate 20, then changes only
`log(z)` by `log(new_on_rate/20)`. For the octahedron, total variation on
sizes `1..24` is `0.138` at on-rate 16 and `0.112` at on-rate 18. Thus the
present contact catalogue captures a broadly similar equilibrium curve but
does **not** quantitatively replicate the kinetic family of curves.
The earlier factorial toy model fits the single baseline curve more closely
(total variation `0.049`), but that is a descriptive fit with different
energies and is not a validated geometry-based explanation.

A faithful equilibrium reproduction would require more candidate scaffold
sizes, properly counted architecture/conformer weights, and shared parameters
that succeed across concentrations. Assigning a separate free energy to each
size could force an exact match, but would only encode the kinetic curve in
new notation and is not the intended test.

### Why the curves are not exactly symmetric; tetrahedral controls

Neither equilibrium nor the kinetic model requires a symmetric size
distribution. For the requested `(20,10,1)` kinetic example, the mean size is
`16.01`, the mode is `16`, and the standardized skewness is `+0.24`: the
envelope is **approximately** centered but has a longer right tail. The
different off-rate denominators at even and odd sizes create alternating
heights. In the scaffold model, the count of *connected* occupied patterns,
the number of intact dimers, and the monomer activity also change unevenly
with size. Binomial symmetry would require independent equivalent sites and
special parameter values; those assumptions do not hold here.

The [structural paper's supplement](https://baldwinlab.chem.ox.ac.uk/publications/2011%20alphaB%203%20supp.pdf)
includes a **tetrahedron with six dimeric edges**, which is a 12-mer under
our one-dimer-per-edge rule. At `(20,10,1)`, `77.8%` of the kinetic
distribution lies above 12, so that tetrahedron alone cannot reproduce it.
As a controlled 24-monomer comparison, `geometry.py` also constructs an
**edge-subdivided tetrahedron**, with two dimeric positions along each
original edge. This latter graph is hypothetical and is not one of the
paper's identified αB-crystallin architectures.

| Parent | Maximum size | Full-curve total variation after fitting | Fitted mode |
| --- | ---: | ---: | ---: |
| Tetrahedron | 12 | 0.778 | 12 |
| Edge-subdivided tetrahedron | 24 | 0.153 | 16 |
| Octahedron | 24 | 0.132 | 16 |

The octahedral fitted curve has skewness `+0.17`, closer to the kinetic
`+0.24` than the subdivided tetrahedron's `−0.65`. Both 24-mer parents
match the modal size, but neither quantitatively reproduces the distribution.
Run `compare_tetrahedra.py` for the distributions and restricted state counts.

The rest of this document retains the exact kinetic reformulation and the
earlier additive-contact examples as checks and contrasts. Their matching
curves are not an independent equilibrium derivation.

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
python -m unittest single.statistical_mechanics.test_equilibrium_scaffold
python -m single.statistical_mechanics.compare_models
python -m single.statistical_mechanics.compare_geometries
python -m single.statistical_mechanics.compare_equilibrium
python -m single.statistical_mechanics.fit_scaffold_to_kinetics
python -m single.statistical_mechanics.compare_tetrahedra
```

The plotting commands save `size_distribution_comparison.png`,
`geometry_comparison.png`, `equilibrium_scaffold_comparison.png`, and
`scaffold_kinetic_fit.png`, and `tetrahedral_comparison.png` here.

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
