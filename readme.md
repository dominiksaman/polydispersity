# Polydispersity

Equilibrium models for polydisperse protein assembly and circular A/B co-assembly.

| Area | What it calculates |
| --- | --- |
| [Single protein](single/readme.md) | The size distribution implied by the modified-helical rate model; a neutral-mass axis conversion |
| [Single-protein statistical mechanics](single/statistical_mechanics/readme.md) | Exact rate reformulation and explicit geometry weights |
| [Single-protein thermodynamics](single/thermodynamics/readme.md) | Formation free energies, monomer/dimer reservoir, and mass balance |
| [Combined proteins](combined/statistical_mechanics/readme.md) | Composition distributions and energy fitting for a fixed-size A/B ring |

**Current reference:** effective on-rate 20, unpaired off-rate 10, paired off-rate 1.
The single-protein distribution peaks at a 16-mer; 3.27% lies above a 24-mer.
These are oligomer number fractions. Geometry populations require additional
shape free energies and statistical weights; they are not inferred from a drawing.

## Run

Python 3.9 or later, from the repository root:

```bash
python -m pip install -r requirements.txt
python -m unittest discover -v
python -m single.example
python -m combined.statistical_mechanics.replicate_counting_model
python -m single.statistical_mechanics.draw_eleven_geometries
python -m single.thermodynamics.simulate_reservoir
```

## Main figures

- [Size and neutral-mass distributions](single/size_distribution.png)
- [Thesis ring co-assembly figure reconstruction](combined/statistical_mechanics/counting_model_replication.png)
- [Two possible 11-monomer contact networks](single/statistical_mechanics/eleven_monomer_geometries.png)
- [Thermodynamic small-species pool and larger oligomers](single/thermodynamics/readme.md)

See the [scientific and software audit](docs/audit.md) for validated results,
corrections, and remaining assumptions. Earlier exploratory fits and sensitivity
plots remain available in Git history; the current tree contains the supported
calculations and the geometry work needed to continue.
