# polydispersity
Repo for dealing with polydisperse proteins and their assembly.

Initially, I will focus on homo-oligomers of small heat-shock proteins, using multiple different approaches: chemical kinetics, thermodynamics, and statistical mechanics.

The `single/statistical_mechanics/` module now includes an independent
equilibrium scaffold model. It counts connected ring and polyhedral contact
configurations, assigns Boltzmann weights, and compares the resulting size
distributions with the existing kinetic model. See
[its README](single/statistical_mechanics/readme.md) for the assumptions,
equations, and figures.

The repository now also includes a statistical-mechanics model for fixed-size,
circular co-assembly of two protein species. A chemical-kinetics model for
polydisperse co-assembly is a possible future step.

The `combined/statistical_mechanics/` module reconstructs the fixed-size ring
co-assembly model described in the thesis. See
[its README](combined/statistical_mechanics/readme.md) for the equations, usage,
and verification.

Then, I want to focus on the description of general complex protein assembly that is at the global energy minimum and limited in maximum size.

This is a hobby project I only have limited time to invest time in every now-and-then. If you want to work on this together, please reach out.

Dom
