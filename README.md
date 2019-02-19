## Installing Rhaco

1) Make sure conda is installed and active on your system (miniconda recommended from https://conda.io/miniconda.html)

2) Make the rhaco conda environment using the specified prerequisites: `conda env create -f environment.yml`

3) Activate the new conda environment with: `source activate rhaco`

4) From the rhaco root directory install rhaco using `pip install -e .`

## Included Command-Line Programs

1) rhaco-create-morph: Used to create the initial conditions for the simulation based on several input parameters (`rhaco-create-morph -h` for more details)

2) rhaco-run-hoomd: Used to interpret the Foyer forcefields and begin a HOOMD molecular dynamics simulation based on several input parameters (`rhaco-run-hoomd -h` for more details)

## Version History

v1.4.1 - Merged in the ad-hoc "create_forcefield" script to the rhaco simulate pipeline so that the required EAM forcefield can be produced at runtime
v1.4.0 - Combinations of positional, rigid, and flexible reactants can now be specified in the system. The --reactant_composition and --reactant_rigid options have been improved to reflect the new functionality.
v1.3.0 - Added functionality for rigid bodies in HOOMD 2 (including orientation quaternions and moments of inertia tensor eigenvalues)
v1.2.0 - Additional functionality added to permit the simulation of metallic nanoparticles on a surface.
v1.1.0 - Rhaco generalized to work with other systems, specifically the interactions of polydimethylsiloxane on Ni-Mn-Ga shape memory alloy
v1.0.0 - Release of Rhaco, set up for small hydrocarbon interactions on an M1 catalyst surface

## Important Citations

Link to paper with M1 coordinates
http://pubs.acs.org/doi/suppl/10.1021/jacs.5b07073/suppl_file/ja5b07073_si_001.pdf

Link to paper with Ethylene coordinates
http://www.emsl.pnl.gov/docs/tms/abinitio/tables/appendixa.pdf

Rhaco has been released under a GPL3 license (please see LICENSE.TXT for conditions).
