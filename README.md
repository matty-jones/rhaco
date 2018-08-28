## Installing Rhaco

1) Make sure conda is installed and active on your system (miniconda recommended from https://conda.io/miniconda.html)

2) Make the rhaco conda environment using the specified prerequisites: `conda env create -f environment.yml`

3) Activate the new conda environment with: `source activate rhaco`

4) From the rhaco root directory install rhaco using `pip install -e .`

## Included Command-Line Programs

1) rhaco-create-morph: Used to create the initial conditions for the simulation based on several input parameters (`rhaco-create-morph -h` for more details)

2) rhaco-run-hoomd: Used to interpret the Foyer forcefields and begin a HOOMD molecular dynamics simulation based on several input parameters (`rhaco-run-hoomd -h` for more details)

## Version History

v1.0 - Release of Rhaco, set up for small hydrocarbon interactions on an M1 catalyst surface
v1.1 - Rhaco generalized to work with other systems, specifically the interactions of polydimethylsiloxane on Ni-Mn-Ga shape memory alloy
v1.2 - Additional functionality added to permit the simulation of metallic nanoparticles on a surface.


## Important Citations

Link to paper with M1 coordinates
http://pubs.acs.org/doi/suppl/10.1021/jacs.5b07073/suppl_file/ja5b07073_si_001.pdf

Link to paper with Ethylene coordinates
http://www.emsl.pnl.gov/docs/tms/abinitio/tables/appendixa.pdf

Rhaco has been released under a GPL3 license (please see LICENSE.TXT for conditions).
