## Installing Lynx

1) Make sure conda is installed and active on your system (miniconda recommended from https://conda.io/miniconda.html)

2) Make the lynx conda environment using the specified prerequisites: `conda env create -f environment.yml`

3) Activate the new conda environment with: `source activate lynx`

4) From the lynx root directory install lynx using `pip install -e .`

## Included Command-Line Programs

1) lynx-create-morph: Used to create the initial conditions for the simulation based on several input parameters (`lynx-create-morph -h` for more details)

2) lynx-run-hoomd: Used to interpret the Foyer forcefields and begin a HOOMD molecular dynamics simulation based on several input parameters (`lynx-run-hoomd -h` for more details)


## Important Citations

Link to paper with M1 coordinates
http://pubs.acs.org/doi/suppl/10.1021/jacs.5b07073/suppl_file/ja5b07073_si_001.pdf

Link to paper with Ethylene coordinates
http://www.emsl.pnl.gov/docs/tms/abinitio/tables/appendixa.pdf
