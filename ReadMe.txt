This repository contains all of the code necessary to recreate the simulations and analysis presented in "Understanding the drivers of dispersal evolution in range expansions and their ecological consequences". All model code is written in R and some scripts have accompanying bash scripts to run them in parallel on the University of Wyoming's high performance computing cluster. 

File descriptions:
---- SimFunctions.R
This script contains all of the necessary functions to run a single simulation of our dispersal evolution model through all stages of each simulation (initial burn-in, unbound range expansion, and a directional range shift). Each function is commented with information on the expected inputs and outputs.

---- CreateAllSimulationFile.R and CreateSesitivitySimulations.R
These scripts are used to create .csv files detailing the different parameter combinations for each simulation to be run for the main analyses and the sensitivity analysis.

---- RunSimulations.R, RunSimulations.sh, RunSensitivitySims.R, RunSensitivitySims.sh
These scripts run all the simulations for the main analysis and the sensitivity analysis in parallel on the University of Wyoming's high performance computing cluster. The R scripts use the .csv files created by the above files to ensure that each simulation is run with appropriate parameters and the bash scripts request sufficient nodes and wall time to run all the simulations.

---- SimulationCheck.R, SimulationCheck.sh, SensSimulationCheck.R, SensSimulationCheck.sh
These scripts check for the existence of each simulation detailed in the .csv files so that if a job is interrupted, any completed simulations may be retained and only incomplete simulations need to be run again.

---- DataExtraction.R, DataExtraction.sh, InitDataExtraction.R, InitDataExtraction.sh, SensDataExtraction.R, SensDataExtraction.sh
These scripts extract data from the individual simulation files for the analysis presented in the paper. The DataExtraction scripts create the data files used for the figures presented in the main manuscript while the InitDataExtraction and SensDataExtraction create the data files used for the supplemental visualizations of initial conditions and the sensitivity analysis respectively.

---- ManuscriptFigures.R, SupplementalFigures.R
Finally, these files create the figures presented in the manuscript.