#!/bin/bash -l

#SBATCH --account=rangeecoevomodels
#SBATCH --time=120:00:00
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cweissle@uwyo.edu
#SBATCH --job-name=DemoSensSims

# Set names of R scripts and log files
Rscript=DemographicSensitivity.R
LogFile=DemographicSensitivity.log

# Change to the relevant working directory
cd /project/rangeecoevomodels/cweissle/DispEv/

# Load R and MPI
module load gcc/7.3.0 r/3.5.3 openmpi/3.1.0 r-rmpi/0.6-9-r353-py27

mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet $Rscript $LogFile
