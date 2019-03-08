#!/bin/bash -l
#PBS -l walltime=02:00:00,nodes=1:ppn=24,mem=62gb
#PBS -m abe
#PBS -M cweissle@umn.edu
#PBS -j oe

# Set the parameter combination to use and generate names of R scripts and log files
Rscript=2_EquilibriumData.R
LogFile=2_EquilibriumData.log

# Change to the relevant working directory
cd ~/DispersalEvolution/

# Load R and MPI
module load R/3.4.4
module load ompi/3.0.0/gnu-7.2.0-centos7

export RMPI_TYPE=OPENMPI
export OMPI_MCA_mpi_warn_on_fork=0

mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet $Rscript $LogFile
