# This script will run range expansion simulations on the MSI compute cluster.

# Set the in file
InFile <- "2019-06-26_DiploidDioSims.csv"
IncreasedSpeed <- 3

# Set the number of processors
nProc <- 24*6

# Set the working directory and load necessary data and libraries
setwd("~/DispersalEvolution/")
library(parallel)
library(Rmpi)

# Read in the data with the SimIDs and corresponding parameter values
SimData <- read.csv(InFile)

# Create the function to be run on the cluster
SimFunc <- function(i){
     SimID <- strsplit(x = as.character(SimData$ID[i]), split = "/")[[1]][4]
     RangeShift(SimDir = SimID, parallel = TRUE, 
                   EquilibriumPrefix = "~/DispersalEvolution/RangeEquilibrium",
                   ShiftPrefix = "~/DispersalEvolution/RangeShift",
                NewSpeed = IncreasedSpeed)
     return(NULL)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- 1:nrow(SimData)

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("SimData", "IncreasedSpeed") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/DispersalEvolution/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


