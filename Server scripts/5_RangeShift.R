# This script will run range expansion simulations on the MSI compute cluster.

SysCommand1 <- "gunzip -r ~/DispersalEvolution/RangeEquilibrium/"
system(SysCommand1)

# Set the in file
InFile1 <- "2019-12-08_DipMono-1-Sims.csv"
#InFile2 <- "2019-10-15_DipDio-Sims.csv"
IncreasedSpeed <- 3

# Set the number of processors
nProc <- 24*6

# Set the working directory and load necessary data and libraries
setwd("~/DispersalEvolution/")
library(parallel)
library(Rmpi)

# Read in the data with the SimIDs and corresponding parameter values
SimData <- read.csv(InFile1)
#SimData2 <- read.csv(InFile2)
#SimData <- rbind(SimData1, SimData2)

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

SysCommand2 <- "gzip -r ~/DispersalEvolution/RangeEquilibrium/"
SysCommand3 <- "gzip -r ~/DispersalEvolution/RangeShift/"
system(SysCommand2)
system(SysCommand3)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


