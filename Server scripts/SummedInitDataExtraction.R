# This script will extract data from the simulations, which can then be used to
#    make the figures for the manuscripts. There are several quantities this script
#    will focus on when extracting the data.
# From the simulations of an unbounded expansion, this script will extract:
#    1) The mean distance travelled (and CV) for each parameter combination
#    2) The change in mean dispersal phenotype (and genetic variation) at the leading
#         edge of expansion compared to the initial population. The script will
#         also record the interquartile ranges for these quantities
# From the simulations of range shifts, this script will extract:
#    1) The change in mean dispersal phenotype (and genetic variation) at the leading
#         edge of the shift compared to the initial population (with IQR)
#    2) The extinction risk for each parameter combination (proportion of simulations
#         to go extinct during the simulation)
# Finally, this script will call the NoEvolShift function to evaluate the extinction
#    risk of each parameter combination in the absence of dispersal evolution

# Set the number of nodes and the number of tasks per node
# NOTE: this should match the .sh file
nodes <- 2
ntasks_per_node <- 32
TotalTasks <- nodes*ntasks_per_node

# Set the working directory and load necessary data and libraries
setwd("/project/rangeecoevomodels/cweissle/DispEv/")
library(parallel)
library(Rmpi)
source("/project/rangeecoevomodels/cweissle/DispEv/SimFunctions.R")

# Load in the matrix with all the simulation information and add columns to hold
#    the needed data
AllSims <- read.csv("SimsWithInitVals.csv")
SummedGenVar <- rep(NA, nrow(AllSims))
AllSims <- cbind(AllSims, SummedGenVar)

# Create the function to be run on the cluster
SimFunc <- function(i){
     # Load in the necessary data objects for each simulation
     CurSim <- paste("/project/rangeecoevomodels/cweissle/DispEv/Sims",
                     AllSims$SimID[i], sep = "/")
     InitialPopMat <- read.csv(paste(CurSim, "EquilibriumPopMat.csv", sep = "/"))
     
     # Read in the relevant paramters to determine the indices of the population matrices
     source(paste(CurSim, "parameters.R", sep = "/"))
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     
     # Next calculate the genetic variance across the whole initial population
     LociSums <- rowSums(InitialPopMat[,PopIndices$DispCols])
     GenVar <- var(LociSums)
     
     # Now return the genetic variance of the initial population
     return(GenVar)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- which((AllSims$Haploid == TRUE & AllSims$L > 1) | AllSims$omega == 1)

# Create the cluster and run the simulations
cl <- makeCluster(TotalTasks - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, "AllSims")

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("/project/rangeecoevomodels/cweissle/DispEv/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# Now process the Sims and save the updated matrix
for(i in 1:length(SimVec)){
        AllSims$SummedGenVar[SimVec[i]] <- Sims[[i]]
}

write.csv(AllSims, file = "/project/rangeecoevomodels/cweissle/DispEv/SimsWithSummedInitVals.csv", 
          row.names = FALSE, quote = FALSE)

