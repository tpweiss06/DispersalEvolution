# This script will extract the within genome variance from all of the simulations.

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
AllSims <- read.csv("SimsWithResults.csv")
WithinExp_1 <- rep(NA, nrow(AllSims))
WithinExp_2 <- rep(NA, nrow(AllSims))
WithinShift <- rep(NA, nrow(AllSims))
WithinInit <- rep(NA, nrow(AllSims))
AllSims <- cbind(AllSims, WithinInit, WithinExp_1, WithinExp_2, WithinShift)

# Define the percentage of available patches constituting the range edge and use
#    that percentage to calculate the equivalent number of patches
EdgeProp <- 0.25
source(paste("/project/rangeecoevomodels/cweissle/DispEv/Sims", AllSims$SimID[1], "parameters.R", sep = "/"))
xSeq <- -100:100
Kvals <- log(R) / PatchAlphas(OccPatches = xSeq, NumPatches = length(xSeq), CurBeta = BetaInit,
                              tau = tau, gamma = gamma, Kmax = Kmax, R = R)
RangeExtent <- sum(Kvals > 0)
if(RangeExtent == length(xSeq)){
     print("Choose a larger range of x values to fully capture the range")
     ErrorOut
}
Edge <- round(EdgeProp * RangeExtent)

# Create the function to be run on the cluster
SimFunc <- function(i){
     # Load in the necessary data objects for each simulation
     CurSim <- paste("/project/rangeecoevomodels/cweissle/DispEv/Sims",
                     AllSims$SimID[i], sep = "/")
     InitialPopMat <- read.csv(paste(CurSim, "EquilibriumPopMat.csv", sep = "/"))
     ExpandPopMat <- read.csv(paste(CurSim, "ExpansionEdgePopMat.csv", sep = "/"))
     ShiftPopMat <- read.csv(paste(CurSim, "ShiftPopMat.csv", sep = "/"))
     
     # Read in the relevant paramters to determine the indices of the population matrices
     source(paste(CurSim, "parameters.R", sep = "/"))
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     
     # Define the edge populations to use for these calculations
     ExpEdgePop_1 <- which(ExpandPopMat$x0 >= max(ExpandPopMat$x0) - Edge)
     ExpEdgePop_2 <- which(ExpandPopMat$x0 <= min(ExpandPopMat$x0) + Edge)
     
     # Calculate the genetic variance for the edge populations in the initial (equilibrium)
     #    population, expansion population, and the range shift population (only if extant)
     InitInds <- nrow(InitialPopMat)
     Exp1Inds <- length(ExpEdgePop_1)
     Exp2Inds <- length(ExpEdgePop_2)
     ShiftInds <- nrow(ShiftPopMat)
     MaxInds <- max(InitInds, Exp2Inds, Exp1Inds, ShiftInds)
     AllWithinInit <- rep(NA, InitInds)
     AllWithinExp_1 <- rep(NA, Exp1Inds)
     AllWithinExp_2 <- rep(NA, Exp2Inds)
     AllWithinShift <- rep(ShiftInds)
     
     for(j in 1:MaxInds){
         if(j <= InitInds){
             AllWithinInit[j] <- var(as.numeric(InitialPopMat[j,PopIndices$DispCols]))
         }
         if(j <= Exp1Inds){
             AllWithinExp_1[j] <- var(as.numeric(ExpandPopMat[ExpEdgePop_1[j],PopIndices$DispCols]))
         }
         if(j <= Exp2Inds){
             AllWithinExp_2[j] <- var(as.numeric(ExpandPopMat[ExpEdgePop_2[j],PopIndices$DispCols]))
         }
         if(j <= ShiftInds){
             AllWithinShift[j] <- var(as.numeric(ShiftPopMat[j,PopIndices$DispCols]))
         }
     }
     
     # Now return a list of each of these and sort them out after the cluster call, saving the final matrix
     Results <- list(WithinInit = mean(AllWithinInit), WithinExp_1 = mean(AllWithinExp_1),
                     WithinExp_2 = mean(AllWithinExp_2), WithinShift = ifelse(ShiftInds == 0, NA, mean(AllWithinShift)))
     return(Results)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- 1001:30000 # Exclude the haploid simulations with 1 locus as within genome variance is undefined for them.

# Create the cluster and run the simulations
cl <- makeCluster(TotalTasks - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("AllSims", "Edge"))

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("/project/rangeecoevomodels/cweissle/DispEv/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# Now process the Sims and save the updated matrix
for(i in 1:length(SimVec)){
    AllSims$WithinInit[SimVec[i]] <- Sims[[i]]$WithinInit
    AllSims$WithinExp_1[SimVec[i]] <- Sims[[i]]$WithinExp_1
    AllSims$WithinExp_2[SimVec[i]] <- Sims[[i]]$WithinExp_2
    AllSims$WithinShift[SimVec[i]] <- Sims[[i]]$WithinShift
}
write.csv(AllSims, file = "/project/rangeecoevomodels/cweissle/DispEv/WithinGenVarResults.csv", 
          row.names = FALSE, quote = FALSE)
