# This script will extract a different version of additive genetic variance from
#   the simulations that are functionally asexual (haploid or obligate selfing).
#   For these simulations, selection acts purely on the sum of loci, rather than
#   on individual loci. As such, a more appropriate metric of variance might be
#   variance among the loci sums. This script will calculate that value and add
#   it to the appropriate columns of the results matrix.

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
SummedDeltaGenExp_1 <- rep(NA, nrow(AllSims))
SummedDeltaGenExp_2 <- rep(NA, nrow(AllSims))
SummedDeltaGenShift <- rep(NA, nrow(AllSims))
AllSims <- cbind(AllSims, SummedDeltaGenExp_1, SummedDeltaGenExp_2, SummedDeltaGenShift)

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
     
     # First determine if the population went extinct during the range shift
     ExtRiskEvol <- ifelse(nrow(ShiftPopMat) == 0, 1, 0)
     
     # Calculate the genetic variance for the edge populations in the initial (equilibrium)
     #    population, expansion population, and the range shift population (only if extant)
     # Initial
     InitLociSums <- rowSums(InitialPopMat[InitEdgePop,PopIndices$DispCols])
     InitGenVar <- var(InitLociSums)
     # Expansion
     ExpLociSums_1 <- rowSums(ExpandPopMat[ExpEdgePop_1,PopIndices$DispCols])
     ExpLociSums_2 <- rowSums(ExpandPopMat[ExpEdgePop_2,PopIndices$DispCols])
     ExpGenVar_1 <- var(ExpLociSums_1)
     ExpGenVar_2 <- var(ExpLociSums_2)
     # Shift
     if(ExtRiskEvol == 0){
          ShiftLociSums <- rowSums(ShiftPopMat[ShiftEdgePop,PopIndices$DispCols])
          ShiftGenVar <- var(ShiftLociSums)
     }
     
     # Now use them to calculate change in genetic variance
     DeltaGenExp_1 <- ExpGenVar_1 - InitGenVar
     DeltaGenExp_2 <- ExpGenVar_2 - InitGenVar
     if(ExtRiskEvol == 0){
         DeltaGenShift <- ShiftGenVar - InitGenVar
     } else{
         DeltaGenShift <- NA
     }
     
     # Now return a list of each of these and sort them out after the cluster call, saving the final matrix
     Results <- list(DeltaGenExp_1 = DeltaGenExp_1, DeltaGenExp_2 = DeltaGenExp_2, DeltaGenShift = DeltaGenShift)
     return(Results)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- which(AllSims$Haploid == TRUE | AllSims$omega == 1)

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
    AllSims$SummedDeltaGenExp_1[SimVec[i]] <- Sims[[i]]$DeltaGenExp_1
    AllSims$SummedDeltaGenExp_2[SimVec[i]] <- Sims[[i]]$DeltaGenExp_2
    AllSims$SummedDeltaGenShift[SimVec[i]] <- Sims[[i]]$DeltaGenShift
}
write.csv(AllSims, file = "/project/rangeecoevomodels/cweissle/DispEv/SimsWithSummedResults.csv", 
          row.names = FALSE, quote = FALSE)
