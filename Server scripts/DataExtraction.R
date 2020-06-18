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
AllSims <- read.csv("AllSimulations.csv")
distance <- rep(NA, nrow(AllSims))
DeltaPhenExp <- rep(NA, nrow(AllSims))
DeltaGenExp <- rep(NA, nrow(AllSims))
DeltaPhenShift <- rep(NA, nrow(AllSims))
DeltaGenShift <- rep(NA, nrow(AllSims))
ExtRiskEvol <- rep(NA, nrow(AllSims))
ExtRiskNoEvol <- rep(NA, nrow(AllSims))
AllSims <- cbind(AllSims, distance, DeltaPhenExp, DeltaGenExp, DeltaPhenShift, DeltaGenShift, ExtRiskEvol, ExtRiskNoEvol)

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
     
     # Next calculate the mean dispersal phenotype and genetic variance of the
     #    edge populations in the initial (equilibrium) population, expansion 
     #    population, and the range shift population (only if extant)
     # Initial (Equilibrium)
     InitEdgePop <- which(InitialPopMat$x0 >= max(InitialPopMat$x0) - Edge)
     InitEdgePhens <- DispPhen(PopMat = InitialPopMat[InitEdgePop,], PopSize = length(InitEdgePop), 
                               PopIndices = PopIndices, Haploid = Haploid, L = L, dmax = dmax, 
                               rho = rho, lambda = lambda)
     InitEdgeMeanPhen <- mean(InitEdgePhens)
     InitGenVarMat <- var(InitialPopMat[InitEdgePop,PopIndices$DispCols])
     InitGenVar <- sum(InitGenVarMat[lower.tri(InitGenVarMat, diag = TRUE)])
     # Expansion
     ExpEdgePop <- which(ExpandPopMat$x0 >= max(ExpandPopMat$x0) - Edge)
     ExpEdgePhens <- DispPhen(PopMat = ExpandPopMat[ExpEdgePop,], PopSize = length(ExpEdgePop), 
                              PopIndices = PopIndices, Haploid = Haploid, L = L, dmax = dmax, 
                              rho = rho, lambda = lambda)
     ExpEdgeMeanPhen <- mean(ExpEdgePhens)
     ExpGenVarMat <- var(ExpandPopMat[ExpEdgePop,PopIndices$DispCols])
     ExpGenVar <- sum(ExpGenVarMat[lower.tri(ExpGenVarMat, diag = TRUE)])
     # Shift
     if(ExtRiskEvol == 0){
         ShiftEdgePop <- which(ShiftPopMat$x0 >= max(ShiftPopMat$x0) - Edge)
         ShiftEdgePhens <- DispPhen(PopMat = ShiftPopMat[ShiftEdgePop,], PopSize = length(ShiftEdgePop), 
                                    PopIndices = PopIndices, Haploid = Haploid, L = L, dmax = dmax, 
                                    rho = rho, lambda = lambda)
         ShiftEdgeMeanPhen <- mean(ShiftEdgePhens)
         ShiftGenVarMat <- var(ShiftPopMat[ShiftEdgePop,PopIndices$DispCols])
         ShiftGenVar <- sum(ShiftGenVarMat[lower.tri(ShiftGenVarMat, diag = TRUE)])
     }
     
     # Now use them to calculate each of the necessary data outputs, starting with
     #    the distance travelled in the expansion
     distance <- max(ExpandPopMat$x0)
     DeltaPhenExp <- ExpEdgeMeanPhen - InitEdgeMeanPhen
     DeltaGenExp <- ExpGenVar - InitGenVar
     if(ExtRiskEvol == 0){
         DeltaPhenShift <- ShiftEdgeMeanPhen - InitEdgeMeanPhen
         DeltaGenShift <- ShiftGenVar - InitGenVar
     } else{
         DeltaPhenShift <- NA
         DeltaGenShift <- NA
     }
     
     # Finally, call the NoEvolShift function to assess the extinction risk under
     #  a no evolution scenario
     ExtRiskNoEvol <- NoEvolShift(SimDir = CurSim, parallel = TRUE)
     
     # Now return a list of each of these and sort them out after the cluster call, saving the final matrix
     Results <- list(distance = distance, DeltaPhenExp = DeltaPhenExp, DeltaGenExp = DeltaGenExp, 
                      DeltaPhenShift = DeltaPhenShift, DeltaGenShift = DeltaGenShift, 
                      ExtRiskEvol = ExtRiskEvol, ExtRiskNoEvol = ExtRiskNoEvol)
     return(Results)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- 1:nrow(AllSims)

# Create the cluster and run the simulations
cl <- makeCluster(TotalTasks - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("AllSims", "Edge"))

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("/project/rangeecoevomodels/cweissle/DispEv/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# Now process the Sims and save the updated matrix
for(i in 1:nrow(AllSims)){
    AllSims$distance <- Sims[[i]]$distance
    AllSims$DeltaPhenExp <- Sims[[i]]$DeltaPhenExp
    AllSims$DeltaGenExp <- Sims[[i]]$DeltaGenExp
    AllSims$DeltaPhenShift <- Sims[[i]]$DeltaPhenShift
    AllSims$DeltaGenShift <- Sims[[i]]$DeltaGenShift
    AllSims$ExtRiskEvol <- Sims[[i]]$ExtRiskEvol
    AllSims$ExtRiskNoEvol <- Sims[[i]]$ExtRiskNoEvol
}
write.csv(AllSims, file = "/project/rangeecoevomodels/cweissle/DispEv/SimsWithResults.csv", 
          row.names = FALSE, quote = FALSE)
