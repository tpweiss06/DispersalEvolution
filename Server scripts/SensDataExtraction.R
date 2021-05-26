# This script will extract data from the sensitivity analysis simulations. Specifically,
#    the script will record the mean change in phenotype, genetic diversity, and distance
#    from unbounded range expansions and the mean change in phenotype, genetic diversity,
#    and extinction risk from climate driven range shifts.

# Set the number of nodes and the number of tasks per node
# NOTE: this should match the .sh file
nodes <- 4
ntasks_per_node <- 32
TotalTasks <- nodes*ntasks_per_node

# Set the working directory and load necessary data and libraries
setwd("/project/rangeecoevomodels/cweissle/DispEv/")
library(parallel)
library(Rmpi)
source("/project/rangeecoevomodels/cweissle/DispEv/SimFunctions.R")

# Load in the matrix with all the simulation information and add columns to hold
#    the needed data
AllSims <- read.csv("SensSimulations.csv")
distance_1 <- rep(NA, nrow(AllSims))
DeltaPhenExp_1 <- rep(NA, nrow(AllSims))
DeltaGenExp_1 <- rep(NA, nrow(AllSims))
distance_2 <- rep(NA, nrow(AllSims))
DeltaPhenExp_2 <- rep(NA, nrow(AllSims))
DeltaGenExp_2 <- rep(NA, nrow(AllSims))
DeltaPhenShift <- rep(NA, nrow(AllSims))
DeltaGenShift <- rep(NA, nrow(AllSims))
ExtRisk <- rep(NA, nrow(AllSims))
AllSims <- cbind(AllSims, distance_1, distance_2, DeltaPhenExp_1, DeltaPhenExp_2, 
                 DeltaGenExp_1, DeltaGenExp_2, DeltaPhenShift, DeltaGenShift, 
                 ExtRisk)

# Define the percentage of available patches constituting the range edge and use
#    that percentage to calculate the equivalent number of patches
EdgeProp <- 0.25
source(paste("/project/rangeecoevomodels/cweissle/DispEv/SensSims", AllSims$SimID[1], "parameters.R", sep = "/"))
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
     CurSim <- paste("/project/rangeecoevomodels/cweissle/DispEv/SensSims",
                     AllSims$SimID[i], sep = "/")
     InitialPopMat <- read.csv(paste(CurSim, "EquilibriumPopMat.csv", sep = "/"))
     ExpandPopMat <- read.csv(paste(CurSim, "ExpansionEdgePopMat.csv", sep = "/"))
     ShiftPopMat <- read.csv(paste(CurSim, "ShiftPopMat.csv", sep = "/"))
     
     # Read in the relevant paramters to determine the indices of the population matrices
     source(paste(CurSim, "parameters.R", sep = "/"))
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     
     # First determine if the population went extinct during the range shift
     ExtRisk <- ifelse(nrow(ShiftPopMat) == 0, 1, 0)
     
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
     ExpEdgePop_1 <- which(ExpandPopMat$x0 >= max(ExpandPopMat$x0) - Edge)
     ExpEdgePop_2 <- which(ExpandPopMat$x0 <= min(ExpandPopMat$x0) + Edge)
     ExpEdgePhens_1 <- DispPhen(PopMat = ExpandPopMat[ExpEdgePop_1,], PopSize = length(ExpEdgePop_1), 
                              PopIndices = PopIndices, Haploid = Haploid, L = L, dmax = dmax, 
                              rho = rho, lambda = lambda)
     ExpEdgePhens_2 <- DispPhen(PopMat = ExpandPopMat[ExpEdgePop_2,], PopSize = length(ExpEdgePop_2), 
                                PopIndices = PopIndices, Haploid = Haploid, L = L, dmax = dmax, 
                                rho = rho, lambda = lambda)
     ExpEdgeMeanPhen_1 <- mean(ExpEdgePhens_1)
     ExpEdgeMeanPhen_2 <- mean(ExpEdgePhens_2)
     ExpGenVarMat_1 <- var(ExpandPopMat[ExpEdgePop_1,PopIndices$DispCols])
     ExpGenVarMat_2 <- var(ExpandPopMat[ExpEdgePop_2,PopIndices$DispCols])
     ExpGenVar_1 <- sum(ExpGenVarMat_1[lower.tri(ExpGenVarMat_1, diag = TRUE)])
     ExpGenVar_2 <- sum(ExpGenVarMat_2[lower.tri(ExpGenVarMat_2, diag = TRUE)])
     # Shift
     if(ExtRisk == 0){
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
     distance_1 <- max(ExpandPopMat$x0)
     distance_2 <- abs(min(ExpandPopMat$x0))
     DeltaPhenExp_1 <- ExpEdgeMeanPhen_1 - InitEdgeMeanPhen
     DeltaPhenExp_2 <- ExpEdgeMeanPhen_2 - InitEdgeMeanPhen
     DeltaGenExp_1 <- ExpGenVar_1 - InitGenVar
     DeltaGenExp_2 <- ExpGenVar_2 - InitGenVar
     if(ExtRisk == 0){
         DeltaPhenShift <- ShiftEdgeMeanPhen - InitEdgeMeanPhen
         DeltaGenShift <- ShiftGenVar - InitGenVar
     } else{
         DeltaPhenShift <- NA
         DeltaGenShift <- NA
     }
     
     # Now return a list of each of these and sort them out after the cluster call, saving the final matrix
     Results <- list(distance_1 = distance_1, distance_2 = distance_2, DeltaPhenExp_1 = DeltaPhenExp_1, 
                     DeltaPhenExp_2 = DeltaPhenExp_2, DeltaGenExp_1 = DeltaGenExp_1, DeltaGenExp_2 = DeltaGenExp_2,
                     DeltaPhenShift = DeltaPhenShift, DeltaGenShift = DeltaGenShift, ExtRisk = ExtRisk)
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
    AllSims$distance_1[i] <- Sims[[i]]$distance_1
    AllSims$DeltaPhenExp_1[i] <- Sims[[i]]$DeltaPhenExp_1
    AllSims$DeltaGenExp_1[i] <- Sims[[i]]$DeltaGenExp_1
    AllSims$distance_2[i] <- Sims[[i]]$distance_2
    AllSims$DeltaPhenExp_2[i] <- Sims[[i]]$DeltaPhenExp_2
    AllSims$DeltaGenExp_2[i] <- Sims[[i]]$DeltaGenExp_2
    AllSims$DeltaPhenShift[i] <- Sims[[i]]$DeltaPhenShift
    AllSims$DeltaGenShift[i] <- Sims[[i]]$DeltaGenShift
    AllSims$ExtRisk[i] <- Sims[[i]]$ExtRisk
}
write.csv(AllSims, file = "/project/rangeecoevomodels/cweissle/DispEv/SensSimsWithResults.csv", 
          row.names = FALSE, quote = FALSE)
