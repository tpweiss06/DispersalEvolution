Init# This script will extract data from the simulations, which can then be used to
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
nodes <- 4
ntasks_per_node <- 32
TotalTasks <- nodes*ntasks_per_node

# Set the working directory and load necessary data and libraries
setwd("/project/rangeecoevomodels/cweissle/DispEv/")
library(parallel)
library(Rmpi)
source("/project/rangeecoevomodels/cweissle/DispEv/SimFunctions.R")

# Load in the AllSimulations.csv file and determine the patches to extract data
#    from.
source(paste("/project/rangeecoevomodels/cweissle/DispEv/Sims", AllSims$SimID[1], "parameters.R", sep = "/"))
xSeq <- -100:100
Kvals <- log(R) / PatchAlphas(OccPatches = xSeq, NumPatches = length(xSeq), CurBeta = BetaInit,
                              tau = tau, gamma = gamma, Kmax = Kmax, R = R)
InitPatches <- xSeq[which(Kvals >= 10)]

# Load in the matrix with all the simulation information and add columns to hold
#    the needed data
AllSims <- read.csv("AllSimulations.csv")
Phen <- matrix(NA, nrow = nrow(AllSims), ncol = length(InitPatches))
Gen <- matrix(NA, nrow = nrow(AllSims), ncol = length(InitPatches))

# Create the function to be run on the cluster
SimFunc <- function(i){
     # Load in the necessary data objects for each simulation
     CurSim <- paste("/project/rangeecoevomodels/cweissle/DispEv/Sims",
                     AllSims$SimID[i], sep = "/")
     InitialPopMat <- read.csv(paste(CurSim, "EquilibriumPopMat.csv", sep = "/"))
     
     # Read in the relevant paramters to determine the indices of the population matrices
     source(paste(CurSim, "parameters.R", sep = "/"))
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     
     # Next calculate the mean dispersal phenotype and genetic variance of each
     #    sufficiently occupied patch (carrying capacity >= 10) of the initial
     #    population
     Phens <- rep(NA, length(InitPatches))
     GenVars <- rep(NA, length(InitPatches))
     for(x in 1:length(InitPatches)){
          CurPop <- which(InitialPopMat$x0 == InitPatches[x])
          # Phenotype
          CurPhens <- DispPhen(PopMat = InitialPopMat[CurPop,], PopSize = length(CurPop), 
                                    PopIndices = PopIndices, Haploid = Haploid, L = L, dmax = dmax, 
                                    rho = rho, lambda = lambda)
          Phens[x] <- mean(CurPhens, na.rm = TRUE)
          # Genetic variance
          CurGenVarMat <- var(InitialPopMat[CurPop,PopIndices$DispCols])
          GenVars[x] <- sum(CurGenVarMat[lower.tri(InitGenVarMat, diag = TRUE)])
     }
     
     # Now return a list of each of these and sort them out after the cluster call, saving the final matrix
     Results <- list(Phens = Phens, GenVars = GenVars)
     return(Results)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- 1:nrow(AllSims)

# Create the cluster and run the simulations
cl <- makeCluster(TotalTasks - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("AllSims", "InitPatches"))

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("/project/rangeecoevomodels/cweissle/DispEv/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# Now process the Sims and save the updated matrices
for(i in 1:nrow(AllSims)){
     Phen[i,] <- Sims[[i]]$Phens
     Gen[i,] <- Sims[[i]]$GenVars
}

# Create two arrays to store the mean values for each simulation
MeanPhenVals <- array(NA, dim = c(6, 5, length(InitPatches)))
MeanGenVals <- array(NA, dim = c(6, 5, length(InitPatches)))
VarPhenVals <- array(NA, dim = c(6, 5, length(InitPatches)))
VarGenVals <- array(NA, dim = c(6, 5, length(InitPatches)))

# First dimension is loci number: 1,2,4,8,16,32
# Second dimension is population type: Asexual, sexual dioecious, obligatory outcrossing, partial selfing, obligatory selfing

Lseq <- c(1, 2, 4, 8, 16, 32)
# NOTE: Haploid is always associated with a monoecious value of TRUE even though
#    it has no meaning for haploid organisms.
Haploid <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
monoecious <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
omega <- c(0, 0, 0, 0.5, 1)

for(l in 1:length(Lseq)){
     for(i in 1:length(Haploid)){
          CurSims <- AllSims$SimID[which(AllSims$L == Lseq[l] & AllSims$Haploid == Haploid[i] & AllSims$monoecious == monoecious[i] &
                                           AllSims$omega == omega[i])]
          CurPhens <- matrix(data = NA, nrow = length(CurSims), ncol = length(InitPatches))
          CurGen <- matrix(data = NA, nrow = length(CurSims), ncol = length(InitPatches))
          for(j in 1:length(CurSims)){
               CurPhens[j,] <- Phen[CurSims[j],]
               CurGen[j,] <- Gen[CurSims[j],]
          }
          MeanPhenVals[l,i,] <- colMeans(CurPhens, na.rm = TRUE)
          MeanGenVals[l,i,] <- colMeans(CurGen, na.rm = TRUE)
          for(j in 1:length(InitPatches)){
               VarPhenVals[l,i,j] <- var(CurPhens[,j], na.rm = TRUE)
               VarGenVals[l,i,j] <- var(CurGen[,j], na.rm = TRUE)
          }
     }
}

save(MeanPhenVals, MeanGenVals, VarPhenVals, VarGenVals, 
     file = "/project/rangeecoevomodels/cweissle/DispEv/InitialConditions.rdata")

