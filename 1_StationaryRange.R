# This script will run simulations on the MSI compute cluster to characterize
#    population structure in a stable range for all parameter combinations.

# Set the number of processors and number of simulations to be run
nProc <- 24*6
NumSims <- 20

# Set the working directory and load necessary data and libraries
setwd("~/DispersalEvolution/")
library(parallel)
library(Rmpi)

# Now set parameters that will be constant across all simulations
BetaInit <- 0
U <- 0.02                # Value taken from Gilbert et al. 2017
sigma <- sqrt(0.02)      # Value taken from Gilbert et al. 2017
R <- 2
Kmax <- 50
kern <- "exp"
LengthShift <- 100
InitPopSize <- Kmax
psi <- 0.5
dmax <- 6
dThresh <- 3.5
EdgeThresh <- 0.95
rho <- 0.1
NumRands <- 1000000
omega <- 0.5
gamma <- 0.02 
tau <- 15 
lambda <- 10
DispVar <- 10
BurnIn <- 5000
monoecious <- FALSE

# Now set parameters that will change across these initial scenarios: haploid
#    and L
HapSeq <- c(TRUE, FALSE)
Lseq <- 1:10
AllParams <- vector(mode = "list", length = length(Lseq) * length(HapSeq))
k <- 1
for(h in HapSeq){
     for(l in Lseq){
          Haploid <- h
          L <- l
          AllParams[[k]] <- list(BetaInit, gamma, tau, omega, U, sigma, L, R, Kmax, Haploid, kern,
                                 monoecious, BurnIn, LengthShift, dThresh, EdgeThresh, InitPopSize, psi,
                                 DispVar, dmax, rho, lambda, NumRands)
          names(AllParams[[k]]) <- c("BetaInit", "gamma", "tau", "omega", "U", "sigma", "L", "R", 
                                     "Kmax", "Haploid", "kern", "monoecious", "BurnIn", "LengthShift", 
                                     "dThresh", "EdgeThresh", "InitPopSize", "psi", "DispVar", "dmax", "rho", 
                                     "lambda", "NumRands")
          k <- k + 1
     }
}

# Create the function to be run on the cluster
SimFunc <- function(i){
     SimDir <- "~/DispersalEvolution/RangeEquilibrium"
     CurParams <- AllParams[[ParamSeq[i]]]
     SimID <- StationarySim(parameters = CurParams, parallel = TRUE, 
                            SimDirectory = SimDir)
     SimInfo <- list(ID = SimID, L = CurParams$L, Haploid = CurParams$Haploid,
                     monoecious = CurParams$monoecious, omega = CurParams$omega)
     # Return a list with the SimID, loci, ploidy, monoecious, omega? (selfing parameter)
     return(SimInfo)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- 1:(NumSims*length(AllParams))
ParamSeq <- rep(1:length(AllParams), each = NumSims)

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("AllParams", "ParamSeq") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/DispersalEvolution/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# Sort all the SimIDs and the related parameters into a single data frame
TotalSims <- length(Sims)
SimIDs <- data.frame(ID = rep(NA, TotalSims), L = rep(NA, TotalSims), 
                     Haploid = rep(NA, TotalSims), monoecious = rep(NA, TotalSims), 
                     omega = rep(NA, TotalSims))
for(i in 1:TotalSims){
     SimIDs$ID[i] <- Sims[[i]]$ID
     SimIDs$L[i] <- Sims[[i]]$L
     SimIDs$Haploid[i] <- Sims[[i]]$Haploid
     SimIDs$monoecious[i] <- Sims[[i]]$monoecious
     SimIDs$omega[i] <- Sims[[i]]$omega
}

OutFile <- paste(Sys.Date(), "Sims.csv", sep = "_")
write.csv(SimIDs, file = OutFile)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


