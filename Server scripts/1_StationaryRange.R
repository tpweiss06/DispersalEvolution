# This script will run simulations on the MSI compute cluster to characterize
#    population structure in a stable range for all parameter combinations.

# Set the type of organisms to be simulated
# NOTE: Haploid is always associated with a monoecious value of TRUE even though
#    it has no meaning for haploid organisms.
Haploid <- FALSE
monoecious <- TRUE
omega <- 1   # start with 0, 0.5, and 1
OutFile <- paste(Sys.Date(), "DipMono-1-Sims.csv", sep = "_")

# Set the number of processors and number of simulations to be run
nProc <- 24*9
NumSims <- 1000

# Set the working directory and load necessary data and libraries
setwd("~/DispersalEvolution/")
library(parallel)
library(Rmpi)

# Now set parameters that will be constant across all simulations
BetaInit <- 0
U <- 0.02                # Value taken from Gilbert et al. 2017
sigma <- sqrt(0.02)      # Value taken from Gilbert et al. 2017
R <- 2
Kmax <- 100
kern <- "exp"
LengthShift <- 200
psi <- 0.5
dmax <- 6
ExpandGens <- 200
rho <- 0.1
NumRands <- 1000000
gamma <- 0.02 
tau <- 15 
lambda <- 10
DispVar <- 10
BurnIn <- 50000
v <- 3

# Now set the L values to be used across simulations
Lseq <- c(1, 2, 4, 8, 16, 32)
AllParams <- vector(mode = "list", length = length(Lseq))
for(l in 1:length(Lseq)){
     L <- Lseq[l]
     AllParams[[l]] <- list(BetaInit, gamma, tau, omega, U, sigma, L, R, Kmax, Haploid, kern,
                            monoecious, BurnIn, LengthShift, ExpandGens, psi,
                            DispVar, dmax, rho, lambda, NumRands, v)
     names(AllParams[[l]]) <- c("BetaInit", "gamma", "tau", "omega", "U", "sigma", "L", "R", 
                                "Kmax", "Haploid", "kern", "monoecious", "BurnIn", "LengthShift", 
                                "ExpandGens", "psi", "DispVar", "dmax", "rho", 
                                "lambda", "NumRands", "v")
}

# Create the function to be run on the cluster
SimFunc <- function(i){
     SimDir <- "~/DispersalEvolution/RangeEquilibrium"
     CurParams <- AllParams[[ParamSeq[i]]]
     SimID <- StationarySim(parameters = CurParams, parallel = TRUE, 
                            SimDirectory = SimDir)
     SimInfo <- list(ID = SimID, L = CurParams$L, Haploid = CurParams$Haploid,
                     monoecious = CurParams$monoecious, omega = CurParams$omega)
     # Return a list with the SimID, loci, ploidy, monoecious, and omega values
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

write.csv(SimIDs, file = OutFile)

system("gzip -r ~/DispersalEvolution/RangeEquilibrium")

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


