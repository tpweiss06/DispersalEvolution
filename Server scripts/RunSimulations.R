# This script will run simulations on the Teton cluster to first equilibriate
#       a population within a stationary range and then allow that population to
#       (1) expand in an unbounded manner and (2) shift in response to shifting
#       conditions (e.g. climate change).

# Set the number of nodes and the number of tasks per node
# NOTE: this should match the .sh file
nodes <- 2
ntasks_per_node <- 32
TotalTasks <- nodes*ntasks_per_node

# Set the number of simulations to run for each parameter combination
NumSims <- 10

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

# Now set the values that will vary across simulations
Lseq <- c(1, 2, 4, 8, 16, 32)
# NOTE: Haploid is always associated with a monoecious value of TRUE even though
#    it has no meaning for haploid organisms.
Haploid <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
monoecious <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
omega <- c(0, 0, 0, 0.5, 1)

# Create and populate a list to hold each paramter combination (30 variations)
ParamCombos <- length(Lseq) * length(Haploid)
AllParams <- vector(mode = "list", length = ParamCombos)
k <- 1
for(i in 1:length(Haploid)){
        CurHap <- Haploid[i]
        CurMon <- monoecious[i]
        CurOmega <- omega[i]
        for(l in 1:length(Lseq)){
                L <- Lseq[l]
                AllParams[[k]] <- list(BetaInit, gamma, tau, CurOmega, U, sigma, L, R, Kmax, CurHap, kern,
                                       CurMon, BurnIn, LengthShift, ExpandGens, psi,
                                       DispVar, dmax, rho, lambda, NumRands, v)
                names(AllParams[[k]]) <- c("BetaInit", "gamma", "tau", "omega", "U", "sigma", "L", "R", 
                                           "Kmax", "Haploid", "kern", "monoecious", "BurnIn", "LengthShift", 
                                           "ExpandGens", "psi", "DispVar", "dmax", "rho", 
                                           "lambda", "NumRands", "v")
                k <- k + 1
        }
}

# Create the function to be run on the cluster
SimFunc <- function(i){
     SimDir <- "/project/rangeecoevomodels/cweissle/DispEv/Sims/"
     CurParams <- AllParams[[ParamSeq[i]]]
     SimID <- StationarySim(parameters = CurParams, parallel = TRUE, 
                            SimDirectory = SimDir)
     RangeExpand(SimDir = SimID, parallel = TRUE)
     RangeShift(SimDir = SimID, parallel = TRUE)
     SimInfo <- list(ID = SimID, L = CurParams$L, Haploid = CurParams$Haploid,
                     monoecious = CurParams$monoecious, omega = CurParams$omega)
     # Return a list with the SimID, loci, ploidy, monoecious, and omega values
     return(SimInfo)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- 1:(NumSims*ParamCombos)
ParamSeq <- rep(1:ParamCombos, each = NumSims)

# Create the cluster and run the simulations
cl <- makeCluster(TotalTasks - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("AllParams", "ParamSeq") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("/project/rangeecoevomodels/cweissle/DispEv/SimFunctions.R") )

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

OutFile <- paste(Sys.Date(), ".csv", sep = "")
write.csv(SimIDs, file = OutFile)

Rmpi::mpi.quit()


