# This script will run simulations on the Teton cluster to first equilibriate
#       a population within a stationary range and then allow that population to
#       (1) expand in an unbounded manner and (2) shift in response to shifting
#       conditions (e.g. climate change).

# Set the number of nodes and the number of tasks per node
# NOTE: this should match the .sh file
nodes <- 20
ntasks_per_node <- 32
TotalTasks <- nodes*ntasks_per_node

# Set the working directory and load necessary data and libraries
setwd("/project/rangeecoevomodels/cweissle/DispEv/")
library(parallel)
library(Rmpi)

# Now set parameters that will be constant across all simulations
BetaInit <- 0
R <- 2
Kmax <- 100
kern <- "exp"
LengthShift <- 200
psi <- 0.5
ExpandGens <- 200
U <- 0.02
sigma <- sqrt(0.02)
NumRands <- 1000000
gamma <- 0.02 
tau <- 15 
lambda <- 10
DispVar <- 10
BurnIn <- 50000
v <- 3

# Make a character vector of the constant parameters to export to the nodes
PartialParams <- c("BetaInit", "gamma", "tau", "R", "Kmax", "kern", 
                   "BurnIn", "LengthShift", "ExpandGens", "psi", "DispVar",
                   "sigma", "U", "lambda", "NumRands", "v")

# Load in the master table of simulations to run
AllSimulations <- read.csv("DispSensSimulations.csv")

# Create the function to be run on the cluster
SimFunc <- function(i){
     # Check if the current simulation has already been completed
     if(AllSimulations$complete[i] == 0){
          # First put together all the parameters in a single list
          CurOmega <- AllSimulations$omega[i]
          CurHap <- AllSimulations$Haploid[i]
          L <- AllSimulations$L[i]
          CurMon <- AllSimulations$monoecious[i]
          CurRho <- AllSimulations$rho[i]
          CurDmax <- AllSimulations$dmax[i]
          CurParams <- list(BetaInit, gamma, tau, CurOmega, U, sigma, L, R, Kmax, 
                            CurHap, kern, CurMon, BurnIn, LengthShift, ExpandGens, 
                            psi, DispVar, CurDmax, CurRho, lambda, NumRands, v)
          names(CurParams) <- c("BetaInit", "gamma", "tau", "omega", "U", "sigma", 
                                "L", "R", "Kmax", "Haploid", "kern", "monoecious", 
                                "BurnIn", "LengthShift", "ExpandGens", "psi", 
                                "DispVar", "dmax", "rho", "lambda", "NumRands", "v")
          SimDir <- "/project/rangeecoevomodels/cweissle/DispEv/DispSensSims"
          SimID <- AllSimulations$SimID[i]
          FullDir <- StationarySim(parameters = CurParams, parallel = TRUE, 
                                   SimID = SimID, SimDirectory = SimDir)
          RangeExpand(SimDir = FullDir, parallel = TRUE)
          RangeShift(SimDir = FullDir, parallel = TRUE)
     }else{
          SimID <- ""
     }
     return(SimID)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- 1:nrow(AllSimulations)

# Create the cluster and run the simulations
cl <- makeCluster(TotalTasks - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c(PartialParams, "AllSimulations"))

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("/project/rangeecoevomodels/cweissle/DispEv/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)




