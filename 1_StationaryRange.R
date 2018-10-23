# This script will run simulations on the MSI compute cluster to characterize
#    population structure in a stable range for all of the 9 range parameter
#    combinations to be explored.

# Set the number of processors and number of simulations to be run
nProc <- 24*1
NumSims <- 23

# Set the working directory and load necessary data and libraries
setwd("~/DispersalEvolution/")
library(parallel)
library(Rmpi)

# Now set parameters that will be constant across all simulations
BetaInit <- 0
gamma <- 0.1
tau <- 20
nu <- 0.001              # Value taken from Gilbert et al. 2017
sigma <- sqrt(0.02)      # Value taken from Gilbert et al. 2017
R <- 2
Kmax <- 50
kern <- "exp"
BurnIn <- 5000
LengthShift <- 0
dThresh <- NA
InitPopSize <- Kmax
psi <- 0.5
DispVar <- 2
dmax <- 10
rho <- 0.1
NumRands <- 1000000

# Now set parameters that will change across scenarios
# This will be done in a more sophisticated way later on
monoecious <- TRUE
Haploid <- TRUE
L <- 1
omega <- 0

Params <- list(BetaInit, gamma, tau, omega, nu, sigma, L, R, Kmax, Haploid, kern,
               monoecious, BurnIn, LengthShift, dThresh, InitPopSize, psi,
               DispVar, dmax, rho, NumRands)
names(Params) <- c("BetaInit", "gamma", "tau", "omega", "nu", "sigma", "L", "R", 
                   "Kmax", "Haploid", "kern", "monoecious", "BurnIn", "LengthShift", 
                   "dThresh", "InitPopSize", "psi", "DispVar", "dmax", "rho", "NumRands")
#AllParams <- vector(mode = "list", length = 9)

#for(i in 1:9){
#     gamma <- RangeParams$gamma[i]
#     lambda <- RangeParams$lambda[i]
#     tau <- RangeParams$tau[i]
#     eta <- RangeParams$eta[i]
     
#     AllParams[[i]] <- list(BetaInit, gamma, tau, lambda, omega, U, Vm, Lf, Ld, Rmax,
#                        Kmax, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
#                        LengthShift, v, InitPopSize, FitInit, FitDiv, DispInit,
#                        DispDiv, eta, NumRands, z, dmax, rho)
#     names(AllParams[[i]]) <- c("BetaInit", "gamma", "tau", "lambda", "omega", "U", "Vm", 
#                            "Lf", "Ld", "Rmax", "Kmax", "width", "kern", "EnvGradType", 
#                            "monoecious", "BurnIn", "BurnOut", "LengthShift", "v", 
#                            "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
#                            "eta", "NumRands", "z", "dmax", "rho")
#}

# Write a function to be passed to various nodes
SimFunc <- function(i){
     StationarySim(parameters = Params, parallel = TRUE, SimID = paste("Trial", i, sep = ""))
     return(i)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- 1:NumSims

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("Params") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/DispersalEvolution/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


