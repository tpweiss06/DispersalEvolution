# This script will run simulations on the MSI compute cluster to characterize
#    population structure in a stable range for all of the 9 range parameter
#    combinations to be explored.

# Set the number of processors and number of simulations to be run
nProc <- 24*18
NumSims <- 200

# Set the working directory and load necessary data and libraries
setwd("~/ShiftingSlopes/MainSim/")
source("~/ShiftingSlopes/SimFunctions.R")
RangeParams <- read.csv("~/ShiftingSlopes/RangeParameters.csv")
library(parallel)
library(Rmpi)

# Now set all parameters
BetaInit <- 0
omega <- 3
U <- c(0.02, 0.02)
Vm <- c(0.0004, 0.0004)
Lf <- 5
Ld <- 5
Rmax <- 2
Kmax <- 100
width <- 10
kern <- "exp"
EnvGradType <- "K"
monoecious <- FALSE
BurnIn <- 2000
BurnOut <- 50
LengthShift <- 100
v <- 1
InitPopSize <- 2500
FitInit <- 0
FitDiv <- 0.025
DispInit <- -1
DispDiv <- 1
NumRands <- 1000000
z <- 0.5
dmax <- 1000
rho <- 0.5

AllParams <- vector(mode = "list", length = 9)

for(i in 1:9){
     gamma <- RangeParams$gamma[i]
     lambda <- RangeParams$lambda[i]
     tau <- RangeParams$tau[i]
     eta <- RangeParams$eta[i]
     
     AllParams[[i]] <- list(BetaInit, gamma, tau, lambda, omega, U, Vm, Lf, Ld, Rmax,
                        Kmax, width, kern, EnvGradType, monoecious, BurnIn, BurnOut, 
                        LengthShift, v, InitPopSize, FitInit, FitDiv, DispInit,
                        DispDiv, eta, NumRands, z, dmax, rho)
     names(AllParams[[i]]) <- c("BetaInit", "gamma", "tau", "lambda", "omega", "U", "Vm", 
                            "Lf", "Ld", "Rmax", "Kmax", "width", "kern", "EnvGradType", 
                            "monoecious", "BurnIn", "BurnOut", "LengthShift", "v", 
                            "InitPopSize", "FitInit", "FitDiv", "DispInit", "DispDiv", 
                            "eta", "NumRands", "z", "dmax", "rho")
}

# Create the population column indices to make an input matrix
PopIndices <- PopMatColNames(Lf = Lf, Ld = Ld, monoecious = monoecious)

# Now initialize the matrix with an entry for each individual
nCol <- ifelse(monoecious, max(PopIndices$DispCols), PopIndices$sex)
InputMat <- matrix(NA, nrow = InitPopSize, ncol = nCol)

# Now step through and populate the matrix
Ylocs <- sample(x = 1:width, size = InitPopSize, replace = TRUE, prob = NULL)
SectorPop <- InitPopSize / 5
CurX <- c(rep(-2, SectorPop), rep(-1, SectorPop), rep(0, SectorPop), 
          rep(1, SectorPop), rep(2, SectorPop))
PatchPop <- SectorPop / width
CurY <- rep(1:width, each = PatchPop)
CurY <- c(CurY, CurY, CurY, CurY, CurY)
for(i in 1:InitPopSize){
     InputMat[i, PopIndices$x1] <- CurX[i]
     InputMat[i, PopIndices$y1] <- CurY[i]
     InputMat[i, PopIndices$FitCols] <- rnorm(n = Lf * 2, mean = FitInit, sd = FitDiv)
     InputMat[i, PopIndices$DispCols] <- rnorm(n = Ld * 2, mean = DispInit,
                                               sd = DispDiv)
     if(!monoecious){
          InputMat[i, PopIndices$sex] <- rbinom(n = 1, size = 1, prob = z)
     }
}

# Write a function to be passed to various nodes
SimFunc <- function(i){
     setwd(paste("~/ShiftingSlopes/MainSim/Params", i, "/", sep = ""))
     FullSim(parameters = AllParams[[i]], parallel = TRUE, PopInit = InputMat)
     return(i)
}

# Create a vector of parameter index values for the parallel computation
SimVec <- rep(1:9, each = NumSims)

# Create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("AllParams", "InputMat") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/ShiftingSlopes/SimFunctions.R") )

# Run the simulations
Sims <- clusterApply(cl, x = SimVec, fun = SimFunc)

# This implementation of openMPI doesn't seem to play nice with Rmpi,
#    so instead of nicely shutting down the cluster within the script,
#    just let the script end here and then MSI will kill the cluster.


