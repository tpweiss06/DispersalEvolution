source("SimFunctions.R")

BetaInit <- 0
U <- 0.02                # Value taken from Gilbert et al. 2017
sigma <- sqrt(0.02)      # Value taken from Gilbert et al. 2017
R <- 2
Kmax <- 100
kern <- "exp"
LengthShift <- 200
psi <- 0.5
dmax <- 6
ExpandGens <- 100
EdgeThresh <- 0.95
rho <- 0.1
NumRands <- 1000000
omega <- 1
gamma <- 0.02 
tau <- 15 
lambda <- 10
DispVar <- 10
BurnIn <- 1000
v <- 1
Haploid <- TRUE
L <- 5
monoecious <- TRUE

AllParams <- list(BetaInit, gamma, tau, omega, U, sigma, L, R, Kmax, Haploid, kern,
                       monoecious, BurnIn, LengthShift, ExpandGens, EdgeThresh, psi,
                       DispVar, dmax, rho, lambda, NumRands, v)
names(AllParams) <- c("BetaInit", "gamma", "tau", "omega", "U", "sigma", "L", "R", 
                           "Kmax", "Haploid", "kern", "monoecious", "BurnIn", "LengthShift", 
                           "ExpandGens", "EdgeThresh", "psi", "DispVar", "dmax", "rho", 
                           "lambda", "NumRands", "v")

EquilibriumDir <- "~/Desktop/TestRun/RangeEquilibrium"
ExpandDir <- "~/Desktop/TestRun/RangeExpand"
ShiftDir <- "~/Desktop/TestRun/RangeShift"

# Test the stationary range function
SimID <- StationarySim(parameters = AllParams, parallel = FALSE, 
                       SimDirectory = EquilibriumDir)
SimID <- "Sim1"
# Test the RangeExpand function
RangeExpand(SimDir = SimID, parallel = FALSE, EquilibriumPrefix = EquilibriumDir,
            ExpandPrefix = ExpandDir)

# Test the RangeShift function
RangeShift(SimDir = SimID, parallel = FALSE, EquilibriumPrefix = EquilibriumDir, 
           ShiftPrefix = ShiftDir)


# Cool, this seems to be working now and ready to go on the server