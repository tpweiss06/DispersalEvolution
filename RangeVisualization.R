# This script will simply make a plot of the realized carrying capacity across
#    space for the stable range according to model parameters

# First, read in the simulation functions and use those
source("~/Desktop/PostdocResearch/DispersalEvolution/GitRepo/SimFunctions.R")

# Set the parameters to use for the visualization
Beta <- 0
Tau <- 15
Gamma <- 0.02
Kmax <- 150
R <- 3
xSeq <- -75:75

# Now calculate the alpha values and the corresponding carrying capacities
AlphaVals <- PatchAlphas(OccPatches = xSeq, NumPatches = length(xSeq), CurBeta = Beta,
                         tau = Tau, gamma = Gamma, Kmax = Kmax, R = R)
Kvals <- log(R) / AlphaVals

# Now plot them
plot(x = xSeq, y = Kvals, main = "Patch carrying capacities", xlab = "Space", 
     ylab = "K")


