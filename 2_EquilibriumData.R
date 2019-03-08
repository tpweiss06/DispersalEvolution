# This script will extract information on the spatial distribution of dispersal
#    phenotypes from the equilibrium ranges to (1) assess equilibrium conditions
#    across different parameerizations and confirm that all have reached similar
#    spatial distributions and (2) determine a reasonable value for dThresh to 
#    represent a value requiring evolution to reach, but not overly difficult.

nProc <- 24*1
NumSims <- 100

# Set the working directory and load necessary data and libraries
setwd("~/DispersalEvolution/")
library(parallel)
library(Rmpi)

# Read in the data with the SimIDs and corresponding parameter values
#    NOTE: this section can be adjusted to combine the SimID data from multiple
#    different simulation runs if necessary by using rbind() and loading in
#    multiple data frames.
InFile <- "2019-02-18_Sims.csv"
SimData <- read.csv(InFile)

# Make a list with the rows in the SimData matrix corresponding to each
#    parameter combination
HapSeq <- unique(SimData$Haploid)
Lseq <- unique(SimData$L)
monoeciousSeq <- unique(SimData$monoecious)
omegaSeq <- unique(SimData$omega)
ParamCombos <- vector(mode = "list", length = length(HapSeq) * length(Lseq) * 
                           length(monoeciousSeq) * length(omegaSeq))
k <- 1
for(h in HapSeq){
     for(l in Lseq){
          for(m in monoeciousSeq){
               for(o in omegaSeq){
                    ParamSims <- which((SimData$Haploid == h) & (SimData$L == l) &
                                   (SimData$monoecious == m) & (SimData$omega == o))
                    ParamCombos[[k]] <- ParamSims
                    k <- k + 1
               }
          }
     }
}


# For each parameter combination, I want a matrix with a row for every x value
#    and 3 columns: mean, upper quartile, lower quartile. 
DispExtract <- function(i){
     # First, read in the different equilibrium matrices and record the lowest
     #    and highest x values
     CurSims <- SimData[ParamCombos[[i]],]
     xLwr <- Inf
     xUpr <- -Inf
     for(j in 1:nrow(CurSims)){
          CurData <- read.csv(paste(CurSims$ID[j], "EquilibriumPopMat.csv", sep = "/"))
          xRange <- range(CurData$x0)
          xLwr <- min(xLwr, xRange[1])
          xUpr <- max(xUpr, xRange[2])
     }
     rm(CurData)
     
     # Now make an object to hold all the individual phenotypes at each x value for each simulation
     xVals <- xLwr:xUpr
     NumPatches <- length(xVals)
     DispPhens <- vector(length = NumPatches, mode = "list")
     
     # Now walk through each simulation and population the DispPhens object, starting
     #    by loading in the necessary parameters
     source(paste(CurSims$ID[1], "parameters.R", sep = "/"))
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     for(j in 1:nrow(CurSims)){
          CurData <- read.csv(paste(CurSims$ID[j], "EquilibriumPopMat.csv", sep = "/"))
          for(x in 1:NumPatches){
               CurX <- subset(CurData, x0 == xVals[x])
               # Only proceed if the current simulation has individuals in that patch
               if(dim(CurX)[1] > 0){
                    # Calculate the genotypes of local individuals
                    if((L == 1) & (Haploid)){
                         genotypes <- CurX[,PopIndices$DispCols]
                    }else{
                         if(nrow(CurX) > 1){
                              genotypes <- rowSums(CurX[,PopIndices$DispCols])
                         } else{
                              genotypes <- sum(CurX[,PopIndices$DispCols])
                         }
                    }
                    # Now calculate the genotypes and add them to the appropriate entry
                    #    of DispPhens
                    phenotypes <- (dmax * exp(rho * genotypes)) / (1 + exp(rho * genotypes))
                    DispPhens[[x]] <- c(DispPhens[[x]], phenotypes)
               }
          }
     }
     rm(CurData)
     
     # Now calculate the summary statistics of the data to return
     DispSummary <- matrix(NA, nrow = NumPatches, ncol = 3)
     for(j in 1:NumPatches){
          DispSummary[j,1] <- mean(DispPhens[[j]])
          DispQuantiles <- quantile(DispPhens[[j]], probs = c(0.25, 0.75))
          DispSummary[j,2:3] <- DispQuantiles
     }
     rm(DispPhens)
     return(DispSummary)
}

# Create the cluster and run the extractions
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("SimData", "ParamCombos") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/DispersalEvolution/SimFunctions.R") )

# Run the extractions
ExtractVec <- 1:length(ParamCombos)
DispData <- clusterApply(cl, x = ExtractVec, fun = DispExtract)

# Save the list returned by clusterApply for graphing
save(DispData, ParamCombos, SimData, file = "EquilibriumDispersal.rdata")

