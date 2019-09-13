# This script will extract information on the spatial distribution of dispersal
#    phenotypes during the range expansion simulations

# Set the in and outfiles for the simulations and data
InFile <- ""
OutFile <- ""

# Set the working directory and load necessary data and libraries
setwd("~/DispersalEvolution/")
library(parallel)
library(Rmpi)

# Read in the data with the SimIDs and corresponding parameter values
SimData <- read.csv(InFile)

# Make a list with the rows in the SimData matrix corresponding to each
#    parameter combination
Lseq <- unique(SimData$L)
#omegaSeq <- unique(SimData$omega)
ParamCombos <- vector(mode = "list", length = length(Lseq))
for(l in 1:length(Lseq)){
     ParamCombos[[l]] <- which( (SimData$L == Lseq[l]) )
}

# Set the number of processors needed for the number of parameter combinations
nProc <- length(ParamCombos) + 1

# Create other objects to be passed to the function
ExampleParams <- paste(SimData$ID[1], "/parameters.R", sep = "")
source(ExampleParams)
Gens <- seq(0, ExpandGens, by = ExpandGens/10)
Ngens <- length(Gens)

# For each parameter combination, the function will return a list with an entry
#    for each time point consisting of a data frame with columns for x coordinate,
#    mean phenotype, lower IQR, and upper IQR.
DispExtract <- function(i){
     # First, subset the simulations corresponding to the current parameter
     #    combination and record the x bounds through time (to construct the
     #    data frames later)
     CurSims <- SimData[ParamCombos[[i]],]
     Nsims <- nrow(CurSims)
     xLwr <- rep(Inf, Ngens)
     xUpr <- rep(-Inf, Ngens)
     for(t in 1:Ngens){
          for(j in 1:Nsims){
               SimID <- strsplit(x = as.character(CurSims$ID[j]), split = "/")[[1]][4]
               CurFile <- paste("~/DispersalEvolution/RangeExpansion", SimID, "SummaryStats.csv", sep = "/")
               CurData <- read.csv(CurFile)
               CurData <- subset(CurData, gen == Gens[t])
               xRange <- range(CurData$x)
               xLwr[t] <- min(xLwr[t], xRange[1])
               xUpr[t] <- max(xUpr[t], xRange[2])
          }
     }
     rm(CurData)
     
     # Construct the list to hold the results through time and load parameters
     DispPhens <- vector(length = Ngens, mode = "list")
     DispGenVar <- vector(length = Ngens, mode = "list")
     source(paste("~/DispersalEvolution/RangeEquilibrium", SimID, "parameters.R", sep = "/"))
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     
     # Now walk through time, filling the list with information from all simulations
     for(t in 1:Ngens){
          CurXseq <- xLwr[t]:xUpr[t]
          N <- length(CurXseq)
          DispPhens[[t]] <- data.frame(x = CurXseq, dBar = rep(NA, N), 
                                       lwr = rep(NA, N), upr = rep(NA, N))
          DispGenVar[[t]] <- data.frame(x = CurXseq, GenVar = rep(NA, N), 
                                        lwr = rep(NA, N), upr = rep(NA, N))
          
          # Create matrices to hold the mean values in each patch from each 
          #    simulation and step through the simulations to fill them up
          All_dBars <- matrix(NA, nrow = N, ncol = Nsims)
          AllGenVar <- matrix(NA, nrow = N, ncol = Nsims)
          for(j in 1:Nsims){
               SimID <- strsplit(x = as.character(CurSims$ID[j]), split = "/")[[1]][4]
               CurFile <- paste("~/DispersalEvolution/RangeExpansion", SimID, "SummaryStats.csv", sep = "/")
               CurData <- read.csv(CurFile)
               CurData <- subset(CurData, gen == Gens[t])
               for(x in 1:N){
                    CurLoc <- which(CurData$x == CurXseq[x])
                    if(length(CurLoc) != 0){
                         All_dBars[x,j] <- CurData$dBar[CurLoc]
                         if(CurData$abund[CurLoc] >= 10){
                              AllGenVar[x,j] <- CurData$GenVar[CurLoc]
                         }
                    }
               }
          }
          rm(CurData)
          
          # Now, calculate the summary statistics across simulations for the
          #    current generation
          for(x in 1:N){
               DispPhens[[t]]$dBar[x] <- mean(All_dBars[x,], na.rm = TRUE)
               dBarIQR <- quantile(All_dBars[x,], probs = c(0.25, 0.75), na.rm = TRUE)
               DispPhens[[t]]$lwr[x] <- dBarIQR[1]
               DispPhens[[t]]$upr[x] <- dBarIQR[2]
               
               DispGenVar[[t]]$GenVar[x] <- mean(AllGenVar[x,], na.rm = TRUE)
               GenVarIQR <- quantile(AllGenVar[x,], probs = c(0.25, 0.75), na.rm = TRUE)
               DispGenVar[[t]]$lwr[x] <- GenVarIQR[1]
               DispGenVar[[t]]$upr[x] <- GenVarIQR[2]
          }
     }
     
     # Put everything together to return
     DispSummary <- list(dBar = DispPhens, GenVar = DispGenVar)
     return(DispSummary)
}

# Create the cluster and run the extractions
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("SimData", "ParamCombos", "Gens", "Ngens") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/DispersalEvolution/SimFunctions.R") )

# Run the extractions
ExtractVec <- 1:length(ParamCombos)
DispData <- clusterApply(cl, x = ExtractVec, fun = DispExtract)

# Save the list returned by clusterApply for graphing
save(DispData, ParamCombos, SimData, file = OutFile)

