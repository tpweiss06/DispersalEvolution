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

# For each parameter combination, I want a matrix with a row for every time point
#    and 7 columns: mean, upper quartile, lower quartile for dbar and GenVar,
#    and a last column with the number of simulations still running at
#    each time point. Additionally, the function will return a vector with the
#    time to evolution of dThresh for each simulation.
ExpandExtract <- function(i){
     # First, read in the different expansion data and record the last time
     #    point of each simulation
     CurSims <- SimData[ParamCombos[[i]],]
     EvolTimes <- rep(NA, nrow(CurSims))
     for(j in 1:nrow(CurSims)){
          CurData <- read.csv(paste(CurSims$ID[j], "SummaryStats.csv", sep = "/"))
          EvolTimes[j] <- nrow(CurData) - 1
     }
     rm(CurData)
     
     # Now make the matrix to hold the values through time
     ExpandVals <- matrix(NA, nrow = max(EvolTimes) + 1, ncol = 7)
     for(t in 1:nrow(ExpandVals)){
          Sim_dBar <- rep(NA, nrow(CurSims))
          Sim_GenVar <- rep(NA, nrow(CurSims))
          for(j in 1:nrow(CurSims)){
               CurData <- read.csv(paste(CurSims$ID[j], "SummaryStats.csv", sep = "/"))
               if(t <= nrow(CurData)){
                    Sim_dBar[j] <- CurData$dBar[t]
                    Siim_GenVar[j] <- CurData$GenVar[t]
               }
          }
          rm(CurData)
          ExpandVals[t,1] <- mean(Sim_dBar, na.rm = TRUE)
          ExpandVals[t,2:3] <- quantile(Sim_dBar, probs = c(0.25,0.75), na.rm = TRUE)
          ExpandVals[t,4] <- mean(Sim_GenVar, na.rm = TRUE)
          ExpandVals[t,5:6] <- quantile(Sim_GenVar, probs = c(0.25,0.75), na.rm = TRUE)
          ExpandVals[t,7] <- sum(!is.na(Sim_dBar))
     }
     # Now add column names to ExpandVals
     colnames(ExpandVals) <- c("dBar_mu","dBar_lwr", "dBar_upr", "GenVar_mu", 
                               "GenVar_lwr", "GenVar_upr", "N")
     
     # Put it all together and return it
     Results <- list(ExpandVals = ExpandVals, EvolTimes = EvolTimes)
     return(Results)
}

# Create the cluster and run the extractions
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("SimData", "ParamCombos") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/DispersalEvolution/SimFunctions.R") )

# Run the extractions
ExtractVec <- 1:length(ParamCombos)
ExpandData <- clusterApply(cl, x = ExtractVec, fun = ExpandExtract)

# Save the list returned by clusterApply for graphing
save(ExpandData, ParamCombos, SimData, file = "ExpansionData.rdata")

