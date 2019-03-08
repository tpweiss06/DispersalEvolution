# This script will extract data from the range shift simulations. Specifically,
#    for each parameter combination, it will extract a data frame with a row for
#    each time point and columns for (1) time, (2) proportion of simulations 
#    extinct, (3) among-simulation mean change in population-wide mean dispersal
#    phenotype from equilibrium, (4 & 5) interquartile range in mean dispersal
#    phenotype changes.

nProc <- 24*1
NumSims <- 20

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

# For each parameter combination, I want a data frame with a row for each time
#    point and columns for (1) time, (2) proportion of simulations extinct,
#    (3) among-simulation mean change in population-wide mean dispersal
#    phenotype from equilibrium, (4 & 5) interquartile range in mean dispersal
#    phenotype changes.
ShiftExtract <- function(i){
     # First read in the parameters for a representative simulation from this
     #    parameter combination
     CurSims <- SimData[ParamCombos[[i]],]
     ParamFile <- paste(as.character(CurSims$ID[1]), "parameters.R", sep = "/")
     source(ParamFile)
     
     # Make an empty data frame of the correct dimensions
     Results <- data.frame(gen = 1:LengthShift, ext = rep(NA, LengthShift),
                           Delta_dBar = rep(NA, LengthShift), 
                           Delta_dBar_lwr = rep(NA, LengthShift),
                           Delta_dBar_upr = rep(NA, LengthShift))
     
     # Loop through each time point, calculating each quantity for each time point
     #    and then deleting the data objects so as to not overload memory
     NumSims <- nrow(CurSims)
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     for(t in 1:LengthShift){
          extinct <- rep(NA, NumSims)
          Delta_dbar <- rep(NA, NumSims)
          for(j in 1:NumSims){
               SimID <- strsplit(x = as.character(CurSims$ID[j]), split = "/")[[1]][4]
               ShiftData <- read.csv(paste("~/DispersalEvolution/RangeShift", 
                                         SimID, "SummaryStats.csv", sep = "/"))
               EquilibriumData <- read.csv(paste(as.character(CurSims$ID[j]), 
                                                 "EquilibriumPopMat.csv", sep = "/"))
               # Calculate the equilibrium genotypes and phenotypes of the population
               if((L == 1) & (Haploid)){
                    genotypes <- EquilibriumData[,PopIndices$DispCols]
               }else{
                    if(nrow(EquilibriumData) > 1){
                         genotypes <- rowSums(EquilibriumData[,PopIndices$DispCols])
                    } else{
                         genotypes <- sum(EquilibriumData[,PopIndices$DispCols])
                    }
               }
               EquilibriumPhenotypes <- (dmax * exp(rho * genotypes)) / (1 + exp(rho * genotypes))
               
               # Determine if the population is extant and if so calculate the current phenotypes
               extinct[j] <- ifelse(max(ShiftData$gen) >= t, 0, 1)
               if(extinct[j] == 0){
                    CurData <- subset(ShiftData, gen == t)
                    ShiftPhenotypes <- CurData$MuPhen
                    Delta_dbar[j] <- mean(ShiftPhenotypes) - mean(EquilibriumPhenotypes)
               }
          }
          Results$ext[t] <- sum(extinct)/NumSims
          Results$Delta_dBar[t] <- mean(Delta_dbar, na.rm = TRUE)
          Results$Delta_dBar_lwr[t] <- quantile(Delta_dbar, probs = 0.25, na.rm = TRUE)
          Results$Delta_dBar_upr[t] <- quantile(Delta_dbar, probs = 0.75, na.rm = TRUE)
     }
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
ShiftData <- clusterApply(cl, x = ExtractVec, fun = ShiftExtract)

# Save the list returned by clusterApply for graphing
save(ShiftData, ParamCombos, SimData, file = "ShiftData.rdata")

