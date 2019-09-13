# This script will extract data from the range shift simulations. Specifically,
#    for each parameter combination, it will extract a data frame with a row for
#    each time point and columns for (1) time, (2) proportion of simulations 
#    extinct, (3) among-simulation mean change in population-wide mean dispersal
#    phenotype from equilibrium, (4 & 5) interquartile range in mean dispersal
#    phenotype changes.

# Set the in and outfiles for the simulations and data
InFile <- "2019-06-26_DiploidDioSims.csv"
OutFile <- "2019-06-26_DiploidDioShift.rdata"

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
Gens <- seq(0, LengthShift, by = LengthShift/10)
Ngens <- length(Gens)

# For each parameter combination, I want a data frame with a row for each time
#    point and columns for (1) time, (2) proportion of simulations extinct,
#    (3) among-simulation mean change in population-wide mean dispersal
#    phenotype from equilibrium, (4 & 5) interquartile range in mean dispersal
#    phenotype changes, and then 6:8 will be the same for genetic variance.
ShiftExtract <- function(i){
     # First read in the parameters for a representative simulation from this
     #    parameter combination
     CurSims <- SimData[ParamCombos[[i]],]
     Nsims <- nrow(CurSims)
     ParamFile <- paste(as.character(CurSims$ID[1]), "parameters.R", sep = "/")
     source(ParamFile)
     
     # Make an empty data frame of the correct dimensions
     Results <- data.frame(gen = Gens, ext = rep(NA, Ngens),
                           Delta_dBar = rep(NA, Ngens), 
                           Delta_dBar_lwr = rep(NA, Ngens),
                           Delta_dBar_upr = rep(NA, Ngens),
                           DeltaGenVar = rep(NA, Ngens),
                           DeltaGenVar_lwr = rep(NA, Ngens),
                           DeltaGenVar_upr = rep(NA, Ngens))
     
     # Loop through each time point, calculating each quantity for each time point
     #    and then deleting the data objects so as to not overload memory
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     for(t in 1:Ngens){
          extinct <- rep(NA, Nsims)
          Delta_dbar <- rep(NA, Nsims)
          DeltaGenVar <- rep(NA, Nsims)
          for(j in 1:Nsims){
               SimID <- strsplit(x = as.character(CurSims$ID[j]), split = "/")[[1]][4]
               ShiftData <- read.csv(paste("~/DispersalEvolution/RangeShift", 
                                         SimID, "SummaryStats.csv", sep = "/"))
               
               # Determine if the population is extant and if so calculate the current phenotypes
               extinct[j] <- ifelse(max(ShiftData$gen) >= Gens[t], 0, 1)
               if(extinct[j] == 0){
                    CurData <- subset(ShiftData, gen == Gens[t])
                    EquilibriumData <- subset(ShiftData, gen == 0)
                    Delta_dbar[j] <- mean(CurData$dBar) - mean(EquilibriumData$dBar)
                    DeltaGenVar[j] <- mean(CurData$GenVar) - mean(EquilibriumData$GenVar)
               }
          }
          Results$ext[t] <- sum(extinct)/Nsims
          Results$Delta_dBar[t] <- mean(Delta_dbar, na.rm = TRUE)
          Results$Delta_dBar_lwr[t] <- quantile(Delta_dbar, probs = 0.25, na.rm = TRUE)
          Results$Delta_dBar_upr[t] <- quantile(Delta_dbar, probs = 0.75, na.rm = TRUE)
          Results$DeltaGenVar[t] <- mean(DeltaGenVar, na.rm = TRUE)
          Results$DeltaGenVar_lwr[t] <- quantile(DeltaGenVar, probs = 0.25, na.rm = TRUE)
          Results$DeltaGenVar_upr[t] <- quantile(DeltaGenVar, probs = 0.75, na.rm = TRUE)
     }
     return(Results)
}

# Create the cluster and run the extractions
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("SimData", "ParamCombos", "Gens", "Ngens") )

# Change the working directory of the worker nodes
temp <- clusterEvalQ(cl, source("~/DispersalEvolution/SimFunctions.R") )

# Run the extractions
ExtractVec <- 1:length(ParamCombos)
ShiftData <- clusterApply(cl, x = ExtractVec, fun = ShiftExtract)

# Save the list returned by clusterApply for graphing
save(ShiftData, ParamCombos, SimData, file = OutFile)

