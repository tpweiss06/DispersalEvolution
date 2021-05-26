# This script will update the AllSimulation.csv file after each partial run to
#    mark down those simulations that are complete and add the correct path to the
#    file

AllSimulations <- read.csv("/project/rangeecoevomodels/cweissle/DispEv/SensSimulations.csv")
SimsToRemove <- NULL

# Check for completeness of all simulations
for(i in 1:nrow(AllSimulations)){
     CurDirectory <- paste("/project/rangeecoevomodels/cweissle/DispEv/SensSims", AllSimulations$SimID[i], sep = "/")
     FinalFile <- paste(CurDirectory, "ShiftPopMat.csv", sep = "/")
     DirectoryCheck <- dir.exists(CurDirectory)
     CompleteCheck <- file.exists(FinalFile)
     if(DirectoryCheck & CompleteCheck){
          AllSimulations$complete[i] <- 1
     } else if(DirectoryCheck){
          SimsToRemove <- c(SimsToRemove, CurDirectory)
     }
}

print(SimsToRemove)

# How many simulations are complete?
print(sum(AllSimulations$complete))

write.csv(AllSimulations, file = "/project/rangeecoevomodels/cweissle/DispEv/SensSimulations.csv", 
          row.names = FALSE, quote = FALSE)



