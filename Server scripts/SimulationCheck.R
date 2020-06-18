# This script will update the AllSimulation.csv file after each partial run to
#    mark down those simulations that are complete and add the correct path to the
#    file

AllSimulations <- read.csv("/project/rangeecoevomodels/cweissle/DispEv/AllSimulations.csv")
SimsToRemove <- NULL

# Check for completeness of all simulations
for(i in 1:nrow(AllSimulations)){
     CurDirectory <- paste("/project/rangeecoevomodels/cweissle/DispEv/Sims", AllSimulations$SimID[i], sep = "/")
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

# Now remove the incomplete simulations
#if(!is.null(SimsToRemove){
#     for(i in 1:length(SimsToRemove)){
#          SysCommand <- paste("rm -r ", SimsToRemove[i], sep = "")
#          system(SysCommand)
#     }
#}

# How many simulations are complete?
print(sum(AllSimulations$complete))

write.csv(AllSimulations, file = "/project/rangeecoevomodels/cweissle/DispEv/AllSimulations.csv", 
          row.names = FALSE, quote = FALSE)



