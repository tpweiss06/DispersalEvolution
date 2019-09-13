# This script will remove all simulations corresponding to the simulation IDs in
#    the input csv file from all subfolders in the DispersalEvolution directory

# First, read in the csv file containing the IDs to be removed
InFile <- "Trial2_Sims.csv"
SimData <- read.csv(InFile)

# Now loop through all simulation IDs in the file and delete them
for(i in 1:nrow(SimData)){
     # Extract just the Simulatin ID so that it can be used on multiple subfolders
     SimID <- strsplit(x = as.character(CurSims$ID[i]), split = "/")[[1]][4]
     # Next, create the file paths to the simulation in all subfolders
     EquilibriumPath <- paste("~/DispersalEvolution/RangeEquilibrium", SimID, sep = "/")
     ExpansionPath <- paste("~/DispersalEvolution/RangeExpansion", SimID, sep = "/")
     ShiftPath <- paste("~/DispersalEvolution/RangeShift", SimID, sep = "/")
     # Finally, check if the simulation exists in each subfolder and delete it
     if(dir.exists(EquilibriumPath)){
          SysCommand <- paste("rm -r ", EquilibriumPath, sep = "")
          system(SysCommand)
     }
     if(dir.exists(ExpansionPath)){
          SysCommand <- paste("rm -r ", ExpansionPath, sep = "")
          system(SysCommand)
     }
     if(dir.exists(ShiftPath)){
          SysCommand <- paste("rm -r ", ShiftPath, sep = "")
          system(SysCommand)
     }
}

# Finally, delete the csv file that now contains useless simulation IDs
SysCommand <- paste("rm ", InFile, sep = "")
system(SysCommand)
