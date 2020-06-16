# This script will create a .csv file with a row for each simulation that
#    needs to be done and columns for the file path, key paramter values,
#    and a column for a completeness indicator (0 = incomplete; 1 = complete).
#    This .csv file will then be used in the RunSimulations.R script to run
#    the simulations that are still incomplete, while saving those that have
#    finished even if a previous run fails.
# Also, this script will delete extra runs for certain parameter values that
#    may have been repeated at this point.

# First make the master simulation data frame, using the parameter combinations
#    to be explored
Lseq <- c(1, 2, 4, 8, 16, 32)
# NOTE: Haploid is always associated with a monoecious value of TRUE even though
#    it has no meaning for haploid organisms.
Haploid <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
monoecious <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
omega <- c(0, 0, 0, 0.5, 1)

# Calculate the number of parameter combinations to use
ParamCombos <- length(Lseq) * length(Haploid)
NumSims <- 1000
TotalSims <- ParamCombos * NumSims

# Create an population the master data frame with the information for each simulation
AllSimulations <- data.frame(SimID = 1:TotalSims, L = rep(NA, TotalSims),
                             Haploid = rep(NA, TotalSims), monoecious = rep(NA, TotalSims),
                             omega = rep(NA, TotalSims), complete = rep(0, TotalSims))
k <- 1
for(l in 1:length(Lseq)){
     for(i in 1:length(Haploid)){
          AllSimulations$L[k:(k+NumSims-1)] <- Lseq[l]
          AllSimulations$Haploid[k:(k+NumSims-1)] <- Haploid[i]
          AllSimulations$monoecious[k:(k+NumSims-1)] <- monoecious[i]
          AllSimulations$omega[k:(k+NumSims-1)] <- omega[i]
          k <- k + 1000
     }
}

# check that the above worked
for(l in 1:length(Lseq)){
     for(i in 1:length(Haploid)){
          CurRows <- nrow(subset(AllSimulations, (L == Lseq[l]) & (Haploid == Haploid[i]) & (monoecious == monoecious[i]) & (omega == omega[i])))
          print(CurRows)
     }
}
sum(is.na(AllSimulations$L))
sum(is.na(AllSimulations$Haploid))
sum(is.na(AllSimulations$monoecious))
sum(is.na(AllSimulations$omega))

write.csv(AllSimulations, file = "~/Desktop/PostdocResearch/DispersalEvolution/GitRepo/AllSimulations.csv", 
          row.names = FALSE, quote = FALSE)


