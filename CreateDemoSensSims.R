# This script will create a .csv file with a row for each simulation that
#    needs to be done and columns for the file path, key paramter values,
#    and a column for a completeness indicator (0 = incomplete; 1 = complete).
#    These simulations are specific to the demographic sensitivity analysis reported in the
#    manuscript, which will be run with the DemographicSensitivity.R script. The
#    structure of the code follows the same logical flow as the simulations used
#    for the main results reported in the manuscript.

library(here)

# First make the master simulation data frame, using the parameter combinations
#    to be explored
L <- 4
# NOTE: Haploid is always associated with a monoecious value of TRUE even though
#    it has no meaning for haploid organisms.
Haploid <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
monoecious <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
omega <- c(0, 0, 0, 0.5, 1)

# Make vectors of the different U and sigma values to use for the sensitivity
#    analysis
R_star <- 2
Kmax_star <- 100
Rseq <- c(1, 1.25, 1.75, 2, 2.5, 3, 4) #The minimum here needs to be at least replacement
KmaxSeq <- c(0.25*Kmax_star, 0.5*Kmax_star, 0.75*Kmax_star, Kmax_star, 1.25*Kmax_star, 2*Kmax_star, 4*Kmax_star)

# Calculate the number of parameter combinations to use
ParamCombos <- length(Rseq) * length(KmaxSeq) * length(Haploid)
NumSims <- 100
TotalSims <- ParamCombos * NumSims

# Create an population the master data frame with the information for each simulation
AllSimulations <- data.frame(SimID = 1:TotalSims, L = rep(L, TotalSims),
                             Haploid = rep(NA, TotalSims), monoecious = rep(NA, TotalSims),
                             omega = rep(NA, TotalSims), complete = rep(0, TotalSims),
                             R = rep(NA, TotalSims), Kmax = rep(NA, TotalSims))
k <- 1
for(u in 1:length(Rseq)){
     for(s in 1:length(KmaxSeq)){
          for(i in 1:length(Haploid)){
               AllSimulations$R[k:(k+NumSims-1)] <- Rseq[u]
               AllSimulations$Kmax[k:(k+NumSims-1)] <- KmaxSeq[s]
               AllSimulations$Haploid[k:(k+NumSims-1)] <- Haploid[i]
               AllSimulations$monoecious[k:(k+NumSims-1)] <- monoecious[i]
               AllSimulations$omega[k:(k+NumSims-1)] <- omega[i]
               k <- k + NumSims
          }
     }
}

write.csv(AllSimulations, file = here("DemoSensSimulations.csv"), 
          row.names = FALSE, quote = FALSE)


