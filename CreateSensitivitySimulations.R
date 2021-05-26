# This script will create a .csv file with a row for each simulation that
#    needs to be done and columns for the file path, key paramter values,
#    and a column for a completeness indicator (0 = incomplete; 1 = complete).
#    These simulations are specific to the sensitivity analysis reported in the
#    manuscript, which will be run with the RunSensitivitySims.R script. The
#    structure of the code follows the same logical flow as the simulations used
#    for the main results reported in the manuscript.

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
U_star <- 0.02
sigma_star <- sqrt(0.02)
Useq <- c(0.25*U_star, 0.5*U_star, 0.75*U_star, U_star, 1.25*U_star, 2*U_star, 4*U_star)
sigmaSeq <- c(0.25*sigma_star, 0.5*sigma_star, 0.75*sigma_star, sigma_star, 1.25*sigma_star, 2*sigma_star, 4*sigma_star)

# Calculate the number of parameter combinations to use
ParamCombos <- length(Useq) * length(sigmaSeq) * length(Haploid)
NumSims <- 100
TotalSims <- ParamCombos * NumSims

# Create an population the master data frame with the information for each simulation
AllSimulations <- data.frame(SimID = 1:TotalSims, L = rep(L, TotalSims),
                             Haploid = rep(NA, TotalSims), monoecious = rep(NA, TotalSims),
                             omega = rep(NA, TotalSims), complete = rep(0, TotalSims),
                             U = rep(NA, TotalSims), sigma = rep(NA, TotalSims))
k <- 1
for(u in 1:length(Useq)){
     for(s in 1:length(sigmaSeq)){
          for(i in 1:length(Haploid)){
               AllSimulations$U[k:(k+NumSims-1)] <- Useq[u]
               AllSimulations$sigma[k:(k+NumSims-1)] <- sigmaSeq[s]
               AllSimulations$Haploid[k:(k+NumSims-1)] <- Haploid[i]
               AllSimulations$monoecious[k:(k+NumSims-1)] <- monoecious[i]
               AllSimulations$omega[k:(k+NumSims-1)] <- omega[i]
               k <- k + NumSims
          }
     }
}

write.csv(AllSimulations, file = "~/Desktop/Wyoming/DispersalEvolution/GitRepo/SensSimulations.csv", 
          row.names = FALSE, quote = FALSE)


