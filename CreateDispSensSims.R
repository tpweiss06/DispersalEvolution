# This script will create a .csv file with a row for each simulation that
#    needs to be done and columns for the file path, key paramter values,
#    and a column for a completeness indicator (0 = incomplete; 1 = complete).
#    These simulations are specific to the dispersal sensitivity analysis reported in the
#    manuscript, which will be run with the DispersalSensitivity.R script. The
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
rho_star <- 0.1
dmax_star <- 6
rhoSeq <- c(0.25*rho_star, 0.5*rho_star, 0.75*rho_star, rho_star, 1.25*rho_star, 2*rho_star, 4*rho_star)
dmaxSeq <- c(0.25*dmax_star, 0.5*dmax_star, 0.75*dmax_star, dmax_star, 1.25*dmax_star, 2*dmax_star, 4*dmax_star)

# Calculate the number of parameter combinations to use
ParamCombos <- length(rhoSeq) * length(dmaxSeq) * length(Haploid)
NumSims <- 100
TotalSims <- ParamCombos * NumSims

# Create an population the master data frame with the information for each simulation
AllSimulations <- data.frame(SimID = 1:TotalSims, L = rep(L, TotalSims),
                             Haploid = rep(NA, TotalSims), monoecious = rep(NA, TotalSims),
                             omega = rep(NA, TotalSims), complete = rep(0, TotalSims),
                             rho = rep(NA, TotalSims), dmax = rep(NA, TotalSims))
k <- 1
for(u in 1:length(Useq)){
     for(s in 1:length(sigmaSeq)){
          for(i in 1:length(Haploid)){
               AllSimulations$rho[k:(k+NumSims-1)] <- rhoSeq[u]
               AllSimulations$dmax[k:(k+NumSims-1)] <- dmaxSeq[s]
               AllSimulations$Haploid[k:(k+NumSims-1)] <- Haploid[i]
               AllSimulations$monoecious[k:(k+NumSims-1)] <- monoecious[i]
               AllSimulations$omega[k:(k+NumSims-1)] <- omega[i]
               k <- k + NumSims
          }
     }
}

write.csv(AllSimulations, file = "~/Desktop/Wyoming/DispersalEvolution/GitRepo/DispSensSimulations.csv", 
          row.names = FALSE, quote = FALSE)


