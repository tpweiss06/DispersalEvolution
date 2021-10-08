# This script will provide a visualization of the derivatives of two possible 
#    functions for linking dispersal genotype to phenotyp. Some type of link
#    is necessary as the genotype can take on negative values, but the 
#    phenotype is strictly positive.

# First define the derivative of the basic link function: d_i = exp(rho*(Genotype-lambda))
BaseDeriv <- function(rho, lambda, genotype){
     return(rho * exp(rho*(genotype - lambda)))
}

# Next, define the derivative of the link function incorporating an upper limit:
#    d_i = (d_hat*exp(rho*(genotype-lambda))) / (1 + exp(rho*(genotype - lambda)))
MaxDerv <- function(rho, lambda, d_hat, genotype){
     numerator <- d_hat * rho * exp(rho * (genotype - lambda))
     denominator <- (1 + exp(rho * (genotype - lambda)))^2
     return(numerator / denominator)
}

Results <- read.csv("SimsWithResults.csv")
InitVals <- read.csv("SimsWithInitVals.csv")

GenotypeVals <- seq(from = -30, to = 50, length.out = 10000)

# Set the parameter values used in the simulations
rho <- 0.1
lambda <- 10
d_hat <- 6


BasicLine <- BaseDeriv(rho = rho, lambda = lambda, genotype = GenotypeVals)
MaxLine <- MaxDerv(rho = rho, lambda = lambda, d_hat = d_hat, genotype = GenotypeVals)

plot(x = GenotypeVals, y = BasicLine, type = "l", xlab = "Genotype", 
     ylab = "Rate of change in phenotype", las = 1, col = "red")
lines(x = GenotypeVals, y = MaxLine, col = "blue", lty = 2)
legend("topleft", legend = c("Without maximum", "With maximum"), lty = c(1,2),
       col = c("red", "blue"), bty = "n")
