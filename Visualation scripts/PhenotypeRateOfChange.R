# This script will provide a visualization of the function transforming dispersal
#    "genotypes" (sum of allele values) to phenotypes. This function is bounded at
#    0 since dispersal distances by definition must be strictly positive and at
#    a maximum achievable dispersal phenotype (d_max). The maximum value has two
#    functions: (1) it mimics biologically realistic scenarios where physiological
#    or physical constraints will produce an upper bound on dispersal (e.g. wind
#    dispersed seeds can only become so light and still be viable) and (2) it
#    allows the transition between the minimum and maximum phenotypes to more
#    closely approximate a linear transition, meaning an allele has similar effects
#    within the range of moderate dispersal phenotypes. The figure will demonstrate
#    these features.

setwd("~/Desktop/GitHubRepos/DispersalEvolution/")
library(RColorBrewer)

# First, create a simplified version of the DispPhen function used in the model
DispPhen <- function(genotype, dmax, rho, lambda){
     # Use the genotype to calculate the expected dispersal distance (phenotype)
     d <- (dmax * exp(rho * (genotype - lambda))) / (1 + exp(rho * (genotype - lambda)))
     return(d)
}

genotypes <- seq(-30, 50, length.out = 1000)
dmax <- 6
rho <- 0.1
lambda <- 10

AltRho <- 1
AltLambda <- 0

NoRhoPhens <- DispPhen(genotype = genotypes, dmax = dmax, rho = AltRho, lambda = AltLambda)
NoRhoLambdaPhens <- DispPhen(genotype = genotypes, dmax = dmax, rho = AltRho, lambda = lambda)
PaperPhens  <- DispPhen(genotype = genotypes, dmax = dmax, rho = rho, lambda = lambda)

PairedCols <- brewer.pal(n = 12, name = "Paired")


# Now make the visualization
pdf(file = "ResultFigures/FigureS1_new.pdf", width = 5, height = 4, onefile = FALSE, paper = "special", useDingbats = FALSE)
     plot(x = genotypes, y = NoRhoPhens, type = "l", lwd = 1.5, ylab = "Phenotype", xlab = "Genotype",
          las = 1, main = "", col = PairedCols[1])
     lines(x = genotypes, y = NoRhoLambdaPhens, lwd = 1.5, col = PairedCols[2])
     lines(x = genotypes, y = PaperPhens, lwd = 1.5, col = PairedCols[3])
     arrows(x0 = 0, y0 = dmax/2, x1 = lambda, y1 = dmax/2, lty = 2, code = 3, angle = 90, length = 0.1)
     text(x = lambda/2, y = dmax/2+0.5, labels = expression(lambda), col = PairedCols[2])
     arrows(x0 = 0.5*lambda, y0 = dmax*(rho*lambda/8 + 0.5 - rho*lambda/4), 
              x1 = 1.5*lambda, y1 = dmax*(3*rho*lambda/8 + 0.5 - rho*lambda/4), 
            lty = 2, col = "black", code = 3, length = 0.15)
     text(x = lambda/2 + 10, y = dmax/2, labels = expression(rho), col = PairedCols[3])
     abline(h = dmax, lty = 2)
     text(x = -1*lambda, y = 5.5, labels = expression(hat(d)))
dev.off()




