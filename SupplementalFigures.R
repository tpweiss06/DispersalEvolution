# This script will make the necessary figures for the supplementary material
#    for the manuscript.

setwd("~/Desktop/GitHubRepos/DispersalEvolution/")

# Set up the same colors and labels as the main manuscript figures
library(RColorBrewer)
ScenCols <- brewer.pal(n = 5, name = "Dark2")
Scenarios <- c("Asexual", "Obligate selfing", "Partial selfing", "Obligate outcrossing", "Sexual (dioecious)")
HapVals <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
monoVals <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
omegaVals <- c(0, 1, 0.5, 0, 0)
Lseq <-c(1,2,4,8,16,32)

############################ Fig. S1 + S2; Initial conditions of the simulations
# a 6 x 2 figure with each row corresponding to a different number of loci, the
#    first column corresponding to the mean phenotype and the second column
#    corresponding to the mean genetic variance. The x-axis of all graphs will
#    correspond to space.
load("InitialConditions.rdata")
xlab <- "Space"
ylab1 <- "Mean phenotype"
ylab2 <- "Mean genetic variance"
PhenRange <- range(MeanPhenVals, na.rm = TRUE)
GenRange <- range(MeanGenVals, na.rm = TRUE)
xSeq <- -59:59
xRange <- c(-60,60)
#Legend_xCoords <- c(-250, -200, -125, -60, 30) - 40
Legend_xCoords <- c(-290, -240, -165, -100, -10)

pdf(file = "ResultFigures/FigureS1.pdf", width = 6, height = 9, onefile = FALSE, paper = "special", useDingbats = FALSE)
     par(mfrow = c(6,2), oma = c(3,3,3,0), mar = c(1.5, 4, 1.25, 1.5))
     for(l in 1:6){
          plot(x = NA, y = NA, xlim = xRange, ylim = PhenRange, main = "", xlab = "",
               ylab = ylab1, las = 1)
          for(p in 1:5){
               lines(x = xSeq, y = MeanPhenVals[l,p,], col = ScenCols[p])
          }
          plot(x = NA, y = NA, xlim = xRange, ylim = GenRange, main = "", xlab = "",
               ylab = ylab2, las = 1)
          for(p in 1:5){
               lines(x = xSeq, y = MeanGenVals[l,p,], col = ScenCols[p])
          }
          if(l == 1){
                  for(i in 1:5){
                          legend(x = Legend_xCoords[i], y = 1.15*GenRange[2], horiz = TRUE, bty = "n", 
                                 legend = Scenarios[i], lty = 1, col = ScenCols[i], xpd = NA)
                  }
          }
     }
     mtext(xlab, side = 1, line = 1, outer = TRUE)
dev.off()


# This figure for variance among simulations might not make it into the supplemental
#    material and instead I might just include a description or some values in the 
#    caption for Figure S1.
ylab1 <- "Phenotype"
ylab2 <- "Genetic variance"
PhenRange <- range(VarPhenVals, na.rm = TRUE)
GenRange <- range(VarGenVals, na.rm = TRUE)
pdf(file = "ResultFigures/FigureS2.pdf", width = 6, height = 9, onefile = FALSE, paper = "special", useDingbats = FALSE)
     par(mfrow = c(6,2), oma = c(3,3,3,0), mar = c(1.5, 4, 1.25, 1.5))
     for(l in 1:6){
          plot(x = NA, y = NA, xlim = xRange, ylim = PhenRange, main = "", xlab = "",
               ylab = "", las = 1)
          mtext(ylab1, side = 2, line = 3.5, cex = 0.75)
          for(p in 1:5){
               lines(x = xSeq, y = VarPhenVals[l,p,], col = ScenCols[p])
          }
          plot(x = NA, y = NA, xlim = xRange, ylim = GenRange, main = "", xlab = "",
               ylab = ylab2, las = 1)
          for(p in 1:5){
               lines(x = xSeq, y = VarGenVals[l,p,], col = ScenCols[p])
          }
          if(l == 1){
                  for(i in 1:5){
                          legend(x = Legend_xCoords[i], y = 1.2*GenRange[2], horiz = TRUE, bty = "n", 
                                 legend = Scenarios[i], lty = 1, col = ScenCols[i], xpd = NA)
                  }
          }
     }
     mtext(xlab, side = 1, line = 1, outer = TRUE)
dev.off()

############### Fig. S3; Sensitivity analysis of mutation parameters (Expansion)
library(plot3D)
# columns are y-axis, rows are x-axis

# Use the manuscriptFigures.R script to process the data (same idea, but IQR widths)
# Make three lists with 5 entries (each a matrix for graphing)
