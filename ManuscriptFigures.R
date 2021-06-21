# This script will make the figures to be used in the manuscript. In particular
#    the script will make the following figures:
# Figure 1: A schematic figure demonstrating the distinction between range
#    expansions and range shifts as we will use them in the paper
# Figure 2: Dispersal evolution in range expansions and shifts (Fig. S1). This will be a
#    two panel figure with the first focusing on dispersal phenotypes and the
#    second focusing on genetic variance. In both panels, the x axis will be 
#    the number of loci in the simulation and different mating systems will be
#    designated by line color (asexual, sexual (dioecious), sexual (obligate selfing),
#    sexual (partial selfing), and sexual (obligate outcrossing)). In the first
#    panel, the y-axis will be the change in the mean dispersal phenotype at the
#    range edge. In the second panel, the y-axis will be the change in genetic variance
#    of the range edge populations.
# Figure 3: Distance traveled in the range expansions. This will be another two
#    panel figure with the first being the mean distance travelled and the second
#    being the CV of distance traveled. The line colors and x-axis will match the 
#    scenarios in Figure 2.
# Figure 4: Extinction risk in range shifts. Another two panel figure with the
#    first panel being a line graph of extinction risk (y-axis) versus number
#    of loci (x-axis). The second panel will be a bar chart of the extinction
#    risk observed in no evolution scenarios. Colors will be consistent with 
#    previous figures
# Figure 5: A single panel showing the change in extinction risk due to evolution
#    across loci.


# Set the working directory and load in the result data
setwd("~/Desktop/GitHubRepos/DispersalEvolution/")
Results <- read.csv("SimsWithResults.csv")

# Set an appropriate color scheme and the scenario names
library(RColorBrewer)
ScenCols <- brewer.pal(n = 5, name = "Dark2")
Scenarios <- c("Asexual", "Obligate selfing", "Partial selfing", "Obligate outcrossing", "Sexual (dioecious)")
#Scenarios <- c("Asexual", "Sexual (dioecious)", "Obligate selfing", "Partial selfing", "Obligate outcrossing")
HapVals <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
#monoVals <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
#omegaVals <- c(0, 0, 1, 0.5, 0)
monoVals <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
omegaVals <- c(0, 1, 0.5, 0, 0)
Lseq <- unique(Results$L)
LseqLocs <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
offset <- 0.025
xLocs <- matrix(NA, nrow = 5, ncol = 6)
xLocs[1,] <- LseqLocs - 2*offset
xLocs[2,] <- LseqLocs - offset
xLocs[3,] <- LseqLocs
xLocs[4,] <- LseqLocs + offset
xLocs[5,] <- LseqLocs + 2*offset
###### Figure 1

###### Figure 2
# First get the data together (mean change in phenotpyes and variance with IQR)
Fig2Data <- vector(mode = "list", length = 5)
for(i in 1:5){
     MeanPhenExp <- rep(NA, 6)
     LwrPhenExp <- rep(NA, 6)
     UprPhenExp <- rep(NA, 6)
     MeanGenExp <- rep(NA, 6)
     LwrGenExp <- rep(NA, 6)
     UprGenExp <- rep(NA, 6)
     MeanPhenShift <- rep(NA, 6)
     LwrPhenShift <- rep(NA, 6)
     UprPhenShift <- rep(NA, 6)
     MeanGenShift <- rep(NA, 6)
     LwrGenShift <- rep(NA, 6)
     UprGenShift <- rep(NA, 6)
     for(l in 1:6){
          CurData <- subset(Results, L == Lseq[l] & Haploid == HapVals[i] & monoecious == monoVals[i] & omega == omegaVals[i])
          MeanPhenExp[l] <- mean(c(CurData$DeltaPhenExp_1, CurData$DeltaPhenExp_2))
          LwrPhenExp[l] <- quantile(c(CurData$DeltaPhenExp_1, CurData$DeltaPhenExp_2), probs = 0.25)
          UprPhenExp[l] <- quantile(c(CurData$DeltaPhenExp_1, CurData$DeltaPhenExp_2), probs = 0.75)
          MeanGenExp[l] <- mean(c(CurData$DeltaGenExp_1, CurData$DeltaGenExp_2))
          LwrGenExp[l] <- quantile(c(CurData$DeltaGenExp_1, CurData$DeltaGenExp_2), probs = 0.25)
          UprGenExp[l] <- quantile(c(CurData$DeltaGenExp_1, CurData$DeltaGenExp_2), probs = 0.75)
          MeanPhenShift[l] <- mean(CurData$DeltaPhenShift, na.rm = TRUE)
          LwrPhenShift[l] <- quantile(CurData$DeltaPhenShift, probs = 0.25, na.rm = TRUE)
          UprPhenShift[l] <- quantile(CurData$DeltaPhenShift, probs = 0.75, na.rm = TRUE)
          MeanGenShift[l] <- mean(CurData$DeltaGenShift, na.rm = TRUE)
          LwrGenShift[l] <- quantile(CurData$DeltaGenShift, probs = 0.25, na.rm = TRUE)
          UprGenShift[l] <- quantile(CurData$DeltaGenShift, probs = 0.75, na.rm = TRUE)
     }
     Fig2Data[[i]] <- data.frame(L = Lseq, MeanPhenExp, LwrPhenExp, UprPhenExp, MeanGenExp,
                                 LwrGenExp, UprGenExp, MeanPhenShift, LwrPhenShift, UprPhenShift,
                                 MeanGenShift, LwrGenShift, UprGenShift)
}

# Now make the main figure
pdf(file = "ResultFigures/Figure1.pdf", width = 5, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,3,1,1), mar = c(1.5,1,1.5,1), family = "serif")
        
        # Make the phenotype figure
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0,1.35), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0, 1.4, by = 0.05), tcl = -0.25, labels = FALSE)
        #mtext(expression(paste(Delta, bar("d"), sep = "")), side = 2, line = 3)
        mtext("Change in dispersal phenotype", side = 2, line = 3)
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig2Data[[i]]$MeanPhenExp, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig2Data[[i]]$LwrPhenExp, x1 = xLocs[i,], y1 = Fig2Data[[i]]$UprPhenExp,
                         col = ScenCols[i], lty = 1)
        }
        mtext("a", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
        
        # Put the legend on the figure
        legend("topleft", legend = Scenarios, pch = 21:25, col = ScenCols, bty = "n", pt.bg = ScenCols, inset = -0.01)
        
        # Make the genetic diversity figure
        plot(NA, NA, xlim = c(-0.05,1.05), ylim = c(-3, 0), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(-3, 0, by = 0.1), tcl = -0.25, labels = FALSE)
        mtext("Change in genetic variance", side = 2, line = 3)
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig2Data[[i]]$MeanGenExp, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig2Data[[i]]$LwrGenExp, x1 = xLocs[i,], y1 = Fig2Data[[i]]$UprGenExp,
                         col = ScenCols[i], lty = 1)
        }
        mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
        mtext("b", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
dev.off()

# Now make the figure with the range shift results
pdf(file = "ResultFigures/Figure2.pdf", width = 5, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,3,1,1), mar = c(1.5,1,1.5,1), family = "serif")

        # Make the phenotype figure
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0,1.65), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0, 1.65, by = 0.05), tcl = -0.25, labels = FALSE)
        #mtext(expression(paste(Delta, bar("d"), sep = "")), side = 2, line = 3)
        mtext("Change in dispersal phenotype", side = 2, line = 3)
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig2Data[[i]]$MeanPhenShift, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig2Data[[i]]$LwrPhenShift, x1 = xLocs[i,], y1 = Fig2Data[[i]]$UprPhenShift,
                         col = ScenCols[i], lty = 1)
        }
        mtext("a", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)

        # Put the legend on the figure
        legend("topleft", legend = Scenarios, pch = 21:25, col = ScenCols, bty = "n", pt.bg = ScenCols, inset = -0.01, ncol = 2)

        # Make the genetic diversity figure
        plot(NA, NA, xlim = c(-0.05,1.05), ylim = c(-6, 0), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(-6, 0, by = 0.25), tcl = -0.25, labels = FALSE)
        mtext("Change in genetic variance", side = 2, line = 3)
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig2Data[[i]]$MeanGenShift, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig2Data[[i]]$LwrGenShift, x1 = xLocs[i,], y1 = Fig2Data[[i]]$UprGenShift,
                         col = ScenCols[i], lty = 1)
        }
        mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
        mtext("b", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
dev.off()

###### Figure 3
# distance with IQR for evolution and no evolution scenarios
# First get the data together (mean distance spread and CV for each scenario)
Fig3Data <- vector(mode = "list", length = 5)
for(i in 1:5){
        DistEvol <- rep(NA, 6)
        LwrDistEvol <- rep(NA, 6)
        UprDistEvol <- rep(NA, 6)
        DistNoEvol <- rep(NA, 6)
        LwrDistNoEvol <- rep(NA, 6)
        UprDistNoEvol <- rep(NA, 6)
        Difference <- rep(NA, 6)
        LwrDiff <- rep(NA, 6)
        UprDiff <- rep(NA, 6)
        for(l in 1:6){
                CurData <- subset(Results, L == Lseq[l] & Haploid == HapVals[i] & monoecious == monoVals[i] & omega == omegaVals[i])
                EvolDistances <- c(CurData$distance_1, CurData$distance_2)
                NoEvolDistances <- c(CurData$DistNoEvol_1, CurData$DistNoEvol_2)
                Diffs <- EvolDistances - NoEvolDistances
                DistEvol[l] <- mean(EvolDistances)
                LwrDistEvol[l] <- quantile(EvolDistances, probs = 0.25)
                UprDistEvol[l] <-  quantile(EvolDistances, probs = 0.75)
                DistNoEvol[l] <- mean(NoEvolDistances)
                LwrDistNoEvol[l] <- quantile(NoEvolDistances, probs = 0.25)
                UprDistNoEvol[l] <- quantile(NoEvolDistances, probs = 0.75)
                Difference[l] <- mean(Diffs)
                LwrDiff[l] <- quantile(Diffs, probs = 0.25)
                UprDiff[l] <- quantile(Diffs, probs = 0.75)
        }
        Fig3Data[[i]] <- data.frame(L = Lseq, DistEvol, LwrDistEvol, UprDistEvol, 
                                    DistNoEvol, LwrDistNoEvol, UprDistNoEvol,
                                    Difference, LwrDiff, UprDiff)
}
# Now make a figure for the mean distance spread with evolution
pdf(file = "ResultFigures/Figure3_EvolDist.pdf", width = 5, height = 3, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(family = "serif", mar = c(3, 4, 1, 1))
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(400, 800), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(400, 800, by = 25), tcl = -0.25, labels = FALSE)
        mtext("Distance spread", side = 2, line = 2.75)
        mtext("Number of loci", side = 1, line = 2)
        # Graph the data
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig3Data[[i]]$DistEvol, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig3Data[[i]]$LwrDistEvol, x1 = xLocs[i,], y1 = Fig3Data[[i]]$UprDistEvol,
                         lty = 1, col = ScenCols[i])
        }
        # Put the legend on the figure
        legend(x = -0.1, y = 750, legend = Scenarios[1:3], pch = 21:23, col = ScenCols[1:3], 
               pt.bg = ScenCols[1:3], bty = "n", yjust = 0.5)
        legend(x = 0.305, y = 710, legend = Scenarios[4:5], pch = 24:25, col = ScenCols[4:5],
               pt.bg = ScenCols[4:5], bty = "n", yjust = 0)
dev.off()

# Now make a figure for the mean distance spread without evolution
pdf(file = "ResultFigures/Figure3_NoEvolDist.pdf", width = 5, height = 3, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(family = "serif", mar = c(3, 4, 1, 1))
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(400, 800), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(400, 800, by = 25), tcl = -0.25, labels = FALSE)
        mtext("Distance spread", side = 2, line = 2.75)
        mtext("Number of loci", side = 1, line = 2)
        # Graph the data
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig3Data[[i]]$DistNoEvol, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig3Data[[i]]$LwrDistNoEvol, x1 = xLocs[i,], y1 = Fig3Data[[i]]$UprDistNoEvol,
                         lty = 1, col = ScenCols[i])
        }
        # Put the legend on the figure
        legend(x = -0.1, y = 750, legend = Scenarios[1:3], pch = 21:23, col = ScenCols[1:3], 
               pt.bg = ScenCols[1:3], bty = "n", yjust = 0.5)
        legend(x = 0.305, y = 710, legend = Scenarios[4:5], pch = 24:25, col = ScenCols[4:5],
               pt.bg = ScenCols[4:5], bty = "n", yjust = 0)
dev.off()

# Now make a figure for the increase in distance spread due to evolution
pdf(file = "ResultFigures/Figure3_DiffDist.pdf", width = 5, height = 3, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(family = "serif", mar = c(3, 4, 1, 1))
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, 250), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0, 300, by = 25), tcl = -0.25, labels = FALSE)
        mtext("Change in distance spread", side = 2, line = 2.75)
        mtext("Number of loci", side = 1, line = 2)
        # Graph the data
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig3Data[[i]]$Difference, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig3Data[[i]]$LwrDiff, x1 = xLocs[i,], y1 = Fig3Data[[i]]$UprDiff,
                         lty = 1, col = ScenCols[i])
        }
        # Put the legend on the figure
        legend(x = -0.1, y = 223, legend = Scenarios[1:3], pch = 21:23, col = ScenCols[1:3], 
               pt.bg = ScenCols[1:3], bty = "n", yjust = 0.5)
        legend(x = 0.305, y = 197, legend = Scenarios[4:5], pch = 24:25, col = ScenCols[4:5],
               pt.bg = ScenCols[4:5], bty = "n", yjust = 0)
        # Finally, add a horizontal line at 0
        abline(h = 0, lty = 2, col = "grey")
dev.off()
###### Figure 4
# First get the extinction risk data together
Fig4Data <- vector(mode = "list", length = 5)
NoEvolExtRisk <- matrix(data = NA, nrow = 5, ncol = 6)
for(i in 1:5){
        EvolExtRisk <- rep(NA, 6)
        ExtRiskLwr <- rep(NA, 6)
        ExtRiskUpr <- rep(NA, 6)
        for(l in 1:6){
                CurData <- subset(Results, L == Lseq[l] & Haploid == HapVals[i] & monoecious == monoVals[i] & omega == omegaVals[i])
                EvolExtRisk[l] <- mean(CurData$ExtRiskEvol)
                NoEvolExtRisk[i,l] <- mean(CurData$ExtRiskNoEvol)
                CI <- 1.96 * sqrt((EvolExtRisk[l] * (1 - EvolExtRisk[l])) / nrow(CurData))
                ExtRiskLwr[l] <- EvolExtRisk[l] - CI
                ExtRiskUpr[l] <- EvolExtRisk[l] + CI
        }
        Fig4Data[[i]] <- data.frame(L = Lseq, EvolExtRisk, ExtRiskLwr, ExtRiskUpr)
}
names(Fig4Data) <- Scenarios

pdf(file = "ResultFigures/NoEvolRisk.pdf", width = 5, height = 3, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mar = c(5, 4, 1, 1), family = "serif")
        plot(x = NA, y = NA, xlim = c(-0.05, 1.05), ylim = c(0.4, 1), xlab = "Number of loci", 
             ylab = "Extinction risk", las =1, xaxt = "n")
        axis(1, at = LseqLocs, labels = Lseq)
        axis(2, at = seq(0.4, 1, by = 0.02), tcl = -0.25, labels = FALSE)
        legend("bottom", legend = Scenarios, pch = 21:25, col = ScenCols, pt.bg = ScenCols, bty = "n", ncol = 2)
        for(i in 1:5){
                points(x = xLocs[i,], y = NoEvolExtRisk[i,], col = ScenCols[i], 
                      pch = 20 + i, bg = ScenCols[i])
                
        }
        
dev.off()
# There's very little variation among loci (as expected), so I will average across them
#       In fact, I think this would be better demonstrated in a table since there's really
#       only variation between mate finding Allee effect vs. not.
NoEvolRisk <- rowMeans(NoEvolExtRisk)
# Asexual               :0.8028
# Obligate selfing      :0.8332
# Partial selfing       :0.8425
# Obligate outcrossing  :0.9998
# Sexual (dioecious)    :1.0000

pdf(file = "ResultFigures/Figure4.pdf", width = 5, height = 3, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(family = "serif", mar = c(3, 4, 1, 1))
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, 1), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0, 1, by = 0.05), tcl = -0.25, labels = FALSE)
        mtext("Extinction probability", side = 2, line = 2.5)
        mtext("Number of loci", side = 1, line = 2)
        # Graph the data
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig4Data[[i]]$EvolExtRisk, pch = 20+i, col = ScenCols[i], bg = ScenCols[i])
                #segments(x0 = xLocs[i,], y0 = Fig4Data[[i]]$ExtRiskLwr, x1 = xLocs[i,], y1 = Fig4Data[[i]]$ExtRiskUpr,
                #         lty = 1, col = ScenCols[i])
        }
        # Put the legend on the figure
        legend(x = -0.1, y = 0.3, legend = Scenarios[1:3], pch = 21:23, col = ScenCols[1:3], 
               pt.bg = ScenCols[1:3], bty = "n", yjust = 1)
        legend(x = 0.295, y = 0.2, legend = Scenarios[4:5], pch = 24:25, col = ScenCols[4:5],
               pt.bg = ScenCols[4:5], bty = "n", yjust = 1)
dev.off()

###### Figure 5
# First get the data together (mean distance spread and CV for each scenario)
DeltaExtRisk <- matrix(NA, nrow = 5, ncol = 6)
for(i in 1:5){
        DeltaExtRisk[i,] <- Fig4Data[[i]]$EvolExtRisk - NoEvolRisk[i] 
}

# Make the figure
pdf(file = "ResultFigures/Figure5.pdf", width = 5, height = 3, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(family = "serif", mar = c(3, 4, 1, 1))
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(-0.9, 0.05), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(-0.9, 0.1, by = 0.05), tcl = -0.25, labels = FALSE)
        mtext("Change in extinction probability", side = 2, line = 2.75)
        mtext("Number of loci", side = 1, line = 2)
        abline(h = 0, lty = 2, col = "grey")
        # Graph the data
        for(i in 1:5){
                points(x = xLocs[i,], y = DeltaExtRisk[i,], pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
        }
        # Put the legend on the figure
        legend(x = -0.1, y = -0.575, legend = Scenarios[1:3], pch = 21:23, col = ScenCols[1:3], 
               pt.bg = ScenCols[1:3], bty = "n", yjust = 1)
        legend(x = 0.305, y = -0.665, legend = Scenarios[4:5], pch = 24:25, col = ScenCols[4:5],
               pt.bg = ScenCols[4:5], bty = "n", yjust = 1)
dev.off()
