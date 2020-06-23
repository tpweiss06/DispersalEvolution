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
setwd("~/Desktop/Wyoming/DispersalEvolution/GitRepo/")
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
pdf(file = "ResultFigures/Figure2.pdf", width = 5, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,3,1,1), mar = c(1.5,1,1.5,1), family = "serif")
        
        # Make the phenotype figure
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0,1.35), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0, 1.4, by = 0.05), tcl = -0.25, labels = FALSE)
        mtext(expression(paste(Delta, bar("d"), sep = "")), side = 2, line = 3)
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
        mtext("Change in genetic diversity", side = 2, line = 3)
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig2Data[[i]]$MeanGenExp, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig2Data[[i]]$LwrGenExp, x1 = xLocs[i,], y1 = Fig2Data[[i]]$UprGenExp,
                         col = ScenCols[i], lty = 1)
        }
        mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
        mtext("b", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
dev.off()

# Now make the supplemental figure with the range shift results
pdf(file = "ResultFigures/FigureS1_Shift.pdf", width = 5, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,3,1,1), mar = c(1.5,1,1.5,1), family = "serif")

        # Make the phenotype figure
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0,1.65), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0, 1.65, by = 0.05), tcl = -0.25, labels = FALSE)
        mtext(expression(paste(Delta, bar("d"), sep = "")), side = 2, line = 3)
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
        mtext("Change in genetic diversity", side = 2, line = 3)
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig2Data[[i]]$MeanGenShift, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig2Data[[i]]$LwrGenShift, x1 = xLocs[i,], y1 = Fig2Data[[i]]$UprGenShift,
                         col = ScenCols[i], lty = 1)
        }
        mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
        mtext("b", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
dev.off()

###### Figure 3
# First get the data together (mean distance spread and CV for each scenario)
R <- 2
r <- log(R)
Fig3Data <- vector(mode = "list", length = 5)
for(i in 1:5){
        Dist <- rep(NA, 6)
        LwrDist <- rep(NA, 6)
        UprDist <- rep(NA, 6)
        MeanDeviation <- rep(NA, 6)
        LwrDeviation <- rep(NA, 6)
        UprDeviation <- rep(NA, 6)
        
        for(l in 1:6){
                CurData <- subset(Results, L == Lseq[l] & Haploid == HapVals[i] & monoecious == monoVals[i] & omega == omegaVals[i])
                Distances <- c(CurData$distance_1, abs(CurData$distance_2))
                Dist[l] <- mean(Distances)
                LwrDist[l] <- quantile(Distances, probs = 0.25)
                UprDist[l] <-  quantile(Distances, probs = 0.75)
                D <- (1/2) * mean(CurData$InitPhenExp)^2
                ExpDist <- sqrt(2 * r * D) * 200
                Deviations <- Distances - ExpDist
                MeanDeviation[l] <- mean(Deviations)
                LwrDeviation[l] <- quantile(Deviations, probs = 0.25)
                UprDeviation[l] <- quantile(Deviations, probs = 0.75)
        }
        Fig3Data[[i]] <- data.frame(L = Lseq, Dist, LwrDist, UprDist, MeanDeviation, LwrDeviation, UprDeviation)
}

# Now make the main figure
pdf(file = "ResultFigures/Figure3.pdf", width = 5, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,4,3,1.5), mar = c(1.5,1,1.5,1), family = "serif")
        
        # Plot the mean distance spread
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(475, 700), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        mtext("Distance spread", side = 2, line = 3.2)
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig3Data[[i]]$Dist, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig3Data[[i]]$LwrDist, x1 = xLocs[i,], y1 = Fig3Data[[i]]$UprDist,
                         col = ScenCols[i], lty = 1)
        }
        
        # Put the legend on the figure
        legend(x = 0.2, y = 785, legend = Scenarios[1:2], inset = -0.25, lty = 1, col = ScenCols[1:2], bty = "n", xpd = NA)
        legend(x = 0.8, y = 775, legend = Scenarios[3:5], inset = -0.25, lty = 1, col = ScenCols[3:5], 
               xpd = NA, horiz = TRUE, box.lty = 2)
        mtext("Sexual (monoecious)", side = 3, line = 1.25, outer = TRUE, adj = 0.65)
        
        # Plot the deviation from the Skellam approximation
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(-100, 100), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        mtext("Deviation from expectation", side = 2, line = 3.2)
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig3Data[[i]]$MeanDeviation, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = Fig3Data[[i]]$LwrDeviation, x1 = xLocs[i,], y1 = Fig3Data[[i]]$UprDeviation,
                         col = ScenCols[i], lty = 1)
        }
        mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
dev.off()

###### Figure 4
# First get the extinction risk data together
Fig4Data <- vector(mode = "list", length = 5)
NoEvolExtRisk <- matrix(data = NA, nrow = 5, ncol = 6)
for(i in 1:5){
        EvolExtRisk <- rep(NA, 6)
        for(l in 1:6){
                CurData <- subset(Results, L == Lseq[l] & Haploid == HapVals[i] & monoecious == monoVals[i] & omega == omegaVals[i])
                EvolExtRisk[l] <- mean(CurData$ExtRiskEvol)
                NoEvolExtRisk[i,l] <- mean(CurData$ExtRiskNoEvol)
        }
        Fig4Data[[i]] <- data.frame(L = Lseq, EvolExtRisk)
}
names(Fig4Data) <- Scenarios

pdf(file = "ResultFigures/NoEvolRisk.pdf", width = 5, height = 3, onefile = FALSE, paper = "special", useDingbats = FALSE)
        plot(x = NA, y = NA, xlim = range(Lseq), ylim = c(0.7, 1), xlab = "Number of loci", ylab = "Extinction risk", las =1)
        for(i in 1:5){
                lines(x = Lseq, y = NoEvolExtRisk[i,], col = ScenCols[i])
        }
dev.off()
# There's very little variation among loci (as expected), so I will average across them
#       In fact, I think this would be better demonstrated in a table since there's really
#       only variation between mate finding Allee effect vs. not.
NoEvolRisk <- rowMeans(NoEvolExtRisk)
# Asexual               :0.803
# Sexual (dioecious)    :1
# Obligate selfing      :0.832
# Partial selfing       :0.842
# Obligate outcrossing  :1

pdf(file = "ResultFigures/Figure4.pdf", width = 5, height = 3, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(family = "serif")
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, 1), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        mtext("Extinction probability", side = 2, line = 3.2)
        mtext("Number of loci", side = 1, line = 2.5)
        # Put the legend on the figure
        legend(x = 0.2, y = 1.2, legend = Scenarios[1:2], inset = -0.25, lty = 1, col = ScenCols[1:2], bty = "n", xpd = NA)
        legend(x = 0.8, y = 1.17, legend = Scenarios[3:5], inset = -0.25, lty = 1, col = ScenCols[3:5], 
               xpd = NA, horiz = TRUE, box.lty = 2)
        mtext("Sexual (monoecious)", side = 3, line = 2.75, adj = 0.65)
        # Graph the data
        for(i in 1:5){
                points(x = xLocs[i,], y = Fig4Data[[i]]$EvolExtRisk, pch = 20+i, col = ScenCols[i], bg = ScenCols[i])
        }
dev.off()

###### Figure 5
# First get the data together (mean distance spread and CV for each scenario)
DeltaExtRisk <- matrix(NA, nrow = 5, ncol = 6)
for(i in 1:5){
        DeltaExtRisk[i,] <- Fig4Data[[i]]$EvolExtRisk - NoEvolRisk[i] 
}

# Make the figure
pdf(file = "ResultFigures/Figure5.pdf", width = 5, height = 3, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(family = "serif")
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(-0.8, 0.1), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        mtext("Change in extinction probability", side = 2, line = 3.2)
        mtext("Number of loci", side = 1, line = 2.5)
        # Put the legend on the figure
        legend(x = 0.2, y = 0.275, legend = Scenarios[1:2], inset = -0.25, lty = 1, col = ScenCols[1:2], bty = "n", xpd = NA)
        legend(x = 0.8, y = 0.25, legend = Scenarios[3:5], inset = -0.25, lty = 1, col = ScenCols[3:5], 
               xpd = NA, horiz = TRUE, box.lty = 2)
        mtext("Sexual (monoecious)", side = 3, line = 2.75, adj = 0.65)
        abline(h = 0, lty = 2, col = "grey")
        # Graph the data
        for(i in 1:5){
                points(x = xLocs[i,], y = DeltaExtRisk[i,], pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
        }
dev.off()
