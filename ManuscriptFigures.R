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
Scenarios <- c("Asexual", "Sexual (dioecious)", "Obligate selfing", "Partial selfing", "Obligate outcrossing")
HapVals <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
monoVals <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
omegaVals <- c(0, 0, 1, 0.5, 0)
Lseq <- unique(Results$L)
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
          MeanPhenExp[l] <- mean(CurData$DeltaPhenExp)
          LwrPhenExp[l] <- quantile(CurData$DeltaPhenExp, probs = 0.25)
          UprPhenExp[l] <- quantile(CurData$DeltaPhenExp, probs = 0.75)
          MeanGenExp[l] <- mean(CurData$DeltaGenExp)
          LwrGenExp[l] <- quantile(CurData$DeltaGenExp, probs = 0.25)
          UprGenExp[l] <- quantile(CurData$DeltaGenExp, probs = 0.75)
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
names(Fig2Data) <- Scenarios

# Now make the main figure
pdf(file = "ResultFigures/Figure2.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,4,3,1.5), mar = c(1.5,1,1.5,1), family = "serif")
        plot(NA, NA, xlim = range(Lseq), ylim = c(0,1.25), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq, labels = Lseq)
        mtext("Change in mean phenotype", side = 2, line = 3)
        
        # Put the legend on the figure
        legend(x = 1, y = 1.8, legend = Scenarios[1:2], inset = -0.25, lty = 1, col = ScenCols[1:2], bty = "n", xpd = NA)
        legend(x = 9, y = 1.7, legend = Scenarios[3:5], inset = -0.25, lty = 1, col = ScenCols[3:5], 
               xpd = NA, horiz = TRUE, box.lty = 2)
        mtext("Sexual (monoecious)", side = 3, line = 1.25, outer = TRUE, adj = 0.65)
        # Graph the data
        for(i in 1:5){
                Ints <- jitter(Lseq, factor = 1)
                TempLseq <- c(Ints[1], Lseq[2:5], Ints[6])
                lines(x = TempLseq, y = Fig2Data[[i]]$MeanPhenExp, lty = 1, col = ScenCols[i])
                segments(x0 = Ints, y0 = Fig2Data[[i]]$LwrPhenExp, x1 = Ints, y1 = Fig2Data[[i]]$UprPhenExp,
                         col = ScenCols[i], lty = 1)
        }
        plot(NA, NA, xlim = range(Lseq), ylim = c(-3, 0), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq)
        mtext("Change in genetic diversity", side = 2, line = 3)
        for(i in 1:5){
                Ints <- jitter(Lseq, factor = 1)
                TempLseq <- c(Ints[1], Lseq[2:5], Ints[6])
                lines(x = TempLseq, y = Fig2Data[[i]]$MeanGenExp, lty = 1, col = ScenCols[i])
                segments(x0 = Ints, y0 = Fig2Data[[i]]$LwrGenExp, x1 = Ints, y1 = Fig2Data[[i]]$UprGenExp,
                         col = ScenCols[i], lty = 1)
        }
        mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
dev.off()


# Now make the supplemental figure with the range shift results
pdf(file = "ResultFigures/FigureS1_Shift.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,4,3,1.5), mar = c(1.5,1,1.5,1), family = "serif")
        plot(NA, NA, xlim = range(Lseq), ylim = c(0,1.5), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq, labels = Lseq)
        mtext("Change in mean phenotype", side = 2, line = 3)

        # Put the legend on the figure
        legend(x = 1, y = 2.1, legend = Scenarios[1:2], inset = -0.25, lty = 1, col = ScenCols[1:2], bty = "n", xpd = NA)
        legend(x = 9, y = 2, legend = Scenarios[3:5], inset = -0.25, lty = 1, col = ScenCols[3:5], 
               xpd = NA, horiz = TRUE, box.lty = 2)
        mtext("Sexual (monoecious)", side = 3, line = 1.25, outer = TRUE, adj = 0.65)
        # Graph the data
        for(i in 1:5){
                Ints <- jitter(Lseq, factor = 1)
                TempLseq <- c(Ints[1], Lseq[2:5], Ints[6])
                lines(x = TempLseq, y = Fig2Data[[i]]$MeanPhenShift, lty = 1, col = ScenCols[i])
                segments(x0 = Ints, y0 = Fig2Data[[i]]$LwrPhenShift, x1 = Ints, y1 = Fig2Data[[i]]$UprPhenShift,
                         col = ScenCols[i], lty = 1)
        }
        plot(NA, NA, xlim = range(Lseq), ylim = c(-6, 0), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq)
        mtext("Change in genetic diversity", side = 2, line = 3)
        for(i in 1:5){
                Ints <- jitter(Lseq, factor = 1)
                TempLseq <- c(Ints[1], Lseq[2:5], Ints[6])
                lines(x = TempLseq, y = Fig2Data[[i]]$MeanGenShift, lty = 1, col = ScenCols[i])
                segments(x0 = Ints, y0 = Fig2Data[[i]]$LwrGenShift, x1 = Ints, y1 = Fig2Data[[i]]$UprGenShift,
                         col = ScenCols[i], lty = 1)
        }
        mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
dev.off()

###### Figure 3
# First get the data together (mean distance spread and CV for each scenario)
Fig3Data <- vector(mode = "list", length = 5)
for(i in 1:5){
        Dist <- rep(NA, 6)
        CV <- rep(NA, 6)
        StdDev <- rep(NA, 6)
        for(l in 1:6){
                CurData <- subset(Results, L == Lseq[l] & Haploid == HapVals[i] & monoecious == monoVals[i] & omega == omegaVals[i])
                Dist[l] <- mean(CurData$distance)
                StdDev[l] <- sd(CurData$distance)
                CV[l] <-  StdDev[l] / Dist[l]
        }
        Fig3Data[[i]] <- data.frame(L = Lseq, Dist, StdDev, CV)
}
names(Fig3Data) <- Scenarios

# Now make the main figure
pdf(file = "ResultFigures/Figure3.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,4,3,1.5), mar = c(1.5,1,1.5,1), family = "serif")
        plot(NA, NA, xlim = range(Lseq), ylim = c(475, 700), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq, labels = Lseq)
        mtext("Mean distance spread", side = 2, line = 3.2)

        # Put the legend on the figure
        legend(x = 1, y = 785, legend = Scenarios[1:2], inset = -0.25, lty = 1, col = ScenCols[1:2], bty = "n", xpd = NA)
        legend(x = 9, y = 775, legend = Scenarios[3:5], inset = -0.25, lty = 1, col = ScenCols[3:5], 
               xpd = NA, horiz = TRUE, box.lty = 2)
        mtext("Sexual (monoecious)", side = 3, line = 1.25, outer = TRUE, adj = 0.65)
        # Graph the data
        for(i in 1:5){
                lines(x = Lseq, y = Fig3Data[[i]]$Dist, lty = 1, col = ScenCols[i])
        }
        plot(NA, NA, xlim = range(Lseq), ylim = c(0.105, 0.135), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq)
        mtext("CV of distance spread", side = 2, line = 3.2)
        for(i in 1:5){
                lines(x = Lseq, y = Fig3Data[[i]]$CV, lty = 1, col = ScenCols[i])
        }
        mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
dev.off()

pdf(file = "ResultFigures/Figure3Alt.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,4,3,1.5), mar = c(1.5,1,1.5,1), family = "serif")
        plot(NA, NA, xlim = range(Lseq), ylim = c(475, 700), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq, labels = Lseq)
        mtext(expression(mu ["spread"]), side = 2, line = 3)
        
        # Put the legend on the figure
        legend(x = 1, y = 785, legend = Scenarios[1:2], inset = -0.25, lty = 1, col = ScenCols[1:2], bty = "n", xpd = NA)
        legend(x = 9, y = 775, legend = Scenarios[3:5], inset = -0.25, lty = 1, col = ScenCols[3:5], 
               xpd = NA, horiz = TRUE, box.lty = 2)
        mtext("Sexual (monoecious)", side = 3, line = 1.25, outer = TRUE, adj = 0.65)
        # Graph the data
        for(i in 1:5){
                lines(x = Lseq, y = Fig3Data[[i]]$Dist, lty = 1, col = ScenCols[i])
        }
        plot(NA, NA, xlim = range(Lseq), ylim = c(50, 85), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq)
        mtext(expression(sigma ["spread"]), side = 2, line = 3)
        for(i in 1:5){
                lines(x = Lseq, y = Fig3Data[[i]]$StdDev, lty = 1, col = ScenCols[i])
        }
        mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
dev.off()

pdf(file = "ResultFigures/Figure3Alt2.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,4,3,1.5), mar = c(1.5,1,1.5,1), family = "serif")
        plot(NA, NA, xlim = range(Lseq), ylim = c(475, 700), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq, labels = Lseq)
        mtext(expression(mu ["spread"]), side = 2, line = 3)

        # Put the legend on the figure
        legend(x = 1, y = 785, legend = Scenarios[1:2], inset = -0.25, lty = 1, col = ScenCols[1:2], bty = "n", xpd = NA)
        legend(x = 9, y = 775, legend = Scenarios[3:5], inset = -0.25, lty = 1, col = ScenCols[3:5], 
               xpd = NA, horiz = TRUE, box.lty = 2)
        mtext("Sexual (monoecious)", side = 3, line = 1.25, outer = TRUE, adj = 0.65)
        # Graph the data
        for(i in 1:5){
                lines(x = Lseq, y = Fig3Data[[i]]$Dist, lty = 1, col = ScenCols[i])
        }
        plot(NA, NA, xlim = c(475, 700), ylim = c(50, 85), main = "", ylab = "", 
             xlab = "", las = 1)
        mtext(expression(sigma ["spread"]), side = 2, line = 3)
        mtext(expression(mu ["spread"]), side = 1, line = 3)
        for(i in 1:5){
                lines(x = Fig3Data[[i]]$Dist, y = Fig3Data[[i]]$StdDev, lty = 1, col = ScenCols[i])
        }
        #mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
dev.off()

# Maybe try one more alt figure where the bottom panel is loci on the x axis and 
#       deviation from pulled wave expectation (root(2rd)) on the y-axis with d
#       being the final dispersal phenotype and v measured over the last 50 
#       generations...
#       To do this, I'll have to start back with data extraction (or do a separate
#       script). I could just add a velocity column to the existing one and then
#       I would have all the info I need (I think?). Double check that R could be
#       accurately measured for both asexual and sexual (what do I do in the reproduce
#       function?). Check the scripts from the AmNat paper to double check.

# Can take ln(R) to be r; need to calculate p = probability of leaving patch (2017 Nat Commun paper)
#       - I think I just use R even in dioecious scenarios since I am tracking the full population
#               and the equation takes account of adjusting R for females
# To calculate p, take the mean dispersal phenotype and use 1 - pexp(q = 1, rate = 1/dBar) to get the
#       probability of dispersing out of the current patch
# So, I will need to adjust the data extraction script to also give me: (1) mean dispersal phenotype at edge,
#       (2) expansion speed over previous 20 generations, (3) expansion speed over previous 10 generations

###### Figure 4
# First get the data together (mean distance spread and CV for each scenario)
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

pdf(file = "ResultFigures/NoEvolRisk.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
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

pdf(file = "ResultFigures/Figure4.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(family = "serif")
        plot(NA, NA, xlim = range(Lseq), ylim = c(0, 1), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq, labels = Lseq)
        mtext("Extinction probability", side = 2, line = 3.2)
        mtext("Number of loci", side = 1, line = 2.5)
        # Put the legend on the figure
        legend(x = 1, y = 1.2, legend = Scenarios[1:2], inset = -0.25, lty = 1, col = ScenCols[1:2], bty = "n", xpd = NA)
        legend(x = 9, y = 1.17, legend = Scenarios[3:5], inset = -0.25, lty = 1, col = ScenCols[3:5], 
               xpd = NA, horiz = TRUE, box.lty = 2)
        mtext("Sexual (monoecious)", side = 3, line = 2.75, adj = 0.65)
        # Graph the data
        for(i in 1:5){
                lines(x = Lseq, y = Fig4Data[[i]]$EvolExtRisk, lty = 1, col = ScenCols[i])
        }
dev.off()

###### Figure 5
# Figure 5: A single panel showing the change in extinction risk due to evolution
#    across loci.
# First get the data together (mean distance spread and CV for each scenario)
DeltaExtRisk <- matrix(NA, nrow = 5, ncol = 6)
for(i in 1:5){
        DeltaExtRisk[i,] <- Fig4Data[[i]]$EvolExtRisk - NoEvolRisk[i] 
}

# Make the figure
pdf(file = "ResultFigures/Figure5.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(family = "serif")
        plot(NA, NA, xlim = range(Lseq), ylim = c(-0.8, 0.1), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = Lseq, labels = Lseq)
        mtext("Change in extinction probability", side = 2, line = 3.2)
        mtext("Number of loci", side = 1, line = 2.5)
        # Put the legend on the figure
        legend(x = 1, y = 0.275, legend = Scenarios[1:2], inset = -0.25, lty = 1, col = ScenCols[1:2], bty = "n", xpd = NA)
        legend(x = 9, y = 0.25, legend = Scenarios[3:5], inset = -0.25, lty = 1, col = ScenCols[3:5], 
               xpd = NA, horiz = TRUE, box.lty = 2)
        mtext("Sexual (monoecious)", side = 3, line = 2.75, adj = 0.65)
        abline(h = 0, lty = 2, col = "grey")
        # Graph the data
        for(i in 1:5){
                lines(x = Lseq, y = DeltaExtRisk[i,], lty = 1, col = ScenCols[i])
        }
dev.off()
