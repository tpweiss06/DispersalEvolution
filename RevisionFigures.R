# This script will create an updated version of Figures 1 & 2 according to the
#    reviewers' comments. In these versions, each figure will have 4 panels (2x2)
#    with the top row corresponding to phenotype values and the bottom row
#    corresponding to genetic variance values. Figure 1 will still focus on
#    "pure" range expansions and Figure 2 will show range shifts. The first
#    column of each will correspond to scenarios with no genetic mixing
#    (asexual and obligately self-fertilizing populations). The second column
#    will focus on the remaining three scenarios. Within each panel, the x-axis
#    will still show the number of loci and each scenario at each loci value will
#    have a filled and hollow point to correspond to initial and final values.

setwd("~/Desktop/GitHubRepos/DispersalEvolution/")

library(RColorBrewer)
ScenCols <- brewer.pal(n = 5, name = "Dark2")
Scenarios <- c("Asexual", "Obligate selfing", "Partial selfing", "Obligate outcrossing", "Sexual (dioecious)")
HapVals <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
monoVals <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
omegaVals <- c(0, 1, 0.5, 0, 0)
Lseq <-c(1,2,4,8,16,32)
LseqLocs <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
offset <- 0.015
xLocs <- matrix(NA, nrow = 6, ncol = 6)
xLocs[1,] <- LseqLocs - 3*offset
xLocs[2,] <- LseqLocs - 2*offset
xLocs[3,] <- LseqLocs - offset
xLocs[4,] <- LseqLocs + offset
xLocs[5,] <- LseqLocs + 2*offset
xLocs[6,] <- LseqLocs + 3*offset
ScenPch <- c(21, 22, 23, 24, 25)

# First, load in the initial values to use for the graph
InitResults <- read.csv("SimsWithInitVals.csv")
InitPlotData <- vector(mode = "list", length = 5)
for(i in 1:5){
     MeanPhen <- rep(NA, 6)
     LwrPhen <- rep(NA, 6)
     UprPhen <- rep(NA, 6)
     MeanGen <- rep(NA, 6)
     LwrGen <- rep(NA, 6)
     UprGen <- rep(NA, 6)
     for(l in 1:6){
          CurRows <- which(InitResults$L == Lseq[l] & InitResults$Haploid == HapVals[i] & InitResults$monoecious == monoVals[i] & InitResults$omega == omegaVals[i])
          CurData <- InitResults[CurRows,]
          MeanPhen[l] <- mean(CurData$MeanPhen)
          LwrPhen[l] <- quantile(CurData$MeanPhen, probs = 0.25)
          UprPhen[l] <- quantile(CurData$MeanPhen, probs = 0.75)
          MeanGen[l] <- mean(CurData$GenVar)
          LwrGen[l] <- quantile(CurData$GenVar, probs = 0.25)
          UprGen[l] <- quantile(CurData$GenVar, probs = 0.75)
     }
     InitPlotData[[i]] <- data.frame(L = Lseq, MeanPhen, LwrPhen, UprPhen, MeanGen,
                                  LwrGen, UprGen)
}

# Next, load in the results with the final values for the expansion
FinalResults <- read.csv("SimsWithNewResults.csv")
FinalPlotData <- vector(mode = "list", length = 5)
for(i in 1:5){
     MeanPhen <- rep(NA, 6)
     LwrPhen <- rep(NA, 6)
     UprPhen <- rep(NA, 6)
     MeanGen <- rep(NA, 6)
     LwrGen <- rep(NA, 6)
     UprGen <- rep(NA, 6)
     for(l in 1:6){
          CurRows <- which(FinalResults$L == Lseq[l] & FinalResults$Haploid == HapVals[i] & FinalResults$monoecious == monoVals[i] & FinalResults$omega == omegaVals[i])
          CurData <- FinalResults[CurRows,]
          FinalPhen <- c(CurData$PhenExp_1, CurData$PhenExp_2)
          FinalGen <- c(CurData$GenExp_1, CurData$GenExp_2)
          MeanPhen[l] <- mean(FinalPhen)
          LwrPhen[l] <- quantile(FinalPhen, probs = 0.25)
          UprPhen[l] <- quantile(FinalPhen, probs = 0.75)
          MeanGen[l] <- mean(FinalGen)
          LwrGen[l] <- quantile(FinalGen, probs = 0.25)
          UprGen[l] <- quantile(FinalGen, probs = 0.75)
     }
     FinalPlotData[[i]] <- data.frame(L = Lseq, MeanPhen, LwrPhen, UprPhen, MeanGen,
                                     LwrGen, UprGen)
}

# Now make the updated Figure 1
PhenMax <- 3
GenMax <- 3.5
pdf(file = "ResultFigures/Figure1a_new.pdf", width = 10, height = 7, onefile = FALSE, paper = "special", useDingbats = FALSE)
     par(mfrow = c(2,2), oma = c(3,3,1,1), mar = c(2,2,2,2), family = "serif")

     # Make the phenotype figures: first non-mixing scenarios
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, PhenMax), main = "", ylab = "", 
          xlab = "", las = 1, xaxt = "n", cex.axis = 1.5)
     axis(side = 1, at = LseqLocs, labels = Lseq, cex.axis = 1.5)
     axis(side = 2, at = seq(0, PhenMax, by = 0.125), tcl = -0.25, labels = FALSE)
     mtext("Dispersal phenotype", side = 2, line = 3, cex = 1.5)
     InitLocs <- c(2, 4)
     FinalLocs <- c(3, 5)
     for(i in 1:2){
          # Initial values
          segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i]]$LwrPhen, 
                   x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i]]$UprPhen,
                   col = ScenCols[i], lty = 1)
          points(x = xLocs[InitLocs[i],], y = InitPlotData[[i]]$MeanPhen, 
                 pch = ScenPch[i], col = ScenCols[i], bg = "white", cex = 1.5)
          # Final values
          points(x = xLocs[FinalLocs[i],], y = FinalPlotData[[i]]$MeanPhen, 
                 pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i], cex = 1.5)
          segments(x0 = xLocs[FinalLocs[i],], y0 = FinalPlotData[[i]]$LwrPhen, 
                   x1 = xLocs[FinalLocs[i],], y1 = FinalPlotData[[i]]$UprPhen,
                   col = ScenCols[i], lty = 1)
     }
     # Put the letter and legend on the figure
     mtext("a", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
     legend("top", legend = c(Scenarios[1:2], "Initial", "Final"), 
            pch = c(ScenPch[1:2], ScenPch[1], ScenPch[1]), 
            col = c(ScenCols[1:2], "black", "black"), bty = "n", 
            pt.bg = c(ScenCols[1:2], "white", "black"), inset = -0.01, ncol = 2, cex = 1.25)
 
     # Make the phenotype figures for the mixing scenarios
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, PhenMax), main = "", ylab = "", 
          xlab = "", las = 1, xaxt = "n", cex.axis = 1.5)
     axis(side = 1, at = LseqLocs, labels = Lseq, cex.axis = 1.5)
     axis(side = 2, at = seq(0, PhenMax, by = 0.125), tcl = -0.25, labels = FALSE)
     InitLocs <- c(1,3,5)
     FinalLocs <- c(2,4,6)
     for(i in 1:3){
          # Initial values
          segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i+2]]$LwrPhen, 
                   x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i+2]]$UprPhen,
                   col = ScenCols[i+2], lty = 1)
          points(x = xLocs[InitLocs[i],], y = InitPlotData[[i+2]]$MeanPhen, 
                 pch = ScenPch[i+2], col = ScenCols[i+2], bg = "white", cex = 1.5)
          # Final values
          points(x = xLocs[FinalLocs[i],], y = FinalPlotData[[i+2]]$MeanPhen, 
                 pch = ScenPch[i+2], col = ScenCols[i+2], bg = ScenCols[i+2], cex = 1.5)
          segments(x0 = xLocs[FinalLocs[i],], y0 = FinalPlotData[[i+2]]$LwrPhen, 
                   x1 = xLocs[FinalLocs[i],], y1 = FinalPlotData[[i+2]]$UprPhen,
                   col = ScenCols[i+2], lty = 1)
     }
     # Put the letter and legend on the figure
     mtext("b", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
     legend("topleft", legend = Scenarios[3:5], pch = ScenPch[3:5], col = ScenCols[3:5], 
            bty = "n", pt.bg = ScenCols[3:5], inset = -0.01, ncol = 2, cex = 1.25)
     
     # Make the genetic variance figures: first non-mixing scenarios
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, GenMax), main = "", ylab = "", 
          xlab = "", las = 1, xaxt = "n", cex.axis = 1.5)
     axis(side = 1, at = LseqLocs, labels = Lseq, cex.axis = 1.5)
     axis(side = 2, at = seq(0, GenMax, by = 0.125), tcl = -0.25, labels = FALSE)
     mtext("Genetic variance", side = 2, line = 3, cex = 1.5)
     InitLocs <- c(2, 4)
     FinalLocs <- c(3, 5)
     for(i in 1:2){
          # Initial values
          segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i]]$LwrGen, 
                   x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i]]$UprGen,
                   col = ScenCols[i], lty = 1)
          points(x = xLocs[InitLocs[i],], y = InitPlotData[[i]]$MeanGen, 
                 pch = ScenPch[i], col = ScenCols[i], bg = "white", cex = 1.5)
          # Final values
          points(x = xLocs[FinalLocs[i],], y = FinalPlotData[[i]]$MeanGen, 
                 pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i], cex = 1.5)
          segments(x0 = xLocs[FinalLocs[i],], y0 = FinalPlotData[[i]]$LwrGen, 
                   x1 = xLocs[FinalLocs[i],], y1 = FinalPlotData[[i]]$UprGen,
                   col = ScenCols[i], lty = 1)
     }
     # Put the letter on the figure
     mtext("c", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
     
     # Make the genetic variance figures for the mixing scenarios
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, GenMax), main = "", ylab = "", 
          xlab = "", las = 1, xaxt = "n", cex.axis = 1.5)
     axis(side = 1, at = LseqLocs, labels = Lseq, cex.axis = 1.5)
     axis(side = 2, at = seq(0, GenMax, by = 0.125), tcl = -0.25, labels = FALSE)
     InitLocs <- c(1,3,5)
     FinalLocs <- c(2,4,6)
     for(i in 1:3){
          # Initial values
          segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i+2]]$LwrGen, 
                   x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i+2]]$UprGen,
                   col = ScenCols[i+2], lty = 1)
          points(x = xLocs[InitLocs[i],], y = InitPlotData[[i+2]]$MeanGen, 
                 pch = ScenPch[i+2], col = ScenCols[i+2], bg = "white", cex = 1.5)
          # Final values
          points(x = xLocs[FinalLocs[i],], y = FinalPlotData[[i+2]]$MeanGen, 
                 pch = ScenPch[i+2], col = ScenCols[i+2], bg = ScenCols[i+2], cex = 1.5)
          segments(x0 = xLocs[FinalLocs[i],], y0 = FinalPlotData[[i+2]]$LwrGen, 
                   x1 = xLocs[FinalLocs[i],], y1 = FinalPlotData[[i+2]]$UprGen,
                   col = ScenCols[i+2], lty = 1)
     }
     # Put the letter and legend on the figure
     mtext("d", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
     
     # Add the x axis label
     mtext("Number of loci", outer = TRUE, line = 1.5, side = 1, cex = 1.5)
dev.off()

######### Next make the updated figure 2 for range shifts
# Remake the final plot data to include the range shift results
FinalPlotData <- vector(mode = "list", length = 5)
for(i in 1:5){
     MeanPhen <- rep(NA, 6)
     LwrPhen <- rep(NA, 6)
     UprPhen <- rep(NA, 6)
     MeanGen <- rep(NA, 6)
     LwrGen <- rep(NA, 6)
     UprGen <- rep(NA, 6)
     for(l in 1:6){
          CurRows <- which(FinalResults$L == Lseq[l] & FinalResults$Haploid == HapVals[i] & FinalResults$monoecious == monoVals[i] & FinalResults$omega == omegaVals[i])
          CurData <- FinalResults[CurRows,]
          FinalPhen <- CurData$PhenShift
          FinalGen <- CurData$GenShift
          MeanPhen[l] <- mean(FinalPhen, na.rm = TRUE)
          LwrPhen[l] <- quantile(FinalPhen, probs = 0.25, na.rm = TRUE)
          UprPhen[l] <- quantile(FinalPhen, probs = 0.75, na.rm = TRUE)
          MeanGen[l] <- mean(FinalGen, na.rm = TRUE)
          LwrGen[l] <- quantile(FinalGen, probs = 0.25, na.rm = TRUE)
          UprGen[l] <- quantile(FinalGen, probs = 0.75, na.rm = TRUE)
     }
     FinalPlotData[[i]] <- data.frame(L = Lseq, MeanPhen, LwrPhen, UprPhen, MeanGen,
                                      LwrGen, UprGen)
}

PhenMax <- 3.25
GenMax <- 3.5
pdf(file = "ResultFigures/Figure2b_new.pdf", width = 10, height = 7, onefile = FALSE, paper = "special", useDingbats = FALSE)
     par(mfrow = c(2,2), oma = c(3,3,1,1), mar = c(2,2,2,2), family = "serif")

     # Make the phenotype figures: first non-mixing scenarios
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(1.2, PhenMax), main = "", ylab = "", 
          xlab = "", las = 1, xaxt = "n", cex.axis = 1.5)
     axis(side = 1, at = LseqLocs, labels = Lseq, cex.axis = 1.5)
     axis(side = 2, at = seq(0, PhenMax, by = 0.125), tcl = -0.25, labels = FALSE)
     mtext("Dispersal phenotype", side = 2, line = 3, cex = 1.5)
     InitLocs <- c(2, 4)
     FinalLocs <- c(3, 5)
     for(i in 1:2){
          # Initial values
          segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i]]$LwrPhen, 
                   x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i]]$UprPhen,
                   col = ScenCols[i], lty = 1)
          points(x = xLocs[InitLocs[i],], y = InitPlotData[[i]]$MeanPhen, 
                 pch = ScenPch[i], col = ScenCols[i], bg = "white", cex = 1.5)
          # Final values
          points(x = xLocs[FinalLocs[i],], y = FinalPlotData[[i]]$MeanPhen, 
                 pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i], cex = 1.5)
          segments(x0 = xLocs[FinalLocs[i],], y0 = FinalPlotData[[i]]$LwrPhen, 
                   x1 = xLocs[FinalLocs[i],], y1 = FinalPlotData[[i]]$UprPhen,
                   col = ScenCols[i], lty = 1)
     }
     # Put the letter and legend on the figure
     mtext("a", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
     legend("top", legend = c(Scenarios[1:2], "Initial", "Final"), 
            pch = c(ScenPch[1:2], ScenPch[1], ScenPch[1]), 
            col = c(ScenCols[1:2], "black", "black"), bty = "n", 
            pt.bg = c(ScenCols[1:2], "white", "black"), inset = -0.01, ncol = 2, cex = 1.25)

     # Make the phenotype figures for the mixing scenarios
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(1.2, PhenMax), main = "", ylab = "", 
          xlab = "", las = 1, xaxt = "n", cex.axis = 1.5)
     axis(side = 1, at = LseqLocs, labels = Lseq, cex.axis = 1.5)
     axis(side = 2, at = seq(0, PhenMax, by = 0.125), tcl = -0.25, labels = FALSE)
     InitLocs <- c(1,3,5)
     FinalLocs <- c(2,4,6)
     for(i in 1:3){
          # Initial values
          segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i+2]]$LwrPhen, 
                   x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i+2]]$UprPhen,
                   col = ScenCols[i+2], lty = 1)
          points(x = xLocs[InitLocs[i],], y = InitPlotData[[i+2]]$MeanPhen, 
                 pch = ScenPch[i+2], col = ScenCols[i+2], bg = "white", cex = 1.5)
          # Final values
          points(x = xLocs[FinalLocs[i],], y = FinalPlotData[[i+2]]$MeanPhen, 
                 pch = ScenPch[i+2], col = ScenCols[i+2], bg = ScenCols[i+2], cex = 1.5)
          segments(x0 = xLocs[FinalLocs[i],], y0 = FinalPlotData[[i+2]]$LwrPhen, 
                   x1 = xLocs[FinalLocs[i],], y1 = FinalPlotData[[i+2]]$UprPhen,
                   col = ScenCols[i+2], lty = 1)
     }
     # Put the letter and legend on the figure
     mtext("b", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
     legend("topleft", legend = Scenarios[3:5], pch = ScenPch[3:5], col = ScenCols[3:5], 
            bty = "n", pt.bg = ScenCols[3:5], inset = -0.01, ncol = 2, cex = 1.25)

     # Make the genetic variance figures: first non-mixing scenarios
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, GenMax), main = "", ylab = "", 
          xlab = "", las = 1, xaxt = "n", cex.axis = 1.5)
     axis(side = 1, at = LseqLocs, labels = Lseq, cex.axis = 1.5)
     axis(side = 2, at = seq(0, GenMax, by = 0.125), tcl = -0.25, labels = FALSE)
     mtext("Genetic variance", side = 2, line = 3, cex = 1.5)
     InitLocs <- c(2, 4)
     FinalLocs <- c(3, 5)
     for(i in 1:2){
          # Initial values
          segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i]]$LwrGen, 
                   x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i]]$UprGen,
                   col = ScenCols[i], lty = 1)
          points(x = xLocs[InitLocs[i],], y = InitPlotData[[i]]$MeanGen, 
                 pch = ScenPch[i], col = ScenCols[i], bg = "white", cex = 1.5)
          # Final values
          points(x = xLocs[FinalLocs[i],], y = FinalPlotData[[i]]$MeanGen, 
                 pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i], cex = 1.5)
          segments(x0 = xLocs[FinalLocs[i],], y0 = FinalPlotData[[i]]$LwrGen, 
                   x1 = xLocs[FinalLocs[i],], y1 = FinalPlotData[[i]]$UprGen,
                   col = ScenCols[i], lty = 1)
     }
     # Put the letter on the figure
     mtext("c", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)

     # Make the genetic variance figures for the mixing scenarios
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, GenMax), main = "", ylab = "", 
          xlab = "", las = 1, xaxt = "n", cex.axis = 1.5)
     axis(side = 1, at = LseqLocs, labels = Lseq, cex.axis = 1.5)
     axis(side = 2, at = seq(0, GenMax, by = 0.125), tcl = -0.25, labels = FALSE)
     InitLocs <- c(1,3,5)
     FinalLocs <- c(2,4,6)
     for(i in 1:3){
          # Initial values
          segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i+2]]$LwrGen, 
                   x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i+2]]$UprGen,
                   col = ScenCols[i+2], lty = 1)
          points(x = xLocs[InitLocs[i],], y = InitPlotData[[i+2]]$MeanGen, 
                 pch = ScenPch[i+2], col = ScenCols[i+2], bg = "white", cex = 1.5)
          # Final values
          points(x = xLocs[FinalLocs[i],], y = FinalPlotData[[i+2]]$MeanGen, 
                 pch = ScenPch[i+2], col = ScenCols[i+2], bg = ScenCols[i+2], cex = 1.5)
          segments(x0 = xLocs[FinalLocs[i],], y0 = FinalPlotData[[i+2]]$LwrGen, 
                   x1 = xLocs[FinalLocs[i],], y1 = FinalPlotData[[i+2]]$UprGen,
                   col = ScenCols[i+2], lty = 1)
     }
     # Put the letter and legend on the figure
     mtext("d", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)

     # Add the x axis label
     mtext("Number of loci", outer = TRUE, line = 1.5, side = 1, cex = 1.5)
dev.off()

######################### Now do the same thing but for the initial conditions
SummedResults <- read.csv("SimsWithSummedInitVals.csv")

offset <- 0.025
xLocs <- matrix(NA, nrow = 4, ncol = 6)
xLocs[1,] <- LseqLocs - 2*offset
xLocs[2,] <- LseqLocs - offset
xLocs[3,] <- LseqLocs + offset
xLocs[4,] <- LseqLocs + 2*offset
ScenPch <- c(21, 22, 23, 24, 25)

NewData <- vector(mode = "list", length = 2)
for(i in 1:2){
     NewMeanGen <- rep(NA, 6)
     NewLwrGen <- rep(NA, 6)
     NewUprGen <- rep(NA, 6)
     for(l in 1:6){
          CurRows <- which(SummedResults$L == Lseq[l] & SummedResults$Haploid == HapVals[i] & SummedResults$monoecious == monoVals[i] & SummedResults$omega == omegaVals[i])
          CurData <- SummedResults[CurRows,]
          if(i == 1 & l == 1){
               NewVar <- CurData$GenVar
          }else{
               NewVar <- CurData$SummedGenVar
          }
          NewMeanGen[l] <- mean(NewVar, na.rm = TRUE)
          NewLwrGen[l] <- quantile(NewVar, probs = 0.25, na.rm = TRUE)
          NewUprGen[l] <- quantile(NewVar, probs = 0.75, na.rm = TRUE)
     }
     NewData[[i]] <- data.frame(L = Lseq, MeanGen = NewMeanGen,
                                LwrGen = NewLwrGen, UprGen = NewUprGen)
}

pdf(file = "ResultFigures/FigureS2_new.pdf", width = 5, height = 4, onefile = FALSE, paper = "special", useDingbats = FALSE)
     # Make the genetic variance figures: first non-mixing scenarios
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, GenMax), main = "", ylab = "Variance in genotype", 
          xlab = "Number of loci", las = 1, xaxt = "n")
     axis(side = 1, at = LseqLocs, labels = Lseq)
     axis(side = 2, at = seq(0, GenMax, by = 0.125), tcl = -0.25, labels = FALSE)
     InitLocs <- c(2, 4)
     FinalLocs <- c(3, 5)
     for(i in 1:2){
          # Initial values
          segments(x0 = xLocs[i+1,], y0 = NewData[[i]]$LwrGen, 
                   x1 = xLocs[i+1,], y1 = NewData[[i]]$UprGen,
                   col = ScenCols[i], lty = 1)
          points(x = xLocs[i+1,], y = NewData[[i]]$MeanGen, 
                 pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i])
     }
     # Put the legend on the figure
     legend("topleft", legend = Scenarios[1:2], pch = ScenPch[1:2], col = ScenCols[1:2], 
            bty = "n", pt.bg = ScenCols[1:2], inset = -0.01, ncol = 1)
dev.off()


######################### Now plot the within genome variance for all scenarios
WithinResults <- read.csv("WithinGenVarResults.csv")

offset <- 0.015
xLocs <- matrix(NA, nrow = 6, ncol = 6)
xLocs[1,] <- LseqLocs - 3*offset
xLocs[2,] <- LseqLocs - 2*offset
xLocs[3,] <- LseqLocs - offset
xLocs[4,] <- LseqLocs + offset
xLocs[5,] <- LseqLocs + 2*offset
xLocs[6,] <- LseqLocs + 3*offset
ScenPch <- c(21, 22, 23, 24, 25)

InitPlotData <- vector(mode = "list", length = 5)
ExpPlotData <- vector(mode = "list", length = 5)
ShiftPlotData <- vector(mode = "list", length = 5)
for(i in 1:5){
        InitWithin <- array(NA, dim = c(3,6))
        ShiftWithin <- array(NA, dim = c(3,6))
        ExpWithin <- array(NA, dim = c(3,6))
        for(l in 1:6){
                WithinRows <- which(WithinResults$L == Lseq[l] & WithinResults$Haploid == HapVals[i] & WithinResults$monoecious == monoVals[i] & WithinResults$omega == omegaVals[i])
                CurData <- WithinResults[WithinRows,]
                InitWithin[1,l] <- mean(CurData$WithinInit)
                InitWithin[2:3,l] <- quantile(CurData$WithinInit, probs = c(0.25,0.75), na.rm = TRUE)
                ShiftWithin[1,l] <- mean(CurData$WithinShift, na.rm = TRUE)
                ShiftWithin[2:3,l] <- quantile(CurData$WithinShift, na.rm = TRUE, probs = c(0.25,0.75))
                # Combine the expansion values
                CurExpData <- c(CurData$WithinExp_1, CurData$WithinExp_2)
                ExpWithin[1,l] <- mean(CurExpData)
                ExpWithin[2:3,l] <- quantile(CurExpData, na.rm = TRUE, probs = c(0.25,0.75))
        }
        InitPlotData[[i]] <- data.frame(L = Lseq, MeanGen = InitWithin[1,],
                                        LwrGen = InitWithin[2,], UprGen = InitWithin[3,])
        ExpPlotData[[i]] <- data.frame(L = Lseq, MeanGen = ExpWithin[1,],
                                        LwrGen = ExpWithin[2,], UprGen = ExpWithin[3,])
        ShiftPlotData[[i]] <- data.frame(L = Lseq, MeanGen = ShiftWithin[1,],
                                        LwrGen = ShiftWithin[2,], UprGen = ShiftWithin[3,])
}

GenMax <- 25
pdf(file = "ResultFigures/FigureS3_new.pdf", width = 5, height = 4, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,2), oma = c(3,3,1,1), mar = c(2,1.5,2,1.5), family = "serif")

        # Make the expansion figures: first non-mixing scenarios
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, GenMax), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0, GenMax, by = 1), tcl = -0.25, labels = FALSE)
        InitLocs <- c(2, 4)
        FinalLocs <- c(3, 5)
        for(i in 1:2){
                # Initial values
                segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i]]$LwrGen, 
                         x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i]]$UprGen,
                         col = ScenCols[i], lty = 1)
                points(x = xLocs[InitLocs[i],], y = InitPlotData[[i]]$MeanGen, 
                       pch = ScenPch[i], col = ScenCols[i], bg = "white")
                # Final values
                points(x = xLocs[FinalLocs[i],], y = ExpPlotData[[i]]$MeanGen, 
                       pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[FinalLocs[i],], y0 = ExpPlotData[[i]]$LwrGen, 
                         x1 = xLocs[FinalLocs[i],], y1 = ExpPlotData[[i]]$UprGen,
                         col = ScenCols[i], lty = 1)
        }
        # Put the letter and legend on the figure
        mtext("a", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
        legend("topleft", legend = c(Scenarios[1:2], "Initial", "Final"), 
               pch = c(ScenPch[1:2], ScenPch[1], ScenPch[1]), 
               col = c(ScenCols[1:2], "black", "black"), bty = "n", 
               pt.bg = c(ScenCols[1:2], "white", "black"), ncol = 2)

        # Make the expansion figures for the mixing scenarios
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, GenMax), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0, GenMax, by = 1), tcl = -0.25, labels = FALSE)
        InitLocs <- c(1,3,5)
        FinalLocs <- c(2,4,6)
        for(i in 1:3){
                # Initial values
                segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i+2]]$LwrGen, 
                         x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i+2]]$UprGen,
                         col = ScenCols[i+2], lty = 1)
                points(x = xLocs[InitLocs[i],], y = InitPlotData[[i+2]]$MeanGen, 
                       pch = ScenPch[i+2], col = ScenCols[i+2], bg = "white")
                # Final values
                points(x = xLocs[FinalLocs[i],], y = ExpPlotData[[i+2]]$MeanGen, 
                       pch = ScenPch[i+2], col = ScenCols[i+2], bg = ScenCols[i+2])
                segments(x0 = xLocs[FinalLocs[i],], y0 = ExpPlotData[[i+2]]$LwrGen, 
                         x1 = xLocs[FinalLocs[i],], y1 = ExpPlotData[[i+2]]$UprGen,
                         col = ScenCols[i+2], lty = 1)
        }
        # Put the letter and legend on the figure
        mtext("b", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
        legend("topleft", legend = Scenarios[3:5], pch = ScenPch[3:5], col = ScenCols[3:5], 
               bty = "n", pt.bg = ScenCols[3:5], inset = -0.01)
        
        # Make the shift figures: first non-mixing scenarios
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, GenMax), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0, GenMax, by = 1), tcl = -0.25, labels = FALSE)
        InitLocs <- c(2, 4)
        FinalLocs <- c(3, 5)
        for(i in 1:2){
                # Initial values
                segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i]]$LwrGen, 
                         x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i]]$UprGen,
                         col = ScenCols[i], lty = 1)
                points(x = xLocs[InitLocs[i],], y = InitPlotData[[i]]$MeanGen, 
                       pch = ScenPch[i], col = ScenCols[i], bg = "white")
                # Final values
                points(x = xLocs[FinalLocs[i],], y = ShiftPlotData[[i]]$MeanGen, 
                       pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[FinalLocs[i],], y0 = ShiftPlotData[[i]]$LwrGen, 
                         x1 = xLocs[FinalLocs[i],], y1 = ShiftPlotData[[i]]$UprGen,
                         col = ScenCols[i], lty = 1)
        }
        # Put the letter on the figure
        mtext("c", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)

        # Make the shift figures for the mixing scenarios
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, GenMax), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0, GenMax, by = 1), tcl = -0.25, labels = FALSE)
        InitLocs <- c(1,3,5)
        FinalLocs <- c(2,4,6)
        for(i in 1:3){
                # Initial values
                segments(x0 = xLocs[InitLocs[i],], y0 = InitPlotData[[i+2]]$LwrGen, 
                         x1 = xLocs[InitLocs[i],], y1 = InitPlotData[[i+2]]$UprGen,
                         col = ScenCols[i+2], lty = 1)
                points(x = xLocs[InitLocs[i],], y = InitPlotData[[i+2]]$MeanGen, 
                       pch = ScenPch[i+2], col = ScenCols[i+2], bg = "white")
                # Final values
                points(x = xLocs[FinalLocs[i],], y = ShiftPlotData[[i+2]]$MeanGen, 
                       pch = ScenPch[i+2], col = ScenCols[i+2], bg = ScenCols[i+2])
                segments(x0 = xLocs[FinalLocs[i],], y0 = ShiftPlotData[[i+2]]$LwrGen, 
                         x1 = xLocs[FinalLocs[i],], y1 = ShiftPlotData[[i+2]]$UprGen,
                         col = ScenCols[i+2], lty = 1)
        }
        # Put the letter and legend on the figure
        mtext("d", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)

        # Add the x and y axes labels
        mtext("Number of loci", outer = TRUE, line = 1.5, side = 1)
        mtext("Within genome variance", outer = TRUE, line = 1.5, side = 2)
dev.off()


