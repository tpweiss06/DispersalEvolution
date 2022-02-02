# This script will explore the difference in genetic variance calculated
#    on a loci basis versus as the summed allele values for asexual scenarios

setwd("~/Desktop/GitHubRepos/DispersalEvolution/")
SummedResults <- read.csv("SimsWithSummedResults.csv")

# Make a two panel graph with points for each loci number for Haploid (right panel)
#    and obligately selfing (left panel)
# Set up the same colors and labels as the main manuscript figures
library(RColorBrewer)
ScenCols <- brewer.pal(n = 5, name = "Dark2")
Scenarios <- c("Asexual", "Obligate selfing", "Partial selfing", "Obligate outcrossing", "Sexual (dioecious)")
HapVals <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
monoVals <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
omegaVals <- c(0, 1, 0.5, 0, 0)
Lseq <-c(1,2,4,8,16,32)
LseqLocs <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
offset <- 0.025
xLocs <- matrix(NA, nrow = 4, ncol = 6)
xLocs[1,] <- LseqLocs - 2*offset
xLocs[2,] <- LseqLocs - offset
xLocs[3,] <- LseqLocs + offset
xLocs[4,] <- LseqLocs + 2*offset
ScenPch <- c(21, 22, 23, 24, 25)

OldExpData <- vector(mode = "list", length = 2)
NewExpData <- vector(mode = "list", length = 2)
OldShiftData <- vector(mode = "list", length = 2)
NewShiftData <- vector(mode = "list", length = 2)
for(i in 1:2){
     ExpOldMeanGen <- rep(NA, 6)
     ExpOldLwrGen <- rep(NA, 6)
     ExpOldUprGen <- rep(NA, 6)
     ExpNewMeanGen <- rep(NA, 6)
     ExpNewLwrGen <- rep(NA, 6)
     ExpNewUprGen <- rep(NA, 6)
     
     ShiftOldMeanGen <- rep(NA, 6)
     ShiftOldLwrGen <- rep(NA, 6)
     ShiftOldUprGen <- rep(NA, 6)
     ShiftNewMeanGen <- rep(NA, 6)
     ShiftNewLwrGen <- rep(NA, 6)
     ShiftNewUprGen <- rep(NA, 6)
     for(l in 1:6){
          CurRows <- which(SummedResults$L == Lseq[l] & SummedResults$Haploid == HapVals[i] & SummedResults$monoecious == monoVals[i] & SummedResults$omega == omegaVals[i])
          CurData <- SummedResults[CurRows,]
          OldExp <- c(CurData$DeltaGenExp_1, CurData$DeltaGenExp_2)
          OldShift <- CurData$DeltaGenShift
          NewExp <- c(CurData$SummedDeltaGenExp_1, CurData$SummedDeltaGenExp_2)
          NewShift <- CurData$SummedDeltaGenShift
          ExpOldMeanGen[l] <- mean(OldExp)
          ExpOldLwrGen[l] <- quantile(OldExp, probs = 0.25)
          ExpOldUprGen[l] <- quantile(OldExp, probs = 0.75)
          ExpNewMeanGen[l] <- mean(NewExp, na.rm = TRUE)
          ExpNewLwrGen[l] <- quantile(NewExp, probs = 0.25, na.rm = TRUE)
          ExpNewUprGen[l] <- quantile(NewExp, probs = 0.75, na.rm = TRUE)
          
          ShiftOldMeanGen[l] <- mean(OldShift, na.rm = TRUE)
          ShiftOldLwrGen[l] <- quantile(OldShift, probs = 0.25, na.rm = TRUE)
          ShiftOldUprGen[l] <- quantile(OldShift, probs = 0.75, na.rm = TRUE)
          ShiftNewMeanGen[l] <- mean(NewShift, na.rm = TRUE)
          ShiftNewLwrGen[l] <- quantile(NewShift, probs = 0.25, na.rm = TRUE)
          ShiftNewUprGen[l] <- quantile(NewShift, probs = 0.75, na.rm = TRUE)
     }
     OldExpData[[i]] <- data.frame(L = Lseq, MeanGen = ExpOldMeanGen,
                                  LwrGen = ExpOldLwrGen, UprGen = ExpOldUprGen)
     NewExpData[[i]] <- data.frame(L = Lseq, MeanGen = ExpNewMeanGen,
                                   LwrGen = ExpNewLwrGen, UprGen = ExpNewUprGen)
     OldShiftData[[i]] <- data.frame(L = Lseq, MeanGen = ShiftOldMeanGen,
                                   LwrGen = ShiftOldLwrGen, UprGen = ShiftOldUprGen)
     NewShiftData[[i]] <- data.frame(L = Lseq, MeanGen = ShiftNewMeanGen,
                                   LwrGen = ShiftNewLwrGen, UprGen = ShiftNewUprGen)
}

pdf(file = "ResultFigures/GenVarExpCompare.pdf", width = 5, height = 4, onefile = FALSE, paper = "special", useDingbats = FALSE)
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(-3, 0), main = "", ylab = "Genetic variance", 
          xlab = "Number of loci", las = 1, xaxt = "n")
     axis(side = 1, at = LseqLocs, labels = Lseq)
     axis(side = 2, at = seq(-3, 0, by = 0.25), tcl = -0.25, labels = FALSE)
     OldLocs <- c(1,3)
     NewLocs <- c(2,4)
     for(i in 1:2){
          # Old
          if(i == 1){
               points(x = xLocs[OldLocs[i],2:6], y = OldExpData[[i]]$MeanGen[2:6], pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i])
               segments(x0 = xLocs[OldLocs[i],2:6], y0 = OldExpData[[i]]$LwrGen[2:6], x1 = xLocs[OldLocs[i],2:6], y1 = OldExpData[[i]]$UprGen[2:6],
                        col = ScenCols[i], lty = 1)
          }else{
               points(x = xLocs[OldLocs[i],], y = OldExpData[[i]]$MeanGen, pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i])
               segments(x0 = xLocs[OldLocs[i],], y0 = OldExpData[[i]]$LwrGen, x1 = xLocs[OldLocs[i],], y1 = OldExpData[[i]]$UprGen,
                        col = ScenCols[i], lty = 1)
          }
          
          # New
          segments(x0 = xLocs[NewLocs[i],], y0 = NewExpData[[i]]$LwrGen, x1 = xLocs[NewLocs[i],], y1 = NewExpData[[i]]$UprGen,
                   col = ScenCols[i], lty = 1)
          points(x = xLocs[NewLocs[i],], y = NewExpData[[i]]$MeanGen, pch = ScenPch[i], col = ScenCols[i], bg = "white")
     }

     # Put the legend on the figure
     legend("top", legend = c(Scenarios[1:2], "Allele variance", "Genotype variance"), 
            pch = c(ScenPch[1:2], ScenPch[1], ScenPch[1]), 
            col = c(ScenCols[1:2], "black", "black"), bty = "n", 
            pt.bg = c(ScenCols[1:2], "black", "white"), inset = -0.01, ncol = 2)
dev.off()

pdf(file = "ResultFigures/GenVarShiftCompare.pdf", width = 5, height = 4, onefile = FALSE, paper = "special", useDingbats = FALSE)
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(-3, 0), main = "", ylab = "Genetic variance", 
          xlab = "Number of loci", las = 1, xaxt = "n")
     axis(side = 1, at = LseqLocs, labels = Lseq)
     axis(side = 2, at = seq(-3, 0, by = 0.25), tcl = -0.25, labels = FALSE)
     OldLocs <- c(1,3)
     NewLocs <- c(2,4)
     for(i in 1:2){
          # Old
          if(i == 1){
               points(x = xLocs[OldLocs[i],2:6], y = OldShiftData[[i]]$MeanGen[2:6], pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i])
               segments(x0 = xLocs[OldLocs[i],2:6], y0 = OldShiftData[[i]]$LwrGen[2:6], x1 = xLocs[OldLocs[i],2:6], y1 = OldShiftData[[i]]$UprGen[2:6],
                        col = ScenCols[i], lty = 1)
          }else{
               points(x = xLocs[OldLocs[i],], y = OldShiftData[[i]]$MeanGen, pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i])
               segments(x0 = xLocs[OldLocs[i],], y0 = OldShiftData[[i]]$LwrGen, x1 = xLocs[OldLocs[i],], y1 = OldShiftData[[i]]$UprGen,
                        col = ScenCols[i], lty = 1)
          }
          
          # New
          segments(x0 = xLocs[NewLocs[i],], y0 = NewShiftData[[i]]$LwrGen, x1 = xLocs[NewLocs[i],], y1 = NewShiftData[[i]]$UprGen,
                   col = ScenCols[i], lty = 1)
          points(x = xLocs[NewLocs[i],], y = NewShiftData[[i]]$MeanGen, pch = ScenPch[i], col = ScenCols[i], bg = "white")
     }

     # Put the legend on the figure
     legend("top", legend = c(Scenarios[1:2], "Allele variance", "Genotype variance"), 
            pch = c(ScenPch[1:2], ScenPch[1], ScenPch[1]), 
            col = c(ScenCols[1:2], "black", "black"), bty = "n", 
            pt.bg = c(ScenCols[1:2], "black", "white"), inset = -0.01, ncol = 2)
dev.off()


######################### Now do the same thing but for the initial conditions
setwd("~/Desktop/GitHubRepos/DispersalEvolution/")
SummedResults <- read.csv("SimsWithSummedInitVals.csv")

# Make a two panel graph with points for each loci number for Haploid (right panel)
#    and obligately selfing (left panel)
# Set up the same colors and labels as the main manuscript figures
library(RColorBrewer)
ScenCols <- brewer.pal(n = 5, name = "Dark2")
Scenarios <- c("Asexual", "Obligate selfing", "Partial selfing", "Obligate outcrossing", "Sexual (dioecious)")
HapVals <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
monoVals <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
omegaVals <- c(0, 1, 0.5, 0, 0)
Lseq <-c(1,2,4,8,16,32)
LseqLocs <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
offset <- 0.025
xLocs <- matrix(NA, nrow = 4, ncol = 6)
xLocs[1,] <- LseqLocs - 2*offset
xLocs[2,] <- LseqLocs - offset
xLocs[3,] <- LseqLocs + offset
xLocs[4,] <- LseqLocs + 2*offset
ScenPch <- c(21, 22, 23, 24, 25)

OldData <- vector(mode = "list", length = 2)
NewData <- vector(mode = "list", length = 2)
for(i in 1:2){
     OldMeanGen <- rep(NA, 6)
     OldLwrGen <- rep(NA, 6)
     OldUprGen <- rep(NA, 6)
     NewMeanGen <- rep(NA, 6)
     NewLwrGen <- rep(NA, 6)
     NewUprGen <- rep(NA, 6)
     
     for(l in 1:6){
          CurRows <- which(SummedResults$L == Lseq[l] & SummedResults$Haploid == HapVals[i] & SummedResults$monoecious == monoVals[i] & SummedResults$omega == omegaVals[i])
          CurData <- SummedResults[CurRows,]
          OldVar <- CurData$GenVar
          NewVar <- CurData$SummedGenVar
          
          OldMeanGen[l] <- mean(OldVar)
          OldLwrGen[l] <- quantile(OldVar, probs = 0.25)
          OldUprGen[l] <- quantile(OldVar, probs = 0.75)
          NewMeanGen[l] <- mean(NewVar, na.rm = TRUE)
          NewLwrGen[l] <- quantile(NewVar, probs = 0.25, na.rm = TRUE)
          NewUprGen[l] <- quantile(NewVar, probs = 0.75, na.rm = TRUE)
     }
     OldData[[i]] <- data.frame(L = Lseq, MeanGen = OldMeanGen,
                                   LwrGen = OldLwrGen, UprGen = OldUprGen)
     NewData[[i]] <- data.frame(L = Lseq, MeanGen = NewMeanGen,
                                   LwrGen = NewLwrGen, UprGen = NewUprGen)
}

pdf(file = "ResultFigures/InitGenVarExpCompare.pdf", width = 5, height = 4, onefile = FALSE, paper = "special", useDingbats = FALSE)
     plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(0, 3.5), main = "", ylab = "Genetic variance", 
          xlab = "Number of loci", las = 1, xaxt = "n")
     axis(side = 1, at = LseqLocs, labels = Lseq)
     axis(side = 2, at = seq(0, 3.5, by = 0.125), tcl = -0.25, labels = FALSE)
     OldLocs <- c(1,3)
     NewLocs <- c(2,4)
     for(i in 1:2){
          # Old
          if(i == 1){
               points(x = xLocs[OldLocs[i],2:6], y = OldData[[i]]$MeanGen[2:6], pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i])
               segments(x0 = xLocs[OldLocs[i],2:6], y0 = OldData[[i]]$LwrGen[2:6], x1 = xLocs[OldLocs[i],2:6], y1 = OldData[[i]]$UprGen[2:6],
                        col = ScenCols[i], lty = 1)
          }else{
               points(x = xLocs[OldLocs[i],], y = OldData[[i]]$MeanGen, pch = ScenPch[i], col = ScenCols[i], bg = ScenCols[i])
               segments(x0 = xLocs[OldLocs[i],], y0 = OldData[[i]]$LwrGen, x1 = xLocs[OldLocs[i],], y1 = OldData[[i]]$UprGen,
                        col = ScenCols[i], lty = 1)
          }
          
          # New
          segments(x0 = xLocs[NewLocs[i],], y0 = NewData[[i]]$LwrGen, x1 = xLocs[NewLocs[i],], y1 = NewData[[i]]$UprGen,
                   col = ScenCols[i], lty = 1)
          points(x = xLocs[NewLocs[i],], y = NewData[[i]]$MeanGen, pch = ScenPch[i], col = ScenCols[i], bg = "white")
     }

     # Put the legend on the figure
     legend("top", legend = c(Scenarios[1:2], "Allele variance", "Genotype variance"), 
            pch = c(ScenPch[1:2], ScenPch[1], ScenPch[1]), 
            col = c(ScenCols[1:2], "black", "black"), bty = "n", 
            pt.bg = c(ScenCols[1:2], "black", "white"), inset = -0.01, ncol = 2)
dev.off()

