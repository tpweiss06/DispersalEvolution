# This script will make the necessary figures for the supplementary material
#    for the manuscript.

setwd("~/Desktop/GitHubRepos/DispersalEvolution/")

# Set up the same colors and labels as the main manuscript figures
library(RColorBrewer)
library(plot3D)
ScenCols <- brewer.pal(n = 5, name = "Dark2")
Scenarios <- c("Asexual", "Obligate selfing", "Partial selfing", "Obligate outcrossing", "Sexual (dioecious)")
HapVals <- c(TRUE, FALSE, FALSE, FALSE, FALSE)
monoVals <- c(TRUE, TRUE, TRUE, TRUE, FALSE)
omegaVals <- c(0, 1, 0.5, 0, 0)
Lseq <-c(1,2,4,8,16,32)
LseqLocs <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
offset <- 0.025
xLocs <- matrix(NA, nrow = 5, ncol = 6)
xLocs[1,] <- LseqLocs - 2*offset
xLocs[2,] <- LseqLocs - offset
xLocs[3,] <- LseqLocs
xLocs[4,] <- LseqLocs + offset
xLocs[5,] <- LseqLocs + 2*offset

################################# Fig. S1; Initial conditions of the simulations
# This figure will have the same format as the main figures in the manuscript
#       with two panels, each with loci number on the x-axis. The top panel
#       will have the mean phenotype with IQR for the simulations and the bottom
#       panel will have the mean (and IQR) genetic variance across simulations
InitResults <- read.csv("SimsWithInitVals.csv")
FigS1Data <- vector(mode = "list", length = 5)
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
        FigS1Data[[i]] <- data.frame(L = Lseq, MeanPhen, LwrPhen, UprPhen, MeanGen,
                                    LwrGen, UprGen)
}

pdf(file = "ResultFigures/FigureS1.pdf", width = 5, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(mfrow = c(2,1), oma = c(3,3,1,1), mar = c(1.5,1,1.5,1), family = "serif")

        # Make the phenotype figure
        plot(NA, NA, xlim = c(-0.05, 1.05), ylim = c(1.2, 1.9), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(1.2, 1.9, by = 0.025), tcl = -0.25, labels = FALSE)
        mtext("Initial dispersal phenotype", side = 2, line = 3)
        for(i in 1:5){
                points(x = xLocs[i,], y = FigS1Data[[i]]$MeanPhen, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = FigS1Data[[i]]$LwrPhen, x1 = xLocs[i,], y1 = FigS1Data[[i]]$UprPhen,
                         col = ScenCols[i], lty = 1)
        }
        mtext("a", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)

        # Put the legend on the figure
        legend("top", legend = Scenarios, pch = 21:25, col = ScenCols, bty = "n", 
               pt.bg = ScenCols, inset = -0.01, ncol = 2)
        
        # Make the genetic diversity figure
        plot(NA, NA, xlim = c(-0.05,1.05), ylim = c(0.5, 3.5), main = "", ylab = "", 
             xlab = "", las = 1, xaxt = "n")
        axis(side = 1, at = LseqLocs, labels = Lseq)
        axis(side = 2, at = seq(0.5, 3.5, by = 0.125), tcl = -0.25, labels = FALSE)
        mtext("Initial genetic diversity", side = 2, line = 3)
        for(i in 1:5){
                points(x = xLocs[i,], y = FigS1Data[[i]]$MeanGen, pch = 20 + i, col = ScenCols[i], bg = ScenCols[i])
                segments(x0 = xLocs[i,], y0 = FigS1Data[[i]]$LwrGen, x1 = xLocs[i,], y1 = FigS1Data[[i]]$UprGen,
                         col = ScenCols[i], lty = 1)
        }
        mtext("Number of loci", side = 1, outer = TRUE, line = 1.5)
        mtext("b", side = 3, line = 0, adj = 0, font = 2, cex = 1.25)
dev.off()

############### Fig. S2-S5; Sensitivity analysis of mutation parameters
# Fig. S2 and S3 will pertain to the changes during unbounded range expansions
# Fig. S4 and S5 will pertain to the changes during range shifts
# Fig. S2 and S4 will show the mean change for each parameter combination
# Fig. S3 and S5 will show the IQR width for changes under each parameter combination

# Before any graphing, the following loop will extract all the necessary data
Results <- read.csv("SensSimsWithResults.csv")
Uvals <- unique(Results$U)
sigVals <- unique(Results$sigma)
NumU <- length(Uvals)
NumSigma <- length(sigVals)
SensFigData <- array(data = NA, dim = c(4, length(HapVals), 3, NumU, NumSigma))
for(i in 1:5){
        for(j in 1:NumU){
                for(k in 1:NumSigma){
                        CurRows <- which(Results$Haploid == HapVals[i] & 
                                        Results$monoecious == monoVals[i] & 
                                        Results$omega == omegaVals[i] & 
                                        Results$U == Uvals[j] &
                                        Results$sigma == sigVals[k])
                        CurData <- Results[CurRows,]
                        # First extract the expansion data
                        DeltaPhenExp <- c(CurData$DeltaPhenExp_1, CurData$DeltaPhenExp_2)
                        DeltaGenExp <- c(CurData$DeltaGenExp_1, CurData$DeltaGenExp_2)
                        Dists <- c(CurData$distance_1, CurData$distance_2)
                        # Figure S2; mean values for unbounded expansions
                        SensFigData[1,i,1,j,k] <- mean(DeltaPhenExp)
                        SensFigData[1,i,2,j,k] <- mean(DeltaGenExp)
                        SensFigData[1,i,3,j,k] <- mean(Dists)
                        # Figure S3; IQR widths for unbounded expansions
                        SensFigData[2,i,1,j,k] <- quantile(DeltaPhenExp, probs = 0.75, na.rm = TRUE) - quantile(DeltaPhenExp, probs = 0.25, na.rm = TRUE)
                        SensFigData[2,i,2,j,k] <- quantile(DeltaGenExp, probs = 0.75, na.rm = TRUE) - quantile(DeltaGenExp, probs = 0.25, na.rm = TRUE)
                        SensFigData[2,i,3,j,k] <- quantile(Dists, probs = 0.75, na.rm = TRUE) - quantile(Dists, probs = 0.25, na.rm = TRUE)
                        # Figure S4; mean values for range shifts
                        SensFigData[3,i,1,j,k] <- mean(CurData$DeltaPhenShift)
                        SensFigData[3,i,2,j,k] <- mean(CurData$DeltaGenShift)
                        SensFigData[3,i,3,j,k] <- mean(CurData$ExtRisk)
                        # Figure S5; IQR widths for range shifts
                        SensFigData[4,i,1,j,k] <- quantile(CurData$DeltaPhenShift, probs = 0.75, na.rm = TRUE) - quantile(CurData$DeltaPhenShift, probs = 0.25, na.rm = TRUE)
                        SensFigData[4,i,2,j,k] <- quantile(CurData$DeltaGenShift, probs = 0.75, na.rm = TRUE) - quantile(CurData$DeltaGenShift, probs = 0.25, na.rm = TRUE)
                }
        }
}

# columns are y-axis, rows are x-axis

