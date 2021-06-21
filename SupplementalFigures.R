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
        mtext("Initial genetic variance", side = 2, line = 3)
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
                        SensFigData[3,i,1,j,k] <- mean(CurData$DeltaPhenShift, na.rm = TRUE)
                        SensFigData[3,i,2,j,k] <- mean(CurData$DeltaGenShift, na.rm = TRUE)
                        SensFigData[3,i,3,j,k] <- mean(CurData$ExtRisk)
                        # Figure S5; IQR widths for range shifts
                        SensFigData[4,i,1,j,k] <- quantile(CurData$DeltaPhenShift, probs = 0.75, na.rm = TRUE) - quantile(CurData$DeltaPhenShift, probs = 0.25, na.rm = TRUE)
                        SensFigData[4,i,2,j,k] <- quantile(CurData$DeltaGenShift, probs = 0.75, na.rm = TRUE) - quantile(CurData$DeltaGenShift, probs = 0.25, na.rm = TRUE)
                }
        }
}

####### For the supplemental figures using heat maps (all sensitivity analysis figures),
#       create a function to identify the appropriate colors to use for a single graph
#       from a larger selection that correspond to the full range of possible values.
FindRange <- function(minimum, maximum, sequence){
        if(class(sequence) != 'numeric'){
                sequence <- as.numeric(sequence)
        }
        LessThan <- ifelse(sequence <= minimum, 1, 0)
        MoreThan <- ifelse(sequence >= maximum, 1, 0)
        MinIndex <- sum(LessThan)
        MaxIndex <- length(sequence) - sum(MoreThan)
        return(MinIndex:MaxIndex)
}

# Create a matrix to use with the layout command for Figures S2-S4
PlotArray <- matrix(c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), 6, rep(7, 4), rep(8, 4), rep(9, 4), rep(10, 4), rep(11, 4), 12, rep(13, 4), rep(14, 4), rep(15, 4), rep(16, 4), rep(17, 4), 18), 
                    nrow=3, ncol=21, byrow=TRUE)
############################ Fig. S2; Mean changes in unbounded range expansions
# Define the ranges for each row of the plot and appropriate color values
ColRanges <- matrix(data = NA, nrow = 3, ncol = 2)
ColVals <- array(data = NA, dim = c(3,2,1000))
for(i in 1:3){
        ColRanges[i,] <- range(SensFigData[1,,i,,], na.rm = TRUE)
        ColVals[i,1,] <- seq(ColRanges[i,1], ColRanges[i,2], length.out=1000)
        ColVals[i,2,] <- jet.col(n = 1000)
}
RowLabels <- c("Phenotype", "Genetic Variance", "Distance")
pdf(file = "ResultFigures/FigureS2.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(oma = c(3,4,1,2), mar = c(1.5,1.5,1.5,1.5), family = "serif")
        layout(PlotArray)
        for(r in 1:3){
                for(i in 1:5){
                        CurRange <- FindRange(minimum = min(SensFigData[1,i,r,,], na.rm = TRUE), 
                                              maximum = max(SensFigData[1,i,r,,], na.rm = TRUE),
                                              sequence = ColVals[r,1,])
                        image2D(z = SensFigData[1,i,r,,], x = 1:length(Uvals), 
                                y = 1:length(sigVals), resfac = 1, colkey=FALSE,
                                col = ColVals[r,2,CurRange], xlab = "", ylab = "",
                                axes = FALSE)
                        box()
                        axis(1, labels = Uvals, at = 1:7, las=1)
                        axis(2, labels = round(sigVals, digits = 2), at = 1:7, las=1)
                        if(r == 1){
                                mtext(Scenarios[i], side = 3, line = 0.5)
                        }
                        if(i == 1){
                                mtext(RowLabels[r], side = 2, line = 4)
                        }
                }
                # Add the color key
                colkey(col = ColVals[r,2,], clim = ColRanges[r,], cex.axis = 1, width = 15)
        }
        mtext("U", side = 1, line = 1.5, outer = TRUE)
        mtext(expression(sigma), side = 2, line = 1.25, outer = TRUE)
dev.off()

############################ Fig. S3; IQR widths in unbounded range expansions
# Define the ranges for each row of the plot and appropriate color values
ColRanges <- matrix(data = NA, nrow = 3, ncol = 2)
ColVals <- array(data = NA, dim = c(3,2,1000))
for(i in 1:3){
        ColRanges[i,] <- range(SensFigData[2,,i,,], na.rm = TRUE)
        ColVals[i,1,] <- seq(ColRanges[i,1], ColRanges[i,2], length.out=1000)
        ColVals[i,2,] <- jet.col(n = 1000)
}
RowLabels <- c("Phenotype", "Genetic Variance", "Distance")
pdf(file = "ResultFigures/FigureS3.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(oma = c(3,4,1,2), mar = c(1.5,1.5,1.5,1.5), family = "serif")
        layout(PlotArray)
        for(r in 1:3){
                for(i in 1:5){
                        CurRange <- FindRange(minimum = min(SensFigData[2,i,r,,], na.rm = TRUE), 
                                              maximum = max(SensFigData[2,i,r,,], na.rm = TRUE),
                                              sequence = ColVals[r,1,])
                        image2D(z = SensFigData[2,i,r,,], x = 1:length(Uvals), 
                                y = 1:length(sigVals), resfac = 1, colkey=FALSE,
                                col = ColVals[r,2,CurRange], xlab = "", ylab = "",
                                axes = FALSE)
                        box()
                        axis(1, labels = Uvals, at = 1:7, las=1)
                        axis(2, labels = round(sigVals, digits = 2), at = 1:7, las=1)
                        if(r == 1){
                                mtext(Scenarios[i], side = 3, line = 0.5)
                        }
                        if(i == 1){
                                mtext(RowLabels[r], side = 2, line = 4)
                        }
                }
                # Add the color key
                colkey(col = ColVals[r,2,], clim = ColRanges[r,], cex.axis = 1, width = 15)
        }
        mtext("U", side = 1, line = 1.5, outer = TRUE)
        mtext(expression(sigma), side = 2, line = 1.25, outer = TRUE)
dev.off()

########################################## Fig. S4; Mean changes in range shifts
# Define the ranges for each row of the plot and appropriate color values
ColRanges <- matrix(data = NA, nrow = 3, ncol = 2)
ColVals <- array(data = NA, dim = c(3,2,1000))
for(i in 1:3){
        ColRanges[i,] <- range(SensFigData[3,,i,,], na.rm = TRUE)
        ColVals[i,1,] <- seq(ColRanges[i,1], ColRanges[i,2], length.out=1000)
        ColVals[i,2,] <- jet.col(n = 1000)
}
RowLabels <- c("Phenotype", "Genetic Variance", "Extinction risk")
pdf(file = "ResultFigures/FigureS4.pdf", width = 10, height = 6, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(oma = c(3,4,1,2), mar = c(1.5,1.5,1.5,1.5), family = "serif")
        layout(PlotArray)
        for(r in 1:3){
                for(i in 1:5){
                        CurRange <- FindRange(minimum = min(SensFigData[3,i,r,,], na.rm = TRUE), 
                                              maximum = max(SensFigData[3,i,r,,], na.rm = TRUE),
                                              sequence = ColVals[r,1,])
                        image2D(z = SensFigData[3,i,r,,], x = 1:length(Uvals), 
                                y = 1:length(sigVals), resfac = 1, colkey=FALSE,
                                col = ColVals[r,2,CurRange], xlab = "", ylab = "",
                                axes = FALSE)
                        box()
                        axis(1, labels = Uvals, at = 1:7, las=1)
                        axis(2, labels = round(sigVals, digits = 2), at = 1:7, las=1)
                        if(r == 1){
                                mtext(Scenarios[i], side = 3, line = 0.5)
                        }
                        if(i == 1){
                                mtext(RowLabels[r], side = 2, line = 4)
                        }
                }
                # Add the color key
                colkey(col = ColVals[r,2,], clim = ColRanges[r,], cex.axis = 1, width = 15)
        }
        mtext("U", side = 1, line = 1.5, outer = TRUE)
        mtext(expression(sigma), side = 2, line = 1.25, outer = TRUE)
dev.off()

########################################## Fig. S5; IQR widths in range shifts
# This figure only has 2 rows since the proportion of scenarios to go extinct 
#       doesn't have an associated IQR
# Create a new matrix to use with the layout command for Figures S5
PlotArray <- matrix(c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), 6, rep(7, 4), rep(8, 4), rep(9, 4), rep(10, 4), rep(11, 4), 12), 
                    nrow=2, ncol=21, byrow=TRUE)
# Define the ranges for each row of the plot and appropriate color values
ColRanges <- matrix(data = NA, nrow = 2, ncol = 2)
ColVals <- array(data = NA, dim = c(2,2,1000))
for(i in 1:2){
        ColRanges[i,] <- range(SensFigData[4,,i,,], na.rm = TRUE)
        ColVals[i,1,] <- seq(ColRanges[i,1], ColRanges[i,2], length.out=1000)
        ColVals[i,2,] <- jet.col(n = 1000)
}
RowLabels <- c("Phenotype", "Genetic Variance")
pdf(file = "ResultFigures/FigureS5.pdf", width = 10, height = 4, onefile = FALSE, paper = "special", useDingbats = FALSE)
        par(oma = c(3,4,1,2), mar = c(1.5,1.5,1.5,1.5), family = "serif")
                layout(PlotArray)
                for(r in 1:2){
                        for(i in 1:5){
                                CurRange <- FindRange(minimum = min(SensFigData[4,i,r,,], na.rm = TRUE), 
                                                      maximum = max(SensFigData[4,i,r,,], na.rm = TRUE),
                                                      sequence = ColVals[r,1,])
                                image2D(z = SensFigData[4,i,r,,], x = 1:length(Uvals), 
                                        y = 1:length(sigVals), resfac = 1, colkey=FALSE,
                                        col = ColVals[r,2,CurRange], xlab = "", ylab = "",
                                        axes = FALSE)
                                box()
                                axis(1, labels = Uvals, at = 1:7, las=1)
                                axis(2, labels = round(sigVals, digits = 2), at = 1:7, las=1)
                                if(r == 1){
                                        mtext(Scenarios[i], side = 3, line = 0.5)
                                }
                                if(i == 1){
                                        mtext(RowLabels[r], side = 2, line = 4)
                                }
                        }
                        # Add the color key
                        colkey(col = ColVals[r,2,], clim = ColRanges[r,], cex.axis = 1, width = 15)
                }
                mtext("U", side = 1, line = 1.5, outer = TRUE)
                mtext(expression(sigma), side = 2, line = 1.25, outer = TRUE)
dev.off()

