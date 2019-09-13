# This script will create figures relating the evolutionary dynamics of dispersal
#    during range expansion between haploid and diploid dioecious populations, each
#    with dispersal defined by a varying number of loci. 

setwd("~/Desktop/PostdocResearch/DispersalEvolution/GitRepo/")
library(RColorBrewer)

#----------------------------------- First, load in both the haploid and diploid
# dioecious data objects, renaming them to avoid confusion.
# Haploid
load("SimData/2019-06-05_HaploidShift.rdata")
Haploid <- ShiftData
HapData <- SimData
# Diploid dioecious
load("SimData/2019-06-26_DiploidDioShift.rdata")
DipDioe <- ShiftData
DipData <- SimData
# Diploid monoecious, omega = 0
load("SimData/2019-07-18_DiploidMono-0-Shift.rdata")
DipMono_0 <- ShiftData
DipMono_0_Data <- SimData
# Diploid monoecious, omega = 0.5
load("SimData/2019-07-19_DiploidMono-05-Shift.rdata")
DipMono_05 <- ShiftData
DipMono_05_Data <- SimData
# Diploid monoecious, omega = 1
load("SimData/2019-07-23_DiploidMono-1-Shift.rdata")
DipMono_1 <- ShiftData
DipMono_1_Data <- SimData

# Remove the superfluous objects to free up space
rm(ShiftData, SimData, ParamCombos)

#------------------------------- Next, create some useful objects for the graphs
# Objects relating to the loci
Lseq <- unique(HapData$L)
LociLines <- 1:length(Lseq)
LociCols <- brewer.pal(n = length(Lseq), name = "Dark2")
LociLabels <- c("1 locus", paste(Lseq[2:6], "loci", sep = " "))
# Objects related to the number of generations
ShiftGens <- 200
Gens <- seq(0, ShiftGens, by = ShiftGens/10)
Ngens <- length(Gens)

#------------------------------- Now plot the extinction probability versus loci
xRange <- range(Lseq)
yRange <- c(0, 1)
AxisSize <- 2
LabSize <- 2
LegSize <- 2
xLabLine <- 3
yLabLine <- 3.5
PointSize <- 1.5
HapPoints <- 17
DipDioPoints <- 16
DipMonoPoints <- 18
AllCols <- brewer.pal(n = 5, name = "Dark2")
HapCol <- AllCols[1]
DipDioCol <- AllCols[2]
DipMono_0_Col <- AllCols[3]
DipMono_05_Col <- AllCols[4]
DipMono_1_Col <- AllCols[5]

# Load the ExtraShift objects to get slightly cleaned up versions of the extinction
#    probabilities
load("SimData/2019-06-26_DipDioExtraShift.rdata")
DipDioExtProbs <- ExtProb
load("SimData/2019-06-05_HaploidExtraShift.rdata")
HapExtProbs <- ExtProb
load("SimData/2019-07-18_DiploidMono-0-ExtraShift.rdata")
DipMono_0_ExtProbs <- ExtProb
load("SimData/2019-07-19_DiploidMono-05-ExtraShift.rdata")
DipMono_05_ExtProbs <- ExtProb
load("SimData/2019-07-23_DiploidMono-1-ExtraShift.rdata")
DipMono_1_ExtProbs <- ExtProb

#### Haploid
FigName <- "ResultFigures/HaploidExtinction.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = c(5,5,2,2) + 0.1, bg = "white")
     plot(NA, NA, xlim = xRange, ylim = yRange, main = "", xlab = "",
          ylab = "", las = 1, cex.axis = AxisSize)
     mtext("Number of Loci", side = 1, cex = LabSize, line = xLabLine)
     mtext("Extinction probability", side = 2, cex = LabSize, line = yLabLine)
     # Plot the results
     points(x = Lseq, y = HapExtProbs, pch = HapPoints, col = HapCol, cex = PointSize)
     # Create a legend
     legend("topright", bty = "n", legend = c("Haploid"), pch = c(HapPoints), 
            col = c(HapCol), cex = LegSize)
dev.off()

#### Plus Diploid dioecious
FigName <- "ResultFigures/PlusDipDioExtinction.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = c(5,5,2,2) + 0.1, bg = "white")
     plot(NA, NA, xlim = xRange, ylim = yRange, main = "", xlab = "",
          ylab = "", las = 1, cex.axis = AxisSize)
     mtext("Number of Loci", side = 1, cex = LabSize, line = xLabLine)
     mtext("Extinction probability", side = 2, cex = LabSize, line = yLabLine)
     # Plot the results
     points(x = Lseq, y = HapExtProbs, pch = HapPoints, col = HapCol, cex = PointSize)
     points(x = Lseq, y = DipDioExtProbs, pch = DipDioPoints, col = DipDioCol, cex = PointSize)
     # Create a legend
     legend("topright", bty = "n", legend = c("Haploid", "Diploid dioecious"), 
            pch = c(HapPoints, DipDioPoints), col = c(HapCol, DipDioCol), cex = LegSize)
dev.off()

#### Plus Diploid monoecious, omega = 1
FigName <- "ResultFigures/PlusDipMono_1_Extinction.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = c(5,5,2,2) + 0.1, bg = "white")
     plot(NA, NA, xlim = xRange, ylim = yRange, main = "", xlab = "",
          ylab = "", las = 1, cex.axis = AxisSize)
     mtext("Number of Loci", side = 1, cex = LabSize, line = xLabLine)
     mtext("Extinction probability", side = 2, cex = LabSize, line = yLabLine)
     # Plot the results
     points(x = Lseq, y = HapExtProbs, pch = HapPoints, col = HapCol, cex = PointSize)
     points(x = Lseq, y = DipDioExtProbs, pch = DipDioPoints, col = DipDioCol, cex = PointSize)
     points(x = Lseq, y = DipMono_1_ExtProbs, pch = DipMonoPoints, col = DipMono_1_Col, cex = PointSize)
     # Create a legend
     legend("topright", bty = "n", legend = c("Haploid", "Diploid dioecious", "Obligate selfing"), 
            pch = c(HapPoints, DipDioPoints, DipMonoPoints), col = c(HapCol, DipDioCol, DipMono_1_Col), 
            cex = LegSize)
dev.off()

#### Plus Diploid monoecious, omega = 0.5
FigName <- "ResultFigures/PlusDipMono_05_Extinction.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = c(5,5,2,2) + 0.1, bg = "white")
     plot(NA, NA, xlim = xRange, ylim = yRange, main = "", xlab = "",
          ylab = "", las = 1, cex.axis = AxisSize)
     mtext("Number of Loci", side = 1, cex = LabSize, line = xLabLine)
     mtext("Extinction probability", side = 2, cex = LabSize, line = yLabLine)
     # Plot the results
     points(x = Lseq, y = HapExtProbs, pch = HapPoints, col = HapCol, cex = PointSize)
     points(x = Lseq, y = DipDioExtProbs, pch = DipDioPoints, col = DipDioCol, cex = PointSize)
     points(x = Lseq, y = DipMono_05_ExtProbs, pch = DipMonoPoints, col = DipMono_05_Col, cex = PointSize)
     points(x = Lseq, y = DipMono_1_ExtProbs, pch = DipMonoPoints, col = DipMono_1_Col, cex = PointSize)
     # Create a legend
     legend("topright", bty = "n", legend = c("Haploid", "Diploid dioecious", 
                                              "Obligate selfing", "Partial selfing"), 
            pch = c(HapPoints, DipDioPoints, rep(DipMonoPoints, 2)), 
            col = c(HapCol, DipDioCol, DipMono_1_Col, DipMono_05_Col), cex = LegSize)
dev.off()

#### Plus Diploid monoecious, omega = 0
FigName <- "ResultFigures/AllExtinction.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = c(5,5,2,2) + 0.1, bg = "white")
     plot(NA, NA, xlim = xRange, ylim = yRange, main = "", xlab = "",
          ylab = "", las = 1, cex.axis = AxisSize)
     mtext("Number of Loci", side = 1, cex = LabSize, line = xLabLine)
     mtext("Extinction probability", side = 2, cex = LabSize, line = yLabLine)
     # Plot the results
     points(x = Lseq, y = HapExtProbs, pch = HapPoints, col = HapCol, cex = PointSize)
     points(x = Lseq, y = DipDioExtProbs, pch = DipDioPoints, col = DipDioCol, cex = PointSize)
     points(x = Lseq, y = DipMono_0_ExtProbs, pch = DipMonoPoints, col = DipMono_0_Col, cex = PointSize)
     points(x = Lseq, y = DipMono_05_ExtProbs, pch = DipMonoPoints, col = DipMono_05_Col, cex = PointSize)
     points(x = Lseq, y = DipMono_1_ExtProbs, pch = DipMonoPoints, col = DipMono_1_Col, cex = PointSize)
     # Create a legend
     legend("topright", bty = "n", legend = c("Haploid", "Diploid dioecious", 
                                              "Obligate selfing", "Partial selfing", "Minimal selfing"), 
            pch = c(HapPoints, DipDioPoints, rep(DipMonoPoints, 3)), 
            col = c(HapCol, DipDioCol, DipMono_1_Col, DipMono_05_Col, DipMono_0_Col), 
            cex = LegSize)
dev.off()

#### Empty
FigName <- "ResultFigures/EmptyExtinction.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = c(5,5,2,2) + 0.1, bg = "white")
     plot(NA, NA, xlim = xRange, ylim = yRange, main = "", xlab = "",
          ylab = "", las = 1, cex.axis = AxisSize)
     mtext("Number of Loci", side = 1, cex = LabSize, line = xLabLine)
     mtext("Extinction probability", side = 2, cex = LabSize, line = yLabLine)
dev.off()

#------------------------------------- Now plot the evolved changes through time
DeltaDbar <- expression(paste(Delta, bar(d)))
DeltaGenVar <- "Change in genetic variance"
xLabel <- "Generation"
xRange <- c(0, 200)
yRange_dBar <- c(0, 1.25)
yRange_GenVar <- c(-4, 0)
InnerMar <- c(5,5.5,2,2) + 0.1
yLabLine <- 3.5
xLabLine <- 3
AxisSize <- 2
LegLineWidth <- 1.5
LegSize <- 1.5
LabSize <- 2

## Haploid
FigName <- "ResultFigures/HaploidShiftDelta_dBar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaDbar, line = yLabLine, side = 2, cex = LabSize)
     
     # Plot the lines for each number of loci
     for(l in 1:length(Lseq)){
          CurData <- Haploid[[l]]
          lines(x = CurData$gen, y = CurData$Delta_dBar, lty = LociLines[l],
                col = LociCols[l])
     }

     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()

## Diploid dioecious
FigName <- "ResultFigures/DipDioShiftDelta_dBar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaDbar, line = yLabLine, side = 2, cex = LabSize)

     # Plot the lines for each number of loci
     for(l in 1:length(Lseq)){
          CurData <- DipDioe[[l]]
          lines(x = CurData$gen, y = CurData$Delta_dBar, lty = LociLines[l],
                col = LociCols[l])
     }

     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
     
## Diploid monoecious, omega = 0
FigName <- "ResultFigures/DipMono_0_ShiftDelta_dBar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaDbar, line = yLabLine, side = 2, cex = LabSize)

     # Plot the lines for each number of loci
     for(l in 1:length(Lseq)){
          CurData <- DipMono_0[[l]]
          lines(x = CurData$gen, y = CurData$Delta_dBar, lty = LociLines[l],
                col = LociCols[l])
     }

     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()

## Diploid monoecious, omega = 0.5
FigName <- "ResultFigures/DipMono_05_ShiftDelta_dBar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaDbar, line = yLabLine, side = 2, cex = LabSize)

     # Plot the lines for each number of loci
     for(l in 1:length(Lseq)){
          CurData <- DipMono_05[[l]]
          lines(x = CurData$gen, y = CurData$Delta_dBar, lty = LociLines[l],
                col = LociCols[l])
     }

     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()

## Diploid monoecious, omega = 1
FigName <- "ResultFigures/DipMono_1_ShiftDelta_dBar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaDbar, line = yLabLine, side = 2, cex = LabSize)

     # Plot the lines for each number of loci
     for(l in 1:length(Lseq)){
          CurData <- DipMono_1[[l]]
          lines(x = CurData$gen, y = CurData$Delta_dBar, lty = LociLines[l],
                col = LociCols[l])
     }

     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()

## Empty
FigName <- "ResultFigures/EmptyShiftDelta_dBar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaDbar, line = yLabLine, side = 2, cex = LabSize)

     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()


############# Now do the same but for genetic variance
## Haploid
FigName <- "ResultFigures/HaploidShiftDelta_GenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaGenVar, line = yLabLine, side = 2, cex = LabSize)

     # Plot the lines for each number of loci
     for(l in 1:length(Lseq)){
          CurData <- Haploid[[l]]
          lines(x = CurData$gen, y = CurData$DeltaGenVar, lty = LociLines[l],
                col = LociCols[l])
     }
     
     # Add a legend
     legend("bottomleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()

## Diploid dioecious
FigName <- "ResultFigures/DipDioShiftDelta_GenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaGenVar, line = yLabLine, side = 2, cex = LabSize)
          
     # Plot the lines for each number of loci
     for(l in 1:length(Lseq)){
          CurData <- DipDioe[[l]]
          lines(x = CurData$gen, y = CurData$DeltaGenVar, lty = LociLines[l],
                col = LociCols[l])
     }
     
     # Add a legend
     legend("bottomleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()

## Diploid monoecious, omega = 0
FigName <- "ResultFigures/DipMono_0_ShiftDelta_GenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaGenVar, line = yLabLine, side = 2, cex = LabSize)

     # Plot the lines for each number of loci
     for(l in 1:length(Lseq)){
          CurData <- DipMono_0[[l]]
          lines(x = CurData$gen, y = CurData$DeltaGenVar, lty = LociLines[l],
                col = LociCols[l])
     }

     # Add a legend
     legend("bottomleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()

## Diploid monoecious, omega = 0.5
FigName <- "ResultFigures/DipMono_05_ShiftDelta_GenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaGenVar, line = yLabLine, side = 2, cex = LabSize)

     # Plot the lines for each number of loci
     for(l in 1:length(Lseq)){
          CurData <- DipMono_05[[l]]
          lines(x = CurData$gen, y = CurData$DeltaGenVar, lty = LociLines[l],
                col = LociCols[l])
     }

     # Add a legend
     legend("bottomleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()

## Diploid monoecious, omega = 1
FigName <- "ResultFigures/DipMono_1_ShiftDelta_GenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaGenVar, line = yLabLine, side = 2, cex = LabSize)

     # Plot the lines for each number of loci
     for(l in 1:length(Lseq)){
          CurData <- DipMono_1[[l]]
          lines(x = CurData$gen, y = CurData$DeltaGenVar, lty = LociLines[l],
                col = LociCols[l])
     }

     # Add a legend
     legend("bottomleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()

## Empty
FigName <- "ResultFigures/EmptyShiftDelta_GenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = xLabel, line = xLabLine, side = 1, cex = LabSize)
     mtext(text = DeltaGenVar, line = yLabLine, side = 2, cex = LabSize)

     # Add a legend
     legend("bottomleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
