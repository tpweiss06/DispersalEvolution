# This script will create figures relating the evolutionary dynamics of dispersal
#    during range expansion between haploid and diploid dioecious populations, each
#    with dispersal defined by a varying number of loci. 

setwd("~/Desktop/PostdocResearch/DispersalEvolution/GitRepo/")
library(RColorBrewer)

#----------------------------------- First, load in both the haploid and diploid
# dioecious data objects, renaming them to avoid confusion.
# Haploid
load("SimData/2019-06-05_HaploidExpand.rdata")
Haploid <- DispData
HapData <- SimData
# Diploid dioecious
load("SimData/2019-06-05_DipDioExpand.rdata")
DipDioe <- DispData
DipData <- SimData
# Diploid monoecious, omega = 0
load("SimData/2019-07-18_DiploidMono-0-Expand.rdata")
DipMono_0 <- DispData
DipMono_0_Data <- SimData
# Diploid monoecious, omega = 0.5
load("SimData/2019-07-19_DiploidMono-05-Expand.rdata")
DipMono_05 <- DispData
DipMono_05_Data <- SimData
# Diploid monoecious, omega = 1
load("SimData/2019-07-23_DiploidMono-1-Expand.rdata")
DipMono_1 <- DispData
DipMono_1_Data <- SimData

# Now remove the superfluous objects to free up space
rm(DispData, SimData, ParamCombos)

#------------------------------- Next, create some useful objects for the graphs
# Objects relating to the loci
Lseq <- unique(HapData$L)
LociLines <- 1:length(Lseq)
LociCols <- brewer.pal(n = length(Lseq), name = "Dark2")
LociLabels <- c("1 locus", paste(Lseq[2:6], "loci", sep = " "))
# Objects related to the number of generations
ExpandGens <- 200
Gens <- seq(0, ExpandGens, by = ExpandGens/10)
Ngens <- length(Gens)
# Objects controlling the axes ranges
xRange <- c(-1000, 1000)
yRange_dBar <- c(1, 3.5)
yRange_GenVar <- c(0, 4)
# Objects for labeling the axes
dBarLab <- expression(bar(d))
GenVarLab <- "Genetic variance"
# Size objects
AxisSize <- 2
LegSize <- 2
LegLineWidth <- 2
LabSize <- 2
InnerMar <- c(5,5,2,2)

#------------------------------ Now create the animated Gif of dBar through time
# # Set the size of the png image
# FigWidth <- 1200
# FigHeight <- 900
# FigPointSize <- 18
# # Create objects to control the sizing of things in these plots
# TitleLine <- 2
# TitleSize <- 1.5
# PlotTitleLine <- 0.5
# PlotTitleSize <- 1.5
# yLab_x <- -1400
# yLab_y <- 2.5
# yLabSize <- 1.5
# xLabLine <- 3
# xLabSize <- 1.5
# AxisSize <- 1.25
# LegSize <- 1.5
# LegLineWidth <- 2
# OuterMar <- c(5, 4, 4, 2)
# InnerMar <- c(0, 1, 0, 1)
# # Create a temporary directory to store the sill images
# TempDirPath <- "ResultFigures/temp"
# dir.create(TempDirPath)
# # Now, loop through each generation making the still image for that generation
# for(t in 1:Ngens){
#      gen <- Gens[t]
#      if(gen == 0){
#           GenNum <- "000"
#      } else if(gen < 100){
#           GenNum <- paste("0", gen, sep = "")
#      } else{
#           GenNum <- gen
#      }
#      FigName <- paste(TempDirPath, "/gen", GenNum, ".png", sep = "")
#      FigTitle <- paste("Generation", gen, sep = " ")
#      png(filename = FigName, width = FigWidth, height = FigHeight, pointsize = FigPointSize)
#           par(mfrow = c(1,2), mar = InnerMar, oma = OuterMar)
#           # Make the Haploid plot first
#           plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = "", 
#                xlab = "", las = 1, cex.axis = AxisSize)
#           # Plot the lines for each loci
#           for(l in 1:length(Lseq)){
#                CurData <- Haploid[[l]]$dBar[[t]]
#                lines(x = CurData$x, y = CurData$dBar, lty = LociLines[l],
#                      col = LociCols[l])
#           }
#           # Add a legend
#           legend("top", bty = "n", legend = LociLabels, 
#                  lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
#           # Add the plot title and y axis label
#           mtext(side = 3, text = "Haploid", line = PlotTitleLine, cex = PlotTitleSize)
#           text(labels = dBarLab, x = yLab_x, y = yLab_y, cex = xLabSize, srt = 0, xpd = NA)
#           
#           # Now add the Diploid dioecious plot
#           plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = "", 
#                xlab = "", las = 1, yaxt = "n", cex.axis = AxisSize)
#           axis(2, labels = FALSE)
#           # Plot the lines for each loci
#           for(l in 1:length(Lseq)){
#                CurData <- DipDioe[[l]]$dBar[[t]]
#                lines(x = CurData$x, y = CurData$dBar, lty = LociLines[l],
#                      col = LociCols[l])
#           }
#           # Add the plot title
#           mtext(side = 3, text = "Diploid dioecious", line = PlotTitleLine, cex = PlotTitleSize)
#           
#           # Add the other figure text (Generation and axes labels)
#           mtext(text = FigTitle, side = 3, outer = TRUE, line = TitleLine, cex = TitleSize)
#           mtext(text = "Space", side = 1, outer = TRUE, line = xLabLine, cex = xLabSize)
#           #mtext(text = dBarLab, side = 2, outer = TRUE, line = yLabLine, cex = yLabSize)
#           
#      dev.off()
# }
# # Now make the animation
# GifName <- "ResultFigures/dBarEvolution.gif"
# SysCommand1 <- paste("convert -delay 50 ", TempDirPath, "/gen*.png ", GifName, sep = "")
# SysCommand2 <- paste("rm -r ", TempDirPath, sep = "")
# system(SysCommand1)
# system(SysCommand2)

#------------------ Now create the animated Gif of genetic variance through time
# # Set the size of the png image
# FigWidth <- 1200
# FigHeight <- 900
# FigPointSize <- 18
# # Create objects to control the sizing of things in these plots
# TitleLine <- 2
# TitleSize <- 1.5
# PlotTitleLine <- 0.5
# PlotTitleSize <- 1.5
# yLab_x <- -1400
# yLab_y <- 2.5
# yLabSize <- 1.5
# xLabLine <- 3
# xLabSize <- 1.5
# AxisSize <- 1.25
# LegSize <- 1.5
# LegLineWidth <- 2
# OuterMar <- c(5, 4, 4, 2)
# InnerMar <- c(0, 1, 0, 1)
# # Create a temporary directory to store the sill images
# TempDirPath <- "ResultFigures/temp"
# dir.create(TempDirPath)
# # Now, loop through each generation making the still image for that generation
# for(t in 1:Ngens){
#      gen <- Gens[t]
#      if(gen == 0){
#           GenNum <- "000"
#      } else if(gen < 100){
#           GenNum <- paste("0", gen, sep = "")
#      } else{
#           GenNum <- gen
#      }
#      FigName <- paste(TempDirPath, "/gen", GenNum, ".png", sep = "")
#      FigTitle <- paste("Generation", gen, sep = " ")
#      png(filename = FigName, width = FigWidth, height = FigHeight, pointsize = FigPointSize)
#           par(mfrow = c(1,2), mar = InnerMar, oma = OuterMar)
#           # Make the Haploid plot first
#           plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "", 
#                xlab = "", las = 1, cex.axis = AxisSize)
#           # Plot the lines for each loci
#           for(l in 1:length(Lseq)){
#                CurData <- Haploid[[l]]$GenVar[[t]]
#                lines(x = CurData$x, y = CurData$GenVar, lty = LociLines[l],
#                      col = LociCols[l])
#           }
#           # Add a legend
#           legend("topright", bty = "n", legend = LociLabels, 
#                  lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
#           # Add the plot title and y axis label
#           mtext(side = 3, text = "Haploid", line = PlotTitleLine, cex = PlotTitleSize)
#           text(labels = dBarLab, x = yLab_x, y = yLab_y, cex = xLabSize, srt = 0, xpd = NA)
#      
#           # Now add the Diploid dioecious plot
#           plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "", 
#                xlab = "", las = 1, yaxt = "n", cex.axis = AxisSize)
#           axis(2, labels = FALSE)
#           # Plot the lines for each loci
#           for(l in 1:length(Lseq)){
#                CurData <- DipDioe[[l]]$GenVar[[t]]
#                lines(x = CurData$x, y = CurData$GenVar, lty = LociLines[l],
#                      col = LociCols[l])
#           }
#           # Add the plot title
#           mtext(side = 3, text = "Diploid dioecious", line = PlotTitleLine, cex = PlotTitleSize)
#      
#           # Add the other figure text (Generation and axes labels)
#           mtext(text = FigTitle, side = 3, outer = TRUE, line = TitleLine, cex = TitleSize)
#           mtext(text = "Space", side = 1, outer = TRUE, line = xLabLine, cex = xLabSize)
#           #mtext(text = dBarLab, side = 2, outer = TRUE, line = yLabLine, cex = yLabSize)
#      
#      dev.off()
# }
# # Now make the animation
# GifName <- "ResultFigures/GenVarEvolution.gif"
# SysCommand1 <- paste("convert -delay 50 ", TempDirPath, "/gen*.png ", GifName, sep = "")
# SysCommand2 <- paste("rm -r ", TempDirPath, sep = "")
# system(SysCommand1)
# system(SysCommand2)

#----------------------------- Now create a figure of the last generation values 
#                             for dBar for all population types
# Haploid
FigName <- "ResultFigures/HaploidExpansion.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = dBarLab,
          xlab = "Space", las = 1, cex.axis = AxisSize, cex.lab = LabSize)
     for(l in 1:length(Lseq)){
          CurData <- Haploid[[l]]$dBar[[11]]
          lines(x = CurData$x, y = CurData$dBar, lty = LociLines[l], 
                col = LociCols[l])
     }
     # Add a legend
     legend("top", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
# Diploid Dioecious
FigName <- "ResultFigures/DipDioExpansion.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = dBarLab,
          xlab = "Space", las = 1, cex.axis = AxisSize, cex.lab = LabSize)
     for(l in 1:length(Lseq)){
          CurData <- DipDioe[[l]]$dBar[[11]]
          lines(x = CurData$x, y = CurData$dBar, lty = LociLines[l], 
                col = LociCols[l])
     }
     # Add a legend
     legend("top", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
# Diploid Monoecious, omega = 0
FigName <- "ResultFigures/DipMono_0_Expansion.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = dBarLab,
          xlab = "Space", las = 1, cex.axis = AxisSize, cex.lab = LabSize)
     for(l in 1:length(Lseq)){
          CurData <- DipMono_0[[l]]$dBar[[11]]
          lines(x = CurData$x, y = CurData$dBar, lty = LociLines[l], 
                col = LociCols[l])
     }
     # Add a legend
     legend("top", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
# Diploid Monoecious, omega = 0.5
FigName <- "ResultFigures/DipMono_05_Expansion.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = dBarLab,
          xlab = "Space", las = 1, cex.axis = AxisSize, cex.lab = LabSize)
     for(l in 1:length(Lseq)){
          CurData <- DipMono_05[[l]]$dBar[[11]]
          lines(x = CurData$x, y = CurData$dBar, lty = LociLines[l], 
                col = LociCols[l])
     }
     # Add a legend
     legend("top", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
# Diploid Monoecious, omega = 1
FigName <- "ResultFigures/DipMono_1_Expansion.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = dBarLab,
          xlab = "Space", las = 1, cex.axis = AxisSize, cex.lab = LabSize)
     for(l in 1:length(Lseq)){
          CurData <- DipMono_1[[l]]$dBar[[11]]
          lines(x = CurData$x, y = CurData$dBar, lty = LociLines[l], 
                col = LociCols[l])
     }
     # Add a legend
     legend("top", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
# Empty plot for presentations
FigName <- "ResultFigures/EmptyExpansion.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_dBar, main = "", ylab = dBarLab,
          xlab = "Space", las = 1, cex.axis = AxisSize, cex.lab = LabSize)
     # Add a legend
     legend("top", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()

#----------------------------- Now create a figure of the last generation values 
#                             for GenVar for all population types
yLabLine <- 3
xLabLine <- 3
LabSize <- 2
InnerMar <- c(5,6,2,2)
# Haploid
FigName <- "ResultFigures/HaploidExpansionGenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = GenVarLab, line = yLabLine, cex = LabSize, side = 2)
     mtext(text = "Space", line = xLabLine, cex = LabSize, side = 1)
     for(l in 1:length(Lseq)){
          CurData <- Haploid[[l]]$GenVar[[11]]
          lines(x = CurData$x, y = CurData$GenVar, lty = LociLines[l], 
                col = LociCols[l])
     }
     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
# Diploid Dioecious
FigName <- "ResultFigures/DipDioExpansionGenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = GenVarLab, line = yLabLine, cex = LabSize, side = 2)
     mtext(text = "Space", line = xLabLine, cex = LabSize, side = 1)
     for(l in 1:length(Lseq)){
          CurData <- DipDioe[[l]]$GenVar[[11]]
          lines(x = CurData$x, y = CurData$GenVar, lty = LociLines[l], 
                col = LociCols[l])
     }
     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
# Diploid Monoecious, omega = 0
FigName <- "ResultFigures/DipMono_0_ExpansionGenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = GenVarLab, line = yLabLine, cex = LabSize, side = 2)
     mtext(text = "Space", line = xLabLine, cex = LabSize, side = 1)
     for(l in 1:length(Lseq)){
          CurData <- DipMono_0[[l]]$GenVar[[11]]
          lines(x = CurData$x, y = CurData$GenVar, lty = LociLines[l], 
                col = LociCols[l])
     }
     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
# Diploid Monoecious, omega = 0.5
FigName <- "ResultFigures/DipMono_05_ExpansionGenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = GenVarLab, line = yLabLine, cex = LabSize, side = 2)
     mtext(text = "Space", line = xLabLine, cex = LabSize, side = 1)
     for(l in 1:length(Lseq)){
          CurData <- DipMono_05[[l]]$GenVar[[11]]
          lines(x = CurData$x, y = CurData$GenVar, lty = LociLines[l], 
                col = LociCols[l])
     }
     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
# Diploid Monoecious, omega = 1
FigName <- "ResultFigures/DipMono_1_ExpansionGenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = GenVarLab, line = yLabLine, cex = LabSize, side = 2)
     mtext(text = "Space", line = xLabLine, cex = LabSize, side = 1)
     for(l in 1:length(Lseq)){
          CurData <- DipMono_1[[l]]$GenVar[[11]]
          lines(x = CurData$x, y = CurData$GenVar, lty = LociLines[l], 
                col = LociCols[l])
     }
     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()
# Empty plot for presentations
FigName <- "ResultFigures/EmptyExpansionGenVar.pdf"
pdf(file = FigName, width = 12, height = 8, onefile = FALSE, paper = "special")
     par(mar = InnerMar, bg = "white")
     plot(x = NA, y = NA, xlim = xRange, ylim = yRange_GenVar, main = "", ylab = "",
          xlab = "", las = 1, cex.axis = AxisSize)
     mtext(text = GenVarLab, line = yLabLine, cex = LabSize, side = 2)
     mtext(text = "Space", line = xLabLine, cex = LabSize, side = 1)
     # Add a legend
     legend("topleft", bty = "n", legend = LociLabels, 
            lty = LociLines, col = LociCols, lwd = LegLineWidth, cex = LegSize)
dev.off()