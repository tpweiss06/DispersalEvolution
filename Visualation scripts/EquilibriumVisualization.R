# Make some basic plots of the summary stats for the equilibrium populations

load("~/Desktop/PostdocResearch/DispersalEvolution/GitRepo/EquilibriumDispersal.rdata")
dmax <- 6

# Dispersal plot
quartz(width = 8, height = 5)
par(mfrow = c(3,10), mar = c(0,0,0,0), oma = c(1,1,1,1))
for(i in 1:length(ParamCombos)){
     plot(NA, NA, xlim = c(0,140), ylim = c(0, dmax), xlab = "", ylab = "", axes = FALSE)
     points(x = 1:nrow(DispData[[i]]), y = DispData[[i]][,1])
     SegSeq <- seq(1, nrow(DispData[[i]]), by = 3)
     segments(x0 = SegSeq, y0 = DispData[[i]][SegSeq,2], 
              x1 = SegSeq, y1 = DispData[[i]][SegSeq,3])
     box()
     abline(h = 1.25)
     abline(h = 1.75, lty = 2)
     abline(h = 0.75, lty = 2)
}
