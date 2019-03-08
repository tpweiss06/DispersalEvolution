# Make some basic plots of the summary stats for the equilibrium populations

load("~/Desktop/PostdocResearch/DispersalEvolution/GitRepo/EquilibriumDispersal.rdata")
dmax <- 6

# Dispersal plot
quartz(width = 12, height = 7)
par(mfrow = c(2,10), mar = c(0,0,0,0))
for(i in 1:20){
     plot(NA, NA, xlim = c(0,140), ylim = c(0, dmax), xlab = "", ylab = "", axes = FALSE)
     points(x = 1:nrow(DispData[[i]]), y = DispData[[i]][,1])
     segments(x0 = 1:nrow(DispData[[i]]), y0 = DispData[[i]][,2], 
              x1 = 1:nrow(DispData[[i]]), y1 = DispData[[i]][,3])
     box()
}
