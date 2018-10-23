disp <- function(rho, z, dmax){
     disp <- (dmax * exp(rho * z)) / (1 + exp(rho * z))
}

alleles <- function(DispVar, A){
     sigma <- sqrt(DispVar)
     AlleleVals <- seq(-6*sigma, 6*sigma, length.out = A)
     return(AlleleVals)
}

range <- function(DispVar, N){
     sigma <- sqrt(DispVar)
     RangeVals <- seq(-48*sigma, 48*sigma, length.out = N)
     return(RangeVals)
     
}

AlleleVals <- alleles(DispVar = 2, A = 255)
RangeVals <- range(DispVar = 2, N = 10000)

AlleleDisps <- disp(rho = 0.1, z = AlleleVals, dmax = 10)
RangeDisps <- disp(rho = 0.1, z = RangeVals, dmax = 10)

plot(x = RangeVals, y = RangeDisps, type = "l", xlab = "Phenotype", ylab = "Dispersal distance",
     main = "")
points(x = AlleleVals, AlleleDisps, col = "red")
