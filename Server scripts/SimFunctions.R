# This script contains functions necessary to simulate a single realization
#    of the Dispersal Evolution IBM. Each function is documented with the 
#    necessary arguments and corresponding output. Functions are sorted into
#    three categories: biological, environmental, and book keeping.

# ------------------------------------------------------------------------------
# ---------------------------- Biological Functions ----------------------------
# ------------------------------------------------------------------------------

###### DispPhen
# This function calculates the dispersal phenotype for a set of individuals in
#    the simulated population.
### INPUTS
# PopMat:      The population matrix to use
# PopSize:     The size of the current population
# PopIndices: Column indices to use for the population matrix
# Haploid: A boolean indicator of haploid (1) or diploid (0) genetics
# L: The number of loci defining dispersal
# dmax: The maximum expected dispersal distance (in units of discrete patches)
# rho: Determines the slope of the transition from 0 to dmax as the loci sum
#         changes
# lambda: Offsets the inflection point of the function from 0
### OUTPUTS
# A vector of length PopSize with the i-th entry corresponding to the dispersal
#    phenotype of the i-th individual in the population input
DispPhen <- function(PopMat, PopSize, PopIndices, Haploid, L, dmax, rho, lambda){
     # Calculate the sum of the quantitative loci defining dispersal for
     #    each individual, accounting for the different possible cases
     if(PopSize == 1){
          LociSums <- sum(PopMat[PopIndices$DispCols])
     } else{
          if((Haploid == 1) & (L == 1)){
               LociSums <- PopMat[,PopIndices$DispCols]
          } else{
               LociSums <- rowSums(PopMat[,PopIndices$DispCols])
          }
     }
     # Now calculate the expected dispersal distance for each individual
     d <- (dmax * exp(rho * (LociSums - lambda))) / (1 + exp(rho * (LociSums - lambda)))
     return(d)
}

###### Dispersal
# Using the dispersal phenotypes from the above function, this function will
#    implement dispersal according to one of three different kernels and return
#    the realized changes in patch number in the x and y dimensions.
### INPUTS:
# d:           The dispersal phenotypes
# kern:        The dispersal kernel to be used. Can be "norm", "exp", or 
#                   "stud_t"
# direction:      A large binary vector of randomly generated 1s and -1s to indicate
#                   dispersal direction (forward or backward)
# DirectIndex:  Index of the where to look within the vector of random numbers
# PopIndices, PopSize, PopMat:  As previously defined
### OUTPUTS:
# The function will return a matrix with 2 columns and as many rows as individuals
#    within the population. Each row will correspond to the post dispersal 
#    change in patch number for X and Y dimensions in that order.
Disperse <- function(d, kern, direction, DirectIndex, PopSize, PopIndices, PopMat){
     # Generate the dispersal distances according to the type of kernel
     if(kern == "norm"){
          sigma <- (d * sqrt(pi)) / sqrt(2)
          dists <- abs(rnorm(n = PopSize, mean = 0, sd = sigma))
     } else if(kern == "exp"){
          dists <- rexp(n = PopSize, rate = 1 / d)
     } else if(kern == "stud_t"){
          dists <- rStudKern(n = PopSize, d = d)
     }
     
     # Calculate the movement of individuals and the subsequent change to their
     #    patch number assuming dispersal from the center of the patch
     NewX <- dists * direction[DirectIndex:(DirectIndex + PopSize - 1)]
     DeltaX <- ifelse(NewX < 0, ceiling(NewX - 1/2), floor(NewX + 1/2))
     return(DeltaX)
}

###### rStudKern
# This function will perform a specified number of random draws from the one 
#    parameter student's t dispersal kernel used in the Shaw 2014 ppaer on the
#    consequences of risky dispersal.
### INPUTS
# n:      The number of random numbers to generate
# d:      The dispersal trait(s) for the individual(s) involved. This will be 
#              converted to another parameter to use in the probability 
#              distribution
### OUTPUTS
# n random numbers generated from the student's t dispersal kernel
rStudKern <- function(n, d){
     # First translate d to the u parameter used in the distribution
     u <- sqrt(d / 1.571)
     # Then generate n uniform random numbers between 0 and 1
     n_unif <- runif(n = n, min = 0, max = 1)
     # Now map them to the Student's t dispersal kernel using the CDF
     StudKernVals <- sqrt((u * n_unif) / (1 - n_unif))
     return(StudKernVals)
}

###### Reproduce
# This function will use the current abundances of the patch and the type of
#    reproduction being modeled (asexual, monoecious, or dioecious) to calculate
#    the expected population size of each patch in the next generation. Then,
#    it will draw and return the realized population sizes from a Poisson 
#    distribution.
### INPUTS
#    R: The intrinsic growth rate for the Ricker model (Melbourne & Hastings 2008)
#    OccPatches: A vector with the unique x identifiers of each occupied patch
#    psi: The expected sex ratio in the population
#    CurBeta: The current center of the range
#    tau: Defines the width of the range
#    gamma: Defines the starkness of the range edge
#    R: The population growth rate
#    Kmax: The maximum achievable carrying capacity in the landscape
#    Expand: A boolean variable indicating whether the current simulation is a
#              strict range expansion or not
#    omega: The probability of sefl-fertilization. Used to check for 0 population growth
#         if only 1 individual and obligatory outcrossing (omega == 0)
#    PopMat, Haploid, PopIndices as previously defined
### OUTPUTS
# A vector of the realized abundances of each occupied patch in the next generation
Reproduce <- function(R, OccPatches, PopMat, Haploid, PopIndices, psi, CurBeta,
                      tau, gamma, Kmax, Expand, omega){
     # Check that Rmax is within the proper range
     if(R < 1){
          write("R values less than 1 will result in negative alpha values, causing unrestricted population growth",
                stderr())
          return(NULL)
     } else if(R == 1){
          write("Setting R to 1 results in alpha values of either 0 or NaN, resulting in unrestricted growth or errors respectively",
                stderr())
          return(NULL)
     }
     
     # Determine the number of occupied patches 
     NumPatches <- length(OccPatches)
     # Get the patch alpha values based on the location within the range
     if(Expand){
          Alphas <- rep(log(R) / Kmax, NumPatches)
     } else{
          Alphas <- PatchAlphas(OccPatches = OccPatches, NumPatches = NumPatches,
                                CurBeta = CurBeta, tau = tau, gamma = gamma, R = R,
                                Kmax = Kmax)
     }
     
     # Create an empty object to hold the expected population sizes for the next generation
     Ntp1 <- rep(NA, NumPatches)
     
     # Loop through each occupied patch and calculate the expected population
     #    size in the next generation
     for(i in 1:NumPatches){
          # Extract the local population
          CurPatchPop <- which( (PopMat[,PopIndices$x1] == OccPatches[i]) )
          PatchAbund <- length(CurPatchPop)
          
          # Now determine the expected population growth according to the stochastic
          #    Ricker models derived by Melbourne & Hastings (2008) and the 
          #    simulation parameters
          if(is.infinite(Alphas[i])){
               Ntp1[i] <- 0
          } else{
               if(Haploid == 1){
                    Ntp1[i] <- PatchAbund * R * exp(-1 * Alphas[i] * PatchAbund)
               } else{
                    if(is.null(PopIndices$sex)){
                         # check for obligatory outcrossing (omega == 0) and only
                         #    1 individual. If true, then 0 population growth, but
                         #    if not, then normal growth
                         if((omega == 0) & (PatchAbund == 1)){
                              Ntp1[i] <- 0
                         }else{
                              Ntp1[i] <- PatchAbund * R * exp(-1 * Alphas[i] * PatchAbund)
                         }
                    } else{
                         CurFemales <- which( (PopMat[,PopIndices$x1] == OccPatches[i]) & 
                                                   (PopMat[,PopIndices$sex] == 1) )
                         nFem <- length(CurFemales)
                         # Calculate expected population growth with a forced 0
                         #    if the whole population is female
                         AllFemale <- nFem == PatchAbund
                         Ntp1[i] <- (nFem * (R / psi) * exp(-1 * Alphas[i] * PatchAbund)) * 
                                        (1 - AllFemale)
                    }
               }
          }
     }
     # Use the expected population sizes to generate the realized population 
     #    sizes for each patch
     RealizedNtp1 <- rpois(n = NumPatches, lambda = Ntp1)
     return(RealizedNtp1)
}

###### NextPopMat
# This function will produce a new PopMat for the next generation.
### INPUTS
# Ntp1:     The output from the Reproduce function
# SumNtp1:  The total number of new offspring making up the next generation
# L:       The number of loci defining the dispersal trait
# LocusVec: An integer vector of 1 to the number of loci
# AlleleVec: An integer vector of 1 to the number of alleles defining dispersal
#              (in the haploid case, this will be identical to LocusVec)
# SegVals:  A binary matrix with two rows used to determine which alleles each
#              parent contributes to the offspring
# SegIndex: The current column index value to be used with the SegVals matrix
# NumMutVec: A large vector with the randomly generated number of mutations to
#              be used with each offspring
# MutIndex:  The index to be used with NumMutVec
# MutStd:    Previously calculated standard deviation for mutational affects
# SelfRands: A large binary vector indicating whether each offspring results
#              from selfing (1) or outcrossing (0). This is calculated previously
#              using the SelfProb parameter
# PopMat, PopIndices, OccPatches, and psi as defined previously
### OUTPUTS
# This function returns a new population matrix with the same columns but with
#    updated row information for the next generation.
NextPopMat <- function(Ntp1, PopIndices, OccPatches, psi, SumNtp1, L, PopMat, LocusVec,
                       SegVals, SegIndex, NumMutVec, OffspringIndex, AlleleVec, MutStd,
                       SelfRands, NumCols, SexRands){
     NewPopMat <- matrix(NA, nrow = SumNtp1, ncol = NumCols)
     # Filter for only the patches that produced offspring and fill in the matrix
     NewOccPatches <- which(Ntp1 > 0)
     Ntp1 <- Ntp1[NewOccPatches]
     StartRow <- 1
     EndRow <- Ntp1[1]
     for(i in 1:length(NewOccPatches)){
          CurVec <- StartRow:EndRow
          # Fill in the location details
          NewPopMat[CurVec, PopIndices$x0] <- OccPatches[NewOccPatches[i]]
          # Identify the potential parents in the current patch
          locals <- which( (PopMat[,PopIndices$x1] == OccPatches[NewOccPatches[i]]) )
          NumLocals <- length(locals)
          # Determine inheritence of alleles and assignment of parents depending
          #    on the mode of reproduction
          if(Haploid){
               # Determine which individuals produce offspring
               if(NumLocals == 1){
                    parents <- rep(locals, Ntp1[i])
               } else{
                    parents <- sample(locals, size = Ntp1[i], replace = TRUE)
               }
               NewPopMat[CurVec, PopIndices$DispCols] <- HapInherit(PatchPop = Ntp1[i], parents = parents, PopMat = PopMat,
                                                                    PopIndices = PopIndices, NumMutVec = NumMutVec,
                                                                    MutIndex = OffspringIndex, MutStd = MutStd, AlleleVec = AlleleVec)
          } else{
               if(is.null(PopIndices$sex)){
                    if(NumLocals == 1){
                         parent1 <- rep(locals, Ntp1[i])
                         parent2 <- rep(locals, Ntp1[i])
                    } else{
                         parent1 <- sample(locals, size = Ntp1[i], replace = TRUE)
                         parent2 <- sample(locals, size = Ntp1[i], replace = TRUE)
                         SameParent <- which(parent1 == parent2)
                         for(j in SameParent){
                              ParentPool <- setdiff(locals, parent1[j])
                              parent2[j] <- sample(ParentPool, size = 1)
                         }
                         SelfFert <- SelfRands[(OffspringIndex + StartRow):(OffspringIndex + EndRow)]
                         parent2 <- ifelse(SelfFert == 1, parent1, parent2)
                    }
                    parents <- cbind(parent1, parent2)
               } else{
                    NewPopMat[CurVec, PopIndices$sex] <- SexRands[(OffspringIndex + StartRow):(OffspringIndex + EndRow)]
                    # Identify the females and males present in the current patch
                    females <- which( (PopMat[,PopIndices$x1] == OccPatches[NewOccPatches[i]]) &
                                           (PopMat[, PopIndices$sex] == 1) )
                    males <- which( (PopMat[,PopIndices$x1] == OccPatches[NewOccPatches[i]]) &
                                         (PopMat[, PopIndices$sex] == 0) )
                    if(length(females) == 1){
                         parent1 <- rep(females, Ntp1[i])
                    }else{
                         parent1 <- sample(females, size = Ntp1[i], replace = TRUE)
                    }
                    if(length(males) == 1){
                         parent2 <- rep(males, Ntp1[i])
                    }else{
                         parent2 <- sample(males, size = Ntp1[i], replace = TRUE)
                    }
                    parents <- cbind(parent1, parent2)
               }
               NewPopMat[CurVec, PopIndices$DispCols] <- DipInherit(PatchPop = Ntp1[i], L = L, PopMat = PopMat,
                                                                    parents = parents, SegVals = SegVals, SegIndex = SegIndex,
                                                                    NumMutVec = NumMutVec, MutIndex = OffspringIndex, 
                                                                    AlleleVec = AlleleVec, MutStd = MutStd, LocusVec = LocusVec,
                                                                    PopIndices = PopIndices)
          }
          # Update the starting and ending row for the next iteration of the loop
          StartRow <- EndRow + 1
          EndRow <- EndRow + Ntp1[i+1]
     }
     return(NewPopMat)
}

###### HapInheritence
# This function handles allele inheritence and mutation for haploid scenarios
### INPUTS
# parents:     A vector of the parent producing each offspring
# PatchPop:    The size of the population in the current patch
# PopMat, PopIndices, NumMutVec, MutIndex, MutStd, and AlleleVec as previously defined
### OUTPUTS
# This function returns a matrix with a row for each offspring and columns corresponding
#    to the inherited allele values
HapInherit <- function(parents, PopMat, PopIndices, NumMutVec, MutIndex, PatchPop, MutStd,
                       AlleleVec){
     OffspringAlleles <- PopMat[parents, PopIndices$DispCols]
     # Determine the number of mutations affecting each offspring
     NumMut <- NumMutVec[MutIndex:(MutIndex + PatchPop-1)]
     MutOffspring <- which(NumMut != 0)
     # Mutate the affected alleles
     if(length(AlleleVec) > 1){
               for(i in MutOffspring){
                    MutLocus <- sample(AlleleVec, size = NumMut[i], replace = TRUE)
                    if(length(parents) > 1){
                         OffspringAlleles[i, MutLocus] <- rnorm(mean = OffspringAlleles[i, MutLocus],
                                                                sd = MutStd, n = NumMut[i])
                    } else{
                         OffspringAlleles[MutLocus] <- rnorm(mean = OffspringAlleles[MutLocus],
                                                                sd = MutStd, n = NumMut)
                    }
               }
     } else{
          OffspringAlleles[MutOffspring] <- rnorm(mean = OffspringAlleles[MutOffspring],
                                                  sd = MutStd, n = length(MutOffspring))
     }
     return(OffspringAlleles)
}

###### DipInheritence
# This function handles allele inheritence and mutation for diploid scenarios
### INPUTS
# All function inputs are as defined previously
### OUTPUTS
# This function returns a matrix with a row for each offspring and columns corresponding
#    to the inherited allele values
DipInherit <- function(PopIndices, parents, PopMat, PatchPop, L, SegVals, AlleleVec,
                        MutStd, NumMutVec, MutIndex, LocusVec, SegIndex){
     OffspringAlleles <- matrix(NA, nrow = PatchPop, ncol = 2*L)
     for(i in 1:PatchPop){
          ParentLoci <- PopMat[parents[i,], PopIndices$DispCols]
          Parent1Alleles <- LocusVec + L * SegVals[1,(SegIndex + (i - 1) * L):(SegIndex + i * L - 1)]
          Parent2Alleles <- LocusVec + L * SegVals[2,(SegIndex + (i - 1) * L):(SegIndex + i * L - 1)]
          OffspringAlleles[i,] <- c(ParentLoci[1,Parent1Alleles], ParentLoci[2,Parent2Alleles])
     }
     # Extract the number of mutations for offspring from the NumMutVec
     NumMut <- NumMutVec[MutIndex:(MutIndex+PatchPop-1)]
     # Now step through each offspring in which a mutation takes place and alter
     #    allele values appropriately
     MutOffspring <- which(NumMut != 0)
     for(i in MutOffspring){
          MutLocus <- sample(AlleleVec, size = NumMut[i], replace = TRUE)
          OffspringAlleles[i,MutLocus] <- rnorm(mean = OffspringAlleles[i,MutLocus], 
                                              sd = MutStd, n = NumMut[i])
     }
     return(OffspringAlleles)
}

###### NextPopMatNoEvol
# This function will produce a new PopMat for the next generation, but with no evolution.
#    Alleles will be drawn from the vectors representing their frequencies at the
#    beginning of the simulation.
### INPUTS
# InitAlleles: A matrix with a column for each allele containing all the allele
#              values from the initial population for the simulation
# Ntp1, PopIndices, OccPatches, SumNtp1, NumCols, OffspringIndex, and SexRands 
#    are all as defined before
### OUTPUTS
# This function returns a new population matrix with the same columns but with
#    updated row information for the next generation.
NextPopMatNoEvol <- function(Ntp1, PopIndices, OccPatches, SumNtp1, NumCols, SexRands,
                             InitAlleles, OffspringIndex){
     NewPopMat <- matrix(NA, nrow = SumNtp1, ncol = NumCols)
     # Filter for only the patches that produced offspring and fill in the matrix
     NewOccPatches <- which(Ntp1 > 0)
     Ntp1 <- Ntp1[NewOccPatches]
     StartRow <- 1
     EndRow <- Ntp1[1]
     for(i in 1:length(NewOccPatches)){
          CurVec <- StartRow:EndRow
          # Fill in the location details
          NewPopMat[CurVec, PopIndices$x0] <- OccPatches[NewOccPatches[i]]
          # Fill in the sex column
          if(!is.null(PopIndices$sex)){
               NewPopMat[CurVec, PopIndices$sex] <- SexRands[(OffspringIndex + StartRow):(OffspringIndex + EndRow)]
          }
          # Update the starting and ending row for the next iteration of the loop
          StartRow <- EndRow + 1
          EndRow <- EndRow + Ntp1[i+1]
     }
     # Fill in the allele columns
     for(j in 1:length(PopIndices$DispCols)){
          NewPopMat[,PopIndices$DispCols[j]] <- sample(InitAlleles[,j], size = SumNtp1, replace = TRUE)
     }
     return(NewPopMat)
}

# ------------------------------------------------------------------------------
# -------------------------- Environmental Functions ---------------------------
# ------------------------------------------------------------------------------

###### PatchAlphas
# This function will generate a vector of alpha values to use in the Reproduce
#    function. The alphas will each correspond to one of the occupied patches of
#    the landscape.
### INPUTS
# NumPatches:  The number of occupied patches in the landscape
# OccPatches, CurBeta, tau, gamma, Kmax, and R as previously defined
### OUTPUT
# A vector of alpha values to be used in the stochastic Ricker model
PatchAlphas <- function(OccPatches, NumPatches, CurBeta, tau, gamma, Kmax, R){
     PatchK <- rep(NA, NumPatches)
     for(i in 1:NumPatches){
          if(OccPatches[i] <= (CurBeta - tau - 1/gamma)){
               PatchK[i] <- 0
          }else if(OccPatches[i] < (CurBeta - tau)){
               PatchK[i] <- (1 - (CurBeta - tau - OccPatches[i])*gamma)*Kmax
          }else if(OccPatches[i] <= (CurBeta + tau)){
               PatchK[i] <- Kmax
          }else if(OccPatches[i] < (CurBeta + tau + 1/gamma)){
               PatchK[i] <- (1 - (OccPatches[i] - CurBeta - tau)*gamma)*Kmax
          }else{
               PatchK[i] <- 0
          }
     }
     PatchAlpha <- log(R) / PatchK
     return(PatchAlpha)
}

###### ChangeClimate
# This function will create a vector of beta values depending on how fast the 
#    climate is set to change, assuming a constant rate of change
### INPUTS
# BetaInit:         The starting value of beta
# LengthShift:      The length of the climate shift
# v:        The speed with which the climate shifts
### OUTPUT
# A one dimensional vector containing the shifted beta values at each time point
#    for the duration of the climate change period
ChangeClimate <- function(BetaInit, LengthShift, v){
     # First check for a negative speed
     if(v < 0){
          write("Negative speed is not supported for climate change", stderr())
          return(NULL)
     }
     # Next check for a 0 value for the length of climate change
     if(LengthShift == 0){
          return(NULL)
     } else{
          return(BetaInit + v * 1:LengthShift)
     }
}

# ------------------------------------------------------------------------------
# --------------------------- Bookkeeping Functions ----------------------------
# ------------------------------------------------------------------------------

###### Initialize
# This function will initialize a population matrix of founders to start a 
#    simulation.
### INPUTS
# PopIndices:  A list of column indices for different parts of the matrix
# NumCols:     The number of columns in the population matrix
# DispCols:    The column indices for the dispersal loci
# PopSize:     The number of founders to create in the intitial population
# BetaInit:    The starting location of the range center where individuals will
#                   initialize
# psi:           The sex ratio of the founding population (if dioecious). Set to
#                   0.5 by default
# DispVar:     The variance (sigma^2) of the initial distribution of dispersal 
#                   genotypes in a haploid scenario. In diploid scenarios, the
#                   standard deviation of initial genotypes will be 2*sigma^2
# A:           The number of alleles at each loci. Set to 255 by default, matching
#                   the commonly used simulation platform QuantiNemo
### OUTPUTS
# A filled in population matrix to start generation 0
Initialize <- function(PopIndices, BetaInit, psi, DispVar, A = 255, NumCols,
                       tau, gamma, Kmax, R){
     # First determine the number of individuals to initialize by calculating
     #    the realized carrying capacity of each patch to be occupied in the
     #    stable range
     xSeq <- -75:75
     AlphaVals <- PatchAlphas(OccPatches = xSeq, NumPatches = length(xSeq), CurBeta = BetaInit,
                              tau = tau, gamma = gamma, Kmax = Kmax, R = R)
     Kvals <- as.integer(log(R) / AlphaVals)
     TotalPop <- sum(Kvals)
     
     # Now make an empty population matrix with the correct names
     PopMat <- matrix(NA, nrow = TotalPop, ncol = NumCols)
     
     # Next, fill in the x1 column for the founders (the founders are
     #    considered post dispersal and will reproduce next)
     PopIndex <- 1
     for(k in 1:length(Kvals)){
          if(Kvals[k] > 0){
               PopMat[PopIndex:(PopIndex + Kvals[k] - 1), PopIndices$x1] <- xSeq[k]
               PopIndex <- PopIndex + Kvals[k]
          }
     }
     
     # Fill in the sex column if it is present
     if( !(is.null(PopIndices$sex)) ){
          PopMat[,PopIndices$sex] <- rbinom(n = TotalPop, size = 1, prob = psi)
     }
     
     # Now fill in the dispersal allele columns and return the matrix
     sigma <- sqrt(DispVar / length(PopIndices$DispCols))
     AlleleVals <- seq(-6*sigma, 6*sigma, length.out = A)
     InitFreqProbs <- dnorm(x = AlleleVals, mean = 0, sd = sigma)
     for(i in PopIndices$DispCols){
          PopMat[,i] <- sample(AlleleVals, size = TotalPop, replace = TRUE, prob = InitFreqProbs)
     }
     return(PopMat)
}

###### PopMatColNames
# This function generates the column names to use for a given simulation
### INPUTS
# L:           the number of loci defining an individual's dispersal ability
# monoecious:  a boolean value indicating whether the individuals in the simulation
#              are monoecious (i.e. possessing both the female and male 
#              reproductive organs in the same individual) or not.
# Haploid:     a boolean value indicating whether individuals are haploid (TRUE),
#                   or diploid (FALSE)
# example:     a boolean value indicating whether the function should simply 
#                   return an example character vector containing the column 
#                   names and order used
### OUTPUTS
# This function creates a list with the indices for different column types in the
#    population matrix. If example is TRUE, it simply returns a character vector.
PopMatColNames <- function(L, monoecious, Haploid, example = FALSE){
     if(example){
          if(Haploid){
               MatNames <- c("x0", "x1", "disp1_1", "disp1_2", "...", "disp1_L", "sex")
          } else{
               MatNames <- c("x0", "x1", "disp1_1", "disp1_2", "...", "disp1_L",
                             "disp2_1", "...", "disp2_L", "sex")
          }
          return(MatNames)
     }
     
     # Create the list and populate it
     if(monoecious | Haploid){
          PopIndices <- vector(mode = "list", length = 3)
          names(PopIndices) <- c("x0", "x1", "DispCols")
     } else{
          PopIndices <- vector(mode = "list", length = 4)
          names(PopIndices) <- c("x0", "x1", "DispCols", "sex")
          PopIndices$sex <- 3 + 2*L
     }
     PopIndices$x0 <- 1
     PopIndices$x1 <- 2
     if(Haploid){
          PopIndices$DispCols <- 3:(3+L-1)
     } else{
          PopIndices$DispCols <- 3:(3+2*L-1)
     }
     #PopIndices$DispCols <- ifelse(Haploid, 3:(3+L-1), 3:(3+2*L-1))
     return(PopIndices)
}

###### GetSafeID
# This function looks through the current working directory and finds a safe
#    name to use that will not overwrite anything else when saving simulation
#    results
### INPUTS
# ParentDirectory:  the file path to the current working directory
# parallel:         A boolean variable indicating whether the current simulations 
#                        are taking place over multiple nodes so that file names 
#                        should include the name of the current node
### OUTPUTS
# The full file path for a new directory that the current simulation results will
#    be saved to
GetSafeID <- function(ParentDirectory, parallel = FALSE){
     # Start by trying out 1 and then increase as necessary
     newID <- 1
     
     # If the simulations are being done in parallel, include the name of the
     #    current node in the new directory name
     if(parallel){
          NodeName <- Sys.getpid()
          DirName <- paste(ParentDirectory, "/", NodeName, "_Sim", newID, sep = "")
          # Until we find a directory name that does not already exist, continue
          #    to increase the ID variable
          while(dir.exists(DirName)){
               newID <- newID + 1
               DirName <- paste(ParentDirectory, "/", NodeName, "_Sim", newID, sep = "")
          }
     } else{
          DirName <- paste(ParentDirectory, "/Sim", newID, sep="")
          while(dir.exists(DirName)){
               newID <- newID + 1
               DirName <- paste(ParentDirectory, "/Sim", newID, sep="")
          }
     }
     
     # Return the identified, safe directory name
     return(DirName)
}

###### SaveParams
# This function will unpack a list of parameters and save them in an output file
#    that can be sourced to load all parameters into the environment. In this 
#    way, the parameters for each simulation are saved and easily accessible and
#    it allows the parameters for the simulations to be passed in to the main 
#    function in an easier format (i.e. a single list).
### INPUTS
# parameters:  a list containing all the parameter values necessary for a 
#                   simulation.
# FilePath:    a character vector for the results directory for the current
#                   simulation.
### OUTPUTS
# This function will unpack the parameter vector and create a .R file which
#    creates objects with the parameter names and values and named for the
#    specific simulation.
SaveParams <- function(parameters, FilePath){
     # First cheack that all necessary parameters are included
     ParamCheck <- names(parameters) == c("BetaInit", "gamma", "tau", "omega",
                                          "U", "sigma", "L", "R", "Kmax", "Haploid",
                                          "kern", "monoecious", "BurnIn",
                                          "LengthShift", "ExpandGens", "psi",
                                          "DispVar", "dmax", "rho", "lambda", "NumRands", "v")
     if(sum(ParamCheck) != length(ParamCheck)){
          write("Incorrect names or number of input parameters", stderr())
          return(NULL)
     }
     
     # Calculate rho so that initial genotypes fall between 
     #    0.25*dmax and 0.75*dmax
     sigma <- sqrt(parameters$DispVar)
     rho <- log(3)/(6*sigma)
     
     OutFile <- paste(FilePath, "parameters.R", sep = "/")
     sink(OutFile)
     cat("BetaInit <- ", parameters$BetaInit, "\n", sep = "")
     cat("gamma <- ", parameters$gamma, "\n", sep = "")
     cat("tau <- ", parameters$tau, "\n", sep = "")
     cat("omega <- ", parameters$omega, "\n", sep = "")
     cat("U <- ", parameters$U, "\n", sep = "")
     cat("sigma <- ", parameters$sigma, "\n", sep = "")
     cat("L <- ", parameters$L, "\n", sep = "")
     cat("R <- ", parameters$R, "\n", sep = "")
     cat("Kmax <- ", parameters$Kmax, "\n", sep = "")
     cat("Haploid <- ", parameters$Haploid, "\n", sep = "")
     cat("kern <- \"", parameters$kern, "\"\n", sep = "")
     cat("monoecious <- ", parameters$monoecious, "\n", sep = "")
     cat("BurnIn <- ", parameters$BurnIn, "\n", sep = "")
     cat("LengthShift <- ", parameters$LengthShift, "\n", sep = "")
     cat("ExpandGens <- ", parameters$ExpandGens, "\n", sep = "")
     cat("psi <- ", parameters$psi, "\n", sep = "")
     cat("DispVar <- ", parameters$DispVar, "\n", sep = "")
     cat("dmax <- ", parameters$dmax, "\n", sep = "")
     cat("rho <- ", parameters$rho, "\n", sep = "")
     cat("lambda <- ", parameters$lambda, "\n", sep = "")
     cat("NumRands <- ", parameters$NumRands, "\n", sep = "")
     cat("v <- ", parameters$v, "\n", sep = "")
     sink()
}


# ------------------------------------------------------------------------------
# ------------------------ Full Simulation Function ----------------------------
# ------------------------------------------------------------------------------

###### StationarySim
# This function will initiate a population in an empty landscape with set range
#    boundaries. The population will be allowed to grow and spread throughout the
#    range for a set number of generations to equilibriate. The only output from
#    this function will be a population matrix from the final generation.
### INPUTS
# parameters: a list containing all of the parameters necessary to run the
#    entire simulation as well as any necessary for further downstream 
#    simulations, like expansion and range shifts
# parallel: A boolean variable indicating whether the simulations are being run
#              on a server or not (which affects how file paths are determined).
# SimID: An optional argument for providing an explicit name for the current simulation
#         for example, this could be used for trouble shooting and creating "test" sims
### OUTPUT
# The population matrix for the last generation of the simulation
StationarySim <- function(parameters, parallel = FALSE, SimID = NA, SimDirectory = NULL){
     # First generate a safe directory name and create it to save all output
     #    from the simulation
     if(is.na(SimID)){
          ResultsDir <- GetSafeID(ParentDirectory = SimDirectory, parallel = parallel)
          
     } else{
          ResultsDir <- paste(SimDirectory, SimID, sep = "/")
     }
     dir.create(ResultsDir)
     
     # Next, save the parameters used for this simulation and source the file
     #    to have access to them within the function
     SaveParams(parameters = parameters, FilePath = ResultsDir)
     source(paste(ResultsDir, "parameters.R", sep = "/"))
     
     # Next, create the population column indices
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     if(Haploid | monoecious){
          NumCols <- 2 + 2^(1-Haploid) * L
     } else{
          NumCols <- 3 + 2 * L
     }
     LocusVec <- 1:L
     AlleleVec <- 1:(2^(1-Haploid) * L)
     
     # Create vectors of random numbers (including sex determination random numbers if
     #    necessary)
     if(Haploid){
          SegVals <- NULL
          SegIndex <- NA
     }else{
          SegVals <- matrix(NA, nrow = 2, ncol = NumRands)
          SegVals[1,] <- sample(c(0,1), replace = TRUE, size = NumRands)
          SegVals[2,] <- sample(c(0,1), replace = TRUE, size = NumRands)
          SegIndex <- 1
     }
     nu <- U / (2^(1-Haploid)*L)
     NumMut <- rbinom(n = NumRands, size = 2^(1-Haploid)*L, prob = nu)
     if( !(is.null(PopIndices$sex)) ){
          SexRands <- rbinom(n = NumRands, size = 1, prob = psi)
     } else{
          SexRands <- NULL
          if(!Haploid){
               SelfRands <- rbinom(n = NumRands, size = 1, prob = omega)
          }
     }
     OffspringIndex <- 1
     direction <- sample(c(-1,1), size = NumRands, replace = TRUE)
     DirectIndex <- 1
     
     # Next initialize the generation 0 founding population and calculate dispersal
     #    phenotypes
     PopMat <- Initialize(PopIndices = PopIndices, NumCols = NumCols,
                          BetaInit = BetaInit, DispVar = DispVar, psi = psi, 
                          tau = tau, gamma = gamma, Kmax = Kmax, R = R)
     PopSize <- nrow(PopMat)
     
     # Generate the vector of occupied patches
     OccPatches <- unique(PopMat[,PopIndices$x1])
     
     # Generate the population sizes for the next generation
     RealizedNtp1 <- Reproduce(CurBeta = BetaInit, gamma = gamma, tau = tau,
                               R = R, Kmax = Kmax, PopMat = PopMat,
                               PopIndices = PopIndices, psi = psi,
                               OccPatches = OccPatches, Haploid = Haploid,
                               Expand = FALSE, omega = omega)
     
     # Make sure the population isn't going immediately extinct
     PopSize <- sum(RealizedNtp1)
     ExtCounter <- 0
     while((PopSize == 0) & (ExtCounter < 10)){
          RealizedNtp1 <- Reproduce(CurBeta = BetaInit, gamma = gamma, tau = tau,
                                    R = R, Kmax = Kmax, PopMat = PopMat,
                                    PopIndices = PopIndices, psi = psi,
                                    OccPatches = OccPatches, Haploid = Haploid,
                                    Expand = FALSE, omega = omega)
          PopSize <- sum(RealizedNtp1)
          ExtCounter <- ExtCounter + 1
     }
     # Now create a new population matrix for the next generation
     PopMat <- NextPopMat(Ntp1 = RealizedNtp1, PopIndices = PopIndices, 
                          OccPatches = OccPatches, psi = psi, SumNtp1 = PopSize, 
                          L = L, PopMat = PopMat, LocusVec = LocusVec,
                          SegVals = SegVals, SegIndex = SegIndex, NumMutVec = NumMut, 
                          OffspringIndex = OffspringIndex, AlleleVec = AlleleVec, 
                          MutStd = sigma, SelfRands = SelfRands, NumCols = NumCols,
                          SexRands = SexRands)
     # Now update the SegIndices
     SegIndex <- SegIndex + PopSize*L
     # And update the OffspringIndex
     OffspringIndex <- OffspringIndex + PopSize

     # Now run through the actual simulation
     InitialPopMat <- PopMat
     g <- 1
     while((g <= BurnIn) & (ExtCounter < 10)){
          # Check for extinction before dispersal and reproduction
          if(PopSize == 0){
               PopMat <- InitialPopMat
               g <- 1
               ExtCounter <- ExtCounter + 1
               PopSize <- nrow(PopMat)
          } else{
               disps <- DispPhen(PopMat = PopMat, PopSize = PopSize, PopIndices = PopIndices, 
                                 Haploid = Haploid, L = L, dmax = dmax, rho = rho, lambda = lambda)
               
               # Check that there are enough entries left in the direction vector and repopulate
               #    it if necessary
               if( (DirectIndex + PopSize) > NumRands ){
                    direction <- sample(c(-1,1), size = NumRands, replace = TRUE)
                    DirectIndex <- 1
               }
               # Get the new post dispersal locations
               Deltas <- Disperse(d = disps, kern = kern, direction = direction,
                                  DirectIndex = DirectIndex, PopSize = PopSize,
                                  PopIndices = PopIndices, PopMat = PopMat)
               PopMat[,PopIndices$x1] <- PopMat[,PopIndices$x0] + Deltas
               
               # Update the angle index
               DirectIndex <- DirectIndex + PopSize
               
               # Generate the matrix of occupied patches and get population numbers
               #    for the next generation
               OccPatches <- unique(PopMat[,PopIndices$x1])
               RealizedNtp1 <- Reproduce(CurBeta = BetaInit, gamma = gamma, tau = tau,
                                         R = R, Kmax = Kmax, PopMat = PopMat,
                                         PopIndices = PopIndices, psi = psi,
                                         OccPatches = OccPatches, Haploid = Haploid,
                                         Expand = FALSE, omega = omega)
               
               # Now create a new population matrix for the next generation
               PopSize <- sum(RealizedNtp1)
               if(PopSize > 0){
                    # Check that the SegVals vectors contain enough values and resample if not
                    if( (!Haploid) & ((SegIndex + PopSize*L) > NumRands) ){
                         SegVals[1,] <- sample(c(0,1), replace = TRUE, size = NumRands)
                         SegVals[2,] <- sample(c(0,1), replace = TRUE, size = NumRands)
                         SegIndex <- 1
                    }
                    # Check the same for the OffspringIndex
                    if( (OffspringIndex + PopSize) > NumRands){
                         NumMut <- rbinom(n = NumRands, size = L, prob = nu)
                         if( !(is.null(PopIndices$sex)) ){
                              SexRands <- rbinom(n = NumRands, size = 1, prob = psi)
                         }
                         OffspringIndex <- 1
                    }
                    PopMat <- NextPopMat(Ntp1 = RealizedNtp1, PopIndices = PopIndices, 
                                         OccPatches = OccPatches, psi = psi, SumNtp1 = PopSize, 
                                         L = L, PopMat = PopMat, LocusVec = LocusVec,
                                         SegVals = SegVals, SegIndex = SegIndex, NumMutVec = NumMut, 
                                         OffspringIndex = OffspringIndex, AlleleVec = AlleleVec, 
                                         MutStd = sigma, SelfRands = SelfRands, NumCols = NumCols,
                                         SexRands = SexRands)
                    # Now update the SegIndices
                    SegIndex <- SegIndex + PopSize*L
                    # And update the OffspringIndex
                    OffspringIndex <- OffspringIndex + PopSize
               } else{
                    PopMat <- matrix(NA, nrow = 0, ncol = NumCols)
               }
               g <- g + 1
          }
     }
     # Check for the number of extinctions that occured before getting a stable population
     if(ExtCounter == 10){
          write("10 or more extinction events while trying to establish a population. Try other parameter values.",
                stderr())
          return(NULL)
     }
     # Finally, save the results here
     PopMatNames <- c("x0", "x1", paste("disp1", 1:L, sep = "_"))
     if(!Haploid){
          PopMatNames <- c(PopMatNames, paste("disp2", 1:L, sep = "_"))
     }
     if( !(is.null(PopIndices$sex)) ){
          PopMatNames <- c(PopMatNames, "sex")
     }
     colnames(PopMat) <- PopMatNames
     write.csv(PopMat, file = paste(ResultsDir, "EquilibriumPopMat.csv", sep = "/"), 
               row.names = FALSE, quote = FALSE)
     # return the results directory for the simulation to allow for easy sorting
     return(ResultsDir)
}

###### RangeExpand
RangeExpand <- function(SimDir, parallel = FALSE, SumMatSize = 5000, EquilibriumPrefix = NULL,
                        ExpandPrefix = NULL){
     # Load the necessary parameters
     source(paste(EquilibriumPrefix, SimDir, "parameters.R", sep = "/"))
     
     # Create the results directory
     ResultsDir <- paste(ExpandPrefix, SimDir, sep = "/")
     dir.create(ResultsDir)
     
     # Next, create the population column indices
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     if(Haploid | monoecious){
          NumCols <- 2 + 2^(1-Haploid) * L
     } else{
          NumCols <- 3 + 2 * L
     }
     LocusVec <- 1:L
     AlleleVec <- 1:(2^(1-Haploid) * L)
     
     # Create vectors of random numbers (including sex determination random numbers if
     #    necessary)
     if(Haploid){
          SegVals <- NULL
          SegIndex <- NA
     }else{
          SegVals <- matrix(NA, nrow = 2, ncol = NumRands)
          SegVals[1,] <- sample(c(0,1), replace = TRUE, size = NumRands)
          SegVals[2,] <- sample(c(0,1), replace = TRUE, size = NumRands)
          SegIndex <- 1
     }
     nu <- U / (2^(1-Haploid)*L)
     NumMut <- rbinom(n = NumRands, size = 2^(1-Haploid)*L, prob = nu)
     if( !(is.null(PopIndices$sex)) ){
          SexRands <- rbinom(n = NumRands, size = 1, prob = psi)
     } else{
          SexRands <- NULL
          if(!Haploid){
               SelfRands <- rbinom(n = NumRands, size = 1, prob = omega)
          }
     }
     OffspringIndex <- 1
     direction <- sample(c(-1,1), size = NumRands, replace = TRUE)
     DirectIndex <- 1
     
     # Set up an object to hold summary statistics 
     SumStatCols <- list(x = 1, gen = 2, abund = 3, dBar = 4, GenVar = 5)
     SumStats <- matrix(NA, nrow = SumMatSize, ncol = 5)
     CurStatsDim <- SumMatSize
     SumStatRow <- 1
     
     # Load in the equilibrium population matrix
     InputMat <- read.csv(paste(EquilibriumPrefix, SimDir, "EquilibriumPopMat.csv", sep = "/"))
     PopMat <- as.matrix(InputMat)
     PopSize <- nrow(PopMat)
     
     # Calculate the mean dispersal phenotypes and genetic variances throughout
     #    the population and store them in the first rows of SumStats
     OccPatches <- unique(PopMat[,PopIndices$x0])
     NumPatches <- length(OccPatches)
     if( (SumStatRow + NumPatches) > CurStatsDim){
          NewMat <- matrix(NA, nrow = SumMatSize, ncol = 5)
          SumStats <- rbind(SumStats, NewMat)
          CurStatsDim <- CurStatsDim + SumMatSize
     }
     for(i in 1:NumPatches){
          PatchPop <- which((PopMat[,PopIndices$x0] == OccPatches[i]))
          PatchPopSize <- length(PatchPop)
          if(PatchPopSize > 0){
               SumStats[SumStatRow, SumStatCols$gen] <- 0
               SumStats[SumStatRow, SumStatCols$x] <- OccPatches[i]
               SumStats[SumStatRow, SumStatCols$abund] <- PatchPopSize
               GenVars <- var(PopMat[PatchPop,PopIndices$DispCols])
               SumStats[SumStatRow, SumStatCols$GenVar] <- sum(GenVars[lower.tri(GenVars, diag = TRUE)])
               Disps <- DispPhen(PopMat = PopMat[PatchPop,], PopSize = PatchPopSize, PopIndices = PopIndices, 
                                 Haploid = Haploid, L = L, dmax = dmax, rho = rho, lambda = lambda)
               SumStats[SumStatRow, SumStatCols$dBar] <- mean(Disps)
               SumStatRow <- SumStatRow + 1
          }
     }
     
     # Now run the simulation
     for(g in 1:ExpandGens){
          disps <- DispPhen(PopMat = PopMat, PopSize = PopSize, PopIndices = PopIndices, 
                            Haploid = Haploid, L = L, dmax = dmax, rho = rho, lambda = lambda)
               
          # Check that there are enough entries left in the direction vector and repopulate
          #    it if necessary
          if( (DirectIndex + PopSize) > NumRands ){
               direction <- sample(c(-1,1), size = NumRands, replace = TRUE)
               DirectIndex <- 1
          }
          # Get the new post dispersal locations
          Deltas <- Disperse(d = disps, kern = kern, direction = direction,
                             DirectIndex = DirectIndex, PopSize = PopSize,
                             PopIndices = PopIndices, PopMat = PopMat)
          PopMat[,PopIndices$x1] <- PopMat[,PopIndices$x0] + Deltas
               
          # Update the angle index
          DirectIndex <- DirectIndex + PopSize
          
          # Generate the matrix of occupied patches and get population numbers
          #    for the next generation
          OccPatches <- unique(PopMat[,PopIndices$x1])
          RealizedNtp1 <- Reproduce(CurBeta = BetaInit, gamma = gamma, tau = tau,
                                    R = R, Kmax = Kmax, PopMat = PopMat,
                                    PopIndices = PopIndices, psi = psi,
                                    OccPatches = OccPatches, Haploid = Haploid,
                                    Expand = TRUE, omega = omega)
          
          # Now create a new population matrix for the next generation
          PopSize <- sum(RealizedNtp1)
          if(PopSize > 0){
               # Check that the SegVals vectors contain enough values and resample if not
               if( (!Haploid) & ((SegIndex + PopSize*L) > NumRands) ){
                    SegVals[1,] <- sample(c(0,1), replace = TRUE, size = NumRands)
                    SegVals[2,] <- sample(c(0,1), replace = TRUE, size = NumRands)
                    SegIndex <- 1
               }
               # Check the same for the OffspringIndex
               if( (OffspringIndex + PopSize) > NumRands){
                    NumMut <- rbinom(n = NumRands, size = L, prob = nu)
                    if( !(is.null(PopIndices$sex)) ){
                         SexRands <- rbinom(n = NumRands, size = 1, prob = psi)
                    }
                    OffspringIndex <- 1
               }
               PopMat <- NextPopMat(Ntp1 = RealizedNtp1, PopIndices = PopIndices, 
                                    OccPatches = OccPatches, psi = psi, SumNtp1 = PopSize, 
                                    L = L, PopMat = PopMat, LocusVec = LocusVec,
                                    SegVals = SegVals, SegIndex = SegIndex, NumMutVec = NumMut, 
                                    OffspringIndex = OffspringIndex, AlleleVec = AlleleVec, 
                                    MutStd = sigma, SelfRands = SelfRands, NumCols = NumCols,
                                    SexRands = SexRands)
               # Now update the SegIndices
               SegIndex <- SegIndex + PopSize*L
               # And update the OffspringIndex
               OffspringIndex <- OffspringIndex + PopSize
          } else{
               PopMat <- matrix(NA, nrow = 0, ncol = NumCols)
          }
          
          # Record the summary statistics for every 10th generation
          if((g %/% 10) == (g / 10)){
               NumPatches <- length(OccPatches)
               if( (SumStatRow + NumPatches) > CurStatsDim){
                    NewMat <- matrix(NA, nrow = SumMatSize, ncol = 5)
                    SumStats <- rbind(SumStats, NewMat)
                    CurStatsDim <- CurStatsDim + SumMatSize
               }
               for(i in 1:NumPatches){
                    PatchPop <- which((PopMat[,PopIndices$x0] == OccPatches[i]))
                    PatchPopSize <- length(PatchPop)
                    if(PatchPopSize > 0){
                         SumStats[SumStatRow, SumStatCols$gen] <- g
                         SumStats[SumStatRow, SumStatCols$x] <- OccPatches[i]
                         SumStats[SumStatRow, SumStatCols$abund] <- PatchPopSize
                         GenVars <- var(PopMat[PatchPop,PopIndices$DispCols])
                         SumStats[SumStatRow, SumStatCols$GenVar] <- sum(GenVars[lower.tri(GenVars, diag = TRUE)])
                         Disps <- DispPhen(PopMat = PopMat[PatchPop,], PopSize = PatchPopSize, PopIndices = PopIndices, 
                                           Haploid = Haploid, L = L, dmax = dmax, rho = rho, lambda = lambda)
                         SumStats[SumStatRow, SumStatCols$dBar] <- mean(Disps)
                         SumStatRow <- SumStatRow + 1
                    }
               }
          }
     }
     # Save the summary statistics
     colnames(SumStats) <- c("x", "gen", "abund", "dBar", "GenVar")
     SumStats <- SumStats[1:(SumStatRow - 1),]
     write.csv(SumStats, file = paste(ResultsDir, "SummaryStats.csv", sep = "/"),
               row.names = FALSE, quote = FALSE)
     return(NULL)
}

###### RangeShift
# This function will take an equilibriated population and simulate climate
#    change. The function will save summary statistics for the shifting
#    population.
### INPUTS
# SimDir: the directory with the simulation files for the current simulation
# parallel: A boolean variable indicating whether the simulations are being run
#              on a server or not (which affects how file paths are determined).
### OUTPUT
# The population matrix for the last generation of the simulation
RangeShift <- function(SimDir, parallel = FALSE, SumMatSize = 5000, EquilibriumPrefix = NULL,
                       ShiftPrefix = NULL, NewSpeed = NA){
     # Load the necessary parameters
     source(paste(EquilibriumPrefix, SimDir, "parameters.R", sep = "/"))
     if(!is.na(NewSpeed)){
          v <- NewSpeed
     }
     # Create the results directory path
     ResultsDir <- paste(ShiftPrefix, SimDir, sep = "/")
     dir.create(ResultsDir)
     
     # Next, create the population column indices
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     if(Haploid | monoecious){
          NumCols <- 2 + 2^(1-Haploid) * L
     } else{
          NumCols <- 3 + 2 * L
     }
     LocusVec <- 1:L
     AlleleVec <- 1:(2^(1-Haploid) * L)
     
     # Create vectors of random numbers (including sex determination random numbers if
     #    necessary)
     if(Haploid){
          SegVals <- NULL
          SegIndex <- NA
     }else{
          SegVals <- matrix(NA, nrow = 2, ncol = NumRands)
          SegVals[1,] <- sample(c(0,1), replace = TRUE, size = NumRands)
          SegVals[2,] <- sample(c(0,1), replace = TRUE, size = NumRands)
          SegIndex <- 1
     }
     nu <- U / (2^(1-Haploid)*L)
     NumMut <- rbinom(n = NumRands, size = 2^(1-Haploid)*L, prob = nu)
     if( !(is.null(PopIndices$sex)) ){
          SexRands <- rbinom(n = NumRands, size = 1, prob = psi)
     } else{
          SexRands <- NULL
          if(!Haploid){
               SelfRands <- rbinom(n = NumRands, size = 1, prob = omega)
          }
     }
     OffspringIndex <- 1
     direction <- sample(c(-1,1), size = NumRands, replace = TRUE)
     DirectIndex <- 1
     
     # Load in the equilibrium population matrix
     PopMat <- read.csv(paste(EquilibriumPrefix, SimDir, "EquilibriumPopMat.csv", sep = "/"))
     PopMat <- as.matrix(PopMat)
     PopSize <- nrow(PopMat)
     
     # Set up an object to hold summary statistics 
     SumStatCols <- list(gen = 1, beta = 2, x = 3, abund = 4, dBar = 5,
                         GenVar = 6)
     SumStatRow <- 1
     SumStats <- matrix(NA, nrow = SumMatSize, ncol = 6)
     CurStatsDim <- SumMatSize
     
     # Record the summary statistics for the starting generation
     OccPatches <- unique(PopMat[,PopIndices$x0])
     NumPatches <- length(OccPatches)
     if( (SumStatRow + NumPatches) > CurStatsDim){
          NewMat <- matrix(NA, nrow = SumMatSize, ncol = 6)
          SumStats <- rbind(SumStats, NewMat)
          CurStatsDim <- CurStatsDim + SumMatSize
     }
     for(i in 1:NumPatches){
          PatchPop <- which((PopMat[,PopIndices$x0] == OccPatches[i]))
          PatchPopSize <- length(PatchPop)
          if(PatchPopSize > 0){
               SumStats[SumStatRow, SumStatCols$gen] <- 0
               SumStats[SumStatRow, SumStatCols$beta] <- BetaInit
               SumStats[SumStatRow, SumStatCols$x] <- OccPatches[i]
               SumStats[SumStatRow, SumStatCols$abund] <- PatchPopSize
               GenVars <- var(PopMat[PatchPop,PopIndices$DispCols])
               SumStats[SumStatRow, SumStatCols$GenVar] <- sum(GenVars[lower.tri(GenVars, diag = TRUE)])
               DispPhens <- DispPhen(PopMat = PopMat[PatchPop,], PopSize = PatchPopSize,
                                     PopIndices = PopIndices, Haploid = Haploid,
                                     L = L, dmax = dmax, rho = rho, lambda = lambda)
               SumStats[SumStatRow, SumStatCols$dBar] <- mean(DispPhens)
               SumStatRow <- SumStatRow + 1
          }
     }
     
     
     # Generate a vector of beta values for the climate change
     BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
                                v = v)
     
     # Now run the simulation
     for(g in 1:LengthShift){
          if(PopSize > 0){
               disps <- DispPhen(PopMat = PopMat, PopSize = PopSize, PopIndices = PopIndices, 
                                 Haploid = Haploid, L = L, dmax = dmax, rho = rho, lambda = lambda)
               
               # Check that there are enough entries left in the direction vector and repopulate
               #    it if necessary
               if( (DirectIndex + PopSize) > NumRands ){
                    direction <- sample(c(-1,1), size = NumRands, replace = TRUE)
                    DirectIndex <- 1
               }
               # Get the new post dispersal locations
               Deltas <- Disperse(d = disps, kern = kern, direction = direction,
                                  DirectIndex = DirectIndex, PopSize = PopSize,
                                  PopIndices = PopIndices, PopMat = PopMat)
               PopMat[,PopIndices$x1] <- PopMat[,PopIndices$x0] + Deltas
               
               # Update the angle index
               DirectIndex <- DirectIndex + PopSize
               
               # Generate the matrix of occupied patches and get population numbers
               #    for the next generation
               OccPatches <- unique(PopMat[,PopIndices$x1])
               RealizedNtp1 <- Reproduce(CurBeta = BetaShift[g], gamma = gamma, tau = tau,
                                         R = R, Kmax = Kmax, PopMat = PopMat,
                                         PopIndices = PopIndices, psi = psi,
                                         OccPatches = OccPatches, Haploid = Haploid,
                                         Expand = FALSE, omega = omega)
               
               # Now create a new population matrix for the next generation
               PopSize <- sum(RealizedNtp1)
               if(PopSize > 0){
                    # Check that the SegVals vectors contain enough values and resample if not
                    if( (!Haploid) & ((SegIndex + PopSize*L) > NumRands) ){
                         SegVals[1,] <- sample(c(0,1), replace = TRUE, size = NumRands)
                         SegVals[2,] <- sample(c(0,1), replace = TRUE, size = NumRands)
                         SegIndex <- 1
                    }
                    # Check the same for the OffspringIndex
                    if( (OffspringIndex + PopSize) > NumRands){
                         NumMut <- rbinom(n = NumRands, size = L, prob = nu)
                         if( !(is.null(PopIndices$sex)) ){
                              SexRands <- rbinom(n = NumRands, size = 1, prob = psi)
                         }
                         OffspringIndex <- 1
                    }
                    PopMat <- NextPopMat(Ntp1 = RealizedNtp1, PopIndices = PopIndices, 
                                         OccPatches = OccPatches, psi = psi, SumNtp1 = PopSize, 
                                         L = L, PopMat = PopMat, LocusVec = LocusVec,
                                         SegVals = SegVals, SegIndex = SegIndex, NumMutVec = NumMut, 
                                         OffspringIndex = OffspringIndex, AlleleVec = AlleleVec, 
                                         MutStd = sigma, SelfRands = SelfRands, NumCols = NumCols,
                                         SexRands = SexRands)
                    # Now update the SegIndices
                    SegIndex <- SegIndex + PopSize*L
                    # And update the OffspringIndex
                    OffspringIndex <- OffspringIndex + PopSize
               } else{
                    PopMat <- matrix(NA, nrow = 0, ncol = NumCols)
               }
               
               # Record the summary statistics for the current generation
               NumPatches <- length(OccPatches)
               if( (SumStatRow + NumPatches) > CurStatsDim){
                    NewMat <- matrix(NA, nrow = SumMatSize, ncol = 6)
                    SumStats <- rbind(SumStats, NewMat)
                    CurStatsDim <- CurStatsDim + SumMatSize
               }
               for(i in 1:NumPatches){
                    PatchPop <- which((PopMat[,PopIndices$x0] == OccPatches[i]))
                    PatchPopSize <- length(PatchPop)
                    if(PatchPopSize > 0){
                         SumStats[SumStatRow, SumStatCols$gen] <- g
                         SumStats[SumStatRow, SumStatCols$beta] <- BetaShift[g]
                         SumStats[SumStatRow, SumStatCols$x] <- OccPatches[i]
                         SumStats[SumStatRow, SumStatCols$abund] <- PatchPopSize
                         GenVars <- var(PopMat[PatchPop,PopIndices$DispCols])
                         SumStats[SumStatRow, SumStatCols$GenVar] <- sum(GenVars[lower.tri(GenVars, diag = TRUE)])
                         DispPhens <- DispPhen(PopMat = PopMat[PatchPop,], PopSize = PatchPopSize,
                                               PopIndices = PopIndices, Haploid = Haploid,
                                               L = L, dmax = dmax, rho = rho, lambda = lambda)
                         SumStats[SumStatRow, SumStatCols$dBar] <- mean(DispPhens)
                         SumStatRow <- SumStatRow + 1
                    }
               }
               
               
          }
     }
     # Save the summary statistics
     colnames(SumStats) <- c("gen", "beta", "x", "abund", "dBar", "GenVar")
     SumStats <- SumStats[1:(SumStatRow - 1),]
     write.csv(SumStats, file = paste(ResultsDir, "SummaryStats.csv", sep = "/"),
               row.names = FALSE, quote = FALSE)
     
     # Finally, save the results here
     PopMatNames <- c("x0", "x1", paste("disp1", 1:L, sep = "_"))
     if(!Haploid){
          PopMatNames <- c(PopMatNames, paste("disp2", 1:L, sep = "_"))
     }
     if( !(is.null(PopIndices$sex)) ){
          PopMatNames <- c(PopMatNames, "sex")
     }
     colnames(PopMat) <- PopMatNames
     write.csv(PopMat, file = paste(ResultsDir, "ShiftedPopMat.csv", sep = "/"), 
               row.names = FALSE, quote = FALSE)
     return(NULL)
}

###### ExtraShift
# This function will perform extra range shift simulations to keep track of 
#    extinction events over a larger sample size. Actual population matrices
#    won't be saved, it will merely return a 1 or 0 to indicate extinction (1)
#    or survival (0) of the simulated population
### INPUTS
# SimDir: the directory with the simulation files for the current simulation
# parallel: A boolean variable indicating whether the simulations are being run
#              on a server or not (which affects how file paths are determined).
### OUTPUT
# A boolean variable indicating extinction
ExtraShift <- function(SimDir, parallel = FALSE, EquilibriumPrefix = NULL, NewSpeed = NA){
     # Load the necessary parameters
     source(paste(EquilibriumPrefix, SimDir, "parameters.R", sep = "/"))
     if(!is.na(NewSpeed)){
          v <- NewSpeed
     }
     
     # Next, create the population column indices
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     if(Haploid | monoecious){
          NumCols <- 2 + 2^(1-Haploid) * L
     } else{
          NumCols <- 3 + 2 * L
     }
     LocusVec <- 1:L
     AlleleVec <- 1:(2^(1-Haploid) * L)
     
     # Create vectors of random numbers (including sex determination random numbers if
     #    necessary)
     if(Haploid){
          SegVals <- NULL
          SegIndex <- NA
     }else{
          SegVals <- matrix(NA, nrow = 2, ncol = NumRands)
          SegVals[1,] <- sample(c(0,1), replace = TRUE, size = NumRands)
          SegVals[2,] <- sample(c(0,1), replace = TRUE, size = NumRands)
          SegIndex <- 1
     }
     nu <- U / (2^(1-Haploid)*L)
     NumMut <- rbinom(n = NumRands, size = 2^(1-Haploid)*L, prob = nu)
     if( !(is.null(PopIndices$sex)) ){
          SexRands <- rbinom(n = NumRands, size = 1, prob = psi)
     } else{
          SexRands <- NULL
          if(!Haploid){
               SelfRands <- rbinom(n = NumRands, size = 1, prob = omega)
          }
     }
     OffspringIndex <- 1
     direction <- sample(c(-1,1), size = NumRands, replace = TRUE)
     DirectIndex <- 1
     
     # Load in the equilibrium population matrix
     PopMat <- read.csv(paste(EquilibriumPrefix, SimDir, "EquilibriumPopMat.csv", sep = "/"))
     PopMat <- as.matrix(PopMat)
     PopSize <- nrow(PopMat)
     
     # Generate a vector of beta values for the climate change
     BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
                                v = v)
     
     # Now run the simulation
     for(g in 1:LengthShift){
          if(PopSize > 0){
               disps <- DispPhen(PopMat = PopMat, PopSize = PopSize, PopIndices = PopIndices, 
                                 Haploid = Haploid, L = L, dmax = dmax, rho = rho, lambda = lambda)
               
               # Check that there are enough entries left in the direction vector and repopulate
               #    it if necessary
               if( (DirectIndex + PopSize) > NumRands ){
                    direction <- sample(c(-1,1), size = NumRands, replace = TRUE)
                    DirectIndex <- 1
               }
               # Get the new post dispersal locations
               Deltas <- Disperse(d = disps, kern = kern, direction = direction,
                                  DirectIndex = DirectIndex, PopSize = PopSize,
                                  PopIndices = PopIndices, PopMat = PopMat)
               PopMat[,PopIndices$x1] <- PopMat[,PopIndices$x0] + Deltas
               
               # Update the angle index
               DirectIndex <- DirectIndex + PopSize
               
               # Generate the matrix of occupied patches and get population numbers
               #    for the next generation
               OccPatches <- unique(PopMat[,PopIndices$x1])
               RealizedNtp1 <- Reproduce(CurBeta = BetaShift[g], gamma = gamma, tau = tau,
                                         R = R, Kmax = Kmax, PopMat = PopMat,
                                         PopIndices = PopIndices, psi = psi,
                                         OccPatches = OccPatches, Haploid = Haploid,
                                         Expand = FALSE, omega = omega)
               
               # Now create a new population matrix for the next generation
               PopSize <- sum(RealizedNtp1)
               if(PopSize > 0){
                    # Check that the SegVals vectors contain enough values and resample if not
                    if( (!Haploid) & ((SegIndex + PopSize*L) > NumRands) ){
                         SegVals[1,] <- sample(c(0,1), replace = TRUE, size = NumRands)
                         SegVals[2,] <- sample(c(0,1), replace = TRUE, size = NumRands)
                         SegIndex <- 1
                    }
                    # Check the same for the OffspringIndex
                    if( (OffspringIndex + PopSize) > NumRands){
                         NumMut <- rbinom(n = NumRands, size = L, prob = nu)
                         if( !(is.null(PopIndices$sex)) ){
                              SexRands <- rbinom(n = NumRands, size = 1, prob = psi)
                         }
                         OffspringIndex <- 1
                    }
                    PopMat <- NextPopMat(Ntp1 = RealizedNtp1, PopIndices = PopIndices, 
                                         OccPatches = OccPatches, psi = psi, SumNtp1 = PopSize, 
                                         L = L, PopMat = PopMat, LocusVec = LocusVec,
                                         SegVals = SegVals, SegIndex = SegIndex, NumMutVec = NumMut, 
                                         OffspringIndex = OffspringIndex, AlleleVec = AlleleVec, 
                                         MutStd = sigma, SelfRands = SelfRands, NumCols = NumCols,
                                         SexRands = SexRands)
                    # Now update the SegIndices
                    SegIndex <- SegIndex + PopSize*L
                    # And update the OffspringIndex
                    OffspringIndex <- OffspringIndex + PopSize
               } else{
                    PopMat <- matrix(NA, nrow = 0, ncol = NumCols)
               }
          }
     }
     Extinct <- ifelse(PopSize > 0, 0, 1)
     return(Extinct)
}

###### NoEvolShift
# This function will perform range shift simulations without evolution to quantify
#    the impact of evolution on extinction risk. Actual population matrices
#    won't be saved, it will merely return a 1 or 0 to indicate extinction (1)
#    or survival (0) of the simulated population
### INPUTS
# SimDir: the directory with the simulation files for the current simulation
# parallel: A boolean variable indicating whether the simulations are being run
#              on a server or not (which affects how file paths are determined).
### OUTPUT
# A boolean variable indicating extinction
NoEvolShift <- function(SimDir, parallel = FALSE, EquilibriumPrefix = NULL, NewSpeed = NA){
     # Load the necessary parameters
     source(paste(EquilibriumPrefix, SimDir, "parameters.R", sep = "/"))
     if(!is.na(NewSpeed)){
          v <- NewSpeed
     }
     
     # Next, create the population column indices
     PopIndices <- PopMatColNames(L = L, monoecious = monoecious, Haploid = Haploid)
     if(Haploid | monoecious){
          NumCols <- 2 + 2^(1-Haploid) * L
     } else{
          NumCols <- 3 + 2 * L
     }
     
     # Create vectors of random numbers for sex determination and dispersal direction
     #    Since evolution is prevented in this simulation, the other vectors
     #         of random numbers governing allele inheritence and mutation are
     #         unnecessary.
     if( !(is.null(PopIndices$sex)) ){
          SexRands <- rbinom(n = NumRands, size = 1, prob = psi)
     }
     direction <- sample(c(-1,1), size = NumRands, replace = TRUE)
     DirectIndex <- 1
     OffspringIndex <- 1
     
     # Load in the equilibrium population matrix
     PopMat <- read.csv(paste(EquilibriumPrefix, SimDir, "EquilibriumPopMat.csv", sep = "/"))
     PopMat <- as.matrix(PopMat)
     PopSize <- nrow(PopMat)
     
     # Store the initial values for each allele to be sampled from later so that
     #    allele frequencies always are drawn from the founding individuals
     InitAlleles <- as.matrix(PopMat[,PopIndices$DispCols])
     
     # Generate a vector of beta values for the climate change
     BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
                                v = v)
     
     # Now run the simulation
     for(g in 1:LengthShift){
          if(PopSize > 0){
               disps <- DispPhen(PopMat = PopMat, PopSize = PopSize, PopIndices = PopIndices, 
                                 Haploid = Haploid, L = L, dmax = dmax, rho = rho, lambda = lambda)
               
               # Check that there are enough entries left in the direction vector and repopulate
               #    it if necessary
               if( (DirectIndex + PopSize) > NumRands ){
                    direction <- sample(c(-1,1), size = NumRands, replace = TRUE)
                    DirectIndex <- 1
               }
               # Get the new post dispersal locations
               Deltas <- Disperse(d = disps, kern = kern, direction = direction,
                                  DirectIndex = DirectIndex, PopSize = PopSize,
                                  PopIndices = PopIndices, PopMat = PopMat)
               PopMat[,PopIndices$x1] <- PopMat[,PopIndices$x0] + Deltas
               
               # Update the angle index
               DirectIndex <- DirectIndex + PopSize
               
               # Generate the matrix of occupied patches and get population numbers
               #    for the next generation
               OccPatches <- unique(PopMat[,PopIndices$x1])
               RealizedNtp1 <- Reproduce(CurBeta = BetaShift[g], gamma = gamma, tau = tau,
                                         R = R, Kmax = Kmax, PopMat = PopMat,
                                         PopIndices = PopIndices, psi = psi,
                                         OccPatches = OccPatches, Haploid = Haploid,
                                         Expand = FALSE, omega = omega)
               
               # Now create a new population matrix for the next generation
               PopSize <- sum(RealizedNtp1)
               if(PopSize > 0){
                    # Check the same for the OffspringIndex
                    if( (OffspringIndex + PopSize) > NumRands){
                         if( !(is.null(PopIndices$sex)) ){
                              SexRands <- rbinom(n = NumRands, size = 1, prob = psi)
                              OffspringIndex <- 1
                         }
                    }
                    PopMat <- NextPopMatNoEvol(Ntp1 = RealizedNtp1, PopIndices = PopIndices,
                                              OccPatches = OccPatches, SumNtp1 = PopSize, 
                                              NumCols = NumCols, SexRands = SexRands,
                                              InitAlleles = InitAlleles, OffspringIndex = OffspringIndex)
                    # And update the OffspringIndex
                    OffspringIndex <- OffspringIndex + PopSize
               } else{
                    PopMat <- matrix(NA, nrow = 0, ncol = NumCols)
               }
          }
     }
     Extinct <- ifelse(PopSize > 0, 0, 1)
     return(Extinct)
}
