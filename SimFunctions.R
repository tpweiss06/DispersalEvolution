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
# NumLoc: The number of loci defining dispersal
# dmax: The maximum expected dispersal distance (in units of discrete patches)
# eta: The dimensions of each discrete patch in Cartesian space
# rho: Determines the slope of the transition from 0 to dmax as the loci sum
#         changes
### OUTPUTS
# A vector of length PopSize with the i-th entry corresponding to the dispersal
#    phenotype of the i-th individual in the population input
DispPhen <- function(PopMat, PopSize, PopIndices, Haploid, NumLoc,
                     dmax, eta, rho){
     # Calculate the sum of the quantitative loci defining dispersal for
     #    each individual, accounting for the different possible cases
     if(PopSize == 1){
          LociSums <- sum(PopMat[,PopIndices$DispCols])
     } else{
          if((Haploid == 1) & (NumLoc == 1)){
               LociSums <- PopMat[,PopIndices$DispCols]
          } else{
               LociSums <- rowSums(PopMat[,PopIndices$DispCols])
          }
     }
     # Now calculate the expected dispersal distance for each individual
     d <- (dmax * eta * exp(rho * LociSums)) / (1 + exp(rho * LociSums))
     return(d)
}

###### Dispersal
# Using the dispersal phenotypes from the above function, this function will
#    implement dispersal according to one of three different kernels and return
#    the realized changes in patch number in the x and y dimensions.
### INPUTS:
# d:           The dispersal phenotypes
# width:       The number of patches spanning the y dimension of the landscape
# kern:        The dispersal kernel to be used. Can be "norm", "exp", or 
#                   "stud_t"
# angles:      A large vector of random numbers previously generated
# AngleIndex:  Index of the where to look within the vector of random numbers
# PopIndices, PopSize, PopMat, eta, rho, and dmax:  As previously defined
### OUTPUTS:
# The function will return a matrix with 2 columns and as many rows as individuals
#    within the population. Each row will correspond to the post dispersal 
#    change in patch number for X and Y dimensions in that order.
Disperse <- function(d, width, kern, eta, angles, AngleIndex, rho,
                     PopSize, PopIndices, PopMat, dmax){
     # Generate the dispersal distances according to the type of kernel
     if(kern == "norm"){
          sigma <- (d^2 * pi) / 2
          dists <- abs(rnorm(n = PopSize, mean = 0, sd = sigma))
     } else if(kern == "exp"){
          dists <- rexp(n = PopSize, rate = 1 / d)
     } else if(kern == "stud_t"){
          dists <- rStudKern(n = PopSize, d = d)
     }
     
     # Calculate the new x and y coordinates of each individual as if the 
     #    center of their current patch is the origin
     NewX <- dists * cos(angles[AngleIndex:(AngleIndex + PopSize - 1)])
     NewY <- dists * sin(angles[AngleIndex:(AngleIndex + PopSize - 1)])
     
     # Use the eta parameter to determine the number of x and y
     #    patches changed, accounting for movement occuring from the
     #    center by adding or subtracting eta/2 as necessary
     DeltaX <- ifelse(NewX < 0, ceiling((NewX - eta/2) / eta), floor((NewX + eta/2) / eta))
     DeltaY <- ifelse(NewY < 0, ceiling((NewY - eta/2) / eta), floor((NewY + eta/2) / eta))
     return(cbind(DeltaX,DeltaY))
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
#    OccPatches: A matrix with the unique x, y identifiers of each occupied patch
#                   organized into columns (column 1 = x; column 2 = y)
#    z: The expected sex ratio in the population
#    RangeParams
#    PopMat, Haploid, PopIndices as previously defined
### OUTPUTS
# A vector of the realized abundances of each occupied patch in the next generation
Reproduce <- function(R, OccPatches, PopMat, Haploid, PopIndices, z, RangeParams){
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
     NumPatches <- nrow(OccPatches)
     # Get the patch alpha values based on the location within the range
     Alphas <- CalcPatchK(patches = OccPatches[,1], RangeParams) #### Need to replace with necessary range params
     # Create an empty object to hold the expected population sizes for the next generation
     Ntp1 <- rep(NA, NumPatches)
     
     # Loop through each occupied patch and calculate the expected population
     #    size in the next generation
     for(i in 1:NumPatches){
          # Extract the local population
          CurPatchPop <- which( (PopMat[,PopIndices$x1] == OccPatches[i,1]) & 
                                     (PopMat[,PopIndices$y1] == OccPatches[i,2]) )
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
                         Ntp1[i] <- PatchAbund * R * exp(-1 * Alphas[i] * PatchAbund)
                    } else{
                         CurFemales <- which( (PopMat[,PopIndices$x1] == OccPatches[i,1]) & 
                                                   (PopMat[,PopIndices$y1] == OccPatches[i,2]) &
                                                   (PopMat[,PopIndices$sex] == 1) )
                         nFem <- length(CurFemales)
                         # Calculate expected population growth with a forced 0
                         #    if the whole population is female
                         AllFemale <- NumFemales == PatchAbund
                         Ntp1[i] <- (nFem * (R / z) * exp(-1 * Alphas[i] * PatchAbund)) * 
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
# Ld:       The number of loci defining the dispersal trait
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
# PopMat, PopIndices, OccPatches, and z as defined previously
### OUTPUTS
# This function returns a new population matrix with the same columns but with
#    updated row information for the next generation.
NextPopMat <- function(Ntp1, PopIndices, OccPatches, z, SumNtp1, Ld, PopMat, LocusVec,
                       SegVals, SegIndex, NumMutVec, MutIndex, AlleleVec, MutStd,
                       SelfRands){
     NewPopMat <- matrix(NA, nrow = SumNtp1, ncol = NumCols)
     # Filter for only the patches that produced offspring and fill in the matrix
     NewOccPatches <- which(Ntp1 > 0)
     Ntp1 <- Ntp1[NewOccPatches]
     StartRow <- 1
     EndRow <- Ntp1[1]
     for(i in 1:length(NewOccPatches)){
          CurVec <- StartRow:EndRow
          # Fill in the location details
          NewPopMat[CurVec, PopIndices$x0] <- OccPatches[NewOccPatches[i], 1]
          NewPopMat[CurVec, PopIndices$y0] <- OccPatches[NewOccPatches[i], 2]
          # Identify the potential parents in the current patch
          locals <- which( (PopMat[,PopIndices$x1] == OccPatches[NewOccPatches[i], 1]) &
                                PopMat[,PopIndices$y1] == OccPatches[NewOccPatches[i], 2])
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
                                                                    MutIndex = MutIndex, MutStd = MutStd, AlleleVec = AlleleVec)
          } else{
               if(is.null(PopIndices$sex)){
                    if(NumLocals == 1){
                         parent1 <- rep(locals, Ntp1[i])
                         parent2 <- rep(locals, Ntp1[i])
                    } else{
                         parent1 <- sample(locals, size = Ntp1[i], replace = TRUE)
                         SelfFert <- SelfRands[(OffspringIndex + StartRow):(OffspringIndex + EndRow)]
                         parent2 <- sample(locals, size = Ntp1[i], replace = TRUE)
                         SameParent <- which(parent1 == parent2)
                         for(j in SameParent){
                              ParentPool <- setdiff(locals, parent1[j])
                              parent2[j] <- sample(ParentPool, size = 1)
                         }
                         parent2 <- ifelse(SelfFert == 1, parent1, parent2)
                    }
                    parents <- cbind(parent1, parent2)
               } else{
                    NewPopMat[CurVec, PopIndices$sex] <- SexRands[(OffspringIndex + StartRow):(OffspringIndex + EndRow)]
                    # Identify the females and males present in the current patch
                    females <- which(PopMat[locals, PopIndices$sex] == 1)
                    males <- which(PopMat[locals, PopIndices$sex] == 0)
                    parent1 <- ifelse(length(females == 1), rep(females, Ntp1[i]),
                                      sample(females, size = Ntp1[i], replace = TRUE))
                    parent2 <- ifelse(length(males) == 1, rep(males, Ntp1[i]),
                                      sample(males, size = Ntp1[i], replace = TRUE))
                    parents <- cbind(parent1, parent2)
               }
               NewPopMat[CurVec, PopIndices$DispCols] <- DipInherit(PatchPop = Ntp1[i], Ld = Ld, PopMat = PopMat,
                                                                    parents = parents, SegVals = SegVals, SegIndex = SegIndex,
                                                                    NumMutVec = NumMutVec, MutIndex = MutIndex, 
                                                                    AlleleVec = AlleleVec, MutStd = MutStd, LocusVec = LocusVec)
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
     # Mutate the affected alleles
     MutOffspring <- which(NumMut != 0)
     for(i in MutOffspring){
          MutLocus <- sample(AlleleVec, size = NumMut[i], replace = TRUE)
          OffspringAlleles[i, MutLocus] <- rnorm(mean = OffspringAlleles[i, MutLocus],
                                                 sd = MutStd, n = NumMut[i])
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
DipInherit <- function(PopIndices, parents, PopMat, PatchPop, Ld, SegVals, AlleleVec,
                        MutStd, NumMutVec, MutIndex, LocusVec, SegIndex){
     OffspringAlleles <- matrix(NA, nrow = PatchPop, ncol = 2*Ld)
     for(i in 1:PatchPop){
          ParentLoci <- PopMat[parents[i,], PopIndices$DispCols]
          Parent1Alleles <- LocusVec + Ld * SegVals[1,(SegIndex + (i - 1) * Ld):(SegIndex + i * Ld - 1)]
          Parent2Alleles <- LocusVec + Ld * SegVals[2,(SegIndex + (i - 1) * Ld):(SegIndex + i * Ld - 1)]
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

# ------------------------------------------------------------------------------
# -------------------------- Environmental Functions ---------------------------
# ------------------------------------------------------------------------------

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
ChangeClimate <- function(BetaInit, LengthShift, v, eta){
     # First check for a negative speed
     if(v < 0){
          write("Negative speed is not supported for climate change", stderr())
          return(NULL)
     }
     # Next check for a 0 value for the length of climate change
     if(LengthShift == 0){
          return(NULL)
     }
     # Perform the climate shift for the range center
     TimeVec <- 1:LengthShift
     BetaVec <- BetaInit + v * eta * TimeVec
     return(BetaVec)
}

###### CalcEnvMean
# This function will use the relevant environmental parameters to calculate the
#    mean of the range capacity function over a given interval using the Mean
#    Value Theorem.
### INPUTS
# beta, gamma, and tau are all parameters in the range capacity function.
#    See the function write up for details on these parameters.
# a:      Lower interval bound(s)
# b:      Upper interval bound(s)
### OUTPUS
# This function will output a single value for the mean value of the function
#    over the given interval.
CalcEnvMean <- function(beta, gamma, tau, a, b){
     # Determine the number of patches to calculate the mean for, then loop
     #    over them
     NumPatches <- length(a)
     EnvMean <- rep(NA, NumPatches)
     for(i in 1:NumPatches){
          # First determine where the interval falls in the range, then calculate
          #    the appropriate mean (See write up)
          if(b[i] <= beta){
               numerator <- exp(gamma * (b[i] - beta + tau)) + 1
               denominator <- exp(gamma * (a[i] - beta + tau)) + 1
               FullIntegral <- (1 / gamma) * log(numerator / denominator)
               EnvMean[i] <- (1 / (b[i] - a[i])) * FullIntegral
          } else if(a[i] >= beta){
               numerator <- exp(-1 * gamma * (b[i] - beta - tau)) + 1
               denominator <- exp(-1 * gamma * (a[i] - beta - tau)) + 1
               FullIntegral <- -1 * (1 / gamma) * log(numerator / denominator)
               EnvMean[i] <- (1 / (b[i] - a[i])) * FullIntegral
          } else{
               PreNum <- exp(gamma * tau) + 1
               PreDen <- exp(gamma * (a[i] - beta + tau)) + 1
               PostNum <- exp(-1 * gamma * (b[i] - beta - tau)) + 1
               PostDen <- exp(gamma * tau) + 1
               FullIntegral <- (1 / gamma) * log(PreNum / PreDen) + 
                    -1 * (1 /gamma) * log(PostNum / PostDen)
               EnvMean[i] <- (1 / (b[i] - a[i])) * FullIntegral
          }
     }
     return(EnvMean)
}

###### GetEnvQual
# This function uses the CalcEnvMean function to extract a vector of
#    environmental means for a given set of patches. 
### INPUTS
# beta, gamma, and tau are all parameters in the range capacity function.
#    See the function write up for details on these parameters.
# patches:          A vector of the patch numbers to get quality information for
# eta:       A single constant value used to modify the size of patches.
#                        By default this will be set to 1, but it could in 
#                        principle be used to explore the consequences of 
#                        discretizing space by letting it become arbitrarily
#                        small or alternatively it could allow a given range
#                        to contain less patches and thus save computational
#                        power or explore the effects of essentially lowering
#                        the population size.
### OUTPUTS
# This function will generate a vector of environmental quality values for each
#    patch as determined by the GetEnvMean function.
GetEnvQual <- function(beta, gamma, tau, patches, eta = 1){
     # Calculate the lower and upper bounds of each patch on the continuous scale
     #    by first calculating their center points
     centers <- patches * eta
     lowers <- centers - eta * 0.5
     uppers <- centers + eta * 0.5
     
     # Now get and return the environmental quality score for each patch
     EnvQuals <- CalcEnvMean(beta = beta, gamma = gamma, tau = tau, 
                             a = lowers, b = uppers)
     return(EnvQuals)
}


# ------------------------------------------------------------------------------
# --------------------------- Bookkeeping Functions ----------------------------
# ------------------------------------------------------------------------------

###### Initialize
# This function will initialize a population matrix of founders to start a 
#    simulation.
### INPUTS
# PopIndices: A list of column indices for different parts of the matrix
# NumCols:          The number of columns in the population matrix
# FitCols:     The column indices for the fitness loci
# DispCols:    The column indices for the dispersal loci
# PopSize:     The number of founders to create in the intitial population
# BetaInit:    The starting location of the range center where individuals will
#                   initialize
# z:    The sex ratio of the founding population (if dioecious). Set to
#                   0.5 by default
# FitInit:     The mean initial value for all fitness loci
# FitDiv:      The standard deviation of the initial distribution of fitness 
#                   loci
# DispInit:    The mean initial value for all dispersal loci
# DispDiv:     The standard deviation of the initial distribution of dispersal 
#                   loci
### OUTPUS
# A filled in population matrix to start generation 0
Initialize <- function(PopIndices, PopSize, BetaInit, width, 
                       z = 0.5, FitInit, FitDiv, DispInit, DispDiv, NumCols){
     # First make an empty population matrix with the correct names
     PopMat <- matrix(NA, nrow = PopSize, ncol = NumCols)
     
     # Next, fill in the x1 and y1 columns for the founders (the founders are
     #    considered post dispersal and will reproduce next)
     PopMat[,PopIndices$x1] <- BetaInit
     PopMat[,PopIndices$y1] <- sample(1:width, size = PopSize, replace = TRUE)
     
     # Fill in the sex column if it is present
     if( !(is.null(PopIndices$sex)) ){
          PopMat[,PopIndices$sex] <- rbinom(n = PopSize, size = 1, prob = z)
     }
     
     # Now fill in the fitness and dispersal allele columns
     PopMat[,PopIndices$FitCols] <- rnorm(n = PopSize * length(PopIndices$FitCols), mean = FitInit,
                               sd = FitDiv)
     PopMat[,PopIndices$DispCols] <- rnorm(n = PopSize * length(PopIndices$DispCols), mean = DispInit,
                                sd = DispDiv)
     
     # Return the filled in initial population matrix
     return(PopMat)
}

###### PopMatColNames
# This function generates the column names to use for a given simulation
### INPUTS
# Lf:        the number of loci defining an individual's fitness
# Ld:       the number of loci defining an individual's dispersal ability
# monoecious:  a boolean value indicating whether the individuals in the simulation
#              are monoecious (i.e. possessing both the female and male 
#              reproductive organs in the same individual) or not.
# example:     a boolean value indicating whether the function should simply 
#                   return an example character vector containing the column 
#                   names and order used
### OUTPUTS
# This function creates a list with the indices for different column types in the
#    population matrix. If example is TRUE, it simply returns a character vector.
PopMatColNames <- function(Lf, Ld, monoecious, example = FALSE){
     if(example){
          MatNames <- c("x0", "y0", "x1", "y1", "sex", "fit1_1", "fit1_2", "...", "fit1_N",
                        "fit2_1", "...", "fit2_N", "disp1_1", "...", "disp2_N")
          return(MatNames)
     }
     
     # Create the list and population it
     PopIndices <- vector(mode = "list", length = 7)
     names(PopIndices) <- c("x0", "y0", "x1", "y1", "sex", "FitCols", "DispCols")
     PopIndices$x0 <- 1
     PopIndices$y0 <- 2
     PopIndices$x1 <- 3
     PopIndices$y1 <- 4
     PopIndices$FitCols <- 5:(5+2*Lf-1)
     PopIndices$DispCols <- (5 + 2*Lf):(5 + 2*Lf + 2*Ld - 1)
     if(monoecious){
          PopIndices$sex <- NULL
     } else{
          PopIndices$sex <- 5 + 2*Lf + 2*Ld
     }
     
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
     # Set the working directory to easily scan for directory names
     setwd(ParentDirectory)
     
     # Start by trying out 1 and then increase as necessary
     newID <- 1
     
     # If the simulations are being done in parallel, include the name of the
     #    current node in the new directory name
     if(parallel){
          NodeName <- Sys.getpid()
          DirName <- paste(NodeName, "_", "Sim", newID, sep = "")
          # Until we find a directory name that does not already exist, continue
          #    to increase the ID variable
          while(dir.exists(DirName)){
               newID <- newID + 1
               DirName <- paste(NodeName, "_", "Sim", newID, sep = "")
          }
     } else{
          DirName <- paste('Sim', newID, sep='')
          while(dir.exists(DirName)){
               newID <- newID + 1
               DirName <- paste('Sim', newID, sep='')
          }
     }
     
     # Now paste together the ParentDirectory with the new directory for the
     #    full file path
     FullPath <- paste(ParentDirectory, DirName, sep = "/")
     return(FullPath)
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
     ParamCheck <- names(parameters) == c("BetaInit", "gamma", "tau", "lambda", "omega",
                                          "U", "Vm", "Lf", "Ld", "Rmax", "Kmax", "width",
                                          "kern", "EnvGradType", "monoecious", "BurnIn",
                                          "BurnOut", "LengthShift", "v", "InitPopSize",
                                          "FitInit", "FitDiv", "DispInit", "DispDiv", "eta",
                                          "NumRands", "z", "dmax", "rho")
     if(sum(ParamCheck) != length(ParamCheck)){
          write("Incorrect names or number of input parameters", stderr())
          return(NULL)
     }
     
     OutFile <- paste(FilePath, "parameters.R", sep = "/")
     sink(OutFile)
     cat("BetaInit <- ", parameters$BetaInit, "\n", sep = "")
     cat("gamma <- ", parameters$gamma, "\n", sep = "")
     cat("tau <- ", parameters$tau, "\n", sep = "")
     cat("lambda <- ", parameters$lambda, "\n", sep = "")
     cat("omega <- ", parameters$omega, "\n", sep = "")
     cat("U <- c(", parameters$U[1], ",", parameters$U[2], ")\n", sep = "")
     cat("Vm <- c(", parameters$Vm[1], ",", parameters$Vm[2], ")\n", sep = "")
     cat("Lf <- ", parameters$Lf, "\n", sep = "")
     cat("Ld <- ", parameters$Ld, "\n", sep = "")
     cat("Rmax <- ", parameters$Rmax, "\n", sep = "")
     cat("Kmax <- ", parameters$Kmax, "\n", sep = "")
     cat("width <- ", parameters$width, "\n", sep = "")
     cat("kern <- \"", parameters$kern, "\"\n", sep = "")
     cat("EnvGradType <- \"", parameters$EnvGradType, "\"\n", sep = "")
     cat("monoecious <- ", parameters$monoecious, "\n", sep = "")
     cat("BurnIn <- ", parameters$BurnIn, "\n", sep = "")
     cat("BurnOut <- ", parameters$BurnOut, "\n", sep = "")
     cat("LengthShift <- ", parameters$LengthShift, "\n", sep = "")
     cat("v <- ", parameters$v, "\n", sep = "")
     cat("InitPopSize <- ", parameters$InitPopSize, "\n", sep = "")
     cat("FitInit <- ", parameters$FitInit, "\n", sep = "")
     cat("FitDiv <- ", parameters$FitDiv, "\n", sep = "")
     cat("DispInit <- ", parameters$DispInit, "\n", sep = "")
     cat("DispDiv <- ", parameters$DispDiv, "\n", sep = "")
     cat("eta <- ", parameters$eta,  "\n", sep = "")
     cat("NumRands <- ", parameters$NumRands, "\n", sep = "")
     cat("z <- ", parameters$z, "\n", sep = "")
     cat("dmax <- ", parameters$dmax, "\n", sep = "")
     cat("rho <- ", parameters$rho, "\n", sep = "")
     sink()
}


# ------------------------------------------------------------------------------
# ------------------------ Full Simulation Function ----------------------------
# ------------------------------------------------------------------------------

###### FullSim
# Currently this is only mean to stand in as a template so that I can think
#    through the rest of the functions I will need.
### INPUTS
# parameters:  A list with all necessary parameter values for a single model
#                   run. These will be unpacked and stored by the SaveParams
#                   function. It will include all the parameters necessary for
#                   the following functions as well as a few additional time
#                   keeping parameters: 
#                   Functions --
#                   PopMatColNames(), Initialize(), ChangeClimate(), Reproduce(), 
#                        Disperse(), and GetSafeID()
#                   Time keeping parameters --
#                   BurnIn:        Length of time pre shift
#                   LengthShift:   Duration of the shift
#                   BurnOut:       Length of time post shift
# parallel: A boolean variable indicating whether the simulations are being run
#              on a server or not (which affects how file paths are determined).
### OUTPUTS
FullSim <- function(parameters, parallel = FALSE, SumMatSize = 5000, PopInit = NULL, SimID = NA){
     # First generate a safe directory name and create it to save all output
     #    from the simulation
     CurDirectory <- getwd()
     if(is.na(SimID)){
          ResultsDir <- GetSafeID(ParentDirectory = CurDirectory, parallel = parallel)
          
     } else{
          ResultsDir <- paste(CurDirectory, SimID, sep = "/")
     }
     dir.create(ResultsDir)
     
     # Next, save the parameters used for this simulation and source the file
     #    to have access to them within the function
     SaveParams(parameters = parameters, FilePath = ResultsDir)
     source(paste(ResultsDir, "parameters.R", sep = "/"))
     
     # Next, create the population column indices
     PopIndices <- PopMatColNames(Lf = Lf, Ld = Ld, monoecious = monoecious)
     NumCols <- 2*Lf + 2*Ld + 4 + (1-monoecious) * 1
     FitLocusVec <- 1:Lf
     DispLocusVec <- 1:Ld
     FitAlleleVec <- 1:(2*Lf)
     DispAlleleVec <- 1:(2*Ld)
     
     # Create vectors of random numbers (including sex determination random numbers if
     #    necessary)
     FitSegVals1 <- sample(c(0,1), replace = TRUE, size = NumRands)
     FitSegVals2 <- sample(c(0,1), replace = TRUE, size = NumRands)
     DispSegVals1 <- sample(c(0,1), replace = TRUE, size = NumRands)
     DispSegVals2 <- sample(c(0,1), replace = TRUE, size = NumRands)
     FitSegIndex <- 1
     DispSegIndex <- 1
     
     FitPerLocusProb <- U[1] / (2*Lf)
     FitMutStd <- sqrt(Vm[1] / U[1])
     DispPerLocusProb <- U[2] / (2*Ld)
     DispMutStd <- ifelse(U[2] == 0, 0, sqrt(Vm[2] / U[2]))
     FitNumMut <- rbinom(n = NumRands, size = Lf, prob = FitPerLocusProb)
     DispNumMut <- rbinom(n = NumRands, size = Ld, prob = DispPerLocusProb)
     if( !(is.null(PopIndices$sex)) ){
          SexRands <- rbinom(n = NumRands, size = 1, prob = z)
     }
     OffspringIndex <- 1
     
     angles <- runif(n = NumRands, min = 0, max = 2*pi)
     AngleIndex <- 1
     
     # Next initialize the generation 0 founding population and allow it to
     #    reproduce
     if(is.null(PopInit)){
          PopMat <- Initialize(PopIndices = PopIndices, PopSize = InitPopSize, NumCols = NumCols,
                               BetaInit = BetaInit, FitInit = FitInit, FitDiv = FitDiv, 
                               DispInit = DispInit, DispDiv = DispDiv, width = width)
     } else{
          PopMat <- PopInit
     }
     CurPop <- 1:nrow(PopMat)
     PopSize <- length(CurPop)
     traits <- CalcTraits(population = CurPop, PopMat = PopMat, PopSize = PopSize,
                          PopIndices = PopIndices)
     
     # Generate the matrix of occupied patches and get the relative fitness of
     #    each individual
     OccPatches <- unique(PopMat[CurPop,c(PopIndices$x1, PopIndices$y1)])
     RelFits <- RelFit(lambda = lambda, beta = BetaInit, traits = traits, 
                       PopMat = PopMat, individuals = CurPop, omega = omega,
                       PopIndices = PopIndices, eta = eta)
     # Allow the population to reproduce
     RealizedNtp1 <- Reproduce(beta = BetaInit, gamma = gamma, tau = tau, omega = omega,
                         Rmax = Rmax, Kmax = Kmax, traits = traits, PopMat = PopMat, EnvGradType = EnvGradType,
                         PopIndices = PopIndices, z = 0.5, eta = eta,
                         SexRands = SexRands, CurPop = CurPop, Ld = Ld, Lf = Lf, 
                         DispSegVals1 = DispSegVals1, DispSegVals2 = DispSegVals2, 
                         DispSegIndex = DispSegIndex, FitSegVals1 = FitSegVals1, FitSegVals2 = FitSegVals2, 
                         FitSegIndex = FitSegIndex, FitMutStd = FitMutStd, DispMutStd = DispMutStd, 
                         FitNumMut = FitNumMut, DispNumMut = DispNumMut, OffspringIndex = OffspringIndex, 
                         FitLocusVec = FitLocusVec, DispLocusVec = DispLocusVec, NumCols = NumCols, 
                         FitPerLocusProb = FitPerLocusProb, DispPerLocusProb = DispPerLocusProb, NumRands = NumRands,
                         OccPatches = OccPatches, RelFits = RelFits)
     
     # Now update the population matrix
     SumNtp1 <- sum(RealizedNtp1)
     if(SumNtp1 > 0){
          PopMat <- MatFill(RealizedNtp1 = RealizedNtp1, PopIndices = PopIndices,
                             OccPatches = OccPatches, z = z, PopMat = PopMat,
                             RelFits = RelFits, NumCols = NumCols, SumNtp1 = SumNtp1, SexRands = SexRands,
                             CurPop = CurPop, Ld = Ld, Lf = Lf, DispSegVals1 = DispSegVals1,
                             DispSegVals2 = DispSegVals2, DispSegIndex = DispSegIndex,
                             FitSegVals1 = FitSegVals1, FitSegVals2 = FitSegVals2,
                             FitSegIndex = FitSegIndex, FitMutStd = FitMutStd,
                             DispMutStd = DispMutStd, FitNumMut = FitNumMut,
                             DispNumMut = DispNumMut, OffspringIndex = OffspringIndex, 
                             FitLocusVec = FitLocusVec, DispLocusVec = DispLocusVec,
                            FitAlleleVec = FitAlleleVec, DispAlleleVec = DispAlleleVec)
          # Now update the SegIndices
          FitSegIndex <- FitSegIndex + SumNtp1*Lf
          DispSegIndex <- DispSegIndex + SumNtp1*Ld
          # And update the OffspringIndex
          OffspringIndex <- OffspringIndex + SumNtp1
     } else{
          PopMat <- matrix(NA, nrow = 0, ncol = NumCols)
     }
     
     # Recalculate CurPop and PopSize
     PopSize <- nrow(PopMat)
     CurPop <- 1:PopSize
     
     # Calculate the time points for the end of the shifting and the total time
     EndShift <- BurnIn + LengthShift
     TotalTime <- BurnIn + LengthShift + BurnOut
     
     # Set up an object to hold summary statistics 
     SumStatCols <- list(gen = 1, beta = 2, x = 3, y = 4, abund = 5, muFit = 6,
                            sigmaFitPhen = 7, sigmaFitGen = 8, muDisp = 9, 
                            sigmaDispPhen = 10, sigmaDispGen = 11)
     SumStatRow <- 1
     SumStats <- matrix(NA, nrow = SumMatSize, ncol = 11)
     CurStatsDim <- SumMatSize
     
     # Calculate the beta values for during the period of climate change
     BetaShift <- ChangeClimate(BetaInit = BetaInit, LengthShift = LengthShift, 
                                v = v, eta = eta)
     
     # Now run through the actual simulation
     for(g in 1:TotalTime){
          if(g <= BurnIn){
               beta <- BetaInit
          } else if( (g > BurnIn) & (g <= EndShift) ){
               beta <- BetaShift[g - BurnIn]
          } else if( (g > EndShift) ){
               beta <- BetaShift[LengthShift]
          }
          
          # Check for extinction before dispersal and reproduction
          if(PopSize != 0){
               traits <- CalcTraits(population = CurPop, PopMat = PopMat, PopSize = PopSize,
                                    PopIndices = PopIndices)
               
               # Check that there are enough entries left in the angles vector and repopulate
               #    it if necessary
               if( (AngleIndex + PopSize) > NumRands ){
                    angles <- runif(n = NumRands, min = 0, max = 2*pi)
                    AngleIndex <- 1
               }
               # Get the new post dispersal locations
               Deltas <- Disperse(traits = traits, width = width, kern = kern, eta = eta,
                                   angles = angles, AngleIndex = AngleIndex, CurPop = CurPop, PopSize = PopSize,
                                   PopIndices = PopIndices, NumRands = NumRands, PopMat = PopMat,
                                   dmax = dmax, rho = rho)
               
               # Update the values in the population matrix
               PopMat[CurPop, PopIndices$x1] <- PopMat[CurPop, PopIndices$x0] + Deltas[,1]
               PopMat[CurPop, PopIndices$y1] <- PopMat[CurPop, PopIndices$y0] + Deltas[,2]
               
               # Finally, check and fix any y values that fall outside of the allowed
               #    width of the landscape and return the updated population matrix
               CurY <- PopMat[CurPop, PopIndices$y1]
               CurY[(CurY > width) | (CurY < 0)] <- CurY[(CurY > width) | (CurY < 0)] %% width
               CurY[CurY == 0] <- width
               PopMat[CurPop, PopIndices$y1] <- CurY
               
               # Update the angle index
               AngleIndex <- AngleIndex + PopSize
               
               # Generate the matrix of occupied patches and get the relative
               #    fitness of each individual
               OccPatches <- unique(PopMat[CurPop,c(PopIndices$x1, PopIndices$y1)])
               if(is.null(nrow(OccPatches))){
                    OccPatches <- matrix(OccPatches, nrow = 1, ncol = 2)
               }
               RelFits <- RelFit(lambda = lambda, beta = beta, traits = traits, 
                                 PopMat = PopMat, individuals = CurPop, omega = omega,
                                 PopIndices = PopIndices, eta = eta)
               
               # Now allow individuals to reproduce
               RealizedNtp1 <- Reproduce(beta = beta, gamma = gamma, tau = tau, omega = omega,
                                      Rmax = Rmax, Kmax = Kmax, traits = traits, PopMat = PopMat, EnvGradType = EnvGradType,
                                      PopIndices = PopIndices, z = z, eta = eta,
                                      SexRands = SexRands, CurPop = CurPop, Ld = Ld, Lf = Lf, 
                                      DispSegVals1 = DispSegVals1, DispSegVals2 = DispSegVals2, 
                                      DispSegIndex = DispSegIndex, FitSegVals1 = FitSegVals1, FitSegVals2 = FitSegVals2, 
                                      FitSegIndex = FitSegIndex, FitMutStd = FitMutStd, DispMutStd = DispMutStd, 
                                      FitNumMut = FitNumMut, DispNumMut = DispNumMut, OffspringIndex = OffspringIndex, 
                                      FitLocusVec = FitLocusVec, DispLocusVec = DispLocusVec, NumCols = NumCols, 
                                      FitPerLocusProb = FitPerLocusProb, DispPerLocusProb = DispPerLocusProb, NumRands = NumRands,
                                      OccPatches = OccPatches, RelFits = RelFits)
               # Now update the population matrix
               SumNtp1 <- sum(RealizedNtp1)
               if(SumNtp1 > 0){
                    # Check that the SegVals vectors contain enough values and resample if not
                    if( (FitSegIndex + SumNtp1*Lf) > NumRands ){
                         FitSegVals1 <- sample(c(0,1), replace = TRUE, size = NumRands)
                         FitSegVals2 <- sample(c(0,1), replace = TRUE, size = NumRands)
                         FitSegIndex <- 1
                    }
                    if( (DispSegIndex + SumNtp1*Ld) > NumRands ){
                         DispSegVals1 <- sample(c(0,1), replace = TRUE, size = NumRands)
                         DispSegVals2 <- sample(c(0,1), replace = TRUE, size = NumRands)
                         DispSegIndex <- 1
                    }
                    # Check the same for the OffspringIndex
                    if( (OffspringIndex + SumNtp1) > NumRands){
                         FitNumMut <- rbinom(n = NumRands, size = Lf, prob = FitPerLocusProb)
                         DispNumMut <- rbinom(n = NumRands, size = Ld, prob = DispPerLocusProb)
                         if( !(is.null(PopIndices$sex)) ){
                              SexRands <- rbinom(n = NumRands, size = 1, prob = z)
                         }
                         OffspringIndex <- 1
                    }
                    PopMat <- MatFill(RealizedNtp1 = RealizedNtp1, PopIndices = PopIndices,
                                      OccPatches = OccPatches, z = z, PopMat = PopMat,
                                      RelFits = RelFits, NumCols = NumCols, SumNtp1 = SumNtp1, SexRands = SexRands,
                                      CurPop = CurPop, Ld = Ld, Lf = Lf, DispSegVals1 = DispSegVals1,
                                      DispSegVals2 = DispSegVals2, DispSegIndex = DispSegIndex,
                                      FitSegVals1 = FitSegVals1, FitSegVals2 = FitSegVals2,
                                      FitSegIndex = FitSegIndex, FitMutStd = FitMutStd,
                                      DispMutStd = DispMutStd, FitNumMut = FitNumMut,
                                      DispNumMut = DispNumMut, OffspringIndex = OffspringIndex, 
                                      FitLocusVec = FitLocusVec, DispLocusVec = DispLocusVec,
                                      FitAlleleVec = FitAlleleVec, DispAlleleVec = DispAlleleVec)
                    # Now update the SegIndices
                    FitSegIndex <- FitSegIndex + SumNtp1*Lf
                    DispSegIndex <- DispSegIndex + SumNtp1*Ld
                    # And update the OffspringIndex
                    OffspringIndex <- OffspringIndex + SumNtp1
               } else{
                    PopMat <- matrix(NA, nrow = 0, ncol = NumCols)
               }
               
               PopSize <- nrow(PopMat)
               CurPop <- 1:PopSize
               
               # Keep track of all summary statistics here
               if(g >= BurnIn){
                    NumPatches <- nrow(OccPatches)
                    if( (SumStatRow + NumPatches) > CurStatsDim){
                         NewMat <- matrix(NA, nrow = SumMatSize, ncol = 11)
                         SumStats <- rbind(SumStats, NewMat)
                         CurStatsDim <- CurStatsDim + SumMatSize
                    }
                    for(i in 1:NumPatches){
                         PatchPop <- which((PopMat[,PopIndices$x0] == OccPatches[i, 1]) &
                                                (PopMat[,PopIndices$y0] == OccPatches[i,2]))
                         PatchPopSize <- length(PatchPop)
                         if(PatchPopSize > 0){
                              SumStats[SumStatRow, SumStatCols$gen] <- g
                              SumStats[SumStatRow, SumStatCols$beta] <- beta
                              SumStats[SumStatRow, SumStatCols$x] <- OccPatches[i,1]
                              SumStats[SumStatRow, SumStatCols$y] <- OccPatches[i,2]
                              SumStats[SumStatRow, SumStatCols$abund] <- PatchPopSize
                              if(PatchPopSize > 1){
                                   PatchFits <- rowSums(PopMat[PatchPop, PopIndices$FitCols])
                                   DispSums <- rowSums(PopMat[PatchPop, PopIndices$DispCols])
                              } else{
                                   PatchFits <- sum(PopMat[PatchPop, PopIndices$FitCols])
                                   DispSums <- sum(PopMat[PatchPop, PopIndices$DispCols])
                              }
                              # Calculate the expected dispersal distances for the patch population
                              ExpDists <- (dmax * eta * exp(rho * DispSums)) / (1 + exp(rho * DispSums))
                              SumStats[SumStatRow, SumStatCols$muFit] <- mean(PatchFits)
                              SumStats[SumStatRow, SumStatCols$sigmaFitGen] <- sd(PopMat[PatchPop, PopIndices$FitCols])
                              SumStats[SumStatRow, SumStatCols$sigmaFitPhen] <- sd(PatchFits)
                              SumStats[SumStatRow, SumStatCols$muDisp] <- mean(ExpDists)
                              SumStats[SumStatRow, SumStatCols$sigmaDispGen] <- sd(PopMat[PatchPop, PopIndices$DispCols])
                              SumStats[SumStatRow, SumStatCols$sigmaDispPhen] <- sd(ExpDists)
                              SumStatRow <- SumStatRow + 1
                         }
                    }
               }
          }
          # Save the population matrix at the time point immediately prior to climate change
          #    to allow for the quantification of initial conditions as well as simulations
          #    with different speeds of climate change.
          if(g == BurnIn){
               PopMatNames <- c("x0", "y0", "x1", "y1", paste("fit", seq(1,(Lf*2)), sep = "_"),
                                paste("disp", seq(1,(Ld*2)), sep = "_"))
               if( !(is.null(PopIndices$sex)) ){
                    PopMatNames <- c(PopMatNames, "sex")
               }
               PreChangePopMat <- PopMat
               colnames(PreChangePopMat) <- PopMatNames
               write.csv(PreChangePopMat, file = paste(ResultsDir, "InitialPopMat.csv", sep = "/"), 
                         row.names = FALSE, quote = FALSE)
          }
     }
     
     # Finally, save the results here
     PopMatNames <- c("x0", "y0", "x1", "y1", paste("fit", seq(1,(Lf*2)), sep = "_"),
                      paste("disp", seq(1,(Ld*2)), sep = "_"))
     if( !(is.null(PopIndices$sex)) ){
          PopMatNames <- c(PopMatNames, "sex")
     }
     colnames(PopMat) <- PopMatNames
     colnames(SumStats) <- c("gen", "beta", "x", "y", "abund", "muFit", "sigmaFitPhen",
                             "sigmaFitGen", "muDisp", "sigmaDispPhen", "sigmaDispGen")
     SumStats <- SumStats[1:(SumStatRow - 1),]
     write.csv(PopMat, file = paste(ResultsDir, "PopMat.csv", sep = "/"), 
               row.names = FALSE, quote = FALSE)
     write.csv(SumStats, file = paste(ResultsDir, "SummaryStats.csv", sep = "/"),
               row.names = FALSE, quote = FALSE)
     return(NULL)
}