# This script contains all of the functions necessary to run a simulation of
#    population evolution and movement in response to climate change. The 
#    functions are organized into three different categories: biological,
#    bookkeeping, and environmental. Each function is documented with a brief
#    description of the purpose of the function and a list and definition of 
#    the necessary function arguments before the function definition itself.
#    Each broad category of functions is delineated by the dashed section
#    headers and each individual function begins with a single line of 6 pound
#    signs (######)


# ------------------------------------------------------------------------------
# ---------------------------- Biological Functions ----------------------------
# ------------------------------------------------------------------------------

###### CalcTraits
# This function will calculate the fitness and dispersal trait values for a set 
#    of individuals based on their alleles
### INPUTS
# population:  A vector of the population (i.e. row numbers) for which to 
#                   calculate a trait values.
# PopMat:      The population matrix to use
# PopSize:     The size of the current population
# PopIndices: Column indices to use for the population matrix
### OUTPUTS
# The function will create a matrix with three columns. The first column will 
#    correspond to the ID (row number) of the individual, the second column
#    will be the fitness trait value, and the third column will be the dispersal
#    trait. If the traits are calculated for the entire population, as will 
#    normally be the case, then the first column will simply correspond to the 
#    row number, but including this column allows for the calculation of 
#    traits for a subset of the current population.
CalcTraits <- function(population, PopMat, PopSize, PopIndices){
     # First create an empty matrix for the trait values and put
     #    the population IDs in the first column
     traits <- matrix(NA, nrow = PopSize, ncol = 3)
     traits[,1] <- population
     
     # Now calculate the sum of each individual's quantitative loci for each 
     #    trait and return those sums in the traits matrix
     if(PopSize > 1){
          traits[,2] <- rowSums(PopMat[population, PopIndices$FitCols])
          traits[,3] <- rowSums(PopMat[population, PopIndices$DispCols])
     } else{
          traits[,2] <- sum(PopMat[population, PopIndices$FitCols])
          traits[,3] <- sum(PopMat[population, PopIndices$DispCols])
     }
     return(traits)
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

###### Disperse
# This function will use the previously calculated dispersal traits to calculate
#    new, post dispersal locations for individuals according to any of three
#    potential dispersal kernels: normal, exponential, or student's t dispersal
#    kernel (Shaw 2014).
### INPUTS:
# traits:      The three column matrix created by the CalcTrait function
# width:       The number of discrete cells defining the width of the word
#                   being simulated
# kern:        The dispersal kernel to be used. Can be "norm", "exp", or 
#                   "stud_t"
# eta:  A constant value determining the width of individual patches on
#                   the Cartesian coordinate system defining both
#                   environmental quality and dispersal kernels.
# angles:      a large vector of random numbers generated before calling the function
# AngleIndex:  the index of the where to look within the vector
# PopIndices:  A list containing all the index values for different columns in
#                   the matrix
# NumRands: An integer denoting how many random numbers to generate at a time for
#              the random number vectors.
# rho: parameter determining the slope of the change from 0 to dmax in expected
#         dispersal distance
# dmax:   the maximum possible expected dispersal distance (in terms of number
#         of patches). Actual dispersal distance can obviously exceed this as it
#         is drawn from a distribution.
### OUTPUTS:
# The function will update objects in the Population matrix in the parent environment
#    and then return the updated value of AngleIndex
Disperse <- function(traits, width, kern, eta, angles, AngleIndex, rho,
                     CurPop, PopSize, PopIndices, NumRands, PopMat, dmax){
     # First calculate the expected dispersal distance for each individual
     d <- (dmax * eta * exp(rho * traits[CurPop,3])) / (1 + exp(rho * traits[CurPop,3]))
     
     # Next generate the dispersal distances according to the type of kernel
     if(kern == "norm"){
          sigma <- (d^2 * pi) / 2
          dists <- abs(rnorm(n = PopSize, mean = 0, sd = sigma))
     } else if(kern == "exp"){
          dists <- rexp(n = PopSize, rate = 1 / d)
     } else if(kern == "stud_t"){
          dists <- rStudKern(n = PopSize, d = d)
     }
     
     # Now calculate the new x and y coordinates of each individual as if their
     #    current location is the origin
     NewX <- dists * cos(angles[AngleIndex:(AngleIndex + PopSize - 1)])
     NewY <- dists * sin(angles[AngleIndex:(AngleIndex + PopSize - 1)])
     
     # Now use the eta parameter to determine the number of x and y
     #    patches changed, accounting for the movement occuring from the patch
     #    center by adding or subtracting eta/2 as necessary
     DeltaX <- ifelse(NewX < 0, ceiling((NewX - eta/2) / eta), floor((NewX + eta/2) / eta))
     DeltaY <- ifelse(NewY < 0, ceiling((NewY - eta/2) / eta), floor((NewY + eta/2) / eta))
     
     return(cbind(DeltaX,DeltaY))
}

###### RelFit
# This function will calculate the relative fitness values for individuals 
#    depending on their position within the range and their respective phenotype
#    values.
### INPUTS
# lambda:    A value determining the strength of local selection across the
#                   range. Higher values will correspond to a steeper gradient
#                   in optimum phenotype values from one end of the range to the
#                   other.
# beta:        The beta value used to define the center of the range
# eta:         The patch size in the model
# traits:      The trait matrix calculated by the CalcTraits function
# PopMat:      The population matrix containing the focal individuals
# individuals: A vector of the individual IDs (row numbers) of the focal 
#                   individuals.
# omega:       This parameter is the inverse of the strength of stabilizing 
#                   selection.
# PopIndices:  A list containing all the index values for different columns in
#                   the matrix
### OUTPUTS
# This function will generate a vector of relative fitness values for each
#    individual in the same order they are presented in the phenotype vector.
RelFit <- function(lambda, beta, traits, PopMat, individuals, omega,
                   PopIndices, eta){
     # Calculate the Zopt values for each individual depending on where they
     #    are in space.
     Zopt <- lambda * (PopMat[individuals, PopIndices$x1] * eta - beta)
     
     # Now, use the equation for stabilizing selection to calculate each
     #    individual's relative fitness and return that as a vector
     RelFits <- exp(-1 * (traits[individuals, 2] - Zopt)^2 / (2*omega^2))
     return(RelFits)
}

###### MatFill
# This function will fill in the appropriate values for a new population matrix
#    in the Reproduce function.
### INPUTS
# RealizedNtp1:     A vector of the stochastically determined population size for
#                        each occupied patch in the next generation
# PopIndices:       A list containing all the index values for different columns in
#                        the matrix
# OccPatches:       A two column matrix with the x and y coordinates of all
#                        occupied patches in the current generation.
# z:         The sex ratio of the population (set to 0.5 by default)
# RelFits:          A vector of the relative fitness values of individuals in the
#                        population.
# FitCols:       A vector of the columns corresponding to fitness traits
# DispCols:      A vector of the columns corresponding to dispersal
# PopMat:           The population matrix for the current generation
# NumCols:          The overall number of columns in population matrices
# SexRands:         A large vector of random numbers for assigning sex
# CurPop:           A vector of the current individuals in the population
# SumNtp1:          The total number of offspring in the next generation
# Ld:            Number of dispersal loci
# Lf:             Number of fitness loci
# SegVals(1,2 and index for both fit and disp; see Inheritence() documentation)
# MutStd for fit and disp
# NumMut(for both fit and disp)
# OffspringIndex: an all purpose index for the random number vectors for which
#    there is a single entry per offpsring (NumMuts and SexRands)
# LocusVec (for both Fit and Disp; see Inheritence() documentation for details)
# AlleleVec (For both fit and disp; see Inheritence() documentation)
### OUTPUTS
# This function will return a new matrix composed of all the offspring and their
#    location, loci, etc. for the next generation.
MatFill <- function(RealizedNtp1, PopIndices, OccPatches, z, RelFits,
                    PopMat, NumCols, SexRands = NULL, CurPop, SumNtp1, Ld, Lf, 
                    DispSegVals1, DispSegVals2, DispSegIndex, FitSegVals1, FitSegVals2, 
                    FitSegIndex, FitMutStd, DispMutStd, FitNumMut, DispNumMut, OffspringIndex,
                    FitLocusVec, DispLocusVec, FitAlleleVec, DispAlleleVec){
     # Create the new matrix for the next population
     NewPopMat <- matrix(NA, nrow = SumNtp1, ncol = NumCols)
     
     # Filter for only the patches that produced offspring
     NewOccPatches <- which(RealizedNtp1 > 0)
     Ntp1 <- RealizedNtp1[NewOccPatches]
     
     # Now step through and fill in the matrix, starting by setting initital row
     #    values
     StartRow <- 1
     EndRow <- Ntp1[1]
     for(i in 1:length(NewOccPatches)){
          CurVec <- StartRow:EndRow
          # Now update the starting and ending row for the next iteration of the
          #    loop.
          CurEnd <- EndRow
          CurStart <- StartRow
          StartRow <- EndRow + 1
          EndRow <- EndRow + Ntp1[i+1]
          
          # Fill in the location details
          NewPopMat[CurVec, PopIndices$x0] <- OccPatches[NewOccPatches[i], 1]
          NewPopMat[CurVec, PopIndices$y0] <- OccPatches[NewOccPatches[i], 2]
          
          
          # Identify the potential parents in the current patch
          locals <- which( (PopMat[,PopIndices$x1] == OccPatches[NewOccPatches[i], 1]) &
                                PopMat[,PopIndices$y1] == OccPatches[NewOccPatches[i], 2])
          NumLocals <- length(locals)
          # Depending on life history being modeled (monoecious vs. dioecious),
          #    add in sex information if relevant and determine parents for each
          #    offspring
          if( !(is.null(PopIndices$sex)) ){
               NewPopMat[CurVec, PopIndices$sex] <- SexRands[(OffspringIndex + CurStart):(OffspringIndex + CurEnd)]
               
               # Identify the females and males present in the current patch
               females <- which(PopMat[locals, PopIndices$sex] == 1)
               males <- which(PopMat[locals, PopIndices$sex] == 0)
               
               # Now extract the relative fitness values for the males and
               #    females
               FemaleFits <- RelFits[females]
               MaleFits <- RelFits[males]
               
               # Sample from the pool of females and males according to their
               #    relative fitnesses to select parents for each offspring, but
               #    only use the sample command if there are more than one male
               #    or female available to avoid undesired behavior from the
               #    sample function
               if(length(females) == 1){
                    parent1 <- rep(females, Ntp1[i])
               } else{
                    parent1 <- sample(females, size = Ntp1[i], replace = TRUE,
                                      prob = FemaleFits)
               }
               if(length(males) == 1){
                    parent2 <- rep(males, Ntp1[i])
               } else{
                    parent2 <- sample(males, size = Ntp1[i], replace = TRUE,
                                      prob = MaleFits)
               }
               parents <- cbind(parent1, parent2)
          } else{
               # Get the relative fitness values for the potential parents
               ParentFits <- RelFits[locals]
               
               # Next select both parents, again avoiding unwanted behavior
               #    from the sample function
               if(NumLocals == 1){
                    parent1 <- rep(locals, Ntp1[i])
                    parent2 <- rep(locals, Ntp1[i])
               } else if(NumLocals == 2){
                    parent1 <- rep(locals[1], Ntp1[i])
                    parent2 <- rep(locals[2], Ntp1[i])
               } else{
                    parent1 <- sample(locals, size = Ntp1[i], replace = TRUE, 
                                      prob = ParentFits)
                    parent2 <- sample(locals, size = Ntp1[i], replace = TRUE, 
                                      prob = ParentFits)
                    SameParent <- which(parent1 == parent2)
                    for(j in SameParent){
                         ParentPool <- setdiff(locals, parent1[j])
                         parent2[j] <- sample(ParentPool, size = 1,
                                              prob = RelFits[ParentPool])
                    }
               }
               parents <- cbind(parent1, parent2)
          }
          # Fill in the loci according to the inheritence function
          PatchNtp1 <- length(CurVec)
          NewPopMat[CurVec, PopIndices$FitCols] <- Inheritence(Cols = PopIndices$FitCols, SumNtp1 = PatchNtp1,
                                                       parents = parents, PopMat = PopMat,
                                                       NumLoci = Lf, SegVals1 = FitSegVals1,
                                                       SegVals2 = FitSegVals2, SegIndex = FitSegIndex,
                                                       MutStd = FitMutStd, NumMutVec = FitNumMut,
                                                       MutIndex = OffspringIndex, LocusVec = FitLocusVec,
                                                       AlleleVec = FitAlleleVec)
          NewPopMat[CurVec, PopIndices$DispCols] <- Inheritence(Cols = PopIndices$DispCols, SumNtp1 = PatchNtp1,
                                                        parents = parents, PopMat = PopMat,
                                                        NumLoci = Ld, SegVals1 = DispSegVals1,
                                                        SegVals2 = DispSegVals2, SegIndex = DispSegIndex,
                                                        MutStd = DispMutStd, NumMutVec = DispNumMut,
                                                        MutIndex = OffspringIndex, LocusVec = DispLocusVec,
                                                        AlleleVec = DispAlleleVec)
     }
     return(NewPopMat)
}

###### Inheritence
# This function performs inheritence of alleles from both parents, assuming
#    independent segregation and no crossover events. It incorporates mutations
#    in alleles via the U and Vm parameters (see below).
# INPUTS
# parents:     A matrix with two columns and a number of rows equal to the 
#                   number of offpsring being created. Each row contains the
#                   IDs (row indices in PopMat) of the parent(s) of the 
#                   offspring corresponding to that row.
# Cols:        A vector containing the focal columns containing the loci for the
#                   trait being inherited.
# PopMat:      The matrix containing all population information for the parental
#                   generation.
# NumLoci:     The number of loci for the given trait
# SegVals1:    A large vector with pre-randomized values for parent1
# SegVals2:    The same for parent2
# SegIndex:    The index for both vectors
# MutStd:      The standard deviation of mutatiosn for that trait
# NumMutVec:   A pre-randomized vector with values of number of mutations for a 
#                   given trait
# MutIndex:    The index for that vector
# SumNtp1:     The total number of offspring in the next generation
# LocusVec:    A vector of 1:the number of loci for a given trait. Made once in
#                   the main function and then re-used over and over here for efficiency.
# AlleleVec:   Similar to LocusVec, but with an entry for each allele (i.e. it is twice
#                   the size of LocusVec)
# OUTPUTS
# This function will return a matrix consisting of a number of rows equal to
#    the number of offspring produced and a number of columns equal to the 
#    number of loci defining the trait under consideration.
Inheritence <- function(Cols, parents, PopMat, SumNtp1, NumLoci, SegVals1, AlleleVec,
                        MutStd, NumMutVec, MutIndex, LocusVec, SegVals2, SegIndex){
     SegregatedLoci <- matrix(NA, nrow = SumNtp1, ncol = 2*NumLoci)
     for(i in 1:SumNtp1){
          ParentLoci <- PopMat[parents[i,], Cols]
          Parent1Alleles <- LocusVec + NumLoci * SegVals1[(SegIndex + (i - 1) * NumLoci):(SegIndex + i * NumLoci - 1)]
          Parent2Alleles <- LocusVec + NumLoci * SegVals2[(SegIndex + (i - 1) * NumLoci):(SegIndex + i * NumLoci - 1)]
          SegregatedLoci[i,] <- c(ParentLoci[1,Parent1Alleles], ParentLoci[2,Parent2Alleles])
     }
     
     # Extract the number of mutations for offspring from the NumMutVec
     NumMut <- NumMutVec[MutIndex:(MutIndex+SumNtp1-1)]
     
     # Now step through each offspring in which a mutation takes place and alter
     #    allele values appropriately
     MutOffspring <- which(NumMut != 0)
     for(i in MutOffspring){
          MutLocus <- sample(AlleleVec, size = NumMut[i], replace = TRUE)
          SegregatedLoci[i,MutLocus] <- rnorm(mean = SegregatedLoci[i,MutLocus], 
                                              sd = MutStd, n = NumMut[i])
     }
     
     # return the segregated and mutated loci
     return(SegregatedLoci)
}

###### Reproduce
# This function will use the previously calculated fitness traits and perform
#    Ricker population growth in all occupied patches with offspring assigned
#    to parental pairs according to individual fitness values.
### INPUTS:
# beta, gamma, tau:    Parameters corresponding to the range capacity fxn
# omega:            The inverse of the strength of stabilizing selection.             
# Rmax, Kmax:           The growth and carrying capacity parameters for logistic
#                        growth. These values are the maximum attainable values
#                        which are then used to compute realized values
#                        throughout the range.
# U:                A vector of 2 entries with the diploid mutation rate for
#                        each trait.
# Vm:               A vector of 2 entries with the mutational variance for each
#                        trait.
# lambda:         A value determining the gradient in local phenotypic optima
#                        across the range.
# traits:           A two column matrix of trait values created by CalcTraits
# PopMat:           The matrix with all information on the current generation
# EnvGradType:      The type of environmental gradient being measured (either "K"
#                        or "R")
# ColumnNames:      A vector of the column names relevant for the simulation.
# z:         The sex ratio of the population if individuals are dioecious.
#                        By default, this is set to 0.5
# FitCols:       A vector of the column indices corresponding to fitness loci
# DispCols:      A vector of the column indices corresponding to dispersal loci
# eta:       A factor determining the mapping of individual patches onto
#                        the landscape quality
# PopIndices:       A list with the index of different columns in the population
#                        matrix
# SexRands(and index):   See the matfill function for details
# CurPop:           A vector of the current population
# Ld & Lf:     Number of loci for each trait
# SegVals and indices for each trait:   See inheritence function
# MutStd for both traits: See inheritence function
# NumMut and indices for each trait:    See inheritence function
# LocusVec for each trait:              See inheritence function
# NumCols:          The number of columns in the population matrix
# PerLocusProb:     The previously calculated per locus probability of a mutation
# NumRands:         The number of random numbers to generate at a time
# OccPatches:       A matrix with the x and y coordinates of all occupied patches
# RelFits:          A vector of relative fitness values
### OUTPUTS:
# The function will return a vector of population sizes for the occupied patches
Reproduce <- function(beta, gamma, tau, omega, Rmax, Kmax, 
                      traits, PopMat, EnvGradType, ColumnNames, z = 0.5,
                      eta, PopIndices, SexRands = NULL, CurPop, Ld, Lf, 
                      DispSegVals1, DispSegVals2, DispSegIndex, FitSegVals1, FitSegVals2, 
                      FitSegIndex, FitMutStd, DispMutStd, FitNumMut, DispNumMut, OffspringIndex, 
                      FitLocusVec, DispLocusVec, NumCols, FitPerLocusProb, DispPerLocusProb, NumRands,
                      OccPatches, RelFits){
     # Check that Rmax is non-negative before anything else
     if(Rmax < 1){
          write("Rmax values less than 1 will result in negative alpha values, causing unrestricted population growth",
                stderr())
          return(NULL)
     } else if(Rmax == 1){
          write("Setting Rmax to 1 results in alpha values of either 0 or NaN, resulting in unrestricted growth or errors respectively",
                stderr())
     }
     
     # Determine the number of occupied patches and get the patch qualities for
     #    each, then create an object to hold expected population sizes
     NumPatches <- nrow(OccPatches)
     PatchEnvQual <- GetEnvQual(beta = beta, gamma = gamma, 
                                     tau = tau, patches = OccPatches[,1], 
                                     eta = eta)
     Ntp1 <- rep(NA, NumPatches)
     
     # Loop through each occupied patch and calculate the expected population
     #    size in the next generation
     for(i in 1:NumPatches){
          # Extract the quality of the current patch and the population size
          PatchQual <- PatchEnvQual[i]
          CurPatchPop <- which( (PopMat[,PopIndices$x1] == OccPatches[i,1]) & 
                                 (PopMat[,PopIndices$y1] == OccPatches[i,2]) )
          PatchPopSize <- length(CurPatchPop)
          
          # Calculate the relevant demographic parameter and return a helpful
          #    error message if needed
          if(EnvGradType == "K"){
               PatchK <- PatchQual * Kmax
               PatchR <- Rmax
          } else if(EnvGradType == "R"){
               PatchR <- PatchQual * Rmax
               PatchK <- Kmax
          } else{
               write("Invalid type of environmental gradient", stderr())
               return(NULL)
          }
          
          # Calculate the patch alpha value based on the R and K values determined
          #    above 
          PatchAlpha <- log(PatchR) / PatchK

          # Now determine the expected population growth according to the stochastic
          #    Ricker models derived by Melbourne & Hastings (2008) and the sex
          #    structure of the population
          if(!is.null(PopIndices$sex)){
               CurFemales <- which( (PopMat[,PopIndices$x1] == OccPatches[i,1]) & 
                                         (PopMat[,PopIndices$y1] == OccPatches[i,2]) &
                                         (PopMat[,PopIndices$sex] == 1) )
               PatchFit <- mean(RelFits[CurFemales])
               NumFemales <- length(CurFemales)
               if((NumFemales == PatchPopSize) | (NumFemales == 0) | (is.infinite(PatchAlpha))){
                    Ntp1[i] <- 0
               } else{
                    Ntp1[i] <- PatchFit * NumFemales * (PatchR / z) * exp(-1 * PatchAlpha * PatchPopSize)
               }
          } else{
               PatchFit <- mean(RelFits[CurPatchPop])
               if(is.infinite(PatchAlpha)){
                    Ntp1[i] <- 0
               } else{
                    Ntp1[i] <- PatchFit * PatchPopSize * PatchR * exp(-1 * PatchAlpha * PatchPopSize)
               }
          }
     }
     
     # Use the expected population sizes to generate the realized population 
     #    sizes for each patch
     RealizedNtp1 <- rpois(n = NumPatches, lambda = Ntp1)

     return(RealizedNtp1)
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