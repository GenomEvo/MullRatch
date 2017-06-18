#DATE: 20/05/17



# SUBFUNCTIONS ----


# Reproductive Step ----

# The following function simulates the reproduction step of the population,
# also involving selection. Individuals are sampled based on their fitness, and
# the function returns a list of the indices of parent individuals that release
# offspring into the next generation.

get.parents <- function(w)
{
  N       <- length(w)
  parents <- sample(N, 
                    N, 
                    replace = TRUE, 
                    prob    = w)
  parents
}


# Distribution Of The Number Of Mutations Per Generation ----

# The following function calculates the mutational state of the population,
# i.e., a vector containing the number of individuals that have 
# 0, 1, 2, ..., kmax mutations.

get.popstate <- function(pop, kmax)
{
  popstate <- rep(0, kmax+1)
  
  for(i in 1:length(pop))
  { 
    popstate[pop[i]+1] <- popstate[pop[i]+1]+1
  }
  popstate
}


# ************Subtitle For fitness.vector() Function (Cap Each Word) ----

# ************The following function... (add description of what fitness.vector
# ************function does, unless subtitle is self-explanatory).

fitness.vector <- function(pop, beta = -2)         # Beta parameterises 
                                                   # epistasis.
{
  # Selection values for the three groups of mutations (Weak, mild and strong).
  
  selection        <- matrix(c(0.001, 0.1, 0.6),
                             nrow=3,
                             ncol=1)
  
  # ************What is "oneMutation"? Maybe add a comment to be consistant...
  
  oneMutation      <- which(pop == 1, arr.in = TRUE)
  
  # ************What is "SumPerIndividual"? Maybe add a comment to be consistant...
  
  SumPerIndividual <- which(rowSums(pop[oneMutation[, 1], ]) ==1)
  
  # ************What is "loc"? Maybe add a comment to be consistant...
  
  loc              <- oneMutation[SumPerIndividual, ]
  
  # ************What are these? Maybe add a comment to be consistant...
  
  inherentWeak     <- (1-selection[1])
  inherentMild     <- (1-selection[2])
  inherentStrong   <- (1-selection[3])
  
  # ************Combines the three elements in the block above into a 
  # ************matrix/vector?
  
  inherent         <- cbind(inherentWeak, 
                            inherentMild, 
                            inherentStrong)
  
  # Vector of the summed seletion values for each individual.
  
  totalselection   <- pop %*% selection
  
  # Vector of fitness values.
  
  fitness          <- (1-beta*(totalselection))^(1/beta)
  
  if (length(loc) != 0) 
  {
    for (i in 1:length(SumPerIndividual))
    {
      if (length(loc) == 2) 
      {
        fitness[loc[1]]   <- inherent[loc[2]]
      }
      else 
      {
        fitness[loc[i,1]] <- inherent[loc[i, 2]]
      }
    }
  }
  
  fitness
}


# ************Subtitle For add.newMutations() Function (Cap Each Word) ----

# ************The following function... (add description of what 
# ************add.newMutations() function does, unless subtitle is 
# ************self-explanatory). Include: The probability that a given new
# ************mutation is weak (v1), mild (v2) or strong (v3). The per genome 
# ************mutation rate (u).

add.newMutations <- function(N, 
                             v1 = 0.49,
                             v2 = 0.46,
                             v3 = 0.05,
                             u  = 0.02)
{which(x %in% c(2,4))
  
  # Mutations are added to the weak, mild and strong classes of mutations, 
  # respectively, according to Poisson (or gamma ?) distribution.
  
  mWeak                   <- rpois(N, u*v1)
  mMild                   <- rpois(N, u*v2)
  mStrong                 <- rpois(N, u*v3)
  
  # ************Combines the three elements in the block above into a 
  # ************matrix/vector? with specified row names...
  
  newmutmatrix            <- cbind(mWeak,    
                                   mMild, 
                                   mStrong)
  row.names(newmutmatrix) <- (1:N)
  
  newmutmatrix
}



# MAIN SIMULATION CODE ----

# Muller's Ratchet Simulation Function ----

# ************Must edit the description below so that it explains the new 
# ************parameters of the nonConstant ratchet.simulator() function
# ************like timepoints and kmax.
# Input parameters are the population size N, the selection coefficient s and
# the per genome mutation rate u. Moreover, a vector of time points can be
# given at which the state of the population is returned. The maximum of that
# vector determines the total number of generations that the simulation runs.
# Finally, the parameter kmax indicates the maximum number of deteterious
# mutations that is kept track of. (Individuals can still have more mutations,
# but those are treated as if the individual had only kmax mutations.)

ratchet.simulator <- function(N = 400,
                              timepoints = seq(0, 1000, by = 10), 
                              kmax=300) 
{
  # ************What is "tmax"? Maybe add a comment to be consistant...
  
  tmax             <- max(timepoints)
  
  # This is an index counting the time points at which the population state is
  # returned.
  
  j                <- 1				 
  
  # The initial population is completely free of any mutations.
  # Each individual in the population is characterized by the number of 
  # deleterious mutations (weak, mild or strong) that the individual carries.
  # These values are formatted into a matrix with specified column names.
  
  pop              <- matrix(0,
                             nrow = N,
                             ncol = 3)
  
  colnames(pop)    <- c("nWeak","nMild","nStrong")
  
  # ************What is "d"? Maybe add a comment to be consistant...
  
  d                <- dim(matrix(0,
                          nrow = length(timepoints),
                          ncol = kmax+1))
  
  # Matrix dimensions for each class (weak, mild or strong).
  
  nWeakpopstates   <- matrix(0, d[1], d[2])
  nMildpopstates   <- matrix(0, d[1], d[2])
  nStrongpopstates <- matrix(0, d[1], d[2])
  
  # The popstates is the distribution of the number of deleterious mutations in
  # the population for each class, ranging from 0 to kmax.
  
  popstates        <- list(nWeakpopstates, 
                           nMildpopstates,
                           nStrongpopstates)
  
  # ************What is "w"? Maybe add a comment to be consistant...
  # Fitness for each individual per generation/time point.
  
  w                <- matrix(0,
                             nrow = tmax,
                             ncol = N)
  
  rownames(w)      <- seq(1, tmax)
  colnames(w)      <- 1:N
  
  for (i in 1:length(popstates)) 
  {
    rownames(popstates[[i]]) <- timepoints
    colnames(popstates[[i]]) <- 0:kmax
  }

  if (0 %in% timepoints)                 # Return popstate at time point 0. 
  {
    for (i in 1:length(popstates)) 
    {
      popstates[[i]][j,] <- get.popstate(pop[, i], kmax)
    }
    j <- j+1		
  }
  
    
  for(t in 1:tmax)
  {
    w[t,]         <- fitness.vector(pop)        # Vector of fitness values for
                                                # each individual.
    
    parents       <- get.parents(w[t, ])        # Vector of indices giving all 
                                                # parental individuals.
    
    pop           <- pop[parents, ]				      # Offspring population, identical
                                                # to their parents.
    
    pop           <- pop + add.newMutations(N)
    
    pop[pop>kmax] <- kmax				                # Individuals are not allowed to 
                                                # have >kmax mutations.
    
    if (t %in% timepoints)			                # Output the population state at 
                                                # the present timepoint.
    {
      for (i in 1:length(popstates))
      {
        popstates[[i]][j,] <- get.popstate(pop[, i], kmax)
      }
      j <- j+1		
    }
  
  }
  
  for(i in 1:length(popstates)) 
  {
    popstates[[i]] <- popstates[[i]][, colSums(popstates[[i]]==0) !=nrow(popstates[[i]])]
  }
  
  output <- list(popstates, pop, w)
  return(output)
}


# Number of Ratchet Clicks ----

# The following function takes the state of the population as input and from 
# this calculates the number of ratchet clicks that have occured for each 
# class. 

ratchet.clicks <- function(pop.states)
{
  threeFitnessClicks           <- matrix(0, 
                                         nrow = length(pop.states[[1]][,1]), 
                                         ncol = length(pop.states))
  
  for (r in 1:length(pop.states))
  {
    clicks                 <- rep(0, length(pop.states[[r]][,1]))
    for(i in 1:length(clicks))
    {
      j         <- 0
      while ((pop.states[[r]][i,j+1]==0)&&(j<length(pop.states[[r]][1,]))) 
      {
        j <- j+1
      }
      clicks[i] <- j
    }
    threeFitnessClicks[, r] <- clicks
  }
  rownames(threeFitnessClicks) <- row.names(pop.states[[1]])
  threeFitnessClicks
}



# PLOTTING THE OUTPUT FROM MULLER'S RATCHET SIMULATION ----


# Run Of The Simulation ----

time.points <- seq(0, 5000)
output      <- ratchet.simulator(timepoints = time.points)
pop.states  <- output[[1]]
fitness     <- output[[3]]
type        <- output[[2]]
clicks      <- ratchet.clicks(pop.states)


#Plotting Outputs ----

# ************What is the below doing? Add comment explaining.

#gen<-5001
#barplot(pop.states[[1]][gen,], xlab="Number of mutations", ylab="Frequency", main=paste("Distribution of mutations at generation",gen))
#box()

# Plotting the progress of Muller's ratchet as the number of ratchet clicks
# over time (generations).

matplot(x    = rownames(clicks),
        y    = clicks, 
        type = "l", 
        xlab = "Generation", 
        ylab = "Ratchet clicks")

# Plotting the average population fitness level over time (generations).

average.loss <- 1-rowMeans(fitness)
plot(x    = seq(1, 5000),
     y    = average.loss, 
     xlab = "Generation", 
     ylab = "Average Loss")

# Creating a table stating the sum of each mutation class.

mType      <- colSums(type)
mTypeTable <- as.table(mType)
mTypeTable
