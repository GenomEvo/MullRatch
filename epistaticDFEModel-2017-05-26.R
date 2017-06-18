#DATE: 26/05/17


# SUBFUNCTIONS ----


# Reproductive Step ----

# The following function simulates the reproduction step of the population,
# also involving selection. Individuals are sampled based on their fitness, and
# the function returns a list of the indices of parent individuals that release
# offspring into the next generation.

get.parents <- function(w)
{
  N         <- length(w)
  parents   <- sample(N, 
                      N, 
                      replace = TRUE, 
                      prob    = w)
  parents
}


# Distribution Of The Number Of Mutations Per Generation ----

# The following function calculates the mutational state of the population,
# i.e., a vector containing the number of individuals that have 
# 0, 1, 2, ..., kmax mutations.

get.popstate           <- function(pop, kmax)
{
  popstate             <- rep(0, kmax+1)
  
  for(i in 1:length(pop)) 
  {
    popstate[pop[i]+1] <- popstate[pop[i]+1]+1
  }
  
  popstate
}


# ************Subtitle For fitness.vector() Function (Cap Each Word) ----

# ************The following function... (add description of what fitness.vector
# ************function does, unless subtitle is self-explanatory). Include an 
# ************explanation of pop (which is a vector containing count of 
# ************mutations each individual carries) and beta (which parameterises
# ************epistasis).

fitness.vector             <- function(pop, beta = -4) 
{
  # Selection values for the three mutation types (weak = 4.45E-05, 
  # mild = 0.0032 and strong = 0.0218).
  
  selection                <- matrix(c(4.45E-05, 0.0032, 0.0218), 
                                     nrow = 3, 
                                     ncol = 1)
  
  # Find all individuals/rows that carry only one mutation in any of the three
  # mutation types (weak, mild and strong).
  
  oneMutation              <- which(pop == 1, arr.in = TRUE)
  
  # Sum across the three columns (mutation types) of each individual/row and 
  # find which individuals/rows has a total of one mutation.
  
  # If there is only one individual with a total of one mutation, use the "if" 
  # statement. If there is two or more individuals with a total of one 
  # mutation, use the "else" statement.
  
  if (length(oneMutation) == 2)
  {
    SumPerIndividual       <- which(sum(pop[oneMutation[, 1], ]) == 1)
  }
  
  else 
  { 
    SumPerIndividual       <- which(rowSums(pop[oneMutation[, 1], ]) ==1)
  }
  
  # ************Add comment
  
  loc                      <- oneMutation[SumPerIndividual, ]
  
  # The fitness of having one mutation, if it is either weak, mild or strong.
  
  inherentWeak             <- (1-selection[1])
  inherentMild             <- (1-selection[2])
  inherentStrong           <- (1-selection[3])
  
  # ************Combines the individual vectors generated above into a single matrix.
  
  inherent                 <- cbind(inherentWeak, 
                                    inherentMild, 
                                    inherentStrong)
  
  # A vector of the summed selection values for each individual.
  
  totalselection           <- pop %*% selection
  
  # The epistasis model.
  
  fitness                  <- (1-beta*(totalselection))^(1/beta)
  
  
  # ************Add comment
  # ************Find individuals with at least ones only 1 mutation.
  # ************Find individuals with one mutation and replace their fitness with the inherent.fitness.
  
  if (length(loc) != 0) 
  { 
    for (i in 1:length(SumPerIndividual))
    {
      # Only one individual have a total of one mutation.
      
      if (length(loc) == 2) 
      {                  
        fitness[loc[1]]    <- inherent[loc[2]]
      }
      
      # Two or more individuals have a total of one mutation.
      
      else 
      { 
        fitness[loc[i, 1]] <- inherent[loc[i, 2]]
      }
    }
  }
  
  fitness
}


# ************Subtitle For add.newMutations() Function (Cap Each Word) ----

# ************The following function... (add description of what 
# ************add.newMutations() function does, unless subtitle is 
# ************self-explanatory). Include an explanation of N (pop size), v1, 
# ************v2, v3 (the probability that a given new mutation is weak, mild
# ************or strong) and u (human per haploid genome mutation rate modified
# ************from Keightley 2012).

add.newMutations          <- function(N, 
                                      v1 = 0.15, 
                                      v2 = 0.45, 
                                      v3 = 0.4, 
                                      u  = 0.175)
{
  # Mutations are added to the weak class, according to the Poisson 
  # distribution.
  
  mWeak                   <- rpois(N, u*v1)
  
  # Mutations are added to the mild class, according to the Poisson 
  # distribution.
  
  mMild                   <- rpois(N, u*v2)
  
  # Mutations are added to the strong class, according to the Poisson 
  # distribution.
  
  mStrong                 <- rpois(N, u*v3)
  
  # ************Combines the individual vectors generated above into a single 
  # ************matrix with specified row names.  This matrix can then be 
  # ************returned.
  
  newmutmatrix            <- cbind(mWeak, 
                                   mMild, 
                                   mStrong)
  
  row.names(newmutmatrix) <- (1:N)
  
  newmutmatrix
}



# MAIN SIMULATION CODE ----


# Muller's Ratchet Simulation Function ----

# Input parameters include the population size (N) and a vector of time
# points (timepoints) can be given at which the state of the population is 
# returned. The maximum of that vector (tmax) determines the total number of 
# generations that the simulation runs. Finally, the parameter kmax indicates
# the maximum number of deteterious mutations that is kept track of. 
# Individuals can still have more mutations, but those are treated as if the 
# individual had only kmax mutations.

ratchet.simulator            <- function(N          = 400,
                                         timepoints = seq(0, 1000, by = 10), 
                                         kmax       = 1000) 
{
  # The set maximum number of generations.
  
  tmax                       <- max(timepoints)
  
  # This is an index counting the time points at which the population state is
  # returned.
  
  j                          <- 1				 
  
  # The initial population is completely free of mutations.
  
  # Each individual in the population is characterized by the number of 
  # deleterious mutations (weak, mild or strong) that the individual carries.
  
  pop                        <- matrix(0, 
                                       nrow = N, 
                                       ncol = 3)
  
  colnames(pop)              <- c("nWeak", 
                                  "nMild", 
                                  "nStrong")
 
  # Matrix dimensions for each class (weak, mild or strong).
  
  d                          <- dim(matrix(0, 
                                    nrow = length(timepoints),
                                    ncol = kmax+1))
  
  nWeakpopstates             <- matrix(0, 
                                       d[1], 
                                       d[2])
  
  nMildpopstates             <- matrix(0, 
                                       d[1], 
                                       d[2])
  
  nStrongpopstates           <- matrix(0, 
                                       d[1], 
                                       d[2])
  
  # The popstates is the distribution of the number of deleterious mutations in
  # the population for each class, ranging from 0 to kmax.
  
  popstates                  <- list(nWeakpopstates, 
                                     nMildpopstates,
                                     nStrongpopstates)
  
  # The fitness of each individual per time point (or generation).
  
  w                          <- matrix(0,
                                       nrow = tmax+1,
                                       ncol = N) 
  
  rownames(w)                <- seq(0, tmax)
  colnames(w)                <- 1:N
  
  
  # **************Add comment
  
  for (i in 1:length(popstates)) 
  {
    rownames(popstates[[i]]) <- timepoints
    colnames(popstates[[i]]) <- 0:kmax
  }

  
  # Return popstate at time point 0.
  
  if (0 %in% timepoints)
  {
    for (i in 1:length(popstates)) 
    {
      popstates[[i]][j,]     <- get.popstate(pop[, i], kmax)
    }
    
    j                        <- j+1		
  }
  
  
  # ****************Add comment
    
  for(t in 1:tmax)
  {
    # Find the fitness of each individual and puts these values into a vector.
    
    w[t,]                    <- fitness.vector(pop)
    
    # A vector of indices giving all parental individuals.
    
    parents                  <- get.parents(w[t, ])
    
    # Offspring population is identical to their parents.
    
    pop                      <- pop[parents, ]
    
    # Mutations are added according to the add.newMutations() funciton.
    
    pop                      <- pop + add.newMutations(N)
    
    # Individuals are not allowed to have > kmax mutations.
    
    pop[pop>kmax]            <- kmax
    
    
    # Output the population state at the present time point.
    
    if (t %in% timepoints) 
    {
      for (i in 1:length(popstates))
      {
        popstates[[i]][j, ]   <- get.popstate(pop[, i], kmax)
      }
      
      j                       <- j+1		
    }
   
    
    # If timepoint is tmax, then find final fitness of population.
    
    if (t %in% tmax) 
    {
      w[tmax+1, ]             <- fitness.vector(pop)
    } 
  }
  
  
  for(i in 1:length(popstates))
  {
    popstates[[i]]            <- popstates[[i]][, colSums(popstates[[i]] == 0) != nrow(popstates[[i]])]
  }
  

  # Create a list of popstates, pop and fitnesses (w) which can then be returned.
  
  output                      <- list(popstates, pop, w)
  
  return(output)
}



# Three Runs Of The Simulation ----

for (i in 1:3) 
{
  # Run the ratchet.simulator() function for specified number of time points 
  # (or generations).
  
  time.points <- seq(0, 5000)
  output      <- ratchet.simulator(timepoints = time.points)
  
  # Naming and saving data output.
  
  s1 = "NCS_7_rep"
  s2 = i
  
  the.name    <- paste(s1, s2, sep = "")
  saveRDS(output, the.name)
}

