#DATE: 30/05/17


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



# MAIN SIMULATION CODE ----


# Muller's Ratchet Simulation Function ----

# Input parameters are the population size (N), the selection coefficient (s) 
# and the human per haploid genome mutation rate (u) modified from Keightley 
# 2012. 

# Moreover, a vector of time points (timepoints) can be given at which the 
# state of the population is returned. The maximum of that vector (tmax) 
# determines the total number of generations that the simulation runs.

# Finally, the parameter kmax indicates the maximum number of deteterious
# mutations that is kept track of. Individuals can still have more mutations,
# but those are treated as if the individual had only kmax mutations.

ratchet.simulator     <- function(N          = 400, 
                                  s          = 0.01, 
                                  u          = 0.175, 
                                  timepoints = seq(0, 1000, by = 10), 
                                  kmax       = 1000)
{
  # The set maximum number of generations.
  
  tmax                <- max(timepoints)
  
  # This is an index counting the time points at which the population state is
  # returned.
  
  j                   <- 1
  
  # The initial population is completely free of mutations. 
  
  # Each individual in the population is characterized by the number of 
  # deleterious mutations that the individual carries.
  
  pop                 <- rep(0, N)
	
  # The popstate is the distribution of the number of deleterious mutations in 
  # the population, ranging from 0 to kmax.
  
  popstates           <- matrix(0,
                                nrow = length(timepoints),
                                ncol = kmax+1)
  
  rownames(popstates) <- timepoints
  colnames(popstates) <- 0:kmax
  
  # A vector of fitness values for each number of mutations.
	
	fitness.vector      <- (1-s)^(0:kmax)
	
	# The fitness of each individual per time point (generation).
	
	w                   <- matrix(0, nrow = tmax, ncol = N)
	
	rownames(w)         <- seq(1, tmax)
	colnames(w)         <- 1:N
	
	
	# Return popstate at time point 0.
	
	if (0 %in% timepoints) 
	{
		popstates[j, ]    <- get.popstate(pop, kmax)
		j                 <- j+1		
	}
	
	
	# For each time point.
		
	for(t in 1:tmax)
	{
	  # Find the fitness of each individual and puts these values into a vector.
	  
	  w[t, ]            <- fitness.vector[pop+1]
		
	  # Vector of indices giving all parental individuals.
	  
	  parents           <- get.parents(w[t, ])
		
		# Offspring population, identical to their parents.
		
		pop               <- pop[parents]
		
		# Mutations are added according to Poisson distribution.
		
		pop               <- pop + rpois(N, u)
		
		# Individuals are not allowed to have > kmax mutations.
		
		pop[pop>kmax]     <- kmax
		
		
		# Output the population state at the present time point.
		
		if (t %in% timepoints) 
		{
			popstates[j, ]   <- get.popstate(pop, kmax)
			j               <- j+1		
		}
		
		
		# If timepoint is tmax, then find final fitness of population.
		
		if (t %in% tmax) 
		{
		  w[tmax+1, ]     <- fitness.vector[pop+1]
		}
	}
	
	
	# Drop all columns consisting of only zeros.
	
	popstates           <- popstates[, colSums(popstates == 0) != nrow(popstates)]
	
	# Create a list of popstates and fitnesses (w) which can then be returned.
	
	output              <- list(popstates, w)
	
	return(output)
}



# Three Runs Of The Simulation ----

for (i in 1:3) 
{
  # Run the ratchet.simulator() function for specified number of time points 
  # (or generations).
  
  time.points <- seq(0, 10000)
  output      <- ratchet.simulator(timepoints = time.points)
  
  # Naming and saving data output.
  
  s1 = "Base_6_rep"
  s2 = i
  
  the.name    <- paste(s1, s2, sep = "")
  saveRDS(output, the.name)
}

