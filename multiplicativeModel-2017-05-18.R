#DATE:18/05/17

# SUBFUNCTIONS ----

# Reproductive Step ----

# The following function simulates the reproduction step of the population,
# also involving selection. Individuals are sampled based on their fitness, and
# the function returns a list of the indices of parent individuals that release
# offspring into the next generation.

get.parents <- function(w)
{
	N       <- length(w)
	parents <- sample(N, N, replace = TRUE, prob = w)
	parents
}

# Distribution Of The Number Of Mutations Per Generation ----

# The following function calculates the mutational state of the population,
# i.e., a vector containing the number of individuals that have 0, 1, 2, ...,
# kmax mutations.

get.popstate <- function(pop, kmax)
{
	popstate <- rep(0, kmax+1)
	for(i in 1:length(pop))
		popstate[pop[i]+1] <- popstate[pop[i]+1]+1
	popstate
}



# MAIN SIMULATION CODE ----

# Muller's Ratchet Simulation Function ----

# Input parameters are the population size N, the selection coefficient s and
# the per genome mutation rate u. Moreover, a vector of time points can be
# given at which the state of the population is returned. The maximum of that
# vector determines the total number of generations that the simulation runs.
# Finally, the parameter kmax indicates the maximum number of deteterious
# mutations that is kept track of. (Individuals can still have more mutations,
# but those are treated as if the individual had only kmax mutations.)

ratchet.simulator <- function(N          = 400, 
                              s          = 0.01, 
                              u          = 0.02, 
                              timepoints = seq(0, 1000, by = 10), 
                              kmax       = 1000)
{
	tmax                <- max(timepoints)
	j                   <- 1		                                               # This is an index counting the time points at which the population state is returned.
	pop                 <- rep(0,N)                                            # Each individual in the population is characterized by the number of deleterious mutations that the individual carries. The initial population is completely free of mutations.
	popstates           <- matrix(0, nrow = length(timepoints), ncol = kmax+1) # The popstate is the distribution of the number of deleterious mutations in the population, ranging from 0 to kmax.
	rownames(popstates) <- timepoints
	colnames(popstates) <- 0:kmax
	fitness.vector      <- (1-s)^(0:kmax)                                      # vector of fitness values for each number of mutations
	w                   <- matrix(0, nrow = tmax, ncol = N)                    # the fitness of each individual per generation/time point
	rownames(w)         <- seq(1,tmax)
	colnames(w)         <- 1:N
	
	if (0 %in% timepoints)                                                     # return popstate at time point 0 
	{
		popstates[j,] <- get.popstate(pop,kmax)
		j             <- j+1		
	}
		
	for(t in 1:tmax)
	{
	  w[t,]         <- fitness.vector[pop+1]                                   # vector of fitness values for each individual
		parents       <- get.parents(w[t,])                                      # vector of indices giving all parental individuals
		pop           <- pop[parents]				                                     # offspring population, identical to their parents
		pop           <- pop+rpois(N,u)				                                   # mutations are added, according to Poisson distribution
		pop[pop>kmax] <- kmax				                                             # individuals are not allowed to have >kmax mutations
		
		if (t %in% timepoints)			                                             # Output the population state at the present time point 
		{
			popstates[j,] <- get.popstate(pop,kmax)
			j             <- j+1		
		}
	}
	
	popstates <- popstates[,colSums(popstates==0) !=nrow(popstates)]           # drop all columns with only zeros in it
	output    <- list(popstates, w)
	return(output)							                                               # return the list of popstates and fitnesses(w)
}

# Number of Ratchet Clicks ----

# The following function takes the state of the population as input and from
# this calculates the number of ratchet clicks that have occured. 

ratchet.clicks <- function(pop.states)
{
	clicks        <- rep(0,length(pop.states[,1]))
	names(clicks) <- row.names(pop.states)
	for(i in 1:length(clicks))
	{
		j <- 0
		while ((pop.states[i,j+1]==0)&&(j<length(pop.states[1,]))) j <- j+1
		clicks[i] <- j
	}
	clicks
}



# PLOTTING THE OUTPUT FROM MULLER'S RATCHET SIMULATION ----

# Run Of The Simulation ----

time.points <- seq(0,5000)                                                   # Parameter value
output      <- ratchet.simulator(timepoints = time.points)
pop.states  <- output[[1]]
fitness     <- output[[2]]
clicks      <- ratchet.clicks(pop.states)

# Plotting the population state, i.e. the distribution of classes of
# individuals with different numbers of mutations.

gen <-5001
barplot(pop.states[gen,], 
        xlab = "Number of mutations", 
        ylab = "Frequency", 
        main = paste("Distribution of mutations at generation",gen))
box()

# Plotting the progress of Muller's ratchet.

plot(x    = names(clicks), 
     y    = clicks, 
     type = "l",
     xlab = "Generation", 
     ylab = "Ratchet clicks")

# Plotting the average fitness loss.

average.loss <- (1-rowMeans(fitness))
plot(x    = seq(1,5000), 
     y    = average.loss, 
     xlab = "Generation", 
     ylab = "Average Loss")
