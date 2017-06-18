#DATE: 11/05/17

# The following function simulates the reproduction step of the population, involving also selection.
# Individuals are sampled based on their fitness, and the function returns a list of the indices of
# parent individuals that release offspring into the next generation.

get.parents<-function(w)
{
  N<-length(w)
  parents<-sample(N,N,replace=TRUE,prob=w)
  parents
}

# The following function calculates the mutational state of the population,
# i.e., a vector containing the number of individuals that have 0, 1, 2, ..., kmax mutations.

get.popstate<-function(pop,kmax)
{
  pop <- rowSums (pop, na.rm = FALSE, dims = 1)
  popstate<-rep(0,kmax+1)
  for(i in 1:length(pop))
    popstate[pop[i]+1]<-popstate[pop[i]+1]+1
  popstate
}

fitness.vector <- function(pop, beta = -0.3) 
{
  s <- matrix(c(0.001, 0.1, 0.6), nrow=3, ncol=1) #Selection values for the three groups
  coeff <- pop %*% s #Summation
  fitness <-exp((log(1-beta*(coeff))/beta)) #Vector of fitness values
  fitness
}

add.newMutations <- function(N, u=0.02)
{
  add.pop <- matrix(0, nrow=N, ncol = 3)
  newMutations <- rpois(N,u)
  loc <- which(newMutations %in% c(1))
  
  for (j in 1:length(loc)) {
    x <- rgamma(1, shape = 0.23, scale = 0.05)
    if (x <= 0.005) {
      add.pop[(loc[j]),1] = add.pop[(loc[j]),1]+1
    }
    else if(x > 0.005 && x <= 0.05) {
      add.pop[(loc[j]),2] = add.pop[(loc[j]),2]+1
    }
    else {
      add.pop[(loc[j]),3] = add.pop[(loc[j]),3]+1
    }
  }
  
  add.pop
}


# Main simulation function. Input parameters are the population size N, 
# the selection coefficient s and the per genome mutation rate u. Moreover, a vector of 
# time points can be given at which the state of the population is returned. 
# The maximum of that vector determines the total number of generations that the simulation runs.
# Finally, the parameter kmax indicates the maximum number of deteterious mutations 
# that is kept track of. (Individuals can still have more mutations, but those are treated
# as if the individual had only kmax mutations.)

ratchet.simulator<-function(N=400,timepoints=seq(0,1000,by=10), kmax=1000) 
{
  tmax<-max(timepoints)
  j<-1				   # This is an index counting the time points at which the 
  # population state is returned.
  pop<-matrix(0,nrow=N, ncol=3)   	   # Each individual in the population is characterized by the number
  colnames(pop)<-c("nWeak","nMild","nStrong")
  # of deleterious mutations (weak, mild or strong) that the individual carries.
  # The initial population is completely free of any mutations.
  
  popstates<-matrix(0,nrow=length(timepoints),ncol=kmax+1)
  ## The popstate is the distribution of the number of deleterious
  ## mutations in the population, ranging from 0 to kmax.
  rownames(popstates)<-timepoints
  colnames(popstates)<-0:kmax

  if (0 %in% timepoints)               # return popstate at time point 0? 
  {
    popstates[j,]<-get.popstate(pop,kmax)
    j<-j+1		
  }
  
    
  for(t in 1:tmax)
  {
    w<-fitness.vector(pop)         # vector of fitness values for each individual
    parents<-get.parents(w)          # vector of indices giving all parental individuals
    
    pop<-pop[parents,]				 # offspring population, identical to their parents 
    
    pop<-pop + add.newMutations(N)
    pop[pop>kmax]<-kmax				            # individuals are not allowed to have >kmax mutations
    
    if (t %in% timepoints)			 # Output the population state at the present time point? 
    {
      popstates[j,]<-get.popstate(pop,kmax)
      j<-j+1		
    }
    
  }
  
  popstates<-popstates[,colSums(popstates==0) !=nrow(popstates)]  # drop all columns with only zeros in it
  
  pop
  return(pop)
}

# The following function takes the state of the population as input and from this 
# calculates the number of ratchet clicks that have occured. 

ratchet.clicks<-function(pop.states)
{
  clicks<-rep(0,length(pop.states[,1]))
  names(clicks)<-row.names(pop.states)
  for(i in 1:length(clicks))
  {
    j<-0
    while ((pop.states[i,j+1]==0)&&(j<length(pop.states[1,]))) j<-j+1
    clicks[i]<-j
  }
  clicks
}


# Example run of the simulation:

time.points<-seq(0,5000)
pop.states<-ratchet.simulator(timepoints=time.points)
clicks<-ratchet.clicks(pop.states)

# plotting the population state, i.e. the distribution of classes of individuals with different numbers of mutations:

gen<-2000
barplot(pop.states[gen,], xlab="Number of mutations", ylab="Frequency", main=paste("Distribution of mutations at generation",gen))
box()

# plotting the progress of Muller's ratchet:

plot(x=names(clicks),y=clicks, type="l", xlab="Generation", ylab="Ratchet clicks")

# plotting the fitness of individuals over


# "Movie" showing the distribution of mutations over time:

for(i in seq(1,length(time.points),by=2))
{
  barplot(pop.states[i,]/pop.states[1,1], ylim=c(0,0.5), xlab="Number of deleterious mutations", ylab="Frequency", main=paste("Distribution of mutations at generation",time.points[i]))
  Sys.sleep(0.10)
}