#DATE: 15/05/17


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
  popstate<-rep(0,kmax+1)
  for(i in 1:length(pop))
    popstate[pop[i]+1]<-popstate[pop[i]+1]+1
  popstate
}

fitness.vector <- function(pop, beta = -2) 
{
  #Beta parameterises epistasis
  
  selection <- matrix(c(0.001, 0.1, 0.6), nrow=3, ncol=1) #Selection values for the three groups
  totalselection <- pop %*% selection #Vector of the summed selection values for each individual
  fitness <-(1-beta*(totalselection))^(1/beta) #Vector of fitness values
  fitness
}

add.newMutations <- function(N, v1=0.49, v2=0.46, v3=0.05, u = 0.02)
{
  #The probability that a given new mutation is weak, mild or strong (v1,v2,v3). 
  #The per genome mutation rate 'u'
  
  mWeak <- rpois(N, u*v1) # mutations are added to the weak class, according to Poisson distribution
  mMild <- rpois(N, u*v2) # mutations are added to the mild, according to Poisson distribution
  mStrong <- rpois(N, u*v3)
  newmutmatrix <- cbind(mWeak, mMild, mStrong)
  row.names(newmutmatrix) <- (1:N)
  newmutmatrix
}


# Main simulation function. Input parameters are the population size N, 
# a vector of time points can be given at which the state of the population is returned. 
# The maximum of that vector determines the total number of generations that the simulation runs.
# Finally, the parameter kmax indicates the maximum number of deteterious mutations 
# that is kept track of. (Individuals can still have more mutations, but those are treated
# as if the individual had only kmax mutations.)

ratchet.simulator<-function(N=400,timepoints=seq(0,1000,by=10), kmax=300) 
{
  tmax<-max(timepoints)
  j<-1				   # This is an index counting the time points at which the 
  # population state is returned.
  pop<-matrix(0,nrow=N, ncol=3)   	   # Each individual in the population is characterized by the number
  colnames(pop)<-c("nWeak","nMild","nStrong")
  # of deleterious mutations (weak, mild or strong) that the individual carries.
  # The initial population is completely free of any mutations. 
  
  d <- dim(matrix(0,nrow=length(timepoints),ncol=kmax+1)) #Matrix dimensions for each class (weak, mild or strong)
  
  nWeakpopstates<-matrix(0,d[1], d[2])
  nMildpopstates<-matrix(0,d[1], d[2])
  nStrongpopstates<-matrix(0,d[1], d[2])
  
  popstates <- list(nWeakpopstates, nMildpopstates,nStrongpopstates)
  ## The popstates is the distribution of the number of deleterious
  ## mutations in the population for each class, ranging from 0 to kmax.
  
  w<-matrix(0,nrow=tmax,ncol=N) #Fitness for each individual per generation/time point
  rownames(w)<-seq(1,tmax)
  colnames(w)<-1:N
  
  for (i in 1:length(popstates)) {
    rownames(popstates[[i]])<- timepoints
    colnames(popstates[[i]])<- 0:kmax
  }

  if (0 %in% timepoints)               # return popstate at time point 0? 
  {
    for (i in 1:length(popstates)) {
      popstates[[i]][j,]<-get.popstate(pop[,i],kmax)
    }
    j<-j+1		
  }
  
    
  for(t in 1:tmax)
  {
    w[t,]<-fitness.vector(pop)         # vector of fitness values for each individual
    parents<-get.parents(w[t,])          # vector of indices giving all parental individuals
    
    pop<-pop[parents,]				 # offspring population, identical to their parents 
    
    pop<-pop + add.newMutations(N)
    pop[pop>kmax]<-kmax				            # individuals are not allowed to have >kmax mutations
    
    if (t %in% timepoints)			 # Output the population state at the present time point? 
    {
      for (i in 1:length(popstates)){
        popstates[[i]][j,]<-get.popstate(pop[,i],kmax)
      }
      j<-j+1		
    }
    
  }
  
  for(i in 1:length(popstates)) {
    popstates[[i]]<- popstates[[i]][,colSums(popstates[[i]]==0) !=nrow(popstates[[i]])]
  }
  
  output <- list(popstates,pop,w)
  return(output)
}

# The following function takes the state of the population as input and from this 
# calculates the number of ratchet clicks that have occured for each class. 

ratchet.clicks<-function(pop.states)
{
  threeFitnessClicks <- matrix(0, nrow = length(pop.states[[1]][,1]), ncol = length(pop.states))
  for (r in 1:length(pop.states)){
    clicks<-rep(0,length(pop.states[[r]][,1]))
    for(i in 1:length(clicks)){
      j<-0
      while ((pop.states[[r]][i,j+1]==0)&&(j<length(pop.states[[r]][1,]))) j<-j+1
      clicks[i]<-j
    }
    threeFitnessClicks[,r] <-clicks
  }
  rownames(threeFitnessClicks)<-row.names(pop.states[[1]])
  threeFitnessClicks
}


# Example run of the simulation:

time.points<-seq(0,5000)
output<-ratchet.simulator(timepoints=time.points)
pop.states <- output[[1]]
fitness <- output[[3]]
type <- output[[2]]
clicks<-ratchet.clicks(pop.states)

# plotting the population state, i.e. the distribution of classes of individuals with different numbers of mutations:

#gen<-5001
#barplot(pop.states[[1]][gen,], xlab="Number of mutations", ylab="Frequency", main=paste("Distribution of mutations at generation",gen))
#box()

# plotting the progress of Muller's ratchet:

matplot(x=rownames(clicks),y=clicks, type="l", xlab="Generation", ylab="Ratchet clicks")

# plotting the average population fitness level
average.loss<- 1-rowMeans(fitness)
plot(x = seq(1,5000),y = average.loss, xlab = "Generation", ylab = "Average Loss")

# creating a table stating the sum of each mutation class

mType <- colSums(type)
mTypeTable <- as.table(mType)
mTypeTable
