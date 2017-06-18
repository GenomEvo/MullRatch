#DATE: 25/05/17


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
# ************function does, unless subtitle is self-explanatory). Include an 
# ************explanation of pop (which is a vector containing count of 
# ************mutations each individual carries) and beta (which parameterises
# ************epistasis).

fitness.vector       <- function(pop, beta = -4) 
{
  # Selection values for the three mutation types (weak, mild and strong).
  
  selection          <- matrix(c(4.45E-05, 0.0032, 0.0218), 
                               nrow = 3, 
                               ncol = 1)
  
  # Find all individuals/rows that carry only one mutation in any of the three
  # mutation types (weak, mild and strong).
  
  oneMutation        <- which(pop == 1, arr.in = TRUE)
  
  # Sum across the three columns (mutation types) of each individual/row and 
  # find which individuals/rows has a total of one mutation.
  
  # If there is only one individual with a total of one mutation, use the "if" 
  # statement. If there is two or more individuals with a total of one 
  # mutation, use the "else" statement.
  
  if (length(oneMutation) == 2)
  {
    SumPerIndividual <- which(sum(pop[oneMutation[, 1], ]) == 1)
  }
  
  else 
  { 
    SumPerIndividual <- which(rowSums(pop[oneMutation[, 1], ]) ==1)
  }
  
  # ************Add comment
  
  loc                <- oneMutation[SumPerIndividual,]
  
  # ************Add comment
  
  inherentWeak       <- (1-selection[1])
  inherentMild       <- (1-selection[2])
  inherentStrong     <- (1-selection[3])
  
  inherent           <- cbind(inherentWeak, inherentMild, inherentStrong)
  
  # A vector of the summed selection values for each individual.
  
  totalselection     <- pop %*% selection
  
  # A vector of fitness values.
  
  fitness            <- (1-beta*(totalselection))^(1/beta)
  
  
  # *****************************************************************************************************************
  
  if (length(loc) != 0) {                      #At least one individual has only 1 mutation
    for (i in 1:length(SumPerIndividual)){
      if (length(loc) == 2) {                  #Only one individual have a total of one mutation
        fitness[loc[1]] <- inherent[loc[2]]
      }
      else {                                   #Two or more individuals have a total of one mutation
        fitness[loc[i,1]] <- inherent[loc[i,2]]
      }
    }
  }
  
  fitness
}

add.newMutations <- function(N, v1=0.15, v2=0.45, v3=0.4, u = 0.175)
{
  #The probability that a given new mutation is weak, mild or strong (v1,v2,v3). 
  #The per genome mutation rate 'u'
  
  mWeak <- rpois(N, u*v1) # mutations are added to the weak class, according to Poisson distribution
  mMild <- rpois(N, u*v2) # mutations are added to the mild class, according to Poisson distribution
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

ratchet.simulator<-function(N=400,timepoints=seq(0,1000,by=10), kmax=1000) 
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
  
  w<-matrix(0,nrow=tmax+1,ncol=N) #Fitness for each individual per generation/time point
  rownames(w)<-seq(0,tmax)
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
   
    if (t %in% tmax) {
      w[tmax+1,] <- fitness.vector(pop)
    } 
  }
  
  for(i in 1:length(popstates)) {
    popstates[[i]]<- popstates[[i]][,colSums(popstates[[i]]==0) !=nrow(popstates[[i]])]
  }
  
  output <- list(popstates,pop,w)
  return(output)
}



# Three runs of the simulation:

for (i in 1:3) {
  time.points<-seq(0,5000)
  output<-ratchet.simulator(timepoints=time.points)
  
  s1 = "NCS_7_rep"
  s2 = i
  
  the.name <- paste(s1,s2, sep = "")
  saveRDS(output, the.name)
}

