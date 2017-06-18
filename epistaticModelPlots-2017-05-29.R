# DATE: 29/05/17

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

#Two different generation lengths (5000 or 10000)
gen_5000 <- matrix(c(1,3,5,7,9,11,13,15,17,19,21,23), nrow = 1)
gen_10000 <- matrix(c(2,4,6,8,10,12,14,16,18,20,22,24), nrow = 1)

for (n in 1:length(gen_5000)){
  #Read each replicate

  rep1 <- readRDS(paste("CS_", gen_5000[n], "_rep1", sep = ""))
  rep2 <- readRDS(paste("CS_", gen_5000[n], "_rep2", sep = ""))
  rep3 <- readRDS(paste("CS_", gen_5000[n], "_rep3", sep = ""))
  
  pop.states <- list(rep1[[1]], rep2[[1]], rep3[[1]])
  fitness <-  list(rep1[[2]], rep2[[2]], rep3[[2]])
  
  for (i in 1:length(pop.states)) {
    clicks <- ratchet.clicks(pop.states[[i]])
    pop.states[[i]] <- clicks
    
    log.fitness <- log(rowMeans(fitness[[i]]))
    fitness[[i]] <- log.fitness
  }
  
  # plotting the progress of Muller's ratchet
  
  png(paste("RClicks_CS_", gen_5000[n], ".png", sep = ""))
  plot(x=names(pop.states[[1]]),y=pop.states[[1]], type="l", xlab= "Time (generations)", ylab="Number of ratchet clicks")
  
  for (j in 2:length(pop.states)) {
    lines(x = names(pop.states[[j]]), y =  pop.states[[j]], col=1, lwd = 1)
  }
  
  dev.off()
  
  # plotting the logarithm of fitness
  
  png(paste("Fitness_CS_", gen_5000[n], ".png", sep = ""))
  plot(x = names(fitness[[1]]), y = fitness[[1]], type = 'l', lwd = 1, xlab = "Time (generations)", ylab = "log w")
  
  for (j in 2:length(pop.states)) {
    lines(x = names(fitness[[j]]) , y =  fitness[[j]], col=1, lwd = 1)
  }
  
  dev.off()
  
  # plotting the logarithm of fitness against number of ratchet clicks
  
  png(paste("RClicks_Fitness_CS_", gen_5000[n], ".png", sep = ""))
  matplot(x=pop.states[[1]], y=fitness[[1]], type="p", pch = 20, xlab="Number of ratchet clicks", ylab="log w")
  
  for (j in 2:length(pop.states)) {
    points(x = pop.states[[j]], y =  fitness[[j]], col=1, pch = 20)
  }
  
  dev.off()
}




