#Date : 22/05/17

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

gen_5000 <- matrix(c(1,3,5), nrow = 1)
gen_10000 <- matrix(c(2,4,6), nrow = 1)

for (n in 1:length(gen_5000)) {
  # Read each replicate
  rep1 <- readRDS(paste("Base_", gen_10000[n], "_rep1", sep = ""))
  rep2 <- readRDS(paste("Base_", gen_10000[n], "_rep2", sep = ""))
  rep3 <- readRDS(paste("Base_", gen_10000[n], "_rep3", sep = ""))
  
  pop.states <- list(rep1[[1]], rep2[[1]], rep3[[1]])
  fitness <-  list(rep1[[2]], rep2[[2]], rep3[[2]])
  
  
  clicks <- rep(0,10001)
  # Sums the number of clicks per replicate
  for (i in 1:length(pop.states)) {
    clicks <- clicks + ratchet.clicks(pop.states[[i]])
  }
  
  # Averages the number of clicks and fitness
  
  clicks <- round(clicks/length(pop.states))
  log.fitness <- log(rowMeans(Reduce('+', fitness)/length(fitness)))
  
  # plotting the progress of Muller's ratchet
  
  png(paste("RClicks_Base_", gen_10000[n], ".png", sep = ""))
  plot(x=names(clicks),y=clicks, type="l", xlab="Time (generations)", ylab="Number of ratchet clicks")
  dev.off()
  
  # plotting the logarithm of fitness
  
  png(paste("Fitness_Base_", gen_10000[n], ".png", sep = ""))
  plot(x = names(log.fitness), y = log.fitness, xlab = "Time (generations)", ylab = "log w")
  dev.off()
  
  # plotting the logarithm of fitness against ratchet clicks
  
  png(paste("RClicks_Fitness_Base_", gen_10000[n], ".png", sep = ""))
  matplot(x=clicks, y=log.fitness, type="l", xlab="Number of ratchet clicks", ylab="log w")
  dev.off()
  
}

