#DATE: 22/05/17


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


#Two different generation lengths (5000 or 10000)

gen_5000 <- matrix(c(1,3,5,7), nrow = 1)
gen_10000 <- matrix(c(2,4,6,8), nrow = 1)

for (n  in 1:length(gen_5000)) {
  rep1 <- readRDS(paste("NCS_", gen_5000[n], "_rep1", sep = ""))
  rep2 <- readRDS(paste("NCS_", gen_5000[n], "_rep2", sep = ""))
  rep3 <- readRDS(paste("NCS_", gen_5000[n], "_rep1", sep = ""))
  
  pop.states <- list(rep1[[1]], rep2[[1]], rep3[[1]])
  fitness <-  list(rep1[[3]], rep2[[3]], rep3[[3]])
  mTypes <- list(rep1[[2]], rep2[[2]], rep3[[2]])
  
  clicks <- matrix(0, nrow = 5001, ncol = 3) 
  for (i in 1:length(pop.states)) {
    clicks <- clicks + ratchet.clicks(pop.states[[i]])
  }
  
  # Average the number of clicks, fitness and mutation types
  
  clicks <- round(clicks/length(pop.states))
  log.fitness <- log(rowMeans(Reduce('+', fitness)/length(fitness)))
  mTypeTable <- as.table(round(colSums(Reduce('+', mTypes)/length(mTypes))))
  
  #Colour specifications
  
  c.specific <- matrix(c(4,1,2), nrow= 1, ncol = 3)
  
  
  # plotting the logarithm of fitness
  
  png(paste("Fitness_NCS_", gen_5000[n], ".png", sep = ""))
  plot(x = names(log.fitness) ,y = log.fitness, xlab = "Time (generations)", ylab = "log w")
  dev.off()
  
  # plotting the progress of Muller's ratchet
  
  png(paste("RClicks_NCS_", gen_5000[n], ".png", sep = ""))
  matplot(x=rownames(clicks),y=clicks, type="l",  col = c.specific,  xlab="Time (generations)", ylab="Number of ratchet clicks")
  legend("topleft", legend=c("Weak", "Mild", "Strong"), col=c.specific, lty=1, bty = "n")
  dev.off()
  
  # plotting the alogarithm of fitness against ratchet clicks
  
  png(paste("RClicks_Fitness_NCS_", gen_5000[n], ".png", sep = ""))
  matplot(x=clicks,y=log.fitness, type="l", col = c.specific, xlab="Number of ratchet clicks", ylab="log w")
  legend("topright", legend=c("Weak", "Mild", "Strong"), col=c.specific, lty=1, bty = "n")
  dev.off()
  
  # creating a table which contains the total sum of each mutation class
  
  table.name <- paste("mType_NCS_", gen_5000[n], sep = "")
  write.table(mTypeTable, table.name, sep="\t")
}


