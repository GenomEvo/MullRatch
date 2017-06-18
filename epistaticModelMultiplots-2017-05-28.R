#DATE: 28/05/17

#install.packages('latex2exp')

library(latex2exp)

# The following function takes the state of the population as input and from this 
# calculates the number of ratchet clicks that have occured for each class. 

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


#Generation length of 10000

gen <- matrix(c(2,4,6,8,10,12,14,16,18,20,22,24), nrow = 1)

# Modifying the plot parameters
par(mfrow=c(4,3))
par(mar = c(0, 0, 0, 0), oma = c(7,5,4,4))

# Beta values
b <- matrix(c(-0.5,-1,-2,-4), nrow = 1)

# Selection values
s <- matrix(c(0.1,0.05,0.01), nrow = 1)

for (n  in 1:length(gen)) {
  # Read each replicate
  rep1 <- readRDS(paste("CS_", gen[n], "_rep1", sep = ""))
  rep2 <- readRDS(paste("CS_", gen[n], "_rep2", sep = ""))
  rep3 <- readRDS(paste("CS_", gen[n], "_rep3", sep = ""))
  
  pop.states <- list(rep1[[1]], rep2[[1]], rep3[[1]])
  fitness <-  list(rep1[[2]], rep2[[2]], rep3[[2]])
  
  # For each replicate calculate the number of clicks and the logarithm of average fitness per time point
  for (i in 1:length(pop.states)) {
    clicks <- ratchet.clicks(pop.states[[i]])
    pop.states[[i]] <- clicks
    
    log.fitness <- log(rowMeans(fitness[[i]]))
    fitness[[i]] <- log.fitness
  }
  
  # plotting the logarithm of fitness against number of ratchet clicks for the three replicates
  matplot(
    x = pop.states[[1]],
    y = fitness[[1]],
    xlim = c(0, 1500),
    ylim = c(-7, 0),
    type = "l",
    lwd = 1.5,
    xlab = "",
    ylab = "",
    xaxt="n",
    yaxt="n"
  )
  
  
  for (j in 2:length(pop.states)) {
    lines(x = pop.states[[j]], y =  fitness[[j]], col=1, lwd = 1.5)
  }
  
  
  if (n %in% c(1,2,3)) {
    mtext(text = TeX(sprintf("$\\s = %g$", s[n])), las = 1, side=3,line=1, cex = 0.7)
  }
  
  if (n %in% c(1,4,7,10)) {
    axis(2, at = seq(-7, 0, 2))
  }
  
  if (n %in% c(3,6,9,12)) {
    loc_b <- which(c(3,6,9,12) == n)
    mtext(text = TeX(sprintf("$\\beta = %g$", b[loc_b])), las = 1, side=4,line=1, cex = 0.7)
  }
  
  if (n %in% c(10,11,12)) {
    axis(1, at = seq(0, 1500, 300))
  }
  
  box(col = "black")
  
}

mtext("Number of ratchet clicks", side = 1, outer = TRUE, cex = 1, line = 4, col = "grey20")
mtext(TeX('$\\log(w)$'), side = 2, outer = TRUE, cex = 1, line = 3, col = "grey20")
