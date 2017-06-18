#DATE: 29/05/17


# ************META TITLE ----


# Ratchet Clicks Function ----

# The following function takes the state of the population as input and from
# this calculates the number of ratchet clicks that have occured. 

ratchet.clicks  <- function(pop.states)
{
  clicks        <- rep(0, length(pop.states[, 1]))
  names(clicks) <- row.names(pop.states)
  
  for(i in 1:length(clicks))
  {
    j           <- 0
    
    while ((pop.states[i, j+1] == 0)&&(j < length(pop.states[1, ]))) 
    {
      j         <- j+1
    }
    
    clicks[i]   <- j
  }
  
  clicks
}


# ************Subheading??? ----

# Two different generation lengths (5000 or 10000).

gen_5000      <- matrix(c(1, 3, 5), 
                        nrow = 1)

gen_10000     <- matrix(c(2, 4, 6), 
                        nrow = 1)


# ************Add comment about what the below block of code does.

for (n in 1:length(gen_10000)) 
{
  # Read each replicate.
  
  rep1        <- readRDS(paste("Base_", gen_10000[n], "_rep1", sep = ""))
  rep2        <- readRDS(paste("Base_", gen_10000[n], "_rep2", sep = ""))
  rep3        <- readRDS(paste("Base_", gen_10000[n], "_rep3", sep = ""))
  
  # Consolidating all replicates into lists.

  pop.states  <- list(rep1[[1]], 
                      rep2[[1]], 
                      rep3[[1]])
  
  fitness     <-  list(rep1[[2]], 
                       rep2[[2]], 
                       rep3[[2]])
  
  # Run ratchet.clicks() function for specified number of time points (or 
  # generations)
  
  # ************Sum the number of clicks per replicate.
  
  # For each replicate calculate the number of clicks and the logarithm of average fitness per time point
  for (i in 1:length(pop.states)) {
    clicks <- ratchet.clicks(pop.states[[i]])
    pop.states[[i]] <- clicks
    
    log.fitness <- log(rowMeans(fitness[[i]]))
    fitness[[i]] <- log.fitness
  }
  
  
  # PLOTS / CREATING FIGURES ----
  
  
  # Plotting The Progress Of Muller's Ratchet ----
  
  # Saving plot as an image.
  
  png(paste("RClicks_Base_", gen_10000[n], ".png", sep = ""))
  
  # Create plot of the number of ratchet clicks over time (generations).
  
  plot(x    = names(pop.states[[1]]), 
       y    = pop.states[[1]], 
       type = "l", 
       xlab = "Time (generations)", 
       ylab = "Number of ratchet clicks")
  
  for (j in 2:length(pop.states)) {
    lines(x = names(pop.states[[j]]), y =  pop.states[[j]])
  }
  
  # ************What is this?
  
  dev.off()
  
  
  # Plotting The Logarithm Of Fitness ----
  
  # Saving plot as an image.
  
  png(paste("Fitness_Base_", gen_10000[n], ".png", sep = ""))
  
  # Create plot of the average fitness loss over time (generations).
  
  plot(x    = names(fitness[[1]]), 
       y    = fitness[[1]], 
       type = 'p',
       xlab = "Time (generations)", 
       ylab = "log w"
       )
  
  # ************What is this?
  for (j in 2:length(fitness)) {
    points(x = names(fitness[[j]]), y =  fitness[[j]])
  }
  dev.off()
  
  
  # Plotting The Logarithm Of Fitness Against Ratchet Clicks ----
  
  # Saving plot as an image.
  
  png(paste("RClicks_Fitness_Base_", gen_10000[n], ".png", sep = ""))
  
  # Create plot of fitness against increasing number of ratchet clicks.
  
  matplot(x    = pop.states[[1]], 
          y    = fitness[[1]], 
          type = "p", 
          xlab = "Number of ratchet clicks", 
          ylab = "log w"
          )
  
  # ************What is this?
  for (j in 2:length(pop.states)) {
    points(x = pop.states[[j]], y =  fitness[[j]])
  }
    dev.off()
}
