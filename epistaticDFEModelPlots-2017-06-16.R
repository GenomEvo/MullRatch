#DATE: 16/06/17


# OUTPUT OF INDIVIDUAL PLOTS ----


# Ratchet Clicks Function ----

# The following function takes the state of the population as input and from
# this calculates the number of ratchet clicks that have occured. 

ratchet.clicks                 <- function(pop.states)
{
  threeFitnessClicks           <- matrix(0, 
                                         nrow = length(pop.states[[1]][, 1]), 
                                         ncol = length(pop.states))
  
  for (r in 1:length(pop.states))
  {
    clicks                     <- rep(0, length(pop.states[[r]][, 1]))
    
    for(i in 1:length(clicks))
    {
      j                        <- 0
      
      while((pop.states[[r]][i, j+1] == 0)&&(j < length(pop.states[[r]][1, ])))
      {
        j                      <- j+1
      }
      
      clicks[i]                <- j
    }
    
    threeFitnessClicks[, r]    <- clicks
  }
  
  rownames(threeFitnessClicks) <- row.names(pop.states[[1]])
  threeFitnessClicks
}


# Two different generation lengths (5000 or 10000).

gen_5000      <- matrix(c(1, 3, 5, 7), 
                        nrow = 1)

gen_10000     <- matrix(c(2, 4, 6, 8), 
                        nrow = 1)


for (n  in 1:length(gen_10000)) 
{
  # Read each replicate.
  
  rep1 <- readRDS(paste("NCS_", gen_10000[n], "_rep1", sep = ""))
  rep2 <- readRDS(paste("NCS_", gen_10000[n], "_rep2", sep = ""))
  rep3 <- readRDS(paste("NCS_", gen_10000[n], "_rep3", sep = ""))
  
  # Consolidating all replicates into lists.
  
  pop.states  <- list(rep1[[1]], 
                      rep2[[1]], 
                      rep3[[1]])
  
  fitness     <-  list(rep1[[3]], 
                       rep2[[3]], 
                       rep3[[3]])
  
  mTypes      <- list(rep1[[2]], 
                      rep2[[2]], 
                      rep3[[2]])
  
  # For each replicate calculate the number of clicks and the logarithm of 
  # average fitness per time point.
  
  for (i in 1:length(pop.states)) 
  {
    clicks          <- ratchet.clicks(pop.states[[i]])
    pop.states[[i]] <- clicks
    
    log.fitness     <- log(rowMeans(fitness[[i]]))
    fitness[[i]]    <- log.fitness
  }
  
  mTypeTable <- as.table(round(colSums(Reduce('+', mTypes)/length(mTypes))))

  
  # Colour Specifications ----

  c.specific <- matrix(c(4, 1, 2), 
                       nrow = 1, 
                       ncol = 3)

  
  # Plotting The Progress Of Muller's Ratchet ----
  
  # Saving plot as an image.
  
  png(paste("RClicks_NCS_", gen_10000[n], ".png", sep = ""))
  
  # Create plot of the number of ratchet clicks over time (generations).
  
  matplot(x    = rownames(pop.states[[1]]), 
          y    = pop.states[[1]], 
          type = "l", 
          col  = c.specific,  
          xlab = "Time (generations)", 
          ylab = "Number of ratchet clicks")
  
  for (j in 2:length(pop.states)) 
  {
    lines(x   = rownames(pop.states[[j]]), 
          y   = pop.states[[j]][, 1], 
          col = plot_colours[1],
          lwd = 1.5)
    
    lines(x   = rownames(pop.states[[j]]), 
          y   = pop.states[[j]][, 2], 
          col = plot_colours[2],
          lwd = 1.5)
    
    lines(x   = rownames(pop.states[[j]]), 
          y   = pop.states[[j]][, 3], 
          col = plot_colours[3],
          lwd = 1.5)
  }
  
  legend("topleft", 
         legend = c("Weak", "Mild", "Strong"), 
         col    = c.specific, 
         lty    = 1, 
         bty    = "n")
  
  dev.off()
  
  
  # Plotting The Logarithm Of Fitness ----
  
  # Saving plot as an image.
  
  png(paste("Fitness_NCS_", gen_10000[n], ".png", sep = ""))
  
  # Create plot of the average fitness loss over time (generations).
  
  matplot(x    = names(fitness[[1]]), 
          y    = fitness[[1]],  
          type = "l", 
          xlab = "Time (generations)", 
          ylab = "log w")
  
  for (j in 2:length(pop.states)) 
  {
    lines(x = names(fitness[[j]]), 
          y =  fitness[[j]])
  }
  
  dev.off()
  
  
  # Plotting The Logarithm Of Fitness Against Ratchet Clicks ----
  
  # Saving plot as an image.
  
  png(paste("RClicks_Fitness_NCS_", gen_10000[n], ".png", sep = ""))
  
  # Create plot of fitness against increasing number of ratchet clicks.
  
  matplot(x    = pop.states[[1]], 
          y    = fitness[[1]], 
          type = "l", 
          col  = c.specific, 
          xlab = "Number of ratchet clicks", 
          ylab = "log w")
  
  for (j in 2:length(pop.states)) 
  {
    lines(x   = pop.states[[j]][, 1], 
          y   = fitness[[j]], 
          col = plot_colours[1], 
          lwd = 1.5)
    
    lines(x   = pop.states[[j]][, 2], 
          y   = fitness[[j]], 
          col = plot_colours[2], 
          lwd = 1.5)
    
    lines(x   = pop.states[[j]][, 3], 
          y   = fitness[[j]], 
          col = plot_colours[3], 
          lwd = 1.5)
  }
  
  legend("topright", 
         legend = c("Weak", "Mild", "Strong"), 
         col    = c.specific, 
         lty    = 1, 
         bty    = "n")
  
  dev.off()

  
  # Creating a table which contains the total sum of each mutation class.

  table.name <- paste("mType_NCS_", 
                      gen_10000[n], 
                      sep = "")
  
  write.table(mTypeTable, 
              table.name, 
              sep="\t")
}
