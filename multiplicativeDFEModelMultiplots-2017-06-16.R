# DATE: 16/06/17


# OUTPUT OF COMPOSITE MULTIPLOTS PLOTS ----


#install.packages('latex2exp') 

library(latex2exp)


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
      
      while ((pop.states[[r]][i, j+1] == 0) && ( j < length(pop.states[[r]][1,])))
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


# A generation length of 10000.

gen                            <- matrix(c(2), nrow = 1)


# Modifying the plot parameters.

par(mfrow = c(1, 1))

par(cex   = 0.6)

par(tcl   = -0.25)

par(mgp   = c(5, 1, 0))

par(mar   = c(0, 0, 2, 0), 
    oma   = c(7, 7, 5, 5))

plot_colours                   <- matrix(c(4, 1, 2), 
                                         nrow = 1, 
                                         ncol = 3)


for (n  in 1:length(gen)) 
  {
  # Read each replicate.
  
  rep1              <- readRDS(paste("NCS_Base_", gen[n], "_rep1", sep = ""))
  rep2              <- readRDS(paste("NCS_Base_", gen[n], "_rep2", sep = ""))
  rep3              <- readRDS(paste("NCS_Base_", gen[n], "_rep3", sep = ""))

  # Consolidating all replicates into lists.
  
  pop.states        <- list(rep1[[1]], 
                            rep2[[1]], 
                            rep3[[1]])
  
  fitness           <-  list(rep1[[3]], 
                             rep2[[3]], 
                             rep3[[3]])

    
  # For each replicate calculate the number of clicks and the logarithm of 
  # average fitness per time point.
  
  for (i in 1:length(pop.states)) 
  {
    clicks          <- ratchet.clicks(pop.states[[i]])
    
    pop.states[[i]] <- clicks
    
    log.fitness     <- log(rowMeans(fitness[[i]]))
    
    fitness[[i]]    <- log.fitness
  }
  
  
  # Plotting The Progress Of Muller's Ratchet ----
  
  # Plotting the logarithm of fitness against number of ratchet clicks for the 
  # three replicates.
  
  matplot(x    = pop.states[[1]],
          y    = fitness[[1]],
          xlim = c(0, 800),
          ylim = c(-4.5, 0),
          col  = plot_colours,
          type = "l",
          lwd  = 1.5, 
          xlab = "Number of ratchet clicks",
          ylab = TeX('$\\log(w)$'),
          cex  = 0.7)
  
  mtext("Number of ratchet clicks", 
        side  = 1, 
        outer = TRUE, 
        cex   = 1, 
        line  = 4, 
        col   = "grey20")
  
  mtext(TeX('$\\log(w)$'), 
        side  = 2, 
        outer = TRUE, 
        cex   = 1, 
        line  = 3, 
        col   = "grey20")
  
  axis(2, at = seq(-4.5, 0, 1))
  axis(2, at = seq(-4, 0, 1))
  
  
  for (j in 2:length(pop.states)) 
  {
    lines(x   = pop.states[[j]][, 1], 
          y   =  fitness[[j]], 
          col = plot_colours[1], 
          lwd = 1.5)
    
    lines(x   = pop.states[[j]][, 2], 
          y   =  fitness[[j]], 
          col = plot_colours[2], 
          lwd = 1.5)
    
    lines(x   = pop.states[[j]][, 3], 
          y   = fitness[[j]], 
          col = plot_colours[3], 
          lwd = 1.5)
  }
  
  legend(x = "topright", 
         inset  = 0, 
         legend = c("Weak", "Mild", "Strong"), 
         col    = plot_colours, 
         lwd    = 5, 
         cex    = 1,  
         bty    = "n", 
         xpd    = NA)
}

