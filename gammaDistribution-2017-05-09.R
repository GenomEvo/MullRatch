#DATE: 09/05/17

rgamma(n=5, shape=0.23, rate=1/0.01)

?rgamma

N=4
u=0.02
tmax=6

#constructed and populated a vector of my population against time steps
populationeffect <- matrix(0, nrow = tmax, ncol = N)
rownames(populationeffect) <- 1:tmax
colnames(populationeffect) <- 1:N

populationeffect


# determine the number of mutations that each member of the population will get

#example code of a for loop
#for (year in 2010:2015){
 # print(paste("The year is", year))
#}

for (i in 1:length(populationeffect) {
  populationeffect[1,]<-rpois(N,u)
}


for i in 1:length(N){
  x=rpois(N,u)
  for i in 1:x
  effect=rgamma(1, shape = 0.23, rate =1/0.01)
  
  print?
  
}

x=rpois(N,u)
rgamma(1, shape = 0.23, rate =1/0.01 )

