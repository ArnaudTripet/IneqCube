setwd("C://Users//tripeta2//switchdrive//balanced_with_inequalities//")

rm(list=ls())

### Load Function

source("R//CubeIneq.R")
source("R//FastFlightCubeIneq.R")
source("R//IneqCstr.R")

###

### Creater dataset

N <- 100 # Size of the dataset
n <- 20 # Size of the sample

X <- cbind(runif(N),runif(N)) # dataset

library(sampling)
pik <- inclusionprobabilities(runif(N),n) # Inclusion probabilities

sum(pik) # should be equal to n

###

### Inequality Constraints

LG1 <- ineq.cstr(X,pik,"two.sided",1) # inequality less and greater than 1

r <- ceiling(t(LG1) %*% pik)

###

### sample

s <- cube.ineq(X,pik,LG1,r) # Draw sample 

plot(X) # plot the population
points(X[s==1,],pch=19) # show selected sample
