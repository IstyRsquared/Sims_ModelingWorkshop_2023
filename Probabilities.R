## title: Simulation models for modeling methods workshop: Probabilities in R
## author: Kristyna Rysava
## date: September 20, 2023

rm(list=ls())

### 1. Distributions and random number generators
## Normal distribution
# Generate 1000 normally distributed numbers with mean=5 and sd=2 and draw a histogram
y <- rnorm(1000, mean=5, sd=2) 
hist(y, nclass=50, freq=FALSE, main="Normal distribution", xlab="Range", col="grey", border="grey60") 

# Plot the probability density function for mean=5 and sd=2
x <- seq(-10, 20, by=.1) 
y <- dnorm(x, mean=5, sd=2) 
plot(x, y, type="l")

# What is the height of the probability density function at point x=4
dnorm(4, mean=5, sd=2)

# Plot the cummulative distribution function mean=5 and sd=2
x <- seq(-10 , 20, by=.1) 
y <- pnorm(x, mean=5, sd=2) 
plot(x,y,type="l")

# Exponential distribution
# Generate 1000 exponentially distributed numbers with rate=3 and draw a histogram
y <- rexp(1000, rate = 3)
hist(y, nclass=50, freq=FALSE, main="Exponential distribution", xlab="Range", col="grey", border="grey60") 

# Uniform distribution
#  Generate n-1000 numbers between 1 and 100, and draw a histogram
y <- sample(1:100, 1000, replace=TRUE)
hist(y, nclass=100, freq=FALSE, pch = 1, main="Uniform distribution", xlab="Range", col="grey", border="grey60") 

### 2. Random walk
# Generate a random walk with probability of moving in a positive direction p=0.75 
p <- 0.75 # probability of moving in a positive direction
n <- 100 # number of points of the grid/steps
t <- seq(1, n, by=1) # initialize time points
x <- numeric(n) # initialize a vector for direction at each step
u <- runif(n) # uniformly distributed numbers to determine direction of the move for n steps
for (i in 2:n) {
  if (u[i] < p)
    x[i] <- x[i-1]+1 # Move in a positive direction
  else
    x[i] <- x[i-1]-1 # Move in a negative direction	
}
plot(t, x, type="l", main="Random walk", xlab="Steps", ylab="Direction")

### 3. Brownian motion & Brownian bridge
#  Generate the Wiener process with number of points on a grid n=100 and length of an interval T=1.
n <- 100 # number of points of the grid 
T <- 1 # length of the interval 
delta <- T/n # time increment 
t <- seq(0,T, length=n+1) 
W <- c(0, cumsum(sqrt(delta) * rnorm(n))) 
plot( t, W, type="l", xlab="Steps", ylab="Direction")

## realize over 1000 simulations
Wsim <- matrix(NA, nrow=1000, ncol=n+1)
for(i in 1:1000){
  Wsim[i,] <- c(0, cumsum(sqrt(delta) * rnorm(n))) 
}
hist(Wsim[,101]) # check if normally distributed

#  Generate the Brownian Bridge based on Wiener process realisation above
B <- W - (t/T)*W[n+1] 
plot( t, B, type="l", xlab="Steps", ylab="Direction")
lines(t, rep(0, n+1), col="gray50", lty=2)
#  check if the process starts and finishes with the same value
B[1]==B[n+1]

###############################################################################################################
### Rates and probabilities
## 1) Probabilities: 30% of the population gets infected for S=10,000
S=100000; P=0.3
P*S # new infecteds calculated deterministically 
rbinom(1, S, P) # new infecteds drawn form Binomial distribution
rpois(1, S*P) # new infecteds drawn form Poisson distribution

## 2) Rates: infection lasts for three weeks
duration = 3
sigma = 1/3
I = 10
sigma*I # new recovereds every week (in the 1st week)

SIRV.model=function(t, x, params){
  S=x[1]
  I=x[2]
  R=x[3]
  V=x[4]
  
  with(as.list(params),{
    N = S + I + V 
    
    ## infection
    dS = beta*S*I - vc*S 
    dI = beta*S*I - sigma*I 
    dR = sigma*I 
    dV = vc*S 
    
    ## output
    res = c(dS, dI, dR, dV)
    list(res)
    
  })
}

require(deSolve)
parameters <- c(sigma=sigma, beta=0, vc=0)
time <- seq(0, 2, 1) # 2 weeks
start <- c(S=S, I=10, R=0, V=0)
as.data.frame(lsoda(start, time, SIRV.model, parameters))

## 3) Conversion
p=0.2
dt=1
x=-log(1-p)/dt
p.new=1-exp(-x*dt)
p; p.new

## 4) Combine and compare: 50% of population gets vaccinated in a year (dt) 
Pvacc=0.5 # probability
Rvacc=-log(1-Pvacc)/53 # rate (weekly for a year)

parameters <- c(sigma=0, beta=0, vc=Rvacc)
time <- seq(0, 53, 1) # simulate weekly
start <- c(S=S, I=0, R=0, V=0)
out <- as.data.frame(lsoda(start, time, SIRV.model, parameters))

# check the difference
S*Pvacc
rbinom(1, S, Pvacc)
out$V[dim(out)[1]]

