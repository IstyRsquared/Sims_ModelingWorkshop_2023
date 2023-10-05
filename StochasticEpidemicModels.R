## title: Simulation models for modeling methods workshop: Stochastic epidemic models in R
## author: Kristyna Rysava
## date: September 26, 2023

rm(list=ls())

### 1. Discrete time
## Simulate population dynamics for 100 years (weekly) 
## for transmission rate parameter beta=2.8, or and recovery rate sigma=1/3
Tmax=100*52
pop=1000
births.per.biweek.per.pop = 15/(52*1000) # ~ range 11-50 per 1,000 per year
deaths.per.biweek.per.pop = 1030/(52*100000) # 1,030 deaths per 100,000 per year

storeS <- rep(0,Tmax) 	
storeR <- rep(0,Tmax) 	

storeS[1] <- pop
storeR[1] <- 0

# loop
for (j in 2:Tmax) {
  deaths <- rbinom(1, storeS[j-1], deaths.per.biweek.per.pop)  # death
  births <- rbinom(1, storeS[j-1], births.per.biweek.per.pop)  # births
  storeS[j] <- storeS[j-1] - deaths + births
  storeR[j] <- storeR[j-1] + deaths
}

plot(1:Tmax, storeS, type="l", col = 'green', pch = 1, ylim=c(0, max(storeS)))
lines(1:Tmax, storeR, type="l", col = 'blue',pch = 3)

### 2. Continuous time 
## a) Transmission
## Simulate infection dynamics for as long as there are S and I present in the population
## for transmission rate parameter beta=2.8, or and recovery rate sigma=1/3
sigma <- 1/3
beta <- 2.8
S0 <- 999 
I0 <- 1
R0 <- 0

t <- 1
S <- S0
I <- I0
R <- R0
N <- S0 + I0 + R0

times <- t
SS <- S
II <- I
RR <- R

while ((I > 0) && ( S > 0)){
  r <- c(beta * S * I/N, sigma * I)
  rtotal <- sum(r)
  delta <- (-1/rtotal)*log(runif(1)) 
  t <- t + delta # Time to the next event
  times <- c(times, t)
  if (runif(1)*rtotal < r[1]) {
    S <- max(0, S - 1)  # Transmission event
    I <- I + 1
  }
  else {
    I <- max(0, I - 1)  # Recovery event
    R <- R + 1
  }
  SS <- c(SS, S)
  II <- c(II, I)
  RR <- c(RR, R)
}

plot(times, SS, type="l", col = 'green',pch = 1, ylim=c(0, 1000))
lines(times, II, type="l", col = 'red',pch = 2)
lines(times, RR, type="l", col = 'blue',pch = 3)
max(II)

# When did the epidemic die out?
paste("Duration of the outbreak: ", round(t,digits=2), "days")

## b) Demography
## Add birth and death into your model,and run over 100 days
# for birth rate lambda=1.5, and death rate mu=0.75 
sigma <- 1/3
beta <- 2.8 
lambda <- 0.005
mu <- 0.075
S0 <- 999 
I0 <- 1
R0<- 0

t <- 0
S <- S0
I <- I0
R <- R0
N <- S0 + I0 + R0

times <- t
SS <- S
II <- I
RR <- R
NN <- N

while (( N > 0) && ( S > 0) && (I > 0)){
  r <- c(lambda* N, beta * S * I/N, sigma * I, mu * S, mu * I, mu * R)
  rtotal <- sum(r)
  cr <- cumsum(r)
  delta <- (-1/rtotal)*log(runif(1))
  t <- t + delta # Time to the next event
  times <- c(times, t)
  P <- runif(1)*rtotal
  if (P < cr[1]) {
    S <-  S + 1 # Birth
  }
  else if (P < cr[2]) {
    S <- max(0, S - 1) # Transmission
    I <- I + 1
  }
  else if (P < cr[3]) {
    I <- max(0, I - 1) # Recovery
    R <- R + 1
  }
  else if (P < cr[4]) {
    S <- max(0, S - 1) # Death of S
  }
  else if (P < cr[5]) {
    I <- max(0, I - 1) # Death of I
  }
  else {
    R <- max(0, R - 1) # Death of R
  }
  
  SS <- c(SS, S)
  II <- c(II, I)
  RR <- c(RR, R)
  N <- S + I + R
  NN <- c(NN, N)
}

par(mfrow=c(2,2))
plot(times, SS, type="l", col = 'green',pch = 1)
plot(times, II, type="l", col = 'red',pch = 2)
plot(times, RR, type="l", col = 'blue',pch = 3)
plot(times, NN, type="l", col = 'grey',pch = 4)

# EX1. Run the model 100 times and plot distribution of the duration of outbreaks and calculate its mean 
outbreak.duration <-c()
for(run in 1:100){
  
}

par(mfrow=c(1,1))
hist(outbreak.duration, main = NULL) 
summary(outbreak.duration) 

# EX2. Modify stochastic the model to have no births, differential death rate 
# and at a rate alpha * R, recovered become susceptible again 
# alpha to 0.02, muS to 0.0001, muI to 0.1, muR to 0.01
# How many infected there were at the peak of the outbreak?
# Is the overall population growing or declining? Why?
muS <- 0.0001
muI <- 0.1
muR <- 0.01
alpha <- 0.02
S0 <- 999 
I0 <- 1
R0<- 0

t <- 0
S <- S0
I <- I0
R <- R0
N <- S0 + I0 + R0

times <- t
SS <- S
II <- I
RR <- R
NN <- N

while (( N > 0) && ( S > 0) && (I > 0)){
 
}

par(mfrow=c(2,2))
plot(times, SS, type="l", col = 'green',pch = 1)
plot(times, II, type="l", col = 'red',pch = 2)
plot(times, RR, type="l", col = 'blue',pch = 3)
plot(times, NN, type="l", col = 'grey',pch = 4)
max(II) 

## c) Vaccination
# EX3. Modify the stochastic SIR model (no demography) to have vaccination with rate gamma=20 
gamma <- 20
S0 <- 999 
I0 <- 1
R0<- 0
V0 <-0  

t <- 0
S <- S0
I <- I0
R <- R0
N <- S0 + I0 + R0
V <- V0

times <- t
SS <- S
II <- I
RR <- R
VV <- V

while (( I > 0) && ( S > 0)){
  
}

par(mfrow=c(2,2))
plot(times, SS, type="l", col = 'green', pch = 1)
plot(times, II, type="l", col = 'red', pch = 2)
plot(times, RR, type="l", col = 'blue', pch = 3)
plot(times, VV, type="l", col = 'gray', pch = 4)

# How does vaccination change the disease dynamics?


