####################################################################################
# Include temperature fluctuations in pred-prey models from Fussmann et al. (2014) #
####################################################################################

arrhenius <- function(Tt,E0,E1,T0=293.15,k=8.617*10^-5){
  E0 * exp( E1*(Tt-T0) / (k*Tt*T0) )
}

cr <- function(t,y,parms){
  
  tpos <- which(tseq>=t)[1]
  # use params from *element* number matching end of time interval
  
  r <- r[tpos]
  K <- K[tpos]
  a <- a[tpos]
  b <- b[tpos]
  x <- x[tpos]
  eps <- 0.85
  
  R <- y[1]
  C <- y[2]
  dR <- R * ( r*(1-R/K) - a*C/(b+R) )
  dC <- C * ( eps*a*R/(b+R) - x )
  
  list(c(dR,dC))
}
  # R = prey biomass (g/m^2)
  # C = predator biomass (g/m^2)
  # r = prey intrinsic rate of increase (low density, absence of predator)
  # K = prey carrying capacity (absence of predators)
  # a = rate at which predator encounters prey
  # b = prey density at which predator is feeding at 1/2 max capacity
  # eps = number of predators produced for each prey eaten
  # x = consumption rate required for predator to sustain itself
  # tau = step length

rE0 <- 8.715*10^-7
KE0 <- 5.623
aE0 <- 8.408*10^-6
bE0 <- 3.664
xE0 <- 2.689*10^-6

# rE1 <- 0.84
# KE1 <- -0.772
# aE1 <- 0.467
# bE1 <- -0.114
# xE1 <- 0.639
  # From species-level averages
  # change slopes later
  # r units are per SECOND; pop more than triples every 24h

rE1 <- 0.84
KE1 <- -0.508
aE1 <- 0.708
bE1 <- -0.678
xE1 <- 0.428
  # From Fig. S1

nt <- 100 # number of timesteps to calculate densities for
tmax <- 10^5 * 60^2  # maximum length of time in seconds
tseq <- seq(0,tmax,length.out=nt)
nP <- 75 # number of perturbations over time series
lP <- nt/nP # perturbation length in nt units

Tmu <- 10
Psd <- 10
Pvals <- rnorm(nP,mean=0,sd=Psd)
Pt <- rep(Pvals,each=lP)
Tt <- rep(Tmu + 293.15,nt) + Pt

r <- arrhenius(Tt,rE0,rE1)
K <- arrhenius(Tt,KE0,KE1)
a <- arrhenius(Tt,aE0,aE1)
b <- arrhenius(Tt,bE0,bE1)
x <- arrhenius(Tt,xE0,xE1)
eps <- 0.85
  # will have to hold each of these constant for multiple timesteps
  # (otherwise numerical approximation of continuous time doesn't work)

R0 <- 1
C0 <- 1

library(deSolve)
ode1 <- rk(y=c(R0=R0,C0=C0),times=tseq,func=cr,parms=NULL)
matplot(tseq,log(ode1[,-1]),type="l")
abline(v=seq(0,tmax,length.out=nP+1),col="blue",lty=3)  
  # +1 accounts for t=0

  # - sometimes transient dynamics persist for ages; other times, snaps straight into 
  # new regime 
  # - the longer the intervals between temperature switches, the lesss transient 
  # dynamics matter
  # - negative perturbations > positive perturbations?
  # - perturbations every timestep breaks the simulation

  # predator dynamics (rate of increase) slower than prey -> sharp vs blunt peaks

# TODO
# - check approximation for errors, and how this varies with n perturbations
# - chemostat resource
# - resource dynamics when generalist consumer unaffected by resource
# = a practical guide to Ecological Modelling Soetaert & Herman



