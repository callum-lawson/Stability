####################################################################################
# Include temperature fluctuations in pred-prey models from Fussmann et al. (2014) #
####################################################################################

# Functions ---------------------------------------------------------------

arrhenius <- function(Tt,E0,E1,T0=293.15,k=8.617*10^-5){
  E0 * exp( E1*(Tt-T0) / (k*Tt*T0) )
}

cr <- function(y,r,K,a,b,x,eps=0.85){
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

cr_dis <- function(t,y,parms=NULL){
  
  tpos <- which(tseq>=t)[1]
  # use params from *element* number matching end of time interval
  
  r <- rseq[tpos]
  K <- Kseq[tpos]
  a <- aseq[tpos]
  b <- bseq[tpos]
  x <- xseq[tpos]

  cr(y,r,K,a,b,x)
  
}

cr_sin <- function(t,y,parms=NULL){
  
  Tt <- Tmu + 293.15 + Psd*sin(2*pi*t/lP) # these params defined externally
  
  r <- arrhenius(Tt,rE0,rE1)
  K <- arrhenius(Tt,KE0,KE1)
  a <- arrhenius(Tt,aE0,aE1)
  b <- arrhenius(Tt,bE0,bE1)
  x <- arrhenius(Tt,xE0,xE1)

  cr(y,r,K,a,b,x)
  
}

# Species parameters ------------------------------------------------------

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

# Simulations with discrete temperatures ----------------------------------

nt <- 1000 # number of timesteps to calculate densities for
tmax <- 10^4 * 60^2  # maximum length of time in seconds
tseq <- seq(0,tmax,length.out=nt)
nP <- 1000 # number of perturbations over time series
lP <- nt/nP # perturbation length in nt units - has to be whole number

Tmu <- 20
Psd <- 5
Pvals <- rnorm(nP,mean=0,sd=Psd)
Pt <- rep(Pvals,each=lP)
Ttseq <- rep(Tmu + 293.15,nt) + Pt

rseq <- arrhenius(Ttseq,rE0,rE1)
Kseq <- arrhenius(Ttseq,KE0,KE1)
aseq <- arrhenius(Ttseq,aE0,aE1)
bseq <- arrhenius(Ttseq,bE0,bE1)
xseq <- arrhenius(Ttseq,xE0,xE1)
eps <- 0.85
  # will have to hold each of these constant for multiple timesteps
  # (otherwise numerical approximation of continuous time doesn't work)

R0 <- 1
C0 <- 1

### Numerical integration

library(deSolve)
ode1 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=cr_dis,parms=NULL)
matplot(tseq,log(ode1[,-1]),type="l")
abline(v=seq(0,tmax,length.out=nP+1),col="blue",lty=3)  
  # +1 accounts for t=0

### Euler's approximation at small timsteps

nf <- 10^3
Nt <- data.frame(R=rep(NA,nt),C=rep(NA,nt))

dt_i <- tseq[2]-tseq[1]
dt_j <- dt_i/nf

for(i in 1:nt){
  for(j in 1:nf){
    tcur <- (i-1)*dt_i + (j-1)*dt_j
    if(j==1){
      if(i==1){
        Ntcur <- c(R0,C0) + unlist(cr_dis(tcur,y=c(R0,C0)))*dt_j
      }
      if(i>1){
        Ntcur <- c(Nt[i-1,1],Nt[i-1,2]) + unlist(cr_dis(tcur,y=c(Nt[i-1,1],Nt[i-1,2])))*dt_j
      }
    }
    if(j>1){
      Ntcur <- c(Ntcur[1],Ntcur[2]) + unlist(cr_dis(tcur,y=c(Ntcur[1],Ntcur[2])))*dt_j
    }
  }
  Nt[i,] <- Ntcur
    # i indexes *end* of time interval
}
      
matplot(tseq,log(Nt),type="l")

# Simulations with sinusoidal temperatures --------------------------------

nt <- 1000 # number of timesteps to calculate densities for
tmax <- 10^4 * 60^2  # maximum length of time in seconds
tseq <- seq(0,tmax,length.out=nt)
lP <- tmax/10000 # wavelength in seconds

Tmu <- 20 # mean temperature in degrees C
Psd <- 5 # wave amplitude

R0 <- 1
C0 <- 1

### Numerical integration

Psd <- 0
ode1 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=cr_sin,parms=NULL)
Psd <- 5
ode2 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=cr_sin,parms=NULL)
matplot(tseq,log(ode1[,-1]),type="l",lty=1)
matplot(tseq,log(ode2[,-1]),type="l",add=T,lty=2)

matplot(tseq,log(cbind(ode1[,2],ode2[,2])),type="l",lty=1)
matplot(tseq,log(cbind(ode1[,3],ode2[,3])),type="l",lty=1)
  # variance -> both predator and prey fluctuate slightly more and 
  # reach lower densities
  # also, cycles slightly slower
  # increasing variance accentuates these effects

### Brief results summary

  # DISCRETE
  # - sometimes transient dynamics persist for ages; other times, snaps straight into 
  # new regime 
  # - the longer the intervals between temperature switches, the lesss transient 
  # dynamics matter
  # - negative perturbations > positive perturbations?
  # - perturbations every timestep breaks the simulation

  # CONTINUOUS
  # - wavelength -> 0 => constant fluctuations
  # (but how does this compare to mean in cons env?)

  # predator dynamics (rate of increase) slower than prey -> sharp vs blunt peaks

# TODO
# - check approximation for errors, and how this varies with n perturbations
# - chemostat resource
# - resource dynamics when generalist consumer unaffected by resource
# = a practical guide to Ecological Modelling Soetaert & Herman



