###################################################
# Numerical simulators for predator-prey dynamics #
###################################################

arrhenius <- function(Tt,E0,E1,T0=293.15,k=8.617*10^-5){
  E0 * exp( E1*(Tt-T0) / (k*Tt*T0) )
}

# Rosenzweig-MacArthur ----------------------------------------------------

romac <- function(y,r,K,a,b,x,eps=0.85){
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

romac_dis <- function(t,y,parms=NULL){
  
  tpos <- which(tseq>=t)[1]
  # use params from *element* number matching end of time interval
  
  r <- rseq[tpos]
  K <- Kseq[tpos]
  a <- aseq[tpos]
  b <- bseq[tpos]
  x <- xseq[tpos]
  
  romac(y,r,K,a,b,x)
  
}

romac_sin <- function(t,y,parms=NULL){
  
  Tt <- Tmu + 293.15 + Psd*sin(2*pi*t/lP) # these params defined externally
  
  r <- arrhenius(Tt,rE0,rE1)
  K <- arrhenius(Tt,KE0,KE1)
  a <- arrhenius(Tt,aE0,aE1)
  b <- arrhenius(Tt,bE0,bE1)
  x <- arrhenius(Tt,xE0,xE1)
  
  romac(y,r,K,a,b,x)
  
}

# Chemostat ---------------------------------------------------------------

chemo <- function(y,a,b,x,eps=0.85,i=5,e=1){
  R <- y[1]
  C <- y[2]
  dR <- i - R * ( e + a*C/(b+R) )
  dC <- C * ( eps*a*R/(b+R) - x )
  list(c(dR,dC))
}
# i = immigration rate
# e = emmigration rate

chemo_sin <- function(t,y,parms=NULL){
  
  Tt <- Tmu + 293.15 + Psd*sin(2*pi*t/lP) 
  # these params defined externally
  # deterministic -> keeps phase the same for different sims
  
  a <- arrhenius(Tt,aE0,aE1)
  b <- arrhenius(Tt,bE0,bE1)
  x <- arrhenius(Tt,xE0,xE1)
  
  chemo(y,a,b,x)
  
}

# Prey only ---------------------------------------------------------------

preyonly <- function(y,r,K){
  R <- y[1]
  dR <- R * ( r*(1-R/K) )
  list(dR)
}

preyonly_sin <- function(t,y,parms=NULL){
  Tt <- Tmu + 293.15 + Psd*sin(2*pi*t/lP) # these params defined externally
  r <- arrhenius(Tt,rE0,rE1)
  K <- arrhenius(Tt,KE0,KE1)
  preyonly(y,r,K)
}

