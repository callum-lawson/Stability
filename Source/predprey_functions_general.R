### Numerical simulators for predator-prey dynamics ###

# TODO:
# - implement my combined consumer-resource dynamics model
# - generalise for all temperature and population model types
# - two resource species
# - discrete births/deaths
# - Gaussian-process temperatures

arrtemp <- function(Tt,T0=293.15,k=8.6173303*10^-5){ # intercept = 20C!
  (Tt-T0) / (k*Tt*T0)
} 

arrrate <- function(Tt,E0,E1){
  E0 * exp(E1 * arrtemp(Tt))
}

Tt_cyclic <- function(t,Tmu,Tsd,Tperiod,T0=273.15){
  ifelse( Tsd==0, T0 + Tmu, T0 + Tmu + Tsd * sin(2*pi*t / Tperiod) )
}

g <- function(R,u,r,K){
  (u + r*R) * (1 - R/K)
}

f <- function(R,C,a,h){
  R * C * a / (1 + a*h*R)
}

d <- function(C,x){
  C * x
}

dCR <- function(y,u,r,K,a,h,x,alpha){
  R <- y[1]
  C <- y[2]
  dR <- g(R,u,r,K) - f(R,C,a,h)
  dC <- alpha * f(R,C,a,h) - d(C,x)
  list(c(dR=dR,dC=dC))
}

dCRt_cyclic <- function(t,y,parms){
  Tt <- with(parms, Tt_cyclic(t,Tmu,Tsd,Tperiod) )
  parmst <- with(parms, as.list( arrrate(Tt,E0,E1) ) )
  with(parmst, dCR(y,u,r,K,a,h,x,alpha=0.85) )
}

# R = prey biomass (g/m^2)
# C = predator biomass (g/m^2)
# u = immigration rate of resource
# r = prey intrinsic rate of increase (low density, absence of predator)
# K = prey carrying capacity (absence of predators)
# a = rate at which predator encounters prey
# h = handling time
# eps = number of predators produced for each prey eaten
# x = consumption rate required for predator to sustain itself