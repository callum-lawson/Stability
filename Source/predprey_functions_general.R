### Numerical simulators for predator-prey dynamics ###

# TODO:
# - generalise for all temperature and population model types
# - time lags
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

g0 <- function(R,u,r,K,phi){
  (u + r*R) * (phi*1 - R/K)
}
  # phi=0 -> resource dynamcis without births

f0 <- function(R,C,E,a,h,psi){
  R * C * a / (1 + a*h*(R+psi*E))
}
  # psi=1 -> consumer "wastes" time on already-parasitised prey {Rhat-R} 

d <- function(C,x){
  C * x
}

dCR_cont <- function(y,u,r,K,a,h,x,alpha){
  R <- y[1]
  C <- y[2]
  dR <- g(R,u,r,K) - f(R,C,a,h)
  dC <- alpha * f(R,C,a,h) - d(C,x)
  list(c(dR=dR,dC=dC))
}

dCR_disc <- function(y,u,r,K,a,h,x,alpha,Rtype){
  if(Rtype=="replenished") phi <- 1; psi <- 0
  if(Rtype=="removed")     phi <- 0; psi <- 1
  if(Rtype=="persist")     phi <- 0; psi <- 0
  R <- y[1]
  C <- y[2]
  E <- y[3]
  dR <- g0(R,u,r,K,phi) - f0(R,C,E,a,h,psi)
  dC <- 0 - d(C,x)
  dE <- alpha * f0(R,C,E,a,h,psi) - d(E,x)
  list(c(dR=dR,dC=dC,dE=dE))
}
  # assume that eggs die at same rate as adult consumers
  # but could also be assumed that die at same rate as unparasitised hosts

dCRt_delay <- function(t,y,parms,alpha=0.85,tau=60*60*24*7){
  
  Tt <- with(parms, Tt_cyclic(t,Tmu,Tsd,Tperiod) )
  parmst <- with(parms, as.list( arrrate(Tt,E0,E1) ) )
  
  R <- y[1]
  C <- y[2]
  
  if(t <= tau){
    f_lag <- 0
  }
  
  if(t > tau){
    Tt_lag <- with(parms, Tt_cyclic(t-tau,Tmu,Tsd,Tperiod) )
    parmst_lag <- with(parms, as.list( arrrate(Tt_lag,E0,E1) ) )
    CR_lag <- lagvalue(t - tau)
    f_lag <- with(parmst_lag, f(CR_lag[1],CR_lag[2],a,h))
  }
  
  with(parmst, {
    dR <- g(R,u,r,K) - f(R,C,a,h)
    dC <- alpha * f_lag - d(C,x)
    return( list(c(dR=dR,dC=dC)) )
  })
  
}
  # never used in conjunction with dCR_disc because two alternative methods
  # of introducing time delays

dCRt_cyclic <- function(t,y,parms){
  Tt <- with(parms, Tt_cyclic(t,Tmu,Tsd,Tperiod) )
  parmst <- with(parms, as.list( arrrate(Tt,E0,E1) ) )
  with(parmst, dCR(y,u,r,K,a,h,x,alpha=0.85) )
}

# Lag distribution - disfunctional ----------------------------------------

dCRt_siglag <- function(t,y,parms,alpha=0.85,taumu=60*60*24*7,tausig=0.001){
  
  Tt <- with(parms, Tt_cyclic(t,Tmu,Tsd,Tperiod) )
  parmst <- with(parms, as.list( arrrate(Tt,E0,E1) ) )
  
  R <- y[1]
  C <- y[2]
  
  if(tausig==0 & t>taumu){
    Tt_lag <- with(parms, Tt_cyclic(t-taumu,Tmu,Tsd,Tperiod) )
    parmst_lag <- with(parms, as.list( arrrate(Tt_lag,E0,E1) ) )
    CR_lag <- lagvalue(t - taumu)
    f_lag <- with(parmst_lag, f(CR_lag[1],CR_lag[2],a,h))
  }
  
  if(tausig==0 & t<=taumu){
    f_lag <- 0
  }
  
  if(tausig>0){
    CR_lagint <- Vectorize(function(tau){
      Tt_lag <- with(parms, Tt_cyclic(t-tau,Tmu,Tsd,Tperiod) )
      parmst_lag <- with(parms, as.list( arrrate(Tt_lag,E0,E1) ) )
      if(t>tau){
        CR_lag_tau <- lagvalue(t - tau)
        return(
          with(parmst_lag, 
               dlnorm(tau,log(taumu),tausig) * f(CR_lag_tau[1],CR_lag_tau[2],a,h)
               )
          )
      }
      else{
        return(0)
      }
    })
    f_lag <- integrate(CR_lagint,lower=taumu-10*tausig,upper=taumu+10*tausig)$value
  }
  
  with(parmst, {
    dR <- g(R,u,r,K) - f(R,C,a,h)
    dC <- alpha * f_lag - d(C,x)
    return( list(c(dR=dR,dC=dC)) )
  })
  
}
  # doesn't seem to work - probably a problem with telling "lagvalue" where to look

# R = prey biomass (g/m^2)
# C = predator biomass (g/m^2)
# E = egg or parasitised host biomass
# u = immigration rate of resource
# r = prey intrinsic rate of increase (low density, absence of predator)
# K = prey carrying capacity (absence of predators)
# a = rate at which predator encounters prey
# h = handling time
# eps = number of predators produced for each prey eaten
# x = consumption rate required for predator to sustain itself