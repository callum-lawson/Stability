### Numerical simulators for predator-prey dynamics ###

# TODO:
# - generalise for all temperature and population model types
# - time lags
# - two resource species
# - discrete births/deaths
# - Gaussian-process temperatures

arrtemp <- function(zt,z0=293.15,k=8.6173303*10^-5){ # intercept = 20C!
  (zt-z0) / (k*zt*z0)
} 

arrrate <- function(zt,e0,e1){
  e0 * exp(e1 * arrtemp(zt))
}

zt_cyclic <- function(t,zmu,zsig,zl,z0=273.15){
  ifelse( zsig==0, z0 + zmu, z0 + zmu + zsig * sin(2*pi*t / zl) )
}

g <- function(R,m,r,K){
  (m + r*R) * (1 - R/K)
}

f <- function(R,C,a,h){
  R * C * a / (1 + a*h*R)
}

d <- function(C,x){
  C * x
}

g0 <- function(R,m,r,K,phi){
  (m + r*R) * (phi*1 - R/K)
}
# phi=0 -> resource dynamcis without births

f0 <- function(R,C,E,a,h,psi){
  R * C * a / (1 + a*h*(R+psi*E))
}
# psi=1 -> consumer "wastes" time on already-parasitised prey {Rhat-R} 

dRC_cont <- function(y,m,r,K,a,h,x,alpha){
  R <- y[1]
  C <- y[2]
  dR <- g(R,m,r,K) - f(R,C,a,h)
  dC <- alpha * f(R,C,a,h) - d(C,x)
  list(c(dR=dR,dC=dC))
}

dRC_disc <- function(y,m,r,K,a,h,x,alpha,Rtype){
  if(Rtype=="replenish") phi <- 1; psi <- 0
  if(Rtype=="remove")    phi <- 0; psi <- 1
  if(Rtype=="persist")   phi <- 0; psi <- 0
  R <- y[1]
  C <- y[2]
  E <- y[3]
  dR <- g0(R,m,r,K,phi) - f0(R,C,E,a,h,psi)
  dC <- 0 - d(C,x)
  dE <- alpha * f0(R,C,E,a,h,psi) - d(E,x)
  list(c(dR=dR,dC=dC,dE=dE))
}
  # assume that eggs die at same rate as adult consumers
  # but could also be assumed that die at same rate as unparasitised hosts

dRCt_cont <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
  with(parmst, dRC_cont(y,m,r,K,a,h,x,alpha=0.85) )
}

dRCt_disc <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
  with(parmst, dRC_disc(y,m,r,K,a,h,x,alpha=0.85,Rtype=parms$Rtype) )
}

DRCt_disc <- function(tseq,sseq,sstart,y0,parms){

  nt <- length(tseq)
  ns <- max(sseq)
  
  yd <- cbind(t=tseq,R=rep(NA,nt),C=rep(NA,nt))
  yd[1,"R"] <- y0[1]
  yd[1,"C"] <- y0[2]
  
  for(s in 1:(ns-1)){
    
    y1 <- ode(y=c(y0,E=0),
              times=c(sstart[s],tseq[sseq==s],sstart[s+1]),
              func=dRCt_disc,
              parms=parms
              )[,c("R","E")]
    
    nts <- nrow(y1)
    droprows <- c(1,nts)
    saverows <- which(sseq==s)
    yd[saverows,"R"] <- y1[-droprows,"R"] # + births
    yd[saverows,"C"] <- y1[-droprows,"E"] # adult consumers die, eggs become adults
    y0[1] <- y1[nts,"R"]
    y0[2] <- y1[nts,"E"]
    if(s==(ns-1)){
      yd[nt,c("R","C")] <- y0
    }
    
  }
  return(yd)
}

dRCt_delay <- function(t,y,parms,alpha=0.85,tau=60*60*24*7){
  
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
  
  R <- y[1]
  C <- y[2]
  
  if(t <= tau){
    f_lag <- 0
  }
  
  if(t > tau){
    zt_lag <- with(parms, zt_cyclic(t-tau,zmu,zsig,zl) )
    parmst_lag <- with(parms, as.list( arrrate(zt_lag,e0,e1) ) )
    RC_lag <- lagvalue(t - tau)
    f_lag <- with(parmst_lag, f(RC_lag[1],RC_lag[2],a,h))
  }
  
  with(parmst, {
    dR <- g(R,m,r,K) - f(R,C,a,h)
    dC <- alpha * f_lag - d(C,x)
    return( list(c(dR=dR,dC=dC)) )
  })
  
}
  # never used in conjunction with dRC_disc because two alternative methods
  # of introducing time delays

# Lag distribution - disfunctional ----------------------------------------

dRCt_siglag <- function(t,y,parms,alpha=0.85,taumu=60*60*24*7,tausig=0.001){
  
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
  
  R <- y[1]
  C <- y[2]
  
  if(tausig==0 & t>taumu){
    zt_lag <- with(parms, zt_cyclic(t-taumu,zmu,zsig,zl) )
    parmst_lag <- with(parms, as.list( arrrate(zt_lag,e0,e1) ) )
    RC_lag <- lagvalue(t - taumu)
    f_lag <- with(parmst_lag, f(RC_lag[1],RC_lag[2],a,h))
  }
  
  if(tausig==0 & t<=taumu){
    f_lag <- 0
  }
  
  if(tausig>0){
    RC_lagint <- Vectorize(function(tau){
      zt_lag <- with(parms, zt_cyclic(t-tau,zmu,zsig,zl) )
      parmst_lag <- with(parms, as.list( arrrate(zt_lag,e0,e1) ) )
      if(t>tau){
        RC_lag_tau <- lagvalue(t - tau)
        return(
          with(parmst_lag, 
               dlnorm(tau,log(taumu),tausig) * f(RC_lag_tau[1],RC_lag_tau[2],a,h)
               )
          )
      }
      else{
        return(0)
      }
    })
    f_lag <- integrate(RC_lagint,lower=taumu-10*tausig,upper=taumu+10*tausig)$value
  }
  
  with(parmst, {
    dR <- g(R,m,r,K) - f(R,C,a,h)
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