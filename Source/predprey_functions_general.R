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

g <- function(R,u,r,K){
  (u + r*R) * (1 - R/K)
}

f <- function(R,C,a,h){
  R * C * a / (1 + a*h*R)
}

d <- function(C,x){
  C * x
}

g0 <- function(R,u,r,K,phi){
  (u + r*R) * (phi*1 - R/K)
}
# phi=0 -> resource dynamcis without births

f0 <- function(R,C,E,a,h,psi){
  R * C * a / (1 + a*h*(R+psi*E))
}
# psi=1 -> consumer "wastes" time on already-parasitised prey {Rhat-R} 

dRC_cont <- function(y,u,r,K,a,h,x,alpha){
  R <- y[1]
  C <- y[2]
  dR <- g(R,u,r,K) - f(R,C,a,h)
  dC <- alpha * f(R,C,a,h) - d(C,x)
  list(c(dR=dR,dC=dC))
}

dRC_disc <- function(y,u,r,K,a,h,x,alpha,Rtype){
  if(Rtype=="replenish") phi <- 1; psi <- 0
  if(Rtype=="remove")    phi <- 0; psi <- 1
  if(Rtype=="persist")   phi <- 0; psi <- 0
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

dRCt_cont <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
  with(parmst, dRC_cont(y,u,r,K,a,h,x,alpha=0.85) )
}

dRCt_disc <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
  with(parmst, dRC_disc(y,u,r,K,a,h,x,alpha=0.85,Rtype=parms$Rtype) )
}

DRCt_disc <- function(tseq,tmax,y0,parms){
  
  smax <- 10 # number of seasons
  tsseq <- seq(0,tsmax,length.out=smax)
  tmax <- tsmax/smax
  
  mP <- 1000
  nP <- round(1000/smax)
  tmat <- sapply(1:smax, function(x) (x-1)*tmax + seq(0,tmax,length.out=nP))
    # densities calculated for approx mP time points in total,
    # including all transitions between seasons
  
  yd <- cbind(t=as.vector(tmat),R=rep(NA,mP),C=rep(NA,mP))
  yd[1,"R"] <- y0[1]
  yd[1,"C"] <- y0[2]
  
  for(s in 1:(smax-1)){
    
    y1 <- ode(y=c(yd[s,-1],E=0),
              times=tmat[,s],
              func=dRCt_disc,
              parms=parms
              )[,c("R","C")]
    
    yd[match(tmat[,s],tmat),"R"] <- y1[,"R"] # + births
    yd[match(tmat[,s],tmat),"C"] <- y1[,"C"] # adult consumers die, eggs become adults
    
  }

  return(yd[!is.na(yd[,"C"]),])
}

E0 <- 0
y <- c(R0,C0,E0)
parms <- list(zmu=zmu,zsig=zsig,zl=zl,e0=e0,e1=e1,Rtype="replenish")
TT <- 100
ode(y=c(R=R0,C=C0,E=E0),times=c(0,TT),func=dRCt_disc,parms=parms)


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
    dR <- g(R,u,r,K) - f(R,C,a,h)
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