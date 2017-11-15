### Numerical simulators for predator-prey dynamics ###

# Data processing ---------------------------------------------------------

sseqgen <- function(x,y){
  dmat <- outer(x,y,"-")
  apply(dmat,1,function(z) max(which(z>=0)))
}

# Temperature -------------------------------------------------------------

zt_cyclic <- function(t,zmu,zsig,zl){
  if(zsig==0) return( rep(zmu, length(t)) )
  if(zsig>0)  return( zmu + zsig * sin(2*pi*t / zl) )
}

# Basic parameters --------------------------------------------------------

arrtemp <- function(zt,z0=293.15,kB=8.6173303*10^-5){ # intercept = 20C!
  zt / (kB*(zt+z0)*z0) # (zt-z0) / (k*zt*z0)
} 

arrrate <- function(zt,e0,e1){
  e0 * exp(e1 * arrtemp(zt))
}

# Flux rates --------------------------------------------------------------

g <- function(R,m,r,k,x,x0=2.689*10^-6){
  (m + r*R) * (1 - R/k) - (x-x0)*R
}
  # x0 is intercept for mortality function 

f <- function(R,C,a,h){
  R * C * a / (1 + a*h*R)
}

d <- function(C,x){
  C * x
}

g0 <- function(R,m,r,k,x,phi){
  (m + r*R) * (phi*1 - R/k) - x
}
# phi=0 -> resource dynamics without births

f0 <- function(R,C,E,a,h,psi){
  R * C * a / (1 + a*h*(R+psi*E))
}
# psi=1 -> consumer "wastes" time on already-parasitised prey {Rhat-R} 


# Non-linear averaging ----------------------------------------------------

arrrate_z <- function(zt,zmu,zsig,e0,e1){
  dnorm(zt,zmu,zsig) * arrrate(zt,e0,e1)
} 

arrint <- function(zmu,zsig,e0,e1){
  if(zsig==0){
    return( arrrate(zmu,e0,e1) )
  }
  if(zsig>0){
    return(
      integrate(arrrate_z,lower=-200,upper=200,zmu=zmu,zsig=zsig,e0=e0,e1=e1)$value
    )
  }
}

parmsvar <- function(e0,e1,zmu,zsig){
  mapply(arrint,e0=e0,e1=e1,MoreArgs=list(zmu=zmu,zsig=zsig))
}
# applies arrint to multiple pars at once


# Population size derivatives ---------------------------------------------

dRC_cont <- function(y,m,r,k,a,h,x,alpha){
  R <- y[1]
  C <- y[2]
  dR <- g(R,m,r,k,x) - f(R,C,a,h)
  dC <- alpha * f(R,C,a,h) - d(C,x)
  list(c(dR=dR,dC=dC))
}

dRC_disc <- function(y,m,r,k,a,h,x,alpha,Rtype){
  if(Rtype=="replenish") phi <- 1; psi <- 0
  if(Rtype=="remove")    phi <- 0; psi <- 1
  if(Rtype=="persist")   phi <- 0; psi <- 0
  R <- y[1]
  C <- y[2]
  E <- y[3]
  dR <- g0(R,m,r,k,x,phi) - f0(R,C,E,a,h,psi)
  dC <- 0 # - d(C,x)
  dE <- alpha * f0(R,C,E,a,h,psi) - d(E,x)
  list(c(dR=dR,dC=dC,dE=dE))
}
  # assume that eggs die at same rate as adult consumers
  # but could also be assumed that die at same rate as unparasitised hosts

dRCt_cont <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
  with(parmst, dRC_cont(y,m,r,k,a,h,x,alpha=0.85) )
}

rR <- Vectorize(
  function(zt,R,parms){
    parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
    with(parmst, dRC_cont(c(R,0),m,r,k,a,h,x,alpha=0.85)[[1]]["dR"]/R)
  },
  vectorize.args=c("zt","R")
)

dRCt_disc <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
  with(parmst, dRC_disc(y,m,r,k,a,h,x,alpha=0.85,Rtype=parms$Rtype) )
}

dRCt_delay <- function(t,y,parms,alpha=0.85){
  
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
  
  R <- y[1]
  C <- y[2]
  
  if(t <= parms$tau){
    f_lag <- 0
  }
  
  if(t > parms$tau){
    zt_lag <- with(parms, zt_cyclic(t-tau,zmu,zsig,zl) )
    parmst_lag <- with(parms, as.list( arrrate(zt_lag,e0,e1) ) )
    RC_lag <- lagvalue(t - parms$tau)
    f_lag <- with(parmst_lag, f(RC_lag[1],RC_lag[2],a,h))
  }
  
  with(parmst, {
    dR <- g(R,m,r,k,x) - f(R,C,a,h)
    dC <- alpha * f_lag - d(C,x)
    return( list(c(dR=dR,dC=dC)) )
  })
  
}
  # never used in conjunction with dRC_disc because two alternative methods
  # of introducing time delays

# Population size integration ---------------------------------------------

DRCt_disc <- function(y0,tseq,sseq,sstart,parms){

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
              )
    
    nts <- nrow(y1)
    droprows <- c(1,nts)
    saverows <- which(sseq==s)
    yd[saverows,"R"] <- y1[-droprows,"R"]
    yd[saverows,"C"] <- y1[-droprows,"C"] 
    y0[1] <- y1[nts,"R"] # + births
    y0[2] <- y1[nts,"E"] # adult consumers die, eggs become adults
    if(s==(ns-1)){
      yd[nt,c("R","C")] <- y0
    }
    
  }
  return(yd)
}

# DOES NOT WORK - Lag distribution model ----------------------------------

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
    dR <- g(R,m,r,k,x) - f(R,C,a,h)
    dC <- alpha * f_lag - d(C,x)
    return( list(c(dR=dR,dC=dC)) )
  })
  
}
  # doesn't seem to work - probably a problem with telling "lagvalue" where to look

# R = prey biomass (g/m^2)
# C = predator biomass (g/m^2)
# E = egg or parasitised host biomass
# m = immigration rate of resource
# r = prey intrinsic rate of increase (low density, absence of predator)
# K = prey carrying capacity (absence of predators)
# a = rate at which predator encounters prey
# h = handling time
# alpha = number of predators produced for each prey eaten
# x = consumer death rate in absence of prey