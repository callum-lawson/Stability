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

# Parameter values --------------------------------------------------------

e0 <- c(
  m = 10^-6, # 10^-5,
  r = 8.715*10^-7,
  k = 5.623,
  a = 6*10^-7, # 3.181989*10^-9, # estimated from data
  h = 0.61, # 1685.586,     # estimated from data
  x = 2.689*10^-6
  # e = 0
)

e1 <- c(
  m = 0,
  r = 0, # from mortality rates # 0.84,
  k = 0, # -0.772,
  a = -0.03, # 0.5091663,   # estimated from data
  h = -0.19, # -0.4660012, # estimated from data
  x = 0.639
  # e = 0.639
)

# r units are per SECOND; pop more than triples every 24h

# Flux rates --------------------------------------------------------------

g <- function(R,m,r,k){
  (m + r*R) * (1 - R/k)
}

f <- function(R,C,a,h){
  R * C * a / (1 + a*h*R)
}

d <- function(C,x){
  C * x
}

g1 <- function(R,m,r){
  m + r*R
}
# births-only resource model

g0 <- function(R,m,r,k,Rhat){
  (m + r*R) * (Rhat/k)
}
# deaths-only resource model

f0 <- function(R,C,a,h,Rhat){
  R * C * a / (1 + a*h*Rhat)
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
  dR <- g(R,m,r,k) - f(R,C,a,h)
  dC <- alpha * f(R,C,a,h) - d(C,x)
  list(c(dR=dR,dC=dC))
}

dRC_disc <- function(y,m,r,k,a,h,x,alpha,omega,kappa,Rtype,Ctype){
  if(Rtype=="C" & Ctype=="C"){
    stop("resource or consumer should have discrete dynamics")
  }
  if(Rtype=="D") { phi <- 0 }
  if(Rtype=="C") { phi <- 1 }
  if(Ctype=="D") { psi <- 0 }
  if(Ctype=="C") { psi <- 1 }
  R <- y[1]
  C <- y[2]
  B <- y[3]
  E <- y[4]
  Ft <- f0(R,C,a,h,R+kappa*E)
  dR <- g1(R,(1-phi)*m,phi*r) - g0(R,(1-phi)*m,r,k,R) - Ft
  dC <- psi * alpha * Ft - d(C,x)
  dB <- (1-phi) * ( g1(R,m,r) - g0(B,m,r=0,k,B) ) # m=?
  dE <- (1-psi) * alpha * Ft - omega * d(E,x) 
  list(c(dR=dR,dC=dC,dB=dB,dE=dE))
}
  # phi included within dR because want constant migration effect regardless
  #   of whether resource births are continuous
  # capacity (K) of host adults and eggs is equivalent 
  # consumer eggs die at same rate as adult consumers

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
  with(parmst, dRC_disc(y,m,r,k,a,h,x,alpha=0.85,
                        omega=parms$omega,
                        kappa=parms$kappa,
                        Rtype=parms$Rtype,
                        Ctype=parms$Ctype
                        ) )
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
    dR <- g(R,m,r,k) - f(R,C,a,h)
    dC <- alpha * f_lag - d(C,x)
    return( list(c(dR=dR,dC=dC)) )
  })
  
}
  # never used in conjunction with dRC_disc because two alternative methods
  # of introducing time delays

# Population size integration ---------------------------------------------

DRCt_disc <- function(y0,tseq,sf,parms){

  tmax <- max(tseq)
  sstart <- seq(0,tmax,length.out=sf+1)
  sseq <- sseqgen(tseq,sstart)
  
  nt <- length(tseq)
  ns <- max(sseq)
  
  yd <- cbind(t=tseq,R=rep(NA,nt),C=rep(NA,nt))
  yd[1,"R"] <- y0[1]
  yd[1,"C"] <- y0[2]
  
  for(s in 1:ns){
    
    if(s<ns)  savetimes <- c(sstart[s],tseq[sseq==s],sstart[s+1])
    if(s==ns) savetimes <- c(sstart[s],tseq[sseq==s])
    
    y1 <- ode(y=c(y0,B=0,E=0),
              times=savetimes,
              func=dRCt_disc,
              parms=parms
              # atol=10^-16
              )
    
    nts <- nrow(y1)
    if(s<ns)  droprows <- c(1,nts)
    if(s==ns) droprows <- 1
    saverows <- which(sseq==s)
    yd[saverows,"R"] <- y1[-droprows,"R"]
    yd[saverows,"C"] <- y1[-droprows,"C"] 
    y0[1] <- y1[nts,"R"] + y1[nts,"B"]
    y0[2] <- y1[nts,"C"] + y1[nts,"E"] 
      # eggs become adults
      # existing resource and consumers carry over
    
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
    dR <- g(R,m,r,k) - f(R,C,a,h)
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