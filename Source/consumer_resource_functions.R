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

e0 <- c(
  m = 10^-5, # very high value relative to consumer
  k = 10,
  a = 6*10^-7, # 3.181989*10^-9, # estimated from data
  h = 12*60^2, # 0.61, # 1685.586,     # estimated from data
  w = 12*60^2,
  mu = 0.1 * 2.689*10^-6, # 2.689*10^-6,
  alpha = 0.1 * 0.85,
  phi = 0.1 # relative death rate of eggs
)

e1 <- c(
  m = 0, # 0.639,
  k = 0, # -0.772,
  a = 0.5091663, # -0.03,   # estimated from data
  h = 0, # -1.9, # -0.19, # -0.4660012, # estimated from data
  w = 0,
  mu = 0, # 0.639
  alpha = 0,
  phi = 0
)

e2 <- c(
  m = 0,
  k = 0, 
  a = 0, 
  h = 0, 
  w = 0,
  mu = -1/4, # metabolic rate per unit mass (could also be -1/3)
  alpha = 0,
  phi = 0
)

# arratel()

# M <- 1
# 
# zparms <- list(zmu=zmu,zsig=zsig,zl=zl)
# eparms <- list(M=M,e0=e0,e1=e1,e2=e2)

arrtempK <- function(z,kB=8.6173303*10^-5){
  - 1/(kB*z)
}

arrtemp <- function(z,z0=293.15,kB=8.6173303*10^-5){ 
  - 1/kB * ( 1/(z+z0) - 1/z0 )
} 
  # first part re-scales temperature so that 0 = 20°C = 293.15 K 
  # second part re-scales intercept so that gives rates at 20°C

arrrate <- function(zt,M,e0,e1,e2){
  e0 * exp(e1 * arrtemp(zt)) * M ^ e2
}

arrratel <- function(zt,M,e0,e1,e2){
  as.list(arrate(zt,M,e0,e1,e2))
}

# Non-linear averaging ----------------------------------------------------

arrrate_z <- function(zt,M,zmu,zsig,e0,e1,e2){
  dnorm(zt,zmu,zsig) * arrrate(zt,M,e0,e1,e2)
} 

arrint <- function(zmu,zsig,e0,e1){
  if(zsig==0){
    return( arrrate(zmu,M,e0,e1,e2) )
  }
  if(zsig>0){
    return(
      integrate(arrrate_z,lower=-200,upper=200,
                zmu=zmu,zsig=zsig,
                M=M,e0=e0,e1=e1,e2=e2)$value
    )
  }
}

parmsvar <- function(e0,e1,zmu,zsig){
  mapply(arrint,e0=e0,e1=e1,MoreArgs=list(zmu=zmu,zsig=zsig))
}
  # applies arrint to multiple pars at once

# Flux rates --------------------------------------------------------------

g <- function(R,m,k){
  m * (k - R)
}
  # For closed nutrients, 
  # inflow proportional to decomposing material 
  # outflow is zero
  # additional parameter = starting total biomass (R0 + C0 ...)
  # model closed resource growth later: 1/alpha * d(C,mu)

f <- function(R,C,a,h,w,Q=C){
  R * C * a / (1 + a*h*R + a*w*Q)
}
  # - effective handling time can be increased by:
  #   1. already-parasitised prey (discrete-time models):
  #     Rhat = R + E, where E are consumer eggs
  #   2. omega is ratio of predator interference time : prey handling time

d <- function(C,mu){
  mu * C
}
  # Multiple prey: individual prey attack rates with same, summed handling time
  # denominator (koen-Alonso)

# Continuous dynamics -----------------------------------------------------

# return:
# - grouped derivatives
# - individual derivatives
# - grouped derivatives for specified time (temperature series)
# - individual per-capita derivatives (vectorised)
# - grouped derivatives with continuous lag
# - grouped derivatives with eggs (combine with above?)
# - grouped derivatives with eggs for specified time (temperature series)
# - grouped derivatives with eggs for single temperature
# - discrete time series combining eggs and adults at intervals
# - equilibrium resource size (two-species model)
# - equilibrium consumer size (three-species model
# - per-capita derivatives (discrete-time model)
# - grouped derivatives (three- and four-species models)
# - discrete time series (three-species model)

# others:
# - multiple resource or consumer species (with which differences?)
# - protected resources

# y <- 1 # vector of states
#   
# zt <- with(parms$z, zt_cyclic(t,zmu,zsig,zl) ) # add stochastic function

d_chain <- function(t=0,y,parms){
  Y <- length(y)
  fN <- dN <- vector(mode="numeric",length=Y)
  parmst <- vector(mode="list",length=Y)
  for(i in 2:Y){
    parmst[[i]] <- with(parms[[i]], arrratel(zt,M,e0,e1,e2))
    fN[i] <- with(parms[[i]], f(y[i-1],y[i],a,h,w) )
    dN[i] <- with(parms[[i]], alpha * fRC[i] - d(C,mu) )
  }
  dN[1] <- g(N[1],m,k) - fN[2]
  as.list(c(dR,dC)) # automatic naming?
}

d_bridge <- function(t=0,y,parms){
  # protection, eggs, generalist?
  # lag, random, DD, discrete
}

# chain, chain, link top

# dRt_cons <- function(t,y,parms){
#   R <- y[1]
#   with(parms, list( dR = g(R,m,k) -  f(R,C,a,h,w) ) )
# }
# 
# dCt_cons <- function(t,y,parms,alpha){
#   C <- y[1]
#   with(parms, alpha * f(R,C,a,h,w) - d(C,mu) )
# }

dRCt_cont <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,M,e0,e1,e2) ) )
  with(parmst, dRC_cont(y,m,k,a,h,w,mu,alpha) )
}

rR <- Vectorize(
  function(zt,R,parms){
    parmst <- with(parms, as.list( arrrate(zt,M,e0,e1,e2) ) )
    with(parmst, dRC_cont(c(R,0),m,k,a,h,w,mu,alpha)[[1]]["dR"]/R)
  },
  vectorize.args=c("zt","R")
)

rC <- Vectorize(
  function(zt,R,C,parms){
    parmst <- with(parms, as.list( arrrate(zt,M,e0,e1,e2) ) )
    with(parmst, dRC_cont(c(R,C),m,k,a,h,w,mu,alpha)[[1]]["dC"]/C)
  },
  vectorize.args=c("zt","R","C")
)

# Random structure --------------------------------------------------------

  # add "p" function later
  # F_wA = q/(2-q)*F_wJ

# Density-dependent lags --------------------------------------------------
  
  # modelled with additional state variable, or with p function?

# Fixed lags --------------------------------------------------------------

dRCt_delay <- function(t,y,parms,alpha){
  
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,M,e0,e1,e2) ) )
  
  R <- y[1]
  C <- y[2]
  
  if(t <= parms$tau){
    f_lag <- 0
  }
  
  if(t > parms$tau){
    zt_lag <- with(parms, zt_cyclic(t-tau,zmu,zsig,zl) )
    parmst_lag <- with(parms, as.list( arrrate(zt_lag,M,e0,e1,e2) ) )
    RC_lag <- lagvalue(t - parms$tau)
    f_lag <- with(parmst_lag, f(RC_lag[1],RC_lag[2],a,h,w))
  }
  
  with(parmst, {
    dR <- g(R,m,k) - f(R,C,a,h,w)
    dC <- alpha * f_lag - d(C,mu)
    return( list(c(dR=dR,dC=dC)) )
  })
  
}
# never used in conjunction with dRC_disc because two alternative methods
# of introducing time delays

# Discrete lags -----------------------------------------------------------

dRC_disc <- function(y,m,k,a,h,w,mu,alpha,phi){
  R <- y[1]
  C <- y[2]
  E <- y[3]
  fRC <- f(R,C,a,h,w,Q=E)
  dR <- g(R,m,k) - fRC
  dC <- - d(C,mu)
  dE <- alpha * fRC - d(E, phi * mu)
  list(c(dR=dR,dC=dC,dE=dE))
}
  # consumer eggs die at same rate as adult consumers
  # egg mortality can be zero, but consumer mortality can't - 
  #   otherwise doesn't reduce to continuous model at very short time interval

dRCt_disc <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,M,e0,e1,e2) ) )
  with(parmst, dRC_disc(y,m,k,a,h,w,mu,alpha,phi) )
}

dRCt_disc_cons <- function(t,y,parms){
  with(parms, dRC_disc(y,m,k,a,h,w,mu,alpha,phi))
}

DRCt_disc <- function(y0,tseq,sf,parms){
  # write this to accept tau parameter as input instead?

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
    
    y1 <- ode(y=c(y0,E=0),
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
    y0[1] <- y1[nts,"R"]
    y0[2] <- y1[nts,"C"] + y1[nts,"E"] 
      # eggs become adults
      # existing resource and consumers carry over
  }
  return(yd)
}

# Timescale separation ----------------------------------------------------

### Continuous

dR_cont <- function(R,C,m,k,a,h,w,mu){
  g(R,m,k) - f(R,C,a,h,w)
}

dC_cont <- function(R,C,m,k,a,h,w,mu,alpha){
  alpha * f(R,C,a,h,w) - d(C,mu)
}

# uniroto:
rC_separate <- Vectorize(
  function(zt,C,parms){
    parmst <- with(parms, as.list( arrrate(zt,M,e0,e1,e2) ) )
    try1 <- try(
      root1 <- with(parmst, uniroot(dR_cont,
                       C=C,m=m,k=k,a=a,h=h,w=w,mu=mu,
                       lower=10^-10,
                       upper=10^10
                       ))
      )
    if(class(try1)=="try-error"){
      Rstar <- NA
    }
    else{
      Rstar <- root1$root
    } 
    with(parmst, dC_cont(Rstar,C,m,k,a,h,w,mu)/C)
  },
  vectorize.args=c("zt","C")
)

# runsteady

# Equilibrium resource size -----------------------------------------------

### R*: Two-species model

Rstarcalc <- Vectorize(
  function(C,z,eparms){
    require(rootSolve)
    parms <- with(eparms, as.list(c( C=C, arrrate(z,M,e0,e1,e2) )) )
    # C is fixed, so enters as parameter instead of state variable
    steady(y=parms$k, # using resource k as starting value
           parms=parms,
           fun=dRt_cons, 
           times=c(0,Inf),
           method="runsteady"
    )$y
  }, 
  vectorize.args=c("C","z")
)

### C*: three-species model

dRC_cons3 <- function(t,y,parms){
  R <- y[1]
  C1 <- y[2]
  # C2 = constant
  fRC1 <- with(parms$parms1, f(R,C1,a,h,w) )
  fRC2 <- with(parms$parms2, f(C1,C2,a,h,w) )
  dR <- with(parms$parms1, g(R,m,k) - fRC1)
  dC1 <- with(parms$parms1, alpha * fRC1 - d(C1,mu) - fRC2 )
  list(c(dR=dR,dC1=dC1))
}

Cstarcalc <- Vectorize(
  function(C2,z,eparms){
    require(rootSolve)
    parms <- with(eparms, list( 
      parms1 = as.list( arrrate(z,M[1],e0,e1,e2) ),
      parms2 = as.list(c( arrrate(z,M[2],e0,e1,e2), C2=C2 ))
    ))
    # C2 is fixed, so enters as parameter instead of state variable
    stars <- steady(y=c(R=1,C1=1),
           parms=parms,
           fun=dRC_cons3, 
           times=c(0,Inf),
           method="runsteady"
    )$y
  }, 
  vectorize.args=c("C2","z")
)

# Growth rate curves ------------------------------------------------------

dCCdt <- Vectorize(
  function(R,C,z,eparms){
    parms <- with(eparms, as.list(c( R=R, arrrate(z,M,e0,e1,e2) )) )
    dCt_cons(t=0,y=C,parms=parms)/C
  }, 
  vectorize.args=c("R","C","z")
)

DCC <- Vectorize(
  function(R,C,z,eparms,sl){
    parms <- with(eparms, as.list( arrrate(z,M,e0,e1,e2) ) )
    N <- ode(y=c(R=R,C=C,E=0),
        times=c(0,sl),
        func=dRCt_disc_cons,
        parms=parms
    )
    N[2,"C"] + N[2,"E"]
  },
  vectorize.args=c("R","C","z")
)

# Three-species food chain ------------------------------------------------

### Continuous time

dRC3 <- function(y,parms1,parms2){
  R <- y[1]
  C1 <- y[2]  
  C2 <- y[3]
  fRC1 <- with(parms1, f(R,C1,a,h,w) )
  fRC2 <- with(parms2, f(C1,C2,a,h,w) ) # add mass-relative attack rates later
  dR <- with(parms1, g(R,m,k) - fRC1 )
  dC1 <- with(parms1, alpha * fRC1 - d(C1,mu) - fRC2 )
  dC2 <- with(parms2, alpha * fRC2 - d(C2,mu) )
  list(c(dR=dR,dC1=dC1,dC2=dC2))
}

dRCt3 <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parms1 <- with(parms, as.list( arrrate(zt,M[1],e0,e1,e2)) )
  parms2 <- with(parms, as.list( arrrate(zt,M[2],e0,e1,e2)) )
  dRC3(y,parms1,parms2)
}

### Discrete time

dRC_disc3 <- function(y,parms1,parms2){
  R <- y[1]
  C1 <- y[2]
  C2 <- y[3]
  E <- y[4]
  fRC1 <- with(parms1, f(R,C1,a,h,w) )
  fRC2 <- with(parms2, f(C1,C2,a,h,w) )
  dR <- with(parms1, g(R,m,k) - fRC1 )
  dC1 <- with(parms1,  alpha * fRC1 - d(C1,mu) - fRC2 )
  dC2 <- with(parms2, - d(C2,mu) )
  dE <- with(parms2, alpha * fRC2 - d(E, phi * mu) )
  list(c(dR=dR,dC1=dC1,dC2=dC2,dE=dE))
}

dRCt_disc3 <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parms1 <- with(parms, as.list( arrrate(zt,M[1],e0,e1,e2)) )
  parms2 <- with(parms, as.list( arrrate(zt,M[2],e0,e1,e2)) )
  dRC_disc3(y,parms1,parms2)
}

DRCt_disc3 <- function(y0,tseq,sf,parms){

  tmax <- max(tseq)
  sstart <- seq(0,tmax,length.out=sf+1)
  sseq <- sseqgen(tseq,sstart)
  
  nt <- length(tseq)
  ns <- max(sseq)
  
  yd <- cbind(t=tseq,R=rep(NA,nt),C1=rep(NA,nt),C2=rep(NA,nt))
  yd[1,"R"] <- y0[1]
  yd[1,"C1"] <- y0[2]
  yd[1,"C2"] <- y0[3]
  
  for(s in 1:ns){
    
    if(s<ns)  savetimes <- c(sstart[s],tseq[sseq==s],sstart[s+1])
    if(s==ns) savetimes <- c(sstart[s],tseq[sseq==s])
    
    y1 <- ode(y=c(y0,E=0),times=savetimes,func=dRCt_disc3,parms=parms)
    
    nts <- nrow(y1)
    if(s<ns)  droprows <- c(1,nts)
    if(s==ns) droprows <- 1
    saverows <- which(sseq==s)
    yd[saverows,"R"] <- y1[-droprows,"R"]
    yd[saverows,"C1"] <- y1[-droprows,"C1"] 
    yd[saverows,"C2"] <- y1[-droprows,"C2"] 
    
    y0[1] <- y1[nts,"R"]
    y0[2] <- y1[nts,"C1"]
    y0[3] <- y1[nts,"C2"] + y1[nts,"E"] 
  }
  return(yd)
}

# Four-species food chain -------------------------------------------------

dRC4 <- function(y,parms1,parms2,parms3){
  R <- y[1]
  C1 <- y[2]  
  C2 <- y[3]
  C3 <- y[4]
  fRC1 <- with(parms1, f(R,C1,a,h,w) )
  fRC2 <- with(parms2, f(C1,C2,a,h,w) ) 
  fRC3 <- with(parms3, f(C2,C3,a,h,w) ) 
  dR <- with(parms1, g(R,m,k) - fRC1 )
  dC1 <- with(parms1, alpha * fRC1 - d(C1,mu) - fRC2 )
  dC2 <- with(parms2, alpha * fRC2 - d(C2,mu) - fRC3 )
  dC3 <- with(parms3, alpha * fRC3 - d(C3,mu) )
  list(c(dR=dR,dC1=dC1,dC2=dC2,dC3=dC3))
}

dRCt4 <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parms1 <- with(parms, as.list( arrrate(zt,M[1],e0,e1,e2)) )
  parms2 <- with(parms, as.list( arrrate(zt,M[2],e0,e1,e2)) )
  parms3 <- with(parms, as.list( arrrate(zt,M[3],e0,e1,e2)) )
  dRC4(y,parms1,parms2,parms3)
}
