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

# Flux rates --------------------------------------------------------------

g <- function(R,m,r,k){
  (m + r*R) * (1 - R/k)
}

f <- function(R,C,a,h){
  R * C * a / (1 + a*h*R)
}

d <- function(C,x){
  x * C
}

# g1 <- function(R,m,r){
#   m + r*R
# }
# # births-only resource model
# 
# g0 <- function(R,m,r,k){
#   (m + r*R) * (R/k)
# }
# # deaths-only resource model

f0 <- function(R,C,a,h,Rhat){
  R * C * a / (1 + a*h*Rhat)
}
# psi=1 -> consumer "wastes" time on already-parasitised prey {Rhat-R} 

# Continuous dynamics -----------------------------------------------------

dRC_cont <- function(y,m,r,k,a,h,x,alpha){
  R <- y[1]
  C <- y[2]
  dR <- g(R,m,r,k) - f(R,C,a,h)
  dC <- alpha * f(R,C,a,h) - d(C,x)
  list(c(dR=dR,dC=dC))
}

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

rC <- Vectorize(
  function(zt,R,C,parms){
    parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
    with(parmst, dRC_cont(c(R,C),m,r,k,a,h,x,alpha=0.85)[[1]]["dC"]/C)
  },
  vectorize.args=c("zt","R","C")
)

# Lagged dynamics ---------------------------------------------------------

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

# Discrete dynamics -------------------------------------------------------

dRC_disc <- function(y,m,r,k,a,h,x,alpha,omega,kappa){
  R <- y[1]
  C <- y[2]
  E <- y[3]
  Ft <- f0(R,C,a,h,R+kappa*E)
  dR <- g1(R,m,r,k) - Ft
  dC <- - d(C,x)
  dE <- alpha * Ft - omega * d(E,x) 
  list(c(dR=dR,dC=dC,dE=dE))
}
  # consumer eggs die at same rate as adult consumers
  # egg mortality can be zero, but consumer mortality can't - 
  #   otherwise doesn't reduce to continuous model at very short time interval

dRCt_disc <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
  with(parmst, dRC_disc(y,m,r,k,a,h,x,alpha=0.85,
                        omega=parms$omega,
                        kappa=parms$kappa
                        ) )
}

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

dR_cont <- function(R,C,m,r,k,a,h,x){
  g(R,m,r,k) - f(R,C,a,h)
}

dC_cont <- function(R,C,m,r,k,a,h,x,alpha){
  alpha * f(R,C,a,h) - d(C,x)
}

rC_separate <- Vectorize(
  function(zt,C,parms){
    parmst <- with(parms, as.list( arrrate(zt,e0,e1) ) )
    try1 <- try(
      root1 <- with(parmst, uniroot(dR_cont,
                       C=C,m=m,r=r,k=k,a=a,h=h,x=x,
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
    with(parmst, dC_cont(Rstar,C,m,r,k,a,h,x,alpha=0.85)/C)
  },
  vectorize.args=c("zt","C")
)

### Discrete

# try1 <- try(
#   root1 <- uniroot(dR,C=Cs,p=p,a=a,h=h,u=q*u1,v=v1,lower=10^-10,upper=10^10) 
# )
# if(class(try1)=="try-error") K1s <- NA
# else K1s <- root1$root

