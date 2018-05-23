### Numerical simulators for predator-prey dynamics ###

# Data processing ---------------------------------------------------------

sseqgen <- function(x,y){
  dmat <- outer(x,y,"-")
  apply(dmat,1,function(z) max(which(z>=0)))
}

# Temperature -------------------------------------------------------------

zt_cyclic <- function(t,zmu,zsig,zl){
  if(zsig==0) return( rep(zmu, length(t)) )
  if(zsig>0)  return( zmu + zsig * sin(2*pi * t/zl) )
}

# Basic parameters --------------------------------------------------------

arrtemp <- function(z,z0=20,T0=273.15,kB=8.6173303*10^-5){ 
  - 1/kB * ( 1/(T0+z) - 1/(T0+z0) )
} 
  # z in C
  # first part re-scales z so that 0 = 20°C = 293.15 K 
  # second part re-scales intercept so that gives rates at 20°C

ratef <- function(z,M,b){
  with(b, exp( b0 + bz * arrtemp(z) + bm * ln(M) ) )
}

# Flux rates --------------------------------------------------------------

g <- function(R,v,k){
  v * (k - R)
}
  # For closed nutrients, 
  # inflow proportional to decomposing material 
  # outflow is zero
  # additional parameter = starting total biomass (R0 + C0 ...)
  # model closed resource growth later: 1/alpha * x(C,mu)

f <- function(R,C,a,h,psi,p=1,Q=0){
  C * a * R / (1 + a*h*(p*R + (1-p)*Q + psi*C))
}
  # - effective handling time can be increased by:
  #   1. already-parasitised prey (discrete-time models):
  #     Q = E, where E are consumer eggs
  #   2. omega is ratio of predator interference time : prey handling time
  # assumes different resource species have
  #   same handling times, nutritional values, assimilation rates
  #   Koen-Alonso 1.12, 1.21 (Royama)

x <- function(C,mu,phi=1){
  phi * mu * C
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

# Derivative functions ----------------------------------------------------

d_chain <- function(y,b,Y,...){
  with(b, {
    dy <- y # convenient way to assign same names as y
    fy <- f(R=y[-Y],C=y[-1],a,h,psi,...)
    dy[1] <-  g(y[1],v,k) - fy[1]
    dy[-1] <- alpha * fy - x(y[-1],mu) - c(fy[-1],0)
    list(dy=dy,fy=fy)
  })
}

# test d_bridge with real params to see if need ratio of dy/(dy-c)

d_bridge <- function(y1,y2,f1,f2=NULL,bc,bridgetype="fraction"){
  with(bc, {
    if(bridgetype=="fraction"){
      u1 <- u * y2
      u2 <- 1 * y1
      # fraction active at equilibrium = w2/(w1+w2)
    }
    if(bridgetype=="waiting"){
      u1 <- u * y2
      u2 <- f1
      # bigger pop -> less food -> more individuals feeding (are non-eggs)
      #   -> lower fraction protected (= constant refuges)
      # generally maladaptive because just delays
    } 
    if(bridgetype=="hiding"){
      u1 <- f1
      u2 <- (1/u) * y1
      # density-dependent germination: higher R -> lower food -> more protected
    }
    if(bridgetype=="switch"){
      u1 <- f1/y1
      u2 <- f2/y2
      # could be reparameterised to include mortality rates
      # can apply to consumer, or resource moving between patches
      # switching should be according to per-capita fitness, so divide by y
    }
    return( m * (u1 - u2) )
  })
}
  # w = estimated fitness 
  # m = flow rate
  # overall idea: fraction can depend on R (via resource density), or not
  # indirectly depends on C (because -> lower R) 
  # ignore case in which unprotected don't eat (no evolutionary justification)
  # when f -> eggs, fs cancel and adults only grow through development
  # bridge rate describes inflow to left-hand chain
  #   = immigration from y2 minus emigration from y1
  # when m is high, *can be subsumed into C or R equations*

zbarf <- function(t,tau,zmu,zsig,zl){
  tau * zmu + (zsig * ( cos(2*pi * t/zl) - cos(2*pi * (t+tau)/zl) ) ) / (2*pi)
}
  # average temperature between t and t + tau
  # tau multiplier because want summed mortality [=exp(-mu*tau)]
  # https://www.youtube.com/watch?v=YF7Ii5dMYIo

d_web <- function(t,y,parms){
  
  with(parms, {
    
    zt <- with(zparms, zt_cyclic(t,zmu,zsig,zl))
      # t-specific parameters
      # deSolve requires that this be done separately for each t
    bt <- c(lapply(bhat,ratef,M=M,z=zt),bc)
    
    if(slevel=="none"){
      
      Y <- length(y) 
        # total chain length
      return( list(d_chain(y,bt,Y)$dy) )
      
    }
    
    if(slevel!="none"){
      
      if(twochain==FALSE) Y <- length(y) - 1 
      if(twochain==TRUE)  Y <- length(y) / 2 
        # total chain length
      
      if(slevel=="consumer") Y1 <- Y
      if(slevel=="resource") Y1 <- Y-1
      
      y1 <- y[1:Y1]
      y2 <- y[-(1:Y1)]
        # all states after y1 chain assigned to y2
      
      if(generalist==FALSE){
        
        d1 <- d_chain(y1,bt,Y)
        dy1 <- d1$dy
          # generalist requires altered feeding rates (see below)
        
      }

      if(twochain==FALSE){
        
        Y2 <- 1 # y2 consists of just one state variable
        YR <- Y1 - 1 # resource level / feeding rate vector
        dy2 <- with(bt, -x(y2,mu[YR],phi))
          # egg mortality = adult mortality scaled by phi
        
          
        if(bridgetype!="lag"){
          dy21 <- d_bridge(y1[Y1],y2[Y2],
                           d1$fy[YR],f2=NULL,
                           bc=bparms$bc,bridgetype
                           )
        }
        
        if(bridgetype=="lag"){
          # lag in {eggs -> adults} only (not vice-versa)
          if(t <= zparms$tau){
            u_lag <- 0
          }
          if(t > zparms$tau){
            # set dy to without food input
            zt_lag <- with(zparms, zt_cyclic(t-tau,zmu,zsig,zl) )
            zsum_lag <- with(zparms, zbarf(t,tau,zmu,zsig,zl) )
            bt_lag <- with(zparms, lapply(bhat,ratef,M=M,z=zt_lag))
            y_lag <- lagvalue(t - zparms$tau)
            f_lag <- with(bt_lag, f(y_lag[YR],y_lag[Y1],a,h,bparms$bc$psi))
            u_lag <- with(bparms$bc, exp(x(f_lag,zsum_lag,phi)) )
            # = food deposit × fraction surviving
          }
          dy21 <- d_bridge(y1[Y1],y2[Y2],
                           u_lag,f2=d1$fy[YR],
                           bc=bparms$bc,bridgetype="switch"
          )
        }

        # if single chain, then other lifestage doesn't feed
      }
      
      if(twochain==TRUE){
        
        Y2 <- Y1
        YR <- Y1 - 1 # resource level / feeding rate vector
        
        if(generalist==FALSE){
          
          d2 <- d_chain(y2,bt,Y2)
          dy2 <- d2$dy
          dy21 <- d_bridge(y1[Y1],y2[Y2],
                           d1$fy[YR],d2$fy[YR],
                           bc=bparms$bc,bridgetype
          )
          # assumes that consumption of R1 produces more R1 hunters
          # (same as switching being proportional to food intake?)
          
        }

        if(generalist==TRUE){
          
          pt1 <- pt2 <- rep(1,YR)
          pt1[YR] <- y1[Y1]/(y1[Y1]+y2[Y2])
          pt2[YR] <- 1 - pt1[YR]
          Q1 <- Q2 <- rep(0,YR)
          Q1[YR] <- y2[YR]
          Q2[YR] <- y1[YR]
            # Q refers to resource in *other* chain
          
          d1 <- d_chain(y1,bt,Y1,p=pt1,Q=Q1)
          dy1 <- d1$dy
          d2 <- d_chain(y2,bt,Y2,p=pt2,Q=Q2)
          dy2 <- d2$dy

          dy21 <- d_bridge(y1[Y1],y2[Y2],
                           d1$fy[YR],d2$fy[YR],
                           bc=bparms$bc,bridgetype
          )
          
        }
        
      }
        # -1 indexing because basal resource has no parameters
        # feeding rates given same body masses, so same index of Y1
      
      dy1[Y1] <- dy1[Y1] + dy21
      dy2[Y2] <- dy2[Y2] - dy21
      # y2 can be eggs
      # build in lag (instead of "store") version of this later
      
      return( list(c(dy1,dy2)) )
      
    } # end structured models
    
  })
    
}

# Trial runs --------------------------------------------------------------

### Inputs (parms)
zmu <- 0 
zsig <- 0
zl <- 24*7
tau <- 0    # time delay for lag model
zparms <- list(zmu=zmu,zsig=zsig,zl=zl,tau=tau)

slevel <- "consumer"
twochain <- TRUE
bridgetype <- "switch"
  # for resource structure, food chain length must be at least 3
  #   (nutrients assumed inactive, so can't change state)
generalist <- TRUE

bc <- list(
  v = 1,     # max flow rate = k grams per m^2 per hour
  k = 10,    # 10g per m^2
  psi = 0,   # relative handling time
  phi = 0,   # relative death rate of eggs
  m = 100,   # rate of population structure adjustment
  u = 1      # odds of individual being in state 1 at equilibrium
)
bhat <- readRDS("Output/rate_parameters_marginal_22May2018.rds")
M <- c(0.01,1)
  # different prey need different temp responses
  # must have length equal to number of consumer stages
  # generalist must be listed *last*
bparms <- list(bc=bc,bhat=bhat)

y0 <- y <- rep(c(R1=1,C1=1,C2=1),2)
t <- 0
parms <- list(zparms=zparms,
              bparms=bparms,
              M=M,
              slevel=slevel,
              twochain=twochain,
              bridgetype=bridgetype
)

tseq <- seq(0,24*60,length.out=100)
require(deSolve)
lvar <- ode(y=y0,times=tseq,func=d_web,parms=parms)
# lvar <- dede(y=y0,times=tseq,func=d_web,parms=parms)

par(mfrow=c(1,1))
matplot(tseq,log10(lvar[,-1]),type="l")

### TODO
# - different climate responses in different chains
# - discrete-time integration (u=0, so that zero backflow)
# - general timescale separation
# - analytically-simplified f functions? (timescale separation for states)
# - include checks and warnings

# 
# 
# d_bridge_lag <- 
#   
#   if(t > parms$tau){
#     # set dy to without food input
#     zt_lag <- with(zparms, zt_cyclic(t-tau,zmu,zsig,zl) )
#     zbar_lag <- with(zparms, zbarf(t,tau,zmu,zsig,zl) )
#     bt_lag <- with(zparms, lapply(bhat,ratef,M=M,z=zt_lag))
#     y_lag <- lagvalue(t - parms$tau)
#     f_lag <- with(bt_lag, f(y_lag[ypos-1],y_lag[ypos],a,h,psi))
#     dy[ypos] <- with(bparms, dy[ypos] +  exp(x(f_lag,zbar_lag*tau,phi)) )
#     # = food deposit × fraction surviving
#   }
# # if t < lag, do nothing
# return( m * (w2 * y1 * exp(x(1,zbar_lag*tau,phi=1)) 
#              - w1 * y2 * exp(x(1,zbar_lag*tau,phi))) )

# *still need true generalist*:
f(p*R1+(1-p)*R2,C,a,h,psi,Q=C)

  
# dRt_cons <- function(t,y,parms){
#   R <- y[1]
#   with(parms, list( dR = g(R,m,k) -  f(R,C,a,h,w) ) )
# }
# 
# dCt_cons <- function(t,y,parms,alpha){
#   C <- y[1]
#   with(parms, alpha * f(R,C,a,h,w) - x(C,mu) )
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
    dC <- alpha * f_lag - x(C,mu)
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
  dC <- - x(C,mu)
  dE <- alpha * fRC - x(E, phi * mu)
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
  alpha * f(R,C,a,h,w) - x(C,mu)
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
  dC1 <- with(parms$parms1, alpha * fRC1 - x(C1,mu) - fRC2 )
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
  dC1 <- with(parms1, alpha * fRC1 - x(C1,mu) - fRC2 )
  dC2 <- with(parms2, alpha * fRC2 - x(C2,mu) )
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
  dC1 <- with(parms1,  alpha * fRC1 - x(C1,mu) - fRC2 )
  dC2 <- with(parms2, - x(C2,mu) )
  dE <- with(parms2, alpha * fRC2 - x(E, phi * mu) )
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
  dC1 <- with(parms1, alpha * fRC1 - x(C1,mu) - fRC2 )
  dC2 <- with(parms2, alpha * fRC2 - x(C2,mu) - fRC3 )
  dC3 <- with(parms3, alpha * fRC3 - x(C3,mu) )
  list(c(dR=dR,dC1=dC1,dC2=dC2,dC3=dC3))
}

dRCt4 <- function(t,y,parms){
  zt <- with(parms, zt_cyclic(t,zmu,zsig,zl) )
  parms1 <- with(parms, as.list( arrrate(zt,M[1],e0,e1,e2)) )
  parms2 <- with(parms, as.list( arrrate(zt,M[2],e0,e1,e2)) )
  parms3 <- with(parms, as.list( arrrate(zt,M[3],e0,e1,e2)) )
  dRC4(y,parms1,parms2,parms3)
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
