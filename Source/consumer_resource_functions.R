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
  with(b, exp( b0 + bz * arrtemp(z) + bm * log(M) ) )
}

ratefs <- function(bdd,re, ...){
  ratef(b=bdd[re,], ...)
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

f <- function(R,C,a,h,psi=0,omega=1,p=1,Q=0){
  C * omega * a * R / (1 + a*h*(p*R + (1-p)*Q + psi*C))
}
  # - effective handling time can be increased by:
  #   1. already-parasitised prey (discrete-time models):
  #     Q = E, where E are consumer eggs
  #   2. omega is ratio of predator interference time : prey handling time
  # assumes different resource species have
  #   same handling times, nutritional values, assimilation rates
  #   Koen-Alonso 1.12, 1.21 (Royama)

  # gamma for feeding juveniles (de Roos)

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

# functions separated by timescale:
# (allow one variable to adjust while holding others constant; repeat)
# - structure (assuming not feeding-induced)
# - resource size
# - climate: 
#   1. NLA function instead of zt
#   (climate-averaging of rates outside of ode, then run in constant env)
#   2. Allow to reach equilibrium densities before changing z 
#   (runsteady for whole ode)

d_chain <- function(y,b,Y,...){
  with(b, {
    dy <- y # convenient way to assign same names as y
    fy <- alpha * f(R=y[-Y],C=y[-1],a,h,psi,...)
    xy <- x(y[-1],mu) - c(fy[-1],0)
    dy[1] <-  g(y[1],v,k) - fy[1] / alpha
    dy[-1] <- fy - xy
    list(dy=dy,fy=fy,xy=xy)
  })
}

# structure types:
# - feeding conduit + migration
# - selection
# - proportion (selection with no response)
# feeding types:
# - one-chain
# - two-chain

d_bridge_base <- function(y1s,y2s,m,u){
  m * (y2s - u * y1s)
}
  # same equation as dilution rate
  # u = odds of y1 at equilibrium
  # feeding-mature (one-way diffusion): u = 0
  # select: m = (dyp1 - dyp2); u = -1
  # fraction: m and u set independently
  # DD: u = u2/y2s (or u2*y2s for reverse?)
  # Timescale separation - produces Hassell?

  # bridge rate = inflow to left-hand chain
  #   = immigration from y2 minus emigration from y1
  # when m is high, *can be subsumed into C or R equations*

  # feeding delays are in *new* biomass allocation
  # switch delays are in *old* biomass allocation

zbarf <- function(t,tau,zmu,zsig,zl){
  tau * zmu + (zsig * ( cos(2*pi * t/zl) - cos(2*pi * (t+tau)/zl) ) ) / (2*pi)
}
  # average temperature between t and t + tau
  # tau multiplier because want summed mortality [=exp(-mu*tau)]
  # https://www.youtube.com/watch?v=YF7Ii5dMYIo

# dede does not include methods to deal with delays that are smaller than the
# stepsize, although in some cases it may be possible to solve such models.
# For these reasons, it can only solve rather simple delay differential equations.

d_bridge_births <- function(parms){
  with(parms,{
    zsum_lag <- zbarf(t,tau,zmu,zsig,zl)

    bt1_lag <- c(lapply(bd,ratefs,M=M[Yr],z=zt_lag,re=Yr),bc)
    # is Yr the right body mass here?
    y1_lag <- lagvalue(t - tau, Ys)
    r1_lag <- lagvalue(t - tau, Yr)
    f1_lag <- with(bt1_lag, f(r1_lag,y1_lag,a,h,psi))
    
    bt2_lag <- c(lapply(bd,ratefs,M=M[YbB],z=zt_lag,re=YbB),bc)
    y2_lag <- lagvalue(t - tau, Ys1B)
    r2_lag <- lagvalue(t - tau, YrB)
    f2_lag <- with(bt2_lag, f(r2_lag,y2_lag,a,h,psi))
    
    dEy <- exp(x(f1_lag+f2_lag,zsum_lag,phi))
    return(dEy)
  })
}

d_bridge_diffuse <- function(Ys11,Ys22){
  with(parms,{
    zsum_lag <- zbarf(t,tau,zmu,zsig,zl)
    y1_lag <- lagvalue(t - tau, Ys11) # Sort out indices
    y2_lag <- lagvalue(t - tau, Ys22) 
    y1_left <- exp(x(y1_lag,zsum_lag,phi=1))
    y2_left <- exp(x(y2_lag,zsum_lag,phi))
    d_bridge_base(y1_left,y2_left,m,u) # m and u defined externally
  })
}

d_bridge_selection <- function(){
  with(parms,{
    zt_lag <- zt_cyclic(t-tau,zmu,zsig,zl)
    zsum_lag <- zbarf(t,tau,zmu,zsig,zl)
    
    bt1_lag <- c(lapply(bd,ratefs,M=M[Yr],z=zt_lag,re=Yr),bc)
      # is Yr the right body mass here?
    y1_lag <- lagvalue(t - tau, Ys)
    r1_lag <- lagvalue(t - tau, Yr)
    W1_lag <- with(bt1_lag, f(r1_lag,y1_lag,a,h,psi)/y1_lag)
  
    bt2_lag <- c(lapply(bd,ratefs,M=M[YbB],z=zt_lag,re=YbB),bc)
    y2_lag <- lagvalue(t - tau, Ys1B)
    r2_lag <- lagvalue(t - tau, YrB)
    W2_lag <- with(bt2_lag, f(r2_lag,y2_lag,a,h,psi)/y2_lag )
    
    deltaW <- W1_lag - W2_lag # assuming equal mortality rates
    
    d_bridge_base(y1s,y2s,m=deltaW,u=-1)
  })
}
  # currently only works for generalist consumers (not e.g. hiding)
  # do we want mortality weighting (i.e. only biomass alive at time of signal?)

d_web <- function(t,y,parms){
  
  with(parms, {
    
    zt <- zt_cyclic(t,zmu,zsig,zl)
      # t-specific parameters
      # deSolve requires that this be done separately for each t
    bdt <- as.data.frame(lapply(bd,ratef,M=M,z=zt))
    bt <- c(bdt,bc)
    
    if(structure==FALSE){
      
      dy <- d_chain(y,bt,Y)$dy
      
    }
    
    if(structure==TRUE){
      
      # q = food given to feeders instead of eggs
      # p = food given to feeders in first patch
      # dE = total food gathered for timestep
      
      y1 <- y[1:Yc]

      ### GATHER FOOD (one- and two-chain)
      
      if(nchain==1){
        
        dy1 <- d_chain(y1,bt,Y)
          # could modify this to include juvenile-only feeding
        dyE <- d1$fy[Yr]
        
      }
      
      if(nchain==2){
        
        y2 <- y[(Yc+1):(Yc*2)]
        
        bt1 <- bdt[1:Yb,]
        bt2 <- bdt[(Yb+1):(Yb*2),]
        
        pt1 <- pt2 <- rep(1,Yb) # or should this be Yr?
          # "eat your own" works for juvenile feeding too
        Q1  <- Q2  <- rep(0,Yb)

        p <- y1[Ys] / (y1[Ys] + y2[Ys])
        
        if(generalist==TRUE){
          
          pt1[Ys] <- p
          pt2[Ys] <- 1 - p
            # can be fraction eating each resource or fraction in each patch
          
          Q1[Yr] <- y2[Yr]
          Q2[Yr] <- y1[Yr]
            # Q refers to resource in *other* chain
        
        }
        
        d1 <- d_chain(y1,bt1,Y,p=pt1,Q=Q1)
        d2 <- d_chain(y2,bt2,Y,p=pt2,Q=Q2,omega=omega)

        dyE <- d1$fy[Yr] + d2$fy[Yr]

      }

      ### DEPOSIT EGGS (one- and two-chain)
      
      if(birthlag==TRUE){
        
        yE <- y[Y] # eggs = last position in chain

        if(tau==0){
          dEy <- d_bridge_base(y1s,y2s,m,u=0)
        }
        if(tau>0){
          dEy <- d_bridge_births(parms)
        }
          # maybe put in its own function?
        
        dE <- dyE - x(yE,bts$mu,phi) - dEy
          # egg mortality = adult mortality scaled by phi
          # using death rate from first chain
        
      }
      
      ### MIGRATE
      
      if(migrate==TRUE){
        # NEED THIS TO WORK FOR BOTH TWO-CHAIN AND ONE-CHAIN MODELS
        # (change in indexing)
        # add tau_migrate term here later
        
        d21 <- d_bridge_wrapper(y1s,y2s,stuff)
        
      }
      
      ### NET GROWTH
      
      if(nchain==1){
        d1$dy <- dEy - d1$xy[Yr1]
      }
      
      if(nchain==2){
        if(juvfeed==FALSE) q <- p
        if(juvfeed==TRUE)  q <- 0
        # need to override earlier p so that all food makes juveniles
        
        d1$dy[Ys1] <- q     * dEy - d1$xy[Yr1] + dy21
        d2$dy[Ys2] <- (1-q) * dEy - d2$xy[Yr2] - dy21
        # y2 can be eggs
        # need to distinguish feeding and mortality
        # q = fraction of resource allocated to consumer type 2
        # = 1 for feeding-controlled and stage-structured
        # = A/(A+B) for generalists OR fraction 
      }
      
      if(birthlag==FALSE) dya <- dy
      if(birthlag==TRUE)  dya <- c(dy,dE)
      
    } # end structured models
  
    return( list(dya) )

  })
    
}

iparmf <- function(y,sparms){
  
  with(sparms, {
      
    Yc <- chainlength
    
    if(structure==FALSE){
      
      iparms <- list(Y=Y)
      
    } 
    
    if(structure==TRUE){
    
      Ya <- Y * nchain

      if(birthlag==FALSE) Y <- Ya
      if(birthlag==TRUE)  Y <- Ya + 1
        # ***CAN BE USED TO LOCATE OFFSPRING VARIABLE***
  
      if(length(nstart)==1) y0 <- rep(nstart,Y)
      if(length(nstart)>1)  y0 <- nstart
      
      if(slevel=="consumer") Ys <- Yc
      if(slevel=="resource") Ys <- Yc - 1
        # position of structure

      Yr <- Ys - 1 # resource level for focal (structured) species
      Yb <- Y  - 1
  
      # YbB   <- Yb + Yr1
      # Ys1B <- Y + Ys1
      # YrB  <- Y + Yr1

      Ybseq <- rep(1:Yb, nchain)
        # no need for birth lag (egg) parms as selected automatically 
        #   from focal species in *first* chain

      M <- 10 ^ Ybseq
        # body masses 2 orders of magnitude apart, starting at 1g

      bd <- bdselect(bhat,Ybseq)
      
      if(birthlag==TRUE | juvfeed==TRUE) q <- 0 # ?
      
      iparms <- list(Y=Y,Ys1=Ys1,Ys2=Ys2,Yr1=Yr1,Yr2=Yr2,
                     Yb=Yb,
                     YbB=YbB,Ys1B=Ys1B,YrB=YrB)
        
      
    }
    
    names(y0) <- # "P"

    return(iparms)

  })
  
}

# Trial runs --------------------------------------------------------------

### Inputs (parms)
zmu <- 0 
zsig <- 0
zl <- 24*7

zparms <- list(zmu=zmu,zsig=zsig,zl=zl)

chainlength <- 2
  # for resource structure, food chain length must be at least 3
  #   (nutrients assumed inactive, so can't change state)
nstart <- 1
  # single number or vector of length Y
  # eggs (storage structure) always start at 0
structure <- TRUE
slevel <- c("consumer","resource")
birthlag <- TRUE
  # if FALSE, births are allocated directly to adults (via q)
lagtype <- c("none","random","fixed","discrete")
adaptive <- FALSE
  # if TRUE, add fitness difference multiplier to m
  # migration rate and direction determined by whichever increases fitness
nchain <- 1
generalist <- FALSE
  # if FALSE, is mixed specialist
juvfeed <- FALSE
   # if true, juveniles feed (so all births re-routed to P)

if(generalist==TRUE & adaptive==TRUE & lagtype!="none"){
  # m applies to offspring production; set generalist switch rate automatically
}

sparms <- list(
  structure = structure,
  slevel = slevel,
  twochain = twochain,
  bridgetype = bridgetype,
  generalist = generalist
)

y0 <- y <- c(R1=1,C1=1,C2=1,C2s=0) # rep(c(R1=1,C1=1,C2=1),2)
yl <- length(y)
t <- 2

bc <- c(
  v = 1,     # max flow rate = k grams per m^2 per hour
  k = 10,    # 10g per m^2
  psi = 0,   # interference:handling time ratio
  phi = 0,   # relative death rate of eggs
  omega = 0, # relative feeding rate of eggs
  u = 0,
  m = 1,     # rate of population structure adjustment
  q = 0.5,   # probability of individual being in state 1 at equilibrium
             # (only used if bridgetype = "constant")
  tau = 1    # continuous lag in bridge function
)
  # phi and omega could instead by controlled by body masses
  #   (in this case, phi can be fraction of adult body mass)
# bhat <- readRDS("Output/rate_parameters_marginal_23May2018.rds")

if(nchain==2){
  omega_new <- rep(1,iparms$Yb)
  omega_new[iparms$Yr] <- bc$omega # Yr because 
  bc$omega <- omega_new
}
bhat <- readRDS("Output/rate_parameters_simulated_27May2018.rds")

bdselect <- function(bhat,bpos){
  lapply(bhat,function(x) x[bpos,])
}

M <- c(M1=0.01,M2=1)
# different prey need different temp responses
# must have length equal to number of consumer stages
# generalist must be listed *last*
  # can have zero climate responses lower in chains too (= no sensitivity)

# change change constant parameters (e.g. w) among species by including as
#   climate-sensitive parameter but then setting climate sensitivity to zero
 
iparms <- iparmf(y,sparms)

vparms <- list(bd=bd,M=M)
fparms <- c(bc,M,zparms,sparms,iparms)
parms <- c(vparms,fparms)
attach(parms)

tseq <- seq(0,24*60,length.out=100)
require(deSolve)
if(parms$tau==0) lvar <- ode(y=y0,times=tseq,func=d_web,parms=parms)
if(parms$tau>0) lvar <- dede(y=y0,times=tseq,func=d_web,parms=parms)

par(mfrow=c(1,1))
matplot(tseq,log10(lvar[,-1]),type="l")

popint <- function(y0,tseq,parms){
  # if(nrow(bhat[])!=length(M)) stop("wrong masses or params")
  # ! generalist + "resource"
}

### TODO
# - discrete-time integration (u=0, so that zero backflow)
# - general timescale separation, 
#   i.e. function to calculate each set of equilibria (P, R, C)
# - analytically-simplified f functions? (timescale separation for states)
# - include checks and warnings
# - test d_bridge with real params to see if need ratio of dy/(dy-c)

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
