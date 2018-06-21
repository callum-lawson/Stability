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

rate_lp <- function(b,z,M){
  with(b, b0 + bz * arrtemp(z) + bm * log(M))
}

rate_abs <- function(bd,bn,z,M){
  lp <- rate_lp(bd,z,M)
  if(bn!="alpha") return( exp(lp) )
  if(bn=="alpha") return( plogis(lp) )
}

rate_df <- function(bd,z,M){
  bn <- names(bd)
  as.data.frame(mapply(rate_abs,bd,bn,MoreArgs=list(z=z,M=M)))
}
  # bd = list of b param dataframes for each rate (a,h,etc.)

rate_dfs <- function(bdd,re, ...){
  rate_df(b=bdd[re,], ...)
}

zbarf <- function(t,tau,zmu,zsig,zl){
  tau * zmu + (zsig * ( cos(2*pi * t/zl) - cos(2*pi * (t + tau)/zl) ) ) / (2*pi)
}
  # average temperature between t and t + tau
  # tau multiplier because want summed mortality [=exp(-mu*tau)]
  # https://www.youtube.com/watch?v=YF7Ii5dMYIo

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
  #   2. omega = predator interference time in units of prey handling time
  # assumes different resource species have
  #   same handling times, nutritional values, assimilation rates
  #   Koen-Alonso 1.12, 1.21 (Royama)

  # gamma for feeding juveniles (de Roos)

x <- function(C,mu,phi=1){
  phi * mu * C
}
  # Multiple prey: individual prey attack rates with same, summed handling time
  # denominator (koen-Alonso)

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
    xy <- x(y[-1],mu) + c(fy[-1],0)
    dy[1] <-  g(y[1],v,k) - fy[1] / alpha[1]
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

migrate_base <- function(A,B,m,u,mtype){
  # m always supplied externally (=maturation rate or W diff)
  # mtype=="diffuse" -> u pre-defined
  if(mtype=="feeding"){
    A <- 0; u <- 1 # or m <- 1
  }
  if(mtype=="selective"){
    # u indicates fitness; migration proceeds at constant speed m
    if(u>0)  return( + m * B )
    if(u<0)  return( - m * A )
    if(u==0) return( 0 )
  }
  if(mtype!="selective")
  return( m * (B - u * A) )
}
  # same equation as dilution rate
  # u = odds of y1 at equilibrium
  # bridge rate = inflow to left-hand chain
  #   = immigration from y2 minus emigration from y1
  # when m is high, can be subsumed into C or R equations

# dede does not include methods to deal with delays that are smaller than the
# stepsize, although in some cases it may be possible to solve such models.
# For these reasons, it can only solve rather simple delay differential equations.

migrate_lag <- function(A,B,t,tau,mtype,parms){
  with(parms,{
    
    A_lag <- lagvalue(t - tau, Ys)
    B_lag <- lagvalue(t - tau, Ys1B)
    
    if(mtype=="diffuse"){
      u <- u * (A_lag/B_lag) / (A/B)
        # overestimate of A density -> higher than expected movement to B
      migrate_base(A,B,m,u) # m and u defined externally
        # how to deal with lags here?
    }
    
    if(mtype %in% c("feeding","selective")){
      
      btA_lag <- c(lapply(bd,ratefs,M=M[Yr],z=zt_lag,re=Yr),bc)
      # is Yr the right body mass here?
      rA_lag <- lagvalue(t - tau, Yr)
      fA_lag <- with(btA_lag, f(rA_lag,A_lag,a,h,psi))
      
      if(nchain==2){
        btB_lag <- c(lapply(bd,ratefs,M=M[YbB],z=zt_lag,re=YbB),bc)
        rB_lag <- lagvalue(t - tau, YrB)
        fB_lag <- with(btB_lag, f(rB_lag,B_lag,a,h,psi))
      }

      ### *need to add in p etc. to account for generalism* ###
      
      if(mtype=="feeding"){
        if(nchain==1) f_lag <- fA_lag
        if(nchain==2) f_lag <- fA_lag + fB_lag
        zsum_lag <- zbarf(t,tau,zmu,zsig,zl)
        m <- exp(x(f_lag,zsum_lag,phi))
      }
      
      if(mtype=="selective"){
        WA_lag <- fA_lag/A_lag
        WB_lag <- fB_lag/B_lag
        m <- WA_lag - WB_lag # assuming equal mortality rates
      }
    }
    
    delta <- migrate_base(A,B,m,u,mtype)
    return(delta)
  })
}
  # do we want mortality weighting (i.e. only biomass alive at time of signal?)

migrate_all <- function(A,B,m,u,t,tau,mtype,parms){
  if(tau==0){
    delta <- with(parms, migrate_base(A,B,m,u,mtype))
  }
  if(tau>0){
      if(t<tau){
      delta <- 0
    }
    if(t>=tau){
      delta <- migrate_lag(A,B,t,tau,mtype,parms)
    }
  }
  return(delta)
}

d_web <- function(t,y,parms){
  
  with(parms, {
    
    zt <- zt_cyclic(t,zmu,zsig,zl)
      # t-specific parameters
      # deSolve requires that this be done separately for each t
    bdt <- rate_df(bd,M=M,z=zt)
    bt <- c(bdt,bc)
    
    if(structure==FALSE){
      
      dya <- d_chain(y,bt,Yc)$dy
      
    }
    
    if(structure==TRUE){
    
      y1 <- y[Ycseq]

      ### GATHER FOOD
      
      # ft = total food gathered for timestep
      
      if(nchain==1){
        
        d1 <- d_chain(y1,bt,Yc)
        ft <- d1$fy[Yr]
        
      }
      
      if(nchain==2){
        
        y2 <- y[Yc2seq]
        
        bt1 <- c(bdt[Ybseq,],bc)
        bt2 <- c(bdt[Yb2seq,],bc)
        
        pt1 <- pt2 <- rep(1,Yb)
          # pt = fraction of feeding on resource 1 (e.g. time spent in patch 1)
          # "eat your own" works for juvenile feeding too
        Q1  <- Q2  <- rep(0,Yb)
          # Qx = resource density in chain x
          # *Double-check whether this correctly accounts for R *and* C changes*
        
        p <- y1[Yr] / (y1[Yr] + y2[Yr])
          # Yr = Ys - 1
        
        if(generalist==TRUE){
          
          pt1[Yr] <- p
          pt2[Yr] <- 1 - p
          
          Q1[Yr-1] <- y1[Yr]
          Q2[Yr-1] <- y2[Yr]
        
        }
        
        y1ss <- y1
        y2ss <- y2
        yss  <- y1[Ys] + y2[Ys]
        y1ss[Ys] <- yss
        y2ss[Ys] <- yss
          # generalist -> whole pop eats both resources
          # yss variables are only used for calculating derivatives
          # (dy and xy terms used later to calculate actual densities)
        
        d1 <- d_chain(y1ss,bt1,Y,p=pt1,Q=Q2)
        d2 <- d_chain(y2ss,bt2,Y,p=pt2,Q=Q1,omega=omega)
          # Q = distraction by food from opposite chain
          # omega should be >0, otherwise use one-chain model

        ft <- d1$fy[Yr] + d2$fy[Yr]

      }

      ### DEPOSIT AND HATCH EGGS
      
      if(store==FALSE){
        dEy <- ft 
          # eggs hatch immediately - mirrored straight back
      }
      
      if(store==TRUE){
        bte <- bdt[Yr,]
        yE <- y[Ya] # eggs = last position in chain
        if(storetype=="diffuse"){
          dEy <- migrate_all(A=y1[Ys],B=yE,
                             m=m_E,u=u_E,
                             t=t,tau=tau_E,
                             mtype="diffuse",
                             parms=parms)
            # only need to use this option if one-chain model with storage
        }
        if(storetype=="feeding"){
          dEy <- migrate_all(A=0,B=yE,
                             m=m_E,u=u_E,
                             t=t,tau=tau_E,
                             mtype="feeding",
                             parms=parms) - ft 
            # do we actually need u to be set to 0 here?
            # food conduit subtracted because handled separately to maturation
        } 
        dE <- with(bte, -dEy - x(yE,mu,phi) )
          # no selective behaviour (e.g. hiding) for one-chain models
        
      } # close egg store operations
      
      ### NET GROWTH - ONE-CHAIN
      
      if(nchain==1){
        d1$dy[Ys] <- d1$dy[Ys] + dEy 
        # dEy = net transfer {E -> ys}
        dy <- d1$dy
      }
      
      ### MIGRATION AND NET GROWTH - TWO-CHAIN
      
      if(nchain==2){
        
        # q = food given to feeders instead of eggs
        
        if(movetype!="feeding"){
          q <- p # food allocated according to abundance
          if(movetype=="selective"){
            u <- d1$fy[Yr]/y1[Ys] - d1$fy[Yr]/y2[Ys]
          } 
            # for diffuse, m supplied as input parameter
          d21 <- migrate_all(A=y1[Ys],B=y2[Ys],m=m_m,u=u_m,
                             t=t,tau=tau_m,mtype=movetype,
                             parms=parms)
        }  
        
        if(movetype=="feeding"){
          q <- 0 # or food?
          d21 <- migrate_all(A=0,B=y2[Ys],m=m_m,u=u_m,
                             t=t,tau=tau_m,mtype=movetype,
                             parms=parms)
            # food transfer part dealt with below
        } 
        
        d1$dy[Ys] <-  q     * dEy - d1$xy[Yr] + d21
        d2$dy[Ys2] <- (1-q) * dEy - d2$xy[Yr] - d21 # account for phi here?
        # y2 can be eggs
        # need to distinguish feeding and mortality
        # q = fraction of resource allocated to consumer type 2
        # = 1 for feeding-controlled and stage-structured
        # = A/(A+B) for generalists OR fraction 
        # *only works if consumer - otherwise missing deaths due to predation*
        
        dy <- c(d1$dy,d2$dy)
      }
      
      if(store==FALSE) dya <- dy
      if(store==TRUE)  dya <- c(dy,dE)
      
    } # end structured models
  
    return( list(dya) )

  })
    
}

# Index variables:
#               | - Ya            
# Yc           | |- Yc2        
# |- Yb   Ys   | |- Ys2 / Yb2  
# |       Yr   | |- Yr2        
# |       |    | |            

iparmf <- function(bhat,sparms){
  
  with(sparms, {
    
    structure <- nchain > 1 | store==TRUE
    
    Yc <- chainlength
    
    Yc2 <- Yc * nchain
    
    # store==TRUE | tau_E>0?
    
    if(store==FALSE)  Ya <- Yc2
    if(store==TRUE)   Ya <- Yc2 + 1
    
    if(slevel=="consumer") Ys <- Yc
    if(slevel=="resource") Ys <- Yc - 1
      # position of structure
    
    Ys2 <- Yc + Ys
    Yr <-  Ys  - 1 # resource level for focal (structured) species
    Yr2 <- Ys2 - 1
    Yb <-  Yc  - 1
    Yb2 <- Yc2 - nchain 
      # basal resources don't have params
      # *should indexing be automated using nchain?* #
    
    Ycseq  <- 1:Yc
    Yc2seq <- (Yc+1):Yc2
    Ybseq  <- 1:Yb
      # egg parms will be selected from focal species in *first* chain
    Yb2seq <- (Yb+1):Yb2
    
    M <- 10 ^ (2*rep(Ybseq-1,nchain))
      # body masses 2 orders of magnitude apart, starting at 1g
      # same body masses used for same chain positions
    
    if(length(nstart)==1) y0 <- rep(nstart,Yc2)
    if(length(nstart)>1)  y0 <- nstart
    
    if(nchain==1){
      names(y0) <- c("R",paste0(rep("C",Yc-1),1:Yb))
      bd <- bdselect(bhat,Ybseq)
    }
    if(nchain==2){
      names(y0) <- paste0(rep(
        c("R",paste0(rep("C",Yc-1),1:Yb)),2),
        rep(c("A","B"),each=Yc)
      )
      bd <- bdselect(bhat,c(Ybseq,Yb2seq))
    }

    if(store==TRUE){
      y0 <- c(y0,E=0) # set eggs to zero
      names(y0)[Ya] <- paste0("E",names(y0)[Ys])
    } 
    
    list(structure=structure,
         Ya=Ya,Yc=Yc,Yc2=Yc2,Ys=Ys,Ys2=Ys2,
         Yr=Yr,Yr2=Yr2,Yb=Yb,Yb2=Yb2,
         Ycseq=Ycseq,Yc2seq=Yc2seq,
         Ybseq=Ybseq,Yb2seq=Yb2seq,
         bd=bd,M=M,y0=y0
         )
  
  })
  
}

bdselect <- function(bhat,bpos){
  lapply(bhat,function(x) x[bpos,])
}

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

### TODO
# - discrete-time integration (u=0, so that zero backflow)
# - general timescale separation, 
#   i.e. function to calculate each set of equilibria (P, R, C)
# - analytically-simplified f functions? (timescale separation for states)
# - include checks and warnings
# - test migrate with real params to see if need ratio of dy/(dy-c)
# - different prey need different temp responses

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
