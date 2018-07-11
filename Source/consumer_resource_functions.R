### Numerical simulators for predator-prey dynamics ###

# Data processing ---------------------------------------------------------

sseqgen <- function(x,y){
  dmat <- outer(x,y,"-")
  apply(dmat,1,function(z) max(which(z>=0)))
}

# Temperature -------------------------------------------------------------

zt_cyclic <- function(t,zparms){
  with(zparms, {
    if(zsig==0) return( rep(zmu, length(t)) )
    if(zsig>0)  return( zmu + zsig * sin(2*pi * t/zl) )
  })
}

zbarf <- function(t,tau,zmu,zsig,zl){
  tau * zmu + (zsig * ( cos(2*pi * t/zl) - cos(2*pi * (t + tau)/zl) ) ) / (2*pi)
}
  # average *temperature* (not mortality) between t and t + tau
  # tau multiplier because want summed mortality [=exp(-mu*tau)]
  # https://www.youtube.com/watch?v=YF7Ii5dMYIo

arrtemp <- function(z,z0=20,T0=273.15,kB=8.6173303*10^-5){ 
  - 1/kB * ( 1/(T0+z) - 1/(T0+z0) )
} 
# z in C
# first part re-scales z so that 0 = 20°C = 293.15 K 
# second part re-scales intercept so that gives rates at 20°C

# Instantaneous rates -----------------------------------------------------

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
  as.data.frame(mapply(rate_abs,bd,bn,MoreArgs=list(z=z,M=M),SIMPLIFY=FALSE))
}
  # bd = list of b param dataframes for each rate (a,h,etc.)

btf <- function(t,bd,M,parms,ztype="whatever"){
  zt <- zt_cyclic(t,parms)
  bdt <- rate_df(bd,z=zt,M)
  return(bdt)
}

subd <- function(bdi,re){
  bdi[re,]
}

bdsubf <- function(bd,re){
  lapply(bd,subd,re=re)
}

# Non-linear averaging ----------------------------------------------------

rate_weight <- function(z,bi,bni,Mi,parms){
  with(parms, dnorm(z,zmu,zsig) * rate_abs(bi,bni,z,Mi) )
} 

rate_int <- function(bi,bni,Mi,parms){
  with(parms, {
    if(zsig==0){
      return( rate_abs(bi,bni,z=zmu,Mi) )
    }
    if(zsig>0){
      return(
        integrate(rate_weight,lower=-Inf,upper=Inf,
                  bi=bi,bni=bni,Mi=Mi,parms=parms)$value
        # if errors, try finite bounds
      )
    }
  })
}
  # works on one "row" of parameter dataframe

rate_int_d <- function(bdi,bni,parms){
  nbdi <- nrow(bdi)
  several <- nbdi > 1
  mapply(rate_int,
         bi=split(bdi,1:nbdi),
         Mi=parms$M,
         MoreArgs=list(bni=bni,parms=parms),
         USE.NAMES=FALSE,
         SIMPLIFY=several # prevents problems in two-species chains
         )
}
  # works over several rows of dataframe for a given parameter

rate_int_l <- function(bd,bn,parms){
  as.data.frame(
    mapply(rate_int_d,bdi=bd,bni=bn,MoreArgs=list(parms=parms))
  )
}
  # works over list of parameter dataframes

# Feeding rates -----------------------------------------------------------

g <- function(R,v,k){
  v * (k - R)
}
  # For closed nutrients, 
  # inflow proportional to decomposing material 
  # outflow is zero
  # additional parameter = starting total biomass (R0 + C0 ...)
  # model closed resource growth later: 1/alpha * x(C,mu)

f <- function(R,C,a,h,alpha,psi=0,omega=1,p=1,Q=0){
  C * alpha * omega * p * a * R / (1 + a*h*(p*R + (1-p)*Q + psi*C/p))
}
  # - effective handling time can be increased by:
  #   1. already-parasitised prey (discrete-time models):
  #     Q = E, where E are consumer eggs
  #   2. omega = predator interference time in units of prey handling time
  # assumes different resource species have
  #   same handling times, nutritional values, assimilation rates
  #   Koen-Alonso 1.12, 1.21 (Royama)
  # C = C1 OR C2 (not both)

  # gamma for feeding juveniles (de Roos)

x <- function(C,mu,phi=1){
  phi * mu * C
}
  # Multiple prey: individual prey attack rates with same, summed handling time
  # denominator (koen-Alonso)

# Chain derivatives -------------------------------------------------------

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
    fy <- f(R=y[-Y],C=y[-1],a,h,alpha,...)
    xy <- x(y[-1],mu) + c(fy[-1]/alpha[-1],0)
    dy[1] <-  g(y[1],v,k) - fy[1]/alpha[1]
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

# Migration rates ---------------------------------------------------------

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
    B_lag <- lagvalue(t - tau, Ya)
      # *** only works for egg structure ***
    
    if(mtype=="diffuse"){
      u <- u * (A_lag/B_lag) / (A/B)
        # overestimate of A density -> higher than expected movement to B
      migrate_base(A,B,m,u) # m and u defined externally
        # how to deal with lags here?
    }
    
    if(mtype %in% c("feeding","selective")){
      
      bdd <- bdsubf(bd,re=Yr)
      btA_lag <- btf(t=t-tau, bdd, M[Yr], parms)
      rA_lag <- lagvalue(t - tau, Yr)
      fA_lag <- with(btA_lag, f(rA_lag,A_lag,a,h,alpha,psi))
      
      if(nchain==2){
        btB_lag <- c(lapply(bd,ratefs,M=M[YbB],z=zt_lag,re=YbB),bc)
          # calculating zt_lag twice because nested within btf
        rB_lag <- lagvalue(t - tau, YrB)
        fB_lag <- with(btB_lag, f(rB_lag,B_lag,a,h,alpha,psi))
      }

        # ***need to add in p etc. to account for generalism***
      
      if(mtype=="feeding"){
        if(nchain==1) f_lag <- fA_lag
        if(nchain==2) f_lag <- fA_lag + fB_lag
        m <- f_lag / B # assuming zero mortality
        # zsum_lag <- zbarf(t,tau,zmu,zsig,zl)
          # MODIFY ZBARF TO INCLUDE MRATE CALCULATION?
        # m <- f_lag * exp(-mu*phi_E*zsum_lag))
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
    if(parms$discrete==FALSE){
      if(t < tau){
        delta <- 0
      }
      if(t >= tau){
        delta <- migrate_lag(A,B,t,tau,mtype,parms)
      }
    }
    if(parms$discrete==TRUE){
      delta <- with(parms, migrate_base(A,B=0,m,u,mtype))
      # eggs don't hatch, so no migration from B
    }
  }
  return(delta)
}

# Continuous derivatives --------------------------------------------------

d_web <- function(t,y,parms,hold=FALSE){
  
  with(parms, {

    if(is.null(bdt)){
      bdt <- btf(t,bd,M,parms) 
    }
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

        p <- y1[Ys] / (y1[Ys] + y2[Ys])
          # Yr = Ys - 1
        
        if(generalist==FALSE){
          
          d1 <- d_chain(y1,bt1,Yc,p=pt1,Q=Q2)
          d2 <- d_chain(y2,bt2,Yc,p=pt2,Q=Q1,omega=omega)
            # rely on bt2 to supply different phi
            # omega should be >0, otherwise use one-chain model
          
          }
        
        if(generalist==TRUE){
          
          pt1[Yr] <- p
          pt2[Yr] <- 1 - p
            # p = fraction of feeding on chain 1
            
          Q1[Yr] <- y1[Yr]
          Q2[Yr] <- y2[Yr]
            # Q = distraction by food from opposite chain
          
          d1 <- d_chain(y1,bt1,Yc,p=pt1,Q=Q2,omega=1)
          d2 <- d_chain(y2,bt2,Yc,p=pt2,Q=Q1,omega=omega)
            # accounts only for within-chain transfers
          
          bt1s <- c(bdt[Yr,], bc)
          bt2s <- c(bdt[Yb+Yr,],bc)
          
          fR2C1 <- with(bt1s,f(y2[Yr],y1[Ys],a,h,alpha,psi,omega=1,1-p,Q1[Yr]))
          fR1C2 <- with(bt2s,f(y1[Yr],y2[Ys],a,h,alpha,psi,omega,  p,  Q2[Yr]))
          
          d1$fy[Yr] <- d1$fy[Yr] + fR1C2
          d2$fy[Yr] <- d2$fy[Yr] + fR2C1
          d1$dy[Ys] <- d1$dy[Ys] + fR1C2
          d2$dy[Ys] <- d2$dy[Ys] + fR2C1
          d1$dy[Yr] <- d1$dy[Yr] - fR1C2/bt1s$alpha
          d2$dy[Yr] <- d2$dy[Yr] - fR2C1/bt2s$alpha
            # params treated as prey (not pred) characteristics
            # add effects of between-chain transfers
          
        }
        
        ft <- d1$fy[Yr] + d2$fy[Yr]

      } # finish food for two-chain

      ### DEPOSIT AND HATCH EGGS
      
      if(store==FALSE){
        dEy <- ft 
          # no storage -> eggs mirrored straight back
      }
      
      if(store==TRUE){
        yE <- y[Ya] # eggs = last position in chain
        if(storetype=="diffuse"){
          if(nchain==1) yss <- y1[Ys]
          if(nchain==2) yss <- y1[Ys] + y2[Ys]
        }
        if(storetype=="feeding"){
          yss <- 0  
        } 
        dEy <- migrate_all(A=yss,B=yE,
                            m=m_E,u=u_E,
                            t=t,tau=tau_E,
                            mtype=storetype,
                            parms=parms)
          # only need to use this option if one-chain model with storage
          # do we actually need u to be set to 0 here?
          # food conduit subtracted because handled separately to maturation
        dE <- with(bdt[Yr,], ft - dEy - x(yE,mu,phi_E) )
          # using mu param from y1
          # no selective behaviour (e.g. hiding) for one-chain models
      } # close egg store operations
      
      ### NET GROWTH - ONE-CHAIN
      
      if(nchain==1){
        if(store==TRUE) d1$dy[Ys] <- d1$dy[Ys] - d1$fy[Yr] + dEy
          # dEy = migration of E to ys
          # if no storage, no need for this operation
        dy <- d1$dy
      }
      
      ### MIGRATION AND NET GROWTH - TWO-CHAIN
      
      if(nchain==2){
        
        # q = food given to feeders instead of eggs
        
        if(movetype!="feeding"){
          q <- p # food allocated according to abundance
          if(movetype=="selective"){
            u_m <- d1$fy[Yr]/y1[Ys] - d2$fy[Yr]/y2[Ys]
          } 
            # for diffuse, m supplied as input parameter
          y1m <- y1[Ys]
        }  
        if(movetype=="feeding"){
          q <- 0
          y1m <- 0
        } 
        
        d21 <- migrate_all(A=y1m,B=y2[Ys],m=m_m,u=u_m,
                           t=t,tau=tau_m,mtype=movetype,
                           parms=parms)
        
        d1$dy[Ys] <-  q    * dEy - x(y1[Ys],bt1$mu[Yr],phi=1)     + d21
        d2$dy[Ys] <- (1-q) * dEy - x(y2[Ys],bt2$mu[Yr],phi=phi_m) - d21
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
  
    if(hold==TRUE)  dya[c(Yc,Yc2)] <- 0
      # adapt for eggs
    
    return( list(dya) )

  })
    
}

# Index variables:
#              | - Ya            
# Yc          | |- Yc2        
# |- Yb   Ys -| |- Ys2 / Yb2  
# |       Yr -| |- Yr2        
# |       Yl -| |            

# Setup -------------------------------------------------------------------

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
    Yl <-  Yr - 1
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
    
    t0 <- 0 # set automatically
    tseq <- seq(t0,tT,length.out=nt)
    
    list(structure=structure,
         Ya=Ya,Yc=Yc,Yc2=Yc2,Ys=Ys,Ys2=Ys2,
         Yr=Yr,Yr2=Yr2,Yl=Yl,Yb=Yb,Yb2=Yb2,
         Ycseq=Ycseq,Yc2seq=Yc2seq,
         Ybseq=Ybseq,Yb2seq=Yb2seq,
         bd=bd,M=M,y0=y0,
         t0=t0,tseq=tseq
         )
  
  })
  
}

bdselect <- function(bhat,bpos){
  lapply(bhat,function(x) x[bpos,])
}

# Discrete lags -----------------------------------------------------------

D_web <- function(parms){
  
  with(parms, {
    
    sstart <- seq(t0,tT,length.out=sS+1)
    sseq <- sseqgen(tseq,sstart)
    ns <- max(sseq)
    
    yd <- matrix(nr=nt,nc=Ya,dimnames=list(NULL,names(y0)))
    yd[1,] <- y0s <- y0
    
    for(s in 1:ns){
      
      if(s<ns)  savetimes <- c(sstart[s],tseq[sseq==s],sstart[s+1])
      if(s==ns) savetimes <- c(sstart[s],tseq[sseq==s])
      
      yds <- ode(y=y0s,times=savetimes,func=d_web,parms=parms)[,-1]
        # -1 removes t variable

      nts <- nrow(yds) # calculate beforehand?
      if(s<ns)  droprows <- c(1,nts)
      if(s==ns) droprows <- 1
      saverows <- which(sseq==s)
      yd[saverows,] <- yds[-droprows,]
      y0s <- yds[nts,]
      y0s[Ys] <- yds[nts,Ys] + yds[nts,Ya]
      y0s[Ya] <- 0
      # eggs become adults
      # existing resource and consumers carry over
      # ***adds eggs onto y1 only***
    }
    
    return(cbind(t=tseq,yd))
    
  })
  
}


# Integration wrapper -----------------------------------------------------

popint <- function(parms){
  require(deSolve)
  with(parms, {
    if(discrete==FALSE){
      if(tau_E==0) return( ode(y=y0,times=tseq,func=d_web,parms=parms) )
      if(tau_E>0)  return( dede(y=y0,times=tseq,func=d_web,parms=parms) )
        # *** currently only works for eggs ***
    } 
    if(discrete==TRUE)  return( D_web(parms) )
  })
}

# Timescale separation ----------------------------------------------------

# vector of top C densities
# resulting vector of equilibrium R densities
# C growth for each of those

rCf <- function(C,parms){
  require(rootSolve)
  with(parms, {
    y0[Yc] <- C
    Rstar <- steady(y=y0,
                    parms=parms,
                    fun=d_web,
                    times=c(0,Inf),
                    method="runsteady",
                    hold=TRUE
    )$y
    return( unlist(d_web(t=0, y=Rstar, parms))[Yc] )
  })
}

rCfv <- Vectorize(rCf,vectorize.args=c("C"))

# Carrying capacities -----------------------------------------------------

Cstarf <- function(parms){
  require(rootSolve)
  with(parms, {
    steady(y=y0,
           parms=parms,
           fun=d_web,
           times=c(0,Inf),
           method="runsteady",
           hold=FALSE
           )$y[Yc]
  })
}

# Cstarfv <- Vectorize(rCf,vectorize.args=c("zmu","zsd"))
  # needs to have zmu as argument


# Population growth - Discrete --------------------------------------------

RCf <- function(C,parms){
  # *no* timescale separation
  parmsD <- parms
  parmsD$sS <- 1
  parmsD$nt <- 2 # very inefficient - no need to do this every time
  parmsD$tseq <- with(parms, c(t0,tT))
  parmsD$y0[parms$Yc] <- C
  RC <- popint(parmsD)[2,parmsD$Yc + 1]
  return(RC)
}
  # *no* timescale separation 
  # but could build in quasi-separation by first setting the resource to its
  #   equilibrium, and then allowing dynamics to proceed normally
  # - would have to calculate changing R* with declining C during interval

RCfv <- Vectorize(RCf,vectorize.args=c("C"))

# Additions:
# - discrete: run for a single season using C~Rstar relationships
#   (or should we just run without timescale sparation? How does timescale
#   separation help here?)
# - structure: may need to write total C in terms of pC + (1-p)C

### TODO
# - general timescale separation, 
#   i.e. function to calculate each set of equilibria (P, R, C)
# - analytically-simplified f functions? (timescale separation for states)
# - include checks and warnings
