#################################################################################
# Phase-space plots when underlying predator or prey paramaters are fluctuating # 
#################################################################################

require(rootSolve)
require(phaseR)

source("Source/predprey_functions.R")

# Species parameters ------------------------------------------------------

rE0 <- 8.715*10^-7
KE0 <- 5.623
cE0 <- 8.408*10^-6
dE0 <- 3.664
xE0 <- 2.689*10^-6

# rE1 <- 0.84
# KE1 <- -0.772
# cE1 <- 0.467
# dE1 <- -0.114
# xE1 <- 0.639
  # From species-level averages
  # change slopes later
  # r units are per SECOND; pop more than triples every 24h

rE1 <- 0.84
KE1 <- -0.508
cE1 <- 0.708
dE1 <- -0.678
xE1 <- 0.428
# From Fig. S1

# Functions ---------------------------------------------------------------

dR_dt_f <- function(R,C,r,K,c,d){
  R * ( r*(1-R/K) - c*C/(d+R) )
  } 

dC_dt_f <- function(R,C,c,d,x,eps=0.85){
  C * ( eps*c*R/(d+R) - x )
} 

Rstar_f <- function(C,...){
  e <- try( 
    d <- uniroot(dR_dt_f, C=C, lower=10^-100, upper=10^100, ...), 
    silent = TRUE 
  ) 
  if(class(e)=="try-error") { 
    return(NA) 
  } 
  else{ 
    return(d$root)	
  } 
}

romac_phaser <- function(t, y, parameters){
  r <- parameters$r
  K <- parameters$K
  a <- parameters$a
  h <- parameters$h
  x <- parameters$x
  eps <- parameters$eps
  dy <- numeric(2)
  dy[1] <- y[1] * ( r*(1-y[1]/K) - a*y[2]/(1+a*h*y[1]) )
  dy[2] <- y[2] * ( eps*a*y[1]/(1+a*h*y[1]) - x )   
  return(list(dy))
}

deang_phaser <- function(t, y, parameters){
  r <- parameters$r
  K <- parameters$K
  a <- parameters$a
  h <- w <- parameters$h # special case where pred and prey time loss are equal
  x <- parameters$x
  eps <- parameters$eps
  dy <- numeric(2)
  dy[1] <- y[1] * ( r*(1-y[1]/K) - a*y[2]/(1+a*(w*y[2] + h*y[1])) )
  dy[2] <- y[2] * ( eps*a*y[1]/(1+a*(w*y[2] + h*y[1])) - x )   
  return(list(dy))
}
  # predators encounter each other at the same rate as prey (=a)

bazy_phaser <- function(t, y, parameters){
  r <- parameters$r
  K <- parameters$K
  a <- parameters$a
  h <- parameters$h
  x <- parameters$x
  E <- parameters$E
  eps <- parameters$eps
  dy <- numeric(2)
  dy[1] <- y[1] * ( r*(1-y[1]/K) - a*y[2]/(1+a*h*y[1]) )
  dy[2] <- y[2] * ( eps*a*y[1]/(1+a*h*y[1]) - x - E*y[2] )   
  return(list(dy))
}
  # E = strength of consumer density-dependence

regrow_phaser <- function(t, y, parameters){
  r <- parameters$r
  K <- parameters$K
  a <- parameters$a
  h <- parameters$h
  x <- parameters$x
  eps <- parameters$eps
  dy <- numeric(2)
  dy[1] <- r*(1-y[1]/K) - y[1] * a*y[2]/(1+a*h*y[1])
  dy[2] <- y[2] * ( eps*a*y[1]/(1+a*h*y[1]) - x )   
  return(list(dy))
}

chemo_phaser <- function(t, y, parameters){
  i <- parameters$i
  e <- parameters$e
  a <- parameters$a
  h <- parameters$h
  x <- parameters$x
  eps <- parameters$eps
  dy <- numeric(2)
  dy[1] <- i - y[1] * ( e + a*y[2]/(1+a*h*y[1]) )
  dy[2] <- y[2] * ( eps*a*y[1]/(1+a*h*y[1]) - x )   
  return(list(dy))
}

gomp_phaser <- function(t, y, parameters){
  r <- parameters$r
  K <- parameters$K
  a <- parameters$a
  h <- parameters$h
  x <- parameters$x
  eps <- parameters$eps
  dy <- numeric(2)
  if(log(y[1])/log(K)>-Inf){
    dy[1] <- y[1] * ( r*(1-log(y[1])/log(K)) - a*y[2]/(1+a*h*y[1]) )
  } 
  else{
    dy[1] <- 0
  }
  dy[2] <- y[2] * ( eps*a*y[1]/(1+a*h*y[1]) - x )   
  return(list(dy))
}

Tseq <- c(15,25) + 293.15
nT <- length(Tseq)
pd <- data.frame(
  T = Tseq,
  r = arrhenius(Tseq,rE0,rE1),
  K = arrhenius(Tseq,KE0,KE1),
  c = arrhenius(Tseq,cE0,cE1),
  d = arrhenius(Tseq,dE0,dE1),
  x = arrhenius(Tseq,xE0,xE1),
  eps = 0.85
  )

pd$h <- with(pd, 1/c) # handling time
pd$a <- with(pd, 1/(d*h)) # attack rate
  # Turchin 2003 p82

pd$i <- 0.0001 # resource inflow
pd$e <- 0.0001 # resource decay
# Reynolds & Brassil 2013

pd0 <- pd
pd0$h <- 0

mult <- 1 # 100
pdm <- pd
pdm$r <- pd$r*mult

#pdm$x <- pd$x/mult

# Turchin: decrease d/k or c/x (stabilise)
# Reynolds: decrease x/r

# Rosenzweig-MacArthur ----------------------------------------------------

cxlim <- c(0,1.5)
cylim <- c(0,0.2)*mult

par(mfrow=c(1,1))
flowField(romac_phaser,x.lim=cxlim,y.lim=cylim,parameters=pdm[2,],points=30,add=FALSE)
clines <- list()
clines[[1]] <- nullclines(romac_phaser,x.lim = cxlim,y.lim = cylim,
                          parameters = pdm[1,], points = 100,
                          colour=rep("blue",2)
                          )
clines[[2]] <- nullclines(romac_phaser,x.lim = cxlim,y.lim = cylim,
                          parameters = pdm[2,], points = 100,
                          colour=rep("red",2)
                          )

tradj <- list()

equil1 <- c(0.275,0.15) # rough case using locator()
tradj[[1]] <- trajectory(romac_phaser, y0=equil1, t.step=60, t.end=60*10^5, parameters = pdm[1,])
tradj[[2]] <- trajectory(romac_phaser, y0=equil1, t.step=60, t.end=60*10^5, parameters = pdm[2,])

# Zero handling time ------------------------------------------------------

cxlim <- c(0,1)
cylim <- c(0,0.2)

flowField(romac_phaser,x.lim=cxlim,y.lim=cylim,parameters=pd0[2,],points=30,add=FALSE)
# plot(1,1,type="n",xlim=cxlim,ylim=cylim,xlab="R",ylab="C")
clines <- list()
clines[[1]] <- nullclines(romac_phaser, x.lim=cxlim, y.lim=cylim,
                          parameters = pd0[1,], points = 100,
                          colour=rep("blue",2)
)
clines[[2]] <- nullclines(romac_phaser, x.lim=cxlim, y.lim=cylim,
                          parameters=pd0[2,], points=100,
                          colour=rep("red",2)
)
equil1 <- c(0.275,0.15) # rough case using locator()
tradj[[1]] <- trajectory(romac_phaser, y0=equil1, t.step=60, t.end=60*10^5, parameters = pd0[1,])
tradj[[2]] <- trajectory(romac_phaser, y0=equil1, t.step=60, t.end=60*10^5, parameters = pd0[2,])

# DeAngelis ---------------------------------------------------------------

cxlim <- c(0,1.5)
cylim <- c(0,0.2)

flowField(deang_phaser, x.lim=cxlim,y.lim=cylim,parameters=pd[2,],points=30,add=FALSE)
clines <- list()
clines[[1]] <- nullclines(deang_phaser, x.lim=cxlim, y.lim=cylim,
                          parameters=pd[1,], points=100,
                          colour=rep("blue",2)
)
clines[[2]] <- nullclines(deang_phaser, x.lim = cxlim, y.lim = cylim,
                          parameters=pd[2,], points=100,
                          colour=rep("red",2)
)

equil1 <- c(0.33,0.17)
tradj <- trajectory(deang_phaser, y0=equil1, t.step=60, t.end=60*100000, parameters = pd[2,])
  # damped oscillations

# Chemostat ---------------------------------------------------------------

cxlim <- c(0,0.5)
cylim <- c(0,20)

flowField(chemo_phaser, x.lim=cxlim, y.lim=cylim, parameters=pd[2,], points=30, add=FALSE)
clines <- list()
clines[[1]] <- nullclines(chemo_phaser, x.lim=cxlim, y.lim=cylim,
                          parameters=pd[1,], points=100,
                          colour=rep("blue",2)
)
clines[[2]] <- nullclines(chemo_phaser, x.lim=cxlim, y.lim=cylim,
                          parameters=pd[2,], points=100,
                          colour=rep("red",2)
)

# Chemostat - zero handling time ------------------------------------------

cxlim <- c(0,0.5)
cylim <- c(0,20)

clines <- list()
flowField(chemo_phaser, x.lim=cxlim, y.lim=cylim, parameters=pd0[2,], points=30, add=FALSE)
clines[[1]] <- nullclines(chemo_phaser,x.lim=cxlim, y.lim=cylim,
                          parameters=pd0[1,], points=100,
                          colour=rep("blue",2)
)
clines[[2]] <- nullclines(chemo_phaser, x.lim=cxlim, y.lim=cylim,
                          parameters=pd0[2,], points=100,
                          colour=rep("red",2)
)

equil1 <- c(0.22,10.9)
tradj <- trajectory(chemo_phaser, y0=equil1, t.step=60, t.end=60*10^5, parameters = pd0[2,])
  # separation of timescales looks apt here

# Bazykin model -----------------------------------------------------------

pd$E <- 10^-4 # strength of consumer DD - arbitrarily chosen

cxlim <- c(0,1.5)
cylim <- c(0,0.2)

flowField(bazy_phaser, x.lim=cxlim,y.lim=cylim,parameters=pd[2,],points=30,add=FALSE)
clines <- list()
clines[[1]] <- nullclines(bazy_phaser, x.lim=cxlim, y.lim=cylim,
                          parameters=pd[1,], points=100,
                          colour=rep("blue",2)
)
clines[[2]] <- nullclines(bazy_phaser, x.lim=cxlim, y.lim=cylim,
                          parameters=pd[2,], points=100,
                          colour=rep("red",2)
)

equil1 <- c(1.44,0.1)
tradj <- trajectory(bazy_phaser, y0=equil1, t.step=60, t.end=60*10^5, parameters = pd[2,])
  # damped oscillations


# Vegetation regrowth model -----------------------------------------------

cxlim <- c(0,0.5)
cylim <- c(0,2)

par(mfrow=c(1,1))
flowField(regrow_phaser,x.lim=cxlim,y.lim=cylim,parameters=pdm[2,],points=30,add=FALSE)
clines <- list()
clines[[1]] <- nullclines(regrow_phaser,x.lim = cxlim,y.lim = cylim,
                          parameters = pdm[1,], points = 100,
                          colour=rep("blue",2)
)
clines[[2]] <- nullclines(regrow_phaser,x.lim = cxlim,y.lim = cylim,
                          parameters = pdm[2,], points = 100,
                          colour=rep("red",2)
)

equil1 <- c(0.275,0.53) # rough case using locator()
tradj <- list()
tradj[[1]] <- trajectory(regrow_phaser, y0=equil1, t.step=60, t.end=60*10^5, parameters = pdm[1,])
tradj[[2]] <- trajectory(regrow_phaser, y0=equil1, t.step=60, t.end=60*10^5, parameters = pdm[2,])


# Gompertz prey -----------------------------------------------------------

cxlim <- c(0,1.5)
cylim <- c(0,2)

par(mfrow=c(1,1))
flowField(gomp_phaser,x.lim=cxlim,y.lim=cylim,parameters=pdm[2,],points=30,add=FALSE)
clines <- list()
clines[[1]] <- nullclines(gomp_phaser,x.lim = cxlim,y.lim = cylim,
                          parameters = pdm[1,], points = 100,
                          colour=rep("blue",2)
)
clines[[2]] <- nullclines(gomp_phaser,x.lim = cxlim,y.lim = cylim,
                          parameters = pdm[2,], points = 100,
                          colour=rep("red",2)
)

tradj <- list()

equil1 <- c(0.275,0.475) # rough case using locator()
tradj[[2]] <- trajectory(gomp_phaser, y0=equil1, t.step=60, t.end=60*10^5, parameters = pdm[2,])

# Constant generalist predator --------------------------------------------

Cslice <- c(0.025,0.05,0.1,0.15)
nCslice <- length(Cslice)
nR <- 100
Rmax <- 2
Rseq <- seq(0,Rmax,length.out=nR)
dR_fixC <- matrix(nr=nR,nc=nT)

par(mfrow=c(2,2))
for(i in 1:nCslice){
  for(j in 1:nT){
    dR_fixC[,j] <- with(pd[j,], dR_dt_f(R=Rseq,C=Cslice[i],r=r,K=K,c=c,d=d)/Rseq)
  }
  matplot(log(Rseq),dR_fixC,type="l")
  abline(h=0,lty=2)
}
  # when type II generalist predator, can have Allee effect
  # and possibly multiple stable equilibria
  # NB: using loop creates conflict with i parameter

# Predator density-dependence ---------------------------------------------

Rmin <- min(c(Rstar1,Rstar2),na.rm=T)
Rmax <- max(c(Rstar1,Rstar2),na.rm=T)
Rseq <- seq(Rmin,Rmax,length.out=nC)

Cstar_f <- function(R,c){
  e <- try( 
    d <- uniroot(dC_dt_f, R=R, c=c, lower=10^-100, upper=10^100), 
    silent = TRUE 
  ) 
  if(class(e)=="try-error") { 
    return(NA) 
  } 
  else{ 
    return(d$root)	
  } 
}

# Manual phase plot construction ------------------------------------------

# Uses uniroot - but doesn't work well when more than one equilibrium resource
# density for a given consumer density

nC <- 1000
Cmin <- 0
Cmax <- 0.25
Cseq <- seq(Cmin,Cmax,length.out=nC)

Rstar <- matrix(nr=nC,nc=nT)
Rdash <- Cstar <- vector(length=nT)
for(i in 1:nT){
  Rstar[,i] <- with(pd[i,], sapply(Cseq,Rstar_f,r=r,K=K,a=a,d=d))
  Rdash[i] <- with(pd[i,], 
                   uniroot(dC_dt_f,C=1,c=c,d=d,x=x,lower=10^-10,upper=10^10)$root
  )
  Cstar[i] <- with(pd[i,],
                   uniroot(dR_dt_f,R=Rdash[i],r=r,K=K,c=c,d=d,lower=10^-100,upper=10^100)$root
  )
}

lcols <- c("blue","red")
matplot(Rstar,Cseq,type="l",xlab="R",ylab="C",lty=1,col=lcols)
abline(v=Rdash,col=lcols,lty=2)
abline(h=Cstar,col=lcols,lty=3)

with(pd[i,], dR_dt_f(R=seq(0,0.25,length.out=100),C=0.1,r=r,K=K,c=c,d=d))
with(pd[i,], dC_dt_f(R=seq(0,0.25,length.out=100),C=0.1,r=r,K=K,c=c,d=d))

