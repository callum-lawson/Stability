#################################################################################
# Phase-space plots when underlying predator or prey paramaters are fluctuating # 
#################################################################################

require(rootSolve)
require(phaseR)

source("Source/predprey_functions.R")

# Species parameters ------------------------------------------------------

rE0 <- 8.715*10^-7
KE0 <- 5.623
aE0 <- 8.408*10^-6
bE0 <- 3.664
xE0 <- 2.689*10^-6

# rE1 <- 0.84
# KE1 <- -0.772
# aE1 <- 0.467
# bE1 <- -0.114
# xE1 <- 0.639
# From species-level averages
# change slopes later
# r units are per SECOND; pop more than triples every 24h

rE1 <- 0.84
KE1 <- -0.508
aE1 <- 0.708
bE1 <- -0.678
xE1 <- 0.428
# From Fig. S1

# Functions ---------------------------------------------------------------

dR_dt_f <- function(R,C,r,K,a,b){
  R * ( r*(1-R/K) - a*C/(b+R) )
  } 

dC_dt_f <- function(R,C,a,b,x,eps=0.85){
  C * ( eps*a*R/(b+R) - x )
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
  b <- parameters$b
  x <- parameters$x
  eps <- parameters$eps
  dy <- numeric(2)
  dy[1] <- y[1] * ( r*(1-y[1]/K) - a*y[2]/(b+y[1]) )
  dy[2] <- y[2] * ( eps*a*y[1]/(b+y[1]) - x )   
  return(list(dy))
}

case_phaser <- function(t, y, parameters){
  r <- parameters$r
  K <- parameters$K
  a <- parameters$a
  x <- parameters$x
  eps <- parameters$eps
  dy <- numeric(2)
  dy[1] <- y[1] * ( r*(1-y[1]/K) - a*y[2] )
  dy[2] <- y[2] * ( eps*a*y[1] - x )   
  return(list(dy))
}

Tseq <- c(15,25) + 293.15
nT <- length(Tseq)
pd <- data.frame(
  T = Tseq,
  r = arrhenius(Tseq,rE0,rE1),
  K = arrhenius(Tseq,KE0,KE1),
  a = arrhenius(Tseq,aE0,aE1),
  b = arrhenius(Tseq,bE0,bE1),
  x = arrhenius(Tseq,xE0,xE1),
  eps = 0.85
  )

parameters <- pd[i,]


# Rosenzweig-MacArthur ----------------------------------------------------

cxlim <- c(0,1.5)
cylim <- c(0,0.2)

flowField(romac_phaser, x.lim=cxlim,y.lim=cylim,parameters=pd[1,],points=15,add=FALSE)
clines <- list()
clines[[1]] <- nullclines(romac_phaser,x.lim = cxlim,y.lim = cylim,
                          parameters = pd[1,], points = 100,
                          colour=rep("blue",2)
                          )
clines[[2]] <- nullclines(romac_phaser,x.lim = cxlim,y.lim = cylim,
                          parameters = pd[2,], points = 100,
                          colour=rep("red",2)
                          )

# Case zero-handling time model -------------------------------------------

cxlim <- c(0,1)
cylim <- c(0,0.2)

flowField(case_phaser, 
          x.lim = cxlim,
          y.lim = cylim,
          parameters = pd[1,],
          points = 15, add = FALSE)
# plot(1,1,type="n",xlim=cxlim,ylim=cylim,xlab="R",ylab="C")
clines <- list()
clines[[1]] <- nullclines(case_phaser,x.lim = cxlim, y.lim = cylim,
                          parameters = pd[1,], points = 100,
                          colour=rep("blue",2)
)
clines[[2]] <- nullclines(case_phaser,x.lim = cxlim, y.lim = cylim,
                          parameters = pd[2,], points = 100,
                          colour=rep("red",2)
)

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
  Rstar[,i] <- with(pd[i,], sapply(Cseq,Rstar_f,r=r,K=K,a=a,b=b))
  Rdash[i] <- with(pd[i,], 
                   uniroot(dC_dt_f,C=1,a=a,b=b,x=x,lower=10^-10,upper=10^10)$root
  )
  Cstar[i] <- with(pd[i,],
                   uniroot(dR_dt_f,R=Rdash[i],r=r,K=K,a=a,b=b,lower=10^-100,upper=10^100)$root
  )
}

lcols <- c("blue","red")
matplot(Rstar,Cseq,type="l",xlab="R",ylab="C",lty=1,col=lcols)
abline(v=Rdash,col=lcols,lty=2)
abline(h=Cstar,col=lcols,lty=3)

with(pd[i,], dR_dt_f(R=seq(0,0.25,length.out=100),C=0.1,r=r,K=K,a=a,b=b))
with(pd[i,], dC_dt_f(R=seq(0,0.25,length.out=100),C=0.1,r=r,K=K,a=a,b=b))

# Predator density-dependence ---------------------------------------------

Rmin <- min(c(Rstar1,Rstar2),na.rm=T)
Rmax <- max(c(Rstar1,Rstar2),na.rm=T)
Rseq <- seq(Rmin,Rmax,length.out=nC)

Cstar_f <- function(R,a){
  e <- try( 
    d <- uniroot(dC_dt_f, R=R, a=a, lower=10^-100, upper=10^100), 
    silent = TRUE 
  ) 
  if(class(e)=="try-error") { 
    return(NA) 
  } 
  else{ 
    return(d$root)	
  } 
}
