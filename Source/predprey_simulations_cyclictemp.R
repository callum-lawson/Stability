####################################################################################
# Include temperature fluctuations in pred-prey models from Fussmann et al. (2014) #
####################################################################################

### TODO
# - check approximation for errors, and how this varies with n perturbations
# - resource dynamics when generalist consumer unaffected by resource (see Turchin)
# - compare different temp fluctuation sizes (instead of different means)
# - two resource species
# - Gaussian-process temperatures

require(deSolve)
require(fields)
source("Source/predprey_functions_general.R")

# Input parameters --------------------------------------------------------

e0 <- c(
  m = 0,
  r = 8.715*10^-7,
  K = 5.623,
  a = 6*10^-7, # 3.181989*10^-9, # estimated from data
  h = 0.61, # 1685.586,     # estimated from data
  x = 2.689*10^-6
)

e1 <- c(
  m = 0,
  r = 0.84,
  K = -0.772,
  a = -0.03, # 0.5091663,   # estimated from data
  h = -0.19, # -0.4660012, # estimated from data
  x = 0.639
)
# r units are per SECOND; pop more than triples every 24h

tmax <- 60^2 * 24 * 7 * 52 * 10 # maximum length of time in seconds
tf <- 1000
tseq <- seq(0,tmax,length.out=tf)

zmu <- 0
zsig <- 5 # wave amplitude
zf <- 10 # wave frequency over whole time series
zl <- tmax/zf 

sf <- 100 + 1 # number of seasons (+1 because new season starts right at end)
sl <- tmax/sf
sstart <- seq(0,tmax,length.out=sf)
sseq <- sseqgen(tseq,sstart)

parms <- list(zmu=zmu,zsig=zsig,zl=zl,e0=e0,e1=e1,Rtype="replenish")

R0 <- 1
C0 <- 10^-2
y0 <- c(R=R0,C=C0)

# Increasing temperature variability --------------------------------------

baseparms <- list(zmu=zmu,zl=zl,e0=e0,e1=e1)

### Consumer present

lvar <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(baseparms,zsig=0))
mvar <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(baseparms,zsig=5))
hvar <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(baseparms,zsig=10))

par(mfrow=c(1,1))
matplot(tseq,log(lvar[,-1]),type="l",col="blue",bty="n")
matplot(tseq,log(mvar[,-1]),type="l",col="black",add=TRUE)
matplot(tseq,log(hvar[,-1]),type="l",col="red",add=TRUE)
  # variance -> both predator and prey fluctuate slightly more and 
  # reach lower densities
  # also, cycles slightly slower (?)
  # increasing variance accentuates these effects

### Consumer absent

## r plots
mycols <- c("blue","black","red")
library(reshape2)
Rseq <- exp(seq(2,4,length.out=100))
zseq <- seq(-5,5,length.out=3)
rd <- expand.grid(R=Rseq,z=zseq)
rd$r <- with(rd, rR(zt=z,R=R,parms=parms))
ra <- acast(melt(rd,id=c("R","z")),R~z)
matplot(log(Rseq),ra,type="l",lty=1,col=mycols)
abline(h=0,lty=3,col="grey")
  # higher temp -> higher growth at low densities but lower K
  # slightly stronger regulation above than below K

parms5 <- mapply(arrint,e0=e0,e1=e1,MoreArgs=list(zmu=0,zsig=0))
parms10 <- mapply(arrint,e0=e0,e1=e1,MoreArgs=list(zmu=0,zsig=5))

arrint(e0=e0[2],e1=e1[2],zmu=0,zsig=1)
arr_z(0,e0=e0[2],e1=e1[2],zmu=0,zsig=1)

lvar <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(baseparms,zsig=0))
mvar <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(baseparms,zsig=5))
hvar <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(baseparms,zsig=10))

## Dynamics
lvar_Ronly <- ode(y=c(y0[1],C=0),times=tseq,func=dRCt_cont,parms=c(baseparms,zsig=0))
mvar_Ronly <- ode(y=c(y0[1],C=0),times=tseq,func=dRCt_cont,parms=c(baseparms,zsig=5))
hvar_Ronly <- ode(y=c(y0[1],C=0),times=tseq,func=dRCt_cont,parms=c(baseparms,zsig=10))

par(mfrow=c(1,1))
matplot(tseq,
        log(cbind(lvar_Ronly[,2],mvar_Ronly[,2],hvar_Ronly[,2])),
        type="l",
        lty=1,
        col=mycols,
        bty="n"
        )

### Brief results summary

  # DISCRETE
  # - sometimes transient dynamics persist for ages; other times, snaps straight into 
  # new regime 
  # - the longer the intervals between temperature switches, the lesss transient 
  # dynamics matter
  # - negative perturbations > positive perturbations?
  # - perturbations every timestep breaks the simulation

  # CONTINUOUS
  # - wavelength -> 0 => constant fluctuations
  # (but how does this compare to mean in cons env?)

  # predator dynamics (rate of increase) slower than prey -> sharp vs blunt peaks

  # resource follows sin fluctuations (driven by sin fluctuations in pred feeding)
  # at low temps, predator fluctuations:
  # - prey maxima are more rounded
  # - prey fluctuations increase in size
  # - pred fluctuations become "double-peaked"
  # - pred fluctuations shrink slightly but then become v. large at lowest temps

# End scenarios -----------------------------------------------------------

### Constant temperature 
## Slow fluctuations
# K of both prey and pred decreases at higher temps (but faster to reach K)

## Fast fluctuations (non-linear averaging)


# Delays and discrete dynamics --------------------------------------------

### Models
standard <- ode(y=y0,times=tseq,func=dRCt_cont,parms=parms)
delayed <- dede(y=y0,times=tseq,func=dRCt_delay,parms=c(parms,tau=60^2*24*7))
discrete <- DRCt_disc(y0,tseq,sseq,sstart,parms)

### Plots
layout(cbind(c(1,1,1,2)))
par(mar=c(4,4,2,2))
matplot(tseq,log(discrete[,-1]),type="l",col="blue",bty="n")
matplot(tseq,log(standard[,-1]),type="l",col="black",add=TRUE)
matplot(tseq,log(delayed[,-1]),type="l",col="red",add=TRUE)
abline(v=seq(0,tmax,length.out=zf+1),col="grey",lty=3)  
# +1 accounts for t=0
par(mar=c(2,4,2,2))
plot(zt_cyclic(tseq,zmu,zsig,zl)~tseq,type="l",bty="n",xaxt="n")

### Autocorrelation

send <- sapply(1:sf, function(x) max(which(sseq==x)))

acn <- discrete[send,3]
acd <- data.frame(a=acn[1:(sf-5)],
                  b=acn[2:(sf-4)],
                  c=acn[3:(sf-3)],
                  d=acn[4:(sf-2)],
                  e=acn[5:(sf-1)],
                  f=acn[6:(sf-0)]
)
acm <- glm(a~offset(b)+b+c+d+e+f,data=acd)
summary(acm)
# fourth-order model
acf(acn)

