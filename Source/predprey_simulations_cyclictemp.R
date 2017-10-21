####################################################################################
# Include temperature fluctuations in pred-prey models from Fussmann et al. (2014) #
####################################################################################

require(deSolve)
require(fields)
source("Source/predprey_functions_general.R")

# Input parameters --------------------------------------------------------

E0 <- c(
  u = 0,
  r = 8.715*10^-7,
  K = 5.623,
  a = 6*10^-7, # 3.181989*10^-9, # estimated from data
  h = 0.61, # 1685.586,     # estimated from data
  x = 2.689*10^-6
)

E1 <- c(
  u = 0,
  r = 0.84,
  K = -0.772,
  a = -0.03, # 0.5091663,   # estimated from data
  h = -0.19, # -0.4660012, # estimated from data
  x = 0.639
)
# r units are per SECOND; pop more than triples every 24h

nt <- 10^3 # number of timesteps to calculate densities for
nyears <- 100
tmax <- nyears * 365 * 24 * 60^2  # maximum length of time in seconds
tseq <- seq(0,tmax,length.out=nt)

Tmu <- 0
Tsd <- 10
nTwaves <- 2
Tperiod <- 1/nTwaves * tmax
parms <- list(Tmu=Tmu,Tsd=Tsd,Tperiod=Tperiod,E0=E0,E1=E1)

R0 <- 1
C0 <- 10^-3 # 0 -> resource-only model
y <- c(R0,C0)

# Simulate and plot results -----------------------------------------------

par(mfrow=c(1,1))
wow <- ode(y=c(R0,C0),times=tseq,func=dCRt_cyclic,parms=parms)
matplot(tseq,log(wow[,-1]),type="l")
abline(v=seq(0,tmax,length.out=nTwaves+1),col="blue",lty=3)  
# +1 accounts for t=0

# Simulations with discrete temperatures ----------------------------------

nt <- 1000 # number of timesteps to calculate densities for
tmax <- 10^4 * 60^2  # maximum length of time in seconds
tseq <- seq(0,tmax,length.out=nt)
nP <- 2 # number of perturbations over time series
lP <- nt/nP # perturbation length in nt units - has to be whole number

Tmu <- 20
Psd <- 5
Pvals <- rnorm(nP,mean=0,sd=Psd)
Pt <- rep(Pvals,each=lP)
Ttseq <- rep(Tmu + 293.15,nt) + Pt

rseq <- arrhenius(Ttseq,rE0,rE1)
Kseq <- arrhenius(Ttseq,KE0,KE1)
cseq <- arrhenius(Ttseq,cE0,cE1)
dseq <- arrhenius(Ttseq,dE0,dE1)
xseq <- arrhenius(Ttseq,xE0,xE1)
eps <- 0.85
  # will have to hold each of these constant for multiple timesteps
  # (otherwise numerical approximation of continuous time doesn't work)

R0 <- 1
C0 <- 1

### Numerical integration

ode1 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=romac_dis,parms=NULL)
matplot(tseq,log(ode1[,-1]),type="l")
abline(v=seq(0,tmax,length.out=nP+1),col="blue",lty=3)  
  # +1 accounts for t=0

acn <- ode1[,3]
acd <- data.frame(a=acn[1:(nt-5)],
                  b=acn[2:(nt-4)],
                  c=acn[3:(nt-3)],
                  d=acn[4:(nt-2)],
                  e=acn[5:(nt-1)],
                  f=acn[6:(nt-0)]
                  )
acm <- glm(a~offset(b)+b+c+d+e+f,data=acd)
summary(acm)
  # fourth-order model
acf(acn)

# Simulations with sinusoidal temperatures --------------------------------

nt <- 1000 # number of timesteps to calculate densities for
tmax <- 10^4 * 60^2  # maximum length of time in seconds
tseq <- seq(0,tmax,length.out=nt)
lP <- tmax/10^1 # wavelength in seconds

Tmu <- 20 # mean temperature in degrees C
Psd <- 5 # wave amplitude

R0 <- 1
C0 <- 1

### Numerical integration

Psd <- 0
ode1 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=romac_sin,parms=NULL)
Psd <- 5
ode2 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=romac_sin,parms=NULL)
matplot(tseq,log(ode1[,-1]),type="l",lty=1)
matplot(tseq,log(ode2[,-1]),type="l",add=T,lty=2)

matplot(tseq,log(cbind(ode1[,2],ode2[,2])),type="l",lty=1)
matplot(tseq,log(cbind(ode1[,3],ode2[,3])),type="l",lty=1)
  # variance -> both predator and prey fluctuate slightly more and 
  # reach lower densities
  # also, cycles slightly slower
  # increasing variance accentuates these effects

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

# Chemostat resource ------------------------------------------------------

nt <- 10^3 # number of timesteps to calculate densities for
tmax <- 10^4 * 60^2  # maximum length of time in seconds
tseq <- seq(0,tmax,length.out=nt)
lP <- tmax/5 # wavelength in seconds
R0 <- 1
C0 <- 1
Psd <- 0

### Constant temperature

Psd <- 0
Tmu <- 0
ode0 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
Tmu <- 10
ode10 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
Tmu <- 20
ode20 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
Tmu <- 30
ode30 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
Tmu <- 40
ode40 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
Tmu <- 50
ode50 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)

par(mfrow=c(1,2))
matplot(tseq,
        log(cbind(ode0[,2],ode10[,2],ode20[,2],ode30[,2],ode40[,2],ode50[,2])),
        type="l",lty=1,col=tim.colors(6),xlab="t",ylab="R"
        )
matplot(tseq,
        log(cbind(ode0[,3],ode10[,3],ode20[,3],ode30[,3],ode40[,3],ode50[,3])),
        type="l",lty=1,col=tim.colors(6),xlab="t",ylab="C"
        )
  # K of both prey and pred decreases at higher temps (but faster to reach K)

### Extinction?

Tmu <- 25
ode25 <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
matplot(tseq,log(ode25[,2:3]),type="l")

### Fluctuating temperature

Psd <- 10
Tmu <- 0
ode0v <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
Tmu <- 10
ode10v <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
Tmu <- 20
ode20v <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
Tmu <- 30
ode30v <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
Tmu <- 40
ode40v <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)
Tmu <- 50
ode50v <- ode(y=c(R0=R0,C0=C0),times=tseq,func=chemo_sin,parms=NULL)

par(mfrow=c(1,2))
matplot(tseq,
        log(cbind(ode0v[,2],ode10v[,2],ode20v[,2],ode30v[,2],ode40v[,2],ode50v[,2])),
        type="l",lty=1,col=tim.colors(6),xlab="t",ylab="R"
        )
matplot(tseq,
        log(cbind(ode0v[,3],ode10v[,3],ode20v[,3],ode30v[,3],ode40v[,3],ode50v[,3])),
        type="l",lty=1,col=tim.colors(6),ylim=c(-2,2.5),xlab="t",ylab="C"
        )
  # resource follows sin fluctuations (driven by sin fluctuations in pred feeding)
  # at low temps, predator fluctuations:
  # - prey maxima are more rounded
  # - prey fluctuations increase in size
  # - pred fluctuations become "double-peaked"
  # - pred fluctuations shrink slightly but then become v. large at lowest temps

# Predator absent ---------------------------------------------------------

# absence of predator -> logistic prey dynamics

nt <- 10^3 # number of timesteps to calculate densities for
tmax <- 10^4 * 60^2  # maximum length of time in seconds
tseq <- seq(0,tmax,length.out=nt)
lP <- tmax/5 # wavelength in seconds
R0 <- 1
Psd <- 0

### Constant temperature

Psd <- 0
Tmu <- 0
ode0 <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)
Tmu <- 10
ode10 <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)
Tmu <- 20
ode20 <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)
Tmu <- 30
ode30 <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)
Tmu <- 40
ode40 <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)
Tmu <- 50
ode50 <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)

par(mfrow=c(1,2))
matplot(tseq,
        log(cbind(ode0[,2],ode10[,2],ode20[,2],ode30[,2],ode40[,2],ode50[,2])),
        type="l",lty=1,col=tim.colors(6)
)
  # lower temp -> lower K (built into assumptions)
  
Psd <- 5
Tmu <- 0
ode0v <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)
Tmu <- 10
ode10v <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)
Tmu <- 20
ode20v <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)
Tmu <- 30
ode30v <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)
Tmu <- 40
ode40v <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)
Tmu <- 50
ode50v <- ode(y=c(R0=R0),times=tseq,func=preyonly_sin,parms=NULL)

par(mfrow=c(1,2))
matplot(tseq,
        log(cbind(ode0v[,2],ode10v[,2],ode20v[,2],ode30v[,2],ode40v[,2],ode50v[,2])),
        type="l",lty=1,col=tim.colors(6)
)
  # fluctuations similar size   
  # but tracking slower at lower temps, so faster fluctuations may make a difference

rt_preyonly <- function(Rt,Tt){
  TtK <- Tt + 293.15
  r <- arrhenius(TtK,rE0,rE1)
  K <- arrhenius(TtK,KE0,KE1)
  rt <- r*(1-Rt/K)
}

nRseq <- 100
nTseq <- 6
Rtseq <- exp(seq(-1,0.5,length.out=nRseq))
Ttseq <- seq(10,40,length.out=nTseq)
Rxdat <- expand.grid(Rt=Rtseq,Tt=Ttseq)

rtmat <- matrix(nr=nRseq,nc=nTseq)
rtmat[] <- with(Rxdat,rt_preyonly(Rt,Tt))
par(mfrow=c(1,1))
matplot(log(Rtseq),rtmat,type="l",col=tim.colors(nTseq),lty=1)
abline(h=0,lty=3)
  # higher temp -> higher growth at low densities but lower K
  # slightly stronger regulation above than below K

# Constant predator -------------------------------------------------------

# generalist predator unaffected by resource (see Turchin)

### TODO
# - check approximation for errors, and how this varies with n perturbations
# - resource dynamics when generalist consumer unaffected by resource
# - compare different temp fluctuation sizes (instead of different means)


