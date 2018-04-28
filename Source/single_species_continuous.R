### Continuous single-species population dynamics with environmental variability ###

require(deSolve)
require(fields)
source("Source/predprey_functions_general.R")

# Input parameters --------------------------------------------------------

x0 <- arrrate(0,e0["x"],e1["x"])

tmax <- 60^2 * 24 * 7 * 52 * 10 # maximum length of time in seconds
tf <- 1001
tseq <- seq(0,tmax,length.out=tf)

zmu <- 0 # =20°C 
zsig <- 5 # wave amplitude
zf <- 10 # wave frequency over whole time series
zl <- tmax/zf 

e0 <- c(
  m = 10^-5, # very high value relative to consumer
  k = 10,
  a = 6*10^-7, # 3.181989*10^-9, # estimated from data
  h = 61000, # 0.61, # 1685.586,     # estimated from data
  mu = 2.689*10^-6 # 2.689*10^-6
)

e1 <- c(
  m = 0, # 0.639,
  k = 0, # -0.772,
  a = -0.03, # 0.5091663,   # estimated from data
  h = 0, # -1.9, # -0.19, # -0.4660012, # estimated from data
  mu = 0 # 0.639
)

zparms <- list(zmu=zmu,zsig=zsig,zl=zl)
eparms <- list(e0=e0,e1=e1)

R0 <- 10^1
C0 <- 10^1
y0 <- c(R=R0,C=C0)

# Increasing temperature variability --------------------------------------

### Consumer present

lvar <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(eparms,zl=zl,zmu=zmu,zsig=0))
mvar <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(eparms,zl=zl,zmu=zmu,zsig=2.5))
hvar <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(eparms,zl=zl,zmu=zmu,zsig=5))

par(mfrow=c(1,1))
matplot(tseq,log(hvar[,-1]),type="l",col="red",bty="n")
matplot(tseq,log(mvar[,-1]),type="l",col="black",add=TRUE)
matplot(tseq,log(lvar[,-1]),type="l",col="blue",add=TRUE)
# variance -> both predator and prey fluctuate slightly more and 
# reach lower densities
# also, cycles slightly slower (?)
# increasing variance accentuates these effects

#with(zparms, curve(zt_cyclic(x,zmu,zsig,zl),xlim=c(0,tmax),n=1001))

# Population growth curves ------------------------------------------------

dCCdt <- Vectorize(
  function(R,C,z,eparms){
    parms <- with(eparms, as.list(c( R=R, arrrate(z,e0,e1) )) )
    dCt_cons(t=0,y=C,parms=parms)/C
  }, 
  vectorize.args=c("R","C","z")
)

library(reshape2)
Cseq <- exp(seq(-5,5,length.out=100))
zseq <- seq(zmu-5,zmu+5,length.out=10)
rd <- expand.grid(C=Cseq,z=zseq)
Rstar <- with(rd, Rstarcalc(C,z,eparms=eparms))
rd$r <- with(rd, dCCdt(R=Rstar,C=C,z=z,eparms=eparms))

ra <- acast(melt(rd,id=c("C","z")),C~z)
par(mfrow=c(1,1))
matplot(log(Cseq),ra,type="l",lty=1,col=heat.colors(10))
abline(h=0,lty=3,col="grey")

# Discrete time -----------------------------------------------------------
# 
# sf <- 100 + 1 # number of seasons (+1 because new season starts right at end)
# sl <- tmax/sf
# sstart <- seq(0,tmax,length.out=sf)
# sseq <- sseqgen(tseq,sstart)

dRCt_disc2 <- function(t,y,parms){
  with(parms, dRC_disc(y,m,k,a,h,mu))
}

DCC <- Vectorize(
  function(R,C,z){
    parms <- with(eparms, as.list(arrrate(z,e0,e1)) )
    ode(y=c(R=R,C=C,E=0),
        times=c(0,10^10),
        func=dRCt_disc2,
        parms=parms
    )[2,"E"]
  },
  vectorize.args=c("R","C","z")
)

rd2 <- rd
rd2$r <- with(rd2, DCC(Rstar,C,z)/C)
  # using Rstar at start, so still making timescale separation assumption

ra2 <- acast(melt(rd2,id=c("C","z")),C~z)
par(mfrow=c(1,1))
matplot(log(Cseq),ra2,type="l",lty=1,col=heat.colors(10))
abline(h=0,lty=3,col="grey")
# could alternatively use equilibrium solver 
# (for this case, where no egg mortality)
