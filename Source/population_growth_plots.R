### Plot r~X plots for consumer-reosurce models ###
### using timescale separation                  ###

require(deSolve)
require(fields)
source("Source/predprey_functions_general.R")

# Input parameters --------------------------------------------------------

x0 <- arrrate(0,e0["x"],e1["x"])

tmax <- 60^2 * 24 * 7 * 52 * 10 # maximum length of time in seconds
tf <- 1000
tseq <- seq(0,tmax,length.out=tf)

zmu <- 0 # =20°C 
zsig <- 5 # wave amplitude
zf <- 10 # wave frequency over whole time series
zl <- tmax/zf 

sf <- 100 + 1 # number of seasons (+1 because new season starts right at end)
sl <- tmax/sf
sstart <- seq(0,tmax,length.out=sf)
sseq <- sseqgen(tseq,sstart)

zparms <- list(zmu=zmu,zsig=zsig,zl=zl)
eparms <- list(e0=e0,e1=e1)

R0 <- 10^1
C0 <- 10^-2
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

# Population growth plots -------------------------------------------------

## r plots
# using timescale separation 

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
  a = 0.3, # -0.03, # 0.5091663,   # estimated from data
  h = 0, # -0.19, # -0.4660012, # estimated from data
  x = 0 # 0.639
  # e = 0.639
)
eparms <- list(e0=e0,e1=e1)

library(reshape2)
zseq <- seq(zmu-5,zmu+5,length.out=10)
Cseq <- exp(seq(-5,-2,length.out=100))
rd <- expand.grid(C=Cseq,z=zseq)
rd$r <- with(rd, rC_separate(z,C,eparms))
ra <- acast(melt(rd,id=c("C","z")),C~z)
par(mfrow=c(1,1))
matplot(log(Cseq),ra,type="l",lty=1,col=heat.colors(10))
abline(h=0,lty=3,col="grey")
