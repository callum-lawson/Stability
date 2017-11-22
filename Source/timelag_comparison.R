### Compare effects of time lags on consumer-resource dynamics ###

require(deSolve)
require(fields)
source("Source/predprey_functions_general.R")

x0 <- arrrate(0,e0["x"],e1["x"])

tmax <- 60^2 * 24 * 7 * 52 * 10 # maximum length of time in seconds
tf <- 1000
tseq <- seq(0,tmax,length.out=tf)

zmu <- 0
zsig <- 0 # wave amplitude
zf <- 10 # wave frequency over whole time series
zl <- tmax/zf 

sf <- 100 + 1 # number of seasons (+1 because new season starts right at end)
sl <- tmax/sf
sstart <- seq(0,tmax,length.out=sf)
sseq <- sseqgen(tseq,sstart)

zparms <- list(zmu=zmu,zsig=zsig,zl=zl)
eparms <- list(e0=e0,e1=e1,Rtype="replenish")

R0 <- 1
C0 <- 10^-2
y0 <- c(R=R0,C=C0)

discrete <- DRCt_disc(y0,tseq,sseq,sstart,parms=c(eparms,zl=zl,zmu=zmu,zsig=0))
matplot(tseq,log(discrete[,-1]),type="l",col="blue",bty="n")
