### Compare effects of time lags on consumer-resource dynamics ###

### TODO
# - add continuous-time lag in resource as well as consumer births

require(deSolve)
require(fields)
source("Source/predprey_functions_general.R")

x0 <- arrrate(0,e0["x"],e1["x"])

tmax <- 60^2 * 24 * 7 * 52 # * 10 # maximum length of time in seconds
tf <- 10^3
tseq <- seq(0,tmax,length.out=tf)

zmu <- 0
zsig <- 0 # wave amplitude
zf <- 10 # wave frequency over whole time series
zl <- tmax/zf 

sf <- tf # number of seasons 
  # must be <= number of recorded timesteps
sl <- tmax/sf
sstart <- seq(0,tmax,length.out=sf)
sseq <- sseqgen(tseq,sstart)

zparms <- list(zmu=zmu,zsig=zsig,zl=zl,tau=sl)
eparms <- list(e0=e0,e1=e1)

R0 <- 10^1
C0 <- 10^-2
y0 <- c(R=R0,C=C0)

standard <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(zparms,eparms))
delayed <- dede(y=y0,times=tseq,func=dRCt_delay,parms=c(zparms,eparms))
discrete1 <- DRCt_disc(y0,tseq,sseq,sstart,parms=c(zparms,eparms,Rtype="replenish"))
discrete2 <- DRCt_disc(y0,tseq,sseq,sstart,parms=c(zparms,eparms,Rtype="cohabit"))
discrete3 <- DRCt_disc(y0,tseq,sseq,sstart,parms=c(zparms,eparms,Rtype="remove"))
discrete4 <- DRCt_disc(y0,tseq,sseq,sstart,parms=c(zparms,eparms,Rtype="persist"))

require(fields)
cols <- tim.colors(4)

par(mar=c(4,4,2,2))
matplot(tseq,log(discrete1[,-1]),type="l",col=cols[1],bty="n")
matplot(tseq,log(discrete2[,-1]),type="l",col=cols[2],add=TRUE)
matplot(tseq,log(discrete3[,-1]),type="l",col=cols[3],add=TRUE)
matplot(tseq,log(discrete4[,-1]),type="l",col=cols[4],add=TRUE)
matplot(tseq,log(standard[,-1]),type="l",col="black",add=TRUE)
matplot(tseq,log(delayed[,-1]),type="l",col="gray",add=TRUE)
# abline(v=sstart,col="grey",lty=3)  
legend("topright",legend=c("replenish","cohabit","remove","persist"),
       col=cols,lty=1,bty="n")

matplot(tseq,log(cbind(discrete1[,2],discrete2[,2],discrete3[,2],discrete4[,2],
                standard[,2],delayed[,2])),
     type="l",
     col=c(cols,"black","gray"),
     bty="n"
     )
legend("topright",legend=c("replenish","cohabit","remove","persist"),
       col=cols,lty=1,bty="n")

