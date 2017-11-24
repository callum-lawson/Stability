### Compare effects of time lags on consumer-resource dynamics ###

### TODO
# - add continuous-time lag in resource as well as consumer births

require(deSolve)
require(fields)
source("Source/predprey_functions_general.R")

x0 <- arrrate(0,e0["x"],e1["x"])

tmax <- 60^2 * 24 * 365 # maximum length of time in seconds
tf <- 10^3
tseq <- seq(0,tmax,length.out=tf)

zmu <- 0
zsig <- 0 # wave amplitude
zf <- 10 # wave frequency over whole time series
zl <- tmax/zf 

sf <- 365
  # number of seasons - must be <= number of recorded timesteps
sl <- tmax/sf

zparms <- list(zmu=zmu,zsig=zsig,zl=zl,tau=sl)
eparms <- list(e0=e0,e1=e1)

R0 <- 10^1
C0 <- 10^-2
y0 <- c(R=R0,C=C0)

standard <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(zparms,eparms))
delayed <- dede(y=y0,times=tseq,func=dRCt_delay,parms=c(zparms,eparms))
discrete1 <- DRCt_disc(y0,tseq,sf,parms=c(zparms,eparms,Rtype="replenish"))
discrete2 <- DRCt_disc(y0,tseq,sf,parms=c(zparms,eparms,Rtype="cohabit"))
discrete3 <- DRCt_disc(y0,tseq,sf,parms=c(zparms,eparms,Rtype="remove"))
discrete4 <- DRCt_disc(y0,tseq,sf,parms=c(zparms,eparms,Rtype="persist"))

par(mfrow=c(1,2),mar=c(4,4,2,2))
require(fields)
cols <- tim.colors(4)
matplot(tseq,log(cbind(discrete1[,3],discrete2[,3],discrete3[,3],discrete4[,3],
                       standard[,3],delayed[,3])),
        type="l",
        col=c(cols,"black","gray"),
        lty=1:6,
        bty="n"
)
matplot(tseq,log(cbind(discrete1[,2],discrete2[,2],discrete3[,2],discrete4[,2],
                standard[,2],delayed[,2])),
     type="l",
     col=c(cols,"black","gray"),
     lty=1:6,
     bty="n"
     )
legend("topright",legend=c("replenish","cohabit","remove","persist"),
       col=cols,lty=1:4,bty="n")

# Effects of season length ------------------------------------------------

aparms <- c(zparms,eparms)

DRCt_disc_one <- function(y0,tmax,Rtype){
  DRCt_disc(y0,tseq=c(0,tmax),sf=1,parms=c(aparms,Rtype=Rtype))[2,2:3]
  }

nT <- 2
Tmin <- -2
Tmax <- 5
Tseq <- 10^seq(Tmin,Tmax,length.out=nT)
da <- array(dim=c(nT,4,2))

for(i in 1:nT){
  da[i,1,] <- DRCt_disc_one(y0,Tseq[i],Rtype="replenish")
  da[i,2,] <- DRCt_disc_one(y0,Tseq[i],Rtype="cohabit")
  da[i,3,] <- DRCt_disc_one(y0,Tseq[i],Rtype="remove")
  da[i,4,] <- DRCt_disc_one(y0,Tseq[i],Rtype="persist")
}

par(mfrow=c(1,2))
matplot(log10(Tseq),da[,,2],type="l",col=cols)
matplot(log10(Tseq),da[,,1],type="l",col=cols)
legend("topright",legend=c("replenish","cohabit","remove","persist"),
       col=cols,lty=1:4,bty="n")


