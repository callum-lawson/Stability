### Compare effects of time lags on consumer-resource dynamics ###

### TODO
# - add continuous-time lag in resource as well as consumer births

require(deSolve)
require(fields)
source("Source/predprey_functions_general.R")

x0 <- arrrate(0,e0["x"],e1["x"])

tmax <- 60^2 * 24 * 365 # maximum length of time in seconds
tf <- 10^4 # 10^3
tseq <- seq(0,tmax,length.out=tf)

zmu <- 0
zsig <- 0 # wave amplitude
zf <- 10 # wave frequency over whole time series
zl <- tmax/zf 

sf <- 10^4 # 365
  # number of seasons - must be <= number of recorded timesteps
sl <- tmax/sf

zparms <- list(zmu=zmu,zsig=zsig,zl=zl,tau=sl)
eparms <- list(e0=e0,e1=e1,omega=1,kappa=0)

R0 <- 10^1
C0 <- 10^-2
y0 <- c(R=R0,C=C0)

standard <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(zparms,eparms))
delayed <- dede(y=y0,times=tseq,func=dRCt_delay,parms=c(zparms,eparms))
discrete <- DRCt_disc(y0,tseq,sf,parms=c(zparms,eparms))

par(mfrow=c(1,2),mar=c(4,4,2,2))
require(fields)
cols <- tim.colors(1)
allcols <- c("black","gray",cols)
matplot(tseq,log(cbind(standard[,2],delayed[,2],discrete[,2])),
        type="l",
        col=allcols,
        lty=1:length(allcols),
        bty="n"
)
matplot(tseq,log(cbind(standard[,3],delayed[,3],discrete[,3])),
        type="l",
        col=allcols,
        lty=1:length(allcols),
        bty="n"
)
legend("topright",legend=c("standard","delayed","discrete"),
       col=allcols,lty=1:length(allcols),bty="n")

eparms2 <- list(e0=replace(e0,5,sl),e1=e1,omega=1)
standard2 <- ode(y=y0,times=tseq,func=dRCt_cont,parms=c(zparms,eparms2))
lines(tseq,log(standard2[,3]),col="pink")

# Effects of season length ------------------------------------------------

DRCt_disc_one <- function(y0,tmax){
  DRCt_disc(y0,tseq=c(0,tmax),sf=1,parms=c(zparms,eparms))[2,2:3]
  }

nT <- 20
Tmin <- -2
Tmax <- 5
Tseq <- 10^seq(Tmin,Tmax,length.out=nT)
da <- array(dim=c(nT,4,2))

for(i in 1:nT){
  da[i,1,] <- DRCt_disc_one(y0,Tseq[i])
  # da[i,2,] <- DRCt_disc_one(y0,Tseq[i])
}

par(mfrow=c(1,2))
matplot(log10(Tseq),da[,,1],type="l",col=cols,ylab="prey growth")
matplot(log10(Tseq),da[,,2],type="l",col=cols,ylab="pred growth")
