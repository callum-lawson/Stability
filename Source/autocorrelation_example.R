# Simple simulation of population fluctuations under autocorrelated climates #

source("Source/gompertz_functions.R")

mycols <- rep(c("blue","red"),each=2)
mylty <- rep(1:2,times=2)

nt <- 1100
zvals <- c(-1,1)
z1 <- rep(rep(zvals,each=1),length.out=nt)
z2 <- rep(rep(zvals,each=5),length.out=nt)
zsim <- cbind(z1,z2)

# Same K (different resistance) -------------------------------------------

xpars <- rbind(
  data.frame(b0=0,b1=3,b2=0,b3=-0.75,b4=0,b5=0),
  data.frame(b0=0,b1=1,b2=0,b3=-0.25,b4=0,b5=0)
)

rplot_3eg(0,1,xpars,xmin=-10,xmax=10,averages=TRUE)
Kcalc(z=zvals[1],pars=xpars)
Kcalc(z=zvals[2],pars=xpars)
# keeps K constant for all values of z

xmat <- xsim(zsim,xpars,nt=nt,outmat=T,warmup=100,lN0=0)
matplot(xmat[990:1000,],type="l",col=mycols,lty=mylty,ylab="ln N")
dplot(xmat,bw=0.001,col=mycols,lty=mylty,xmin=-3,xmax=3)

# Controlling to keep Ks the same (stronger DD -> weaker resistance):
# - Stronger DD creates bigger range of population sizes regardless of 
# climate autocorrelation

# Different K (same resistance) -------------------------------------------

xpars <- rbind(
  data.frame(b0=0,b1=1,b2=0,b3=-1,b4=0,b5=0),
  data.frame(b0=0,b1=1,b2=0,b3=-0.1,b4=0,b5=0)
)

rplot_3eg(0,1,xpars,xmin=-1,xmax=1)
Kcalc(z=zvals[1],pars=xpars)
Kcalc(z=zvals[2],pars=xpars)

xmat <- xsim(zsim,xpars,nt=nt,outmat=T,warmup=100,lN0=0)
matplot(xmat[990:1000,],type="l",col=mycols,lty=mylty,ylab="ln N")
dplot(xmat,bw=0.001,col=mycols,lty=mylty,xmin=-3,xmax=3)

# Different Ks (same resistance in both cases):
# - When clim is negatively autocorrelated, stronger DD creates bigger range of
# population sizes (pop catches up with new K immediately, so reflects full 
# variation in clim)
# - But when clim is positively autocorrelated, can have bigger variation in 
# population sizes with weak DD, because K varies more in that case

# Simulations without autocorrelation -------------------------------------

xpars <- rbind(
  data.frame(b0=0,b1=1,b2=0,b3=-1,b4=0,b5=0),
  data.frame(b0=0,b1=1,b2=0,b3=-0.1,b4=0,b5=0)
)

set.seed(1)
zmu <- 0
zsd <- 1
nz <- length(zmu)
nt <- 10^4
zsim <- matrix(NA,nr=nt,nc=nz)
zsim[] <- rnorm(nt*nz,rep(zmu,each=nt),sd=rep(zsd,each=nt))

xmat <- xsim(zsim,xpars,lN0=0,warmup=1000,nt=nt,outmat=T)
matplot(xmat[990:1000,],type="l",col=c("blue","red"),ylab="ln N")
dplot(xmat,bw=0.1,col=c("blue","red"),xmin=-10,xmax=10)

# Over-compensation -------------------------------------------------------

xpars <- rbind(
  data.frame(b0=0,b1=1,b2=0,b3=-1.9,b4=0,b5=0),
  data.frame(b0=0,b1=1,b2=0,b3=-1.1,b4=0,b5=0)
)

set.seed(1)
zmu <- 0
zsd <- 1
nz <- length(zmu)
nt <- 10^4
zsim <- matrix(NA,nr=nt,nc=nz)
zsim[] <- rnorm(nt*nz,rep(zmu,each=nt),sd=rep(zsd,each=nt))

xmat <- xsim(zsim,xpars,nt=nt,outmat=T)
matplot(xmat[990:1000,],type="l",col=c("blue","red"),ylab="ln N")
dplot(xmat,bw=0.1,col=c("blue","red"),xmin=-10,xmax=10)

# Autocorrelation vs DD strength plots ------------------------------------

acts <- function(alpha,nt){
  y <- arima.sim(model=list(ar=alpha),n=nt)
  sy <- scale(y)
  return(sy)
}

nac <- 11
acseq <- seq(-0.99,0.99,length.out=nac)

np <- 50 
dmin <- -0.1
dmax <- -1.99
dseq <- seq(dmin,dmax,length.out=np)
xpars <- data.frame(
  b0=rep(0,np),b1=rep(1,np),b2=rep(0,np),b3=dseq,b4=rep(0,np),b5=rep(0,np)
)

nt <- 10^5
zsim <- matrix(NA,nr=nt,nc=nac)
for(i in 1:nac){
  if(acseq[i]!=0) zsim[,i] <- acts(acseq[i],nt)
  else zsim[,i] <- rnorm(nt,0,1)
}
xmat <- xsim(zsim,xpars,nt=nt,lN0=0,warmup=10^4)

xmed <- t(apply(xmat,c(2,3),median))
xsd <- t(apply(xmat,c(2,3),sd))

library(fields)
matplot(dseq,xmed,type="l",col=tim.colors(nac),lty=1,xlab="DD")
matplot(dseq,log10(xsd),type="l",col=tim.colors(nac),lty=1,xlab="DD")
plot(dseq,log10(xsd[,6]),type="l",lty=1,xlab="DD")

# Autocorrelation vs DD strength: same K ----------------------------------

xpars <- data.frame(
  b0=rep(0,np),b1=-dseq,b2=rep(0,np),b3=dseq,b4=rep(0,np),b5=rep(0,np)
  )
# -b1 = a*b31, where a=K/z

xmat <- xsim(zsim,xpars,nt=nt,lN0=0,warmup=10^4)

xmed <- t(apply(xmat,c(2,3),median))
xsd <- t(apply(xmat,c(2,3),sd))

matplot(dseq,xmed,type="l",col=tim.colors(nac),lty=1,xlab="DD")
matplot(dseq,log10(xsd),type="l",col=tim.colors(nac),lty=1,xlab="DD")
plot(dseq,log10(xsd[,6]),type="l",lty=1,xlab="DD")
  # better to be resistant than resilient, especially in neg-autocor envs

# Nonlinear climate effects -----------------------------------------------

xpars <- rbind(
  data.frame(b0=0,b1=-2,b2=1,b3=-1,b4=0,b5=0),
  data.frame(b0=0,b1=-2,b2=1,b3=-0.5,b4=0,b5=0)
)
rplot_3eg(0,1,xpars,xmin=-5,xmax=5,averages=TRUE)
  # linear DD -> 
  # - effects of clim var on mean X are the same regardless of fluctuation speed
  # - but the slower the fluctuation speed, the bigger the fluctuations in X 
  # (not shown)

xpars <- rbind(
  data.frame(b0=0,b1=2,b2=1,b3=-1,b4=0,b5=0.25),
  data.frame(b0=0,b1=2,b2=1,b3=-0.5,b4=0,b5=0.25)
)
rplot_3eg(0,1,xpars,xmin=-5,xmax=5,averages=TRUE)
  # increasing DD strength with clim -> 
  # clim var has stronger effects (pos or neg) when fluctuations are fast
  # (opposite is true for decreasing DD strength with clim)

xpars <- rbind(
  data.frame(b0=0,b1=1,b2=0,b3=-1,b4=-0.1,b5=0)
)
rplot_3eg(0,1,xpars,xmin=-5,xmax=5,averages=TRUE)
  # with non-linear DD, can't calculate K distribution 
  # because there are always some climates where pop goes extinct
