##############################################################################
# Simple simulation of population fluctuations under autocorrelated climates #
##############################################################################

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
  data.frame(b0=0,b1=3,b2=0,b3=-0.75,b4=0),
  data.frame(b0=0,b1=1,b2=0,b3=-0.25,b4=0)
)

rplot_3eg(0,1,xpars,xmin=-1,xmax=1)
Kcalc(z=zvals[1],pars=xpars)
Kcalc(z=zvals[2],pars=xpars)

xmat <- xsim(zsim,xpars,nt=nt,outmat=T,warmup=100,lN0=0)
matplot(xmat[990:1000,],type="l",col=mycols,lty=mylty,ylab="ln N")
dplot(xmat,bw=0.001,col=mycols,lty=mylty,xmin=-3,xmax=3)

# Controlling to keep Ks the same (stronger DD -> weaker resistance):
# - Stronger DD creates bigger range of population sizes regardless of 
# climate autocorrelation

# Different K (same resistance) -------------------------------------------

xpars <- rbind(
  data.frame(b0=0,b1=1,b2=0,b3=-0.75,b4=0),
  data.frame(b0=0,b1=1,b2=0,b3=-0.25,b4=0)
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

