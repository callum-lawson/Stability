#####################################################################################
# How do nonlinear population growth responses to climate affect carrying capacity? #
#####################################################################################

source("Source/gompertz_functions.R")

# b0 = intercept
# b1 = linear climate 	
# b2 = squared climate 	
# b3 = density 
# b4 = climate*density

myxmin <- -2
myxmax <- 2
  # for plotting

set.seed(5)
zmu <- c(0,0)
zsd <- c(0.25,0.50)
zdem <- zmu # demonstrated relationships
nz <- length(zmu)
nt <- 10^4
zsim <- matrix(NA,nr=nt,nc=nz)
zsim[] <- rnorm(nt*nz,rep(zmu,each=nt),sd=rep(zsd,each=nt))

# 0.5,0.75; -0.25,-1

# Two-climate comparisons - linear climate --------------------------------

xpars <- rbind(
  data.frame(b0=0,b1=0.5,b2=0,b3=-0.1,b4=0),
  data.frame(b0=0,b1=0.75,b2=0,b3=-0.75,b4=0)
)
# try intercept 0.1

xmat <- xsim(zsim,xpars,nt=nt,outmat=T)

rplot_3eg(zmu,zsd,xpars,xmin=-1,xmax=1)
abline(v=0,lty=3)
dplot(xmat,col=rep(c("blue","red"),each=2),lty=rep(1:2,times=2),bw=0.1)
# shows that same population fluctuations can be made with different
# resistance / resilience combos
# (because stronger climate response can subsitute for more variability)

# Two-climate comparisons - squared ---------------------------------------

source("Source/gompertz_functions.R")

xpars <- rbind(
  #data.frame(b0=0,b1=0.5,b2=0.75,b3=-0.1,b4=0),
  data.frame(b0=0,b1=0.5,b2=0.75,b3=-0.1,b4=0),
  #data.frame(b0=0,b1=0.5,b2=0.75,b3=-1.5,b4=0),
  data.frame(b0=0,b1=0.5,b2=0.75,b3=-0.5,b4=0)
)
# try intercept 0.1

xmat <- xsim(zsim,xpars,nt=nt,outmat=T)

rplot_3eg(zmu,zsd,xpars,xmin=-1,xmax=1)
abline(v=0,lty=3)
dplot(xmat,bw=0.1,col=rep(c("blue","red"),each=2),lty=rep(1:2,times=2),xmin=-1,xmax=3)
  # stronger effects of nonlinearity when DD weaker
  # but DD strength > 1, may help less to have nonlinear effects
  # when DD strength > 2, continually climbs to higher population sizes and never comes back

apply(xmat,2,median)
apply(xmat,2,mean)

Kcalc(z=zmu[1],pars=xpars)

# K1 <- Kcalc(z=zmu,pars=pars1)
# K2 <- Kcalc(z=zmu,pars=pars2)

xsdc <- apply(xmat,2,sd)
xsdc[2]/xsdc[1]
xsdc[4]/xsdc[3]
xsdc[6]/xsdc[5]
  # ratio between sd change seems identical

# Continuous range - squared climate --------------------------------------

nz <- 50
nt <- 10000
np <- 100

zsd_min <- 0
zsd_max <- 1
set.seed(5)
zmu <- rep(0,nz)
zsd <- seq(zsd_min,zsd_max,length.out=nz)
zsim <- matrix(NA,nr=nt,nc=nz)
# zsim[] <- rnorm(nt*nz,rep(zmu,each=nt),sd=rep(zsd,each=nt))
zsim[] <- outer(rnorm(nt,mean=zmu,sd=1), zsd, "*")
  # simulate same climate series, then scale according to sd

b3_min <- -0.05
b3_max <- -1.5
b3seq <- seq(b3_min,b3_max,length.out=np)
xpars <- data.frame(
  b0=rep(0,np),b1=rep(1,np),b2=rep(0.1,np),b3=b3seq,b4=rep(0,np)
)
xmat <- xsim(zsim,xpars,nt=nt)
xmed <- apply(xmat,c(2,3),median)

library(fields)
matplot(zsd,xmed,type="l",col=tim.colors(np),lty=1,xlab="zsd")
matplot(zsd,xmed,type="l",col=tim.colors(np),lty=1,ylim=c(0,0.1),xlab="zsd")
matplot(zsd,xmed,type="l",col=tim.colors(np),lty=1,ylim=c(-0.015,0.015),xlab="zsd")

  # weaker DD ->
  # - stronger influence of environment
  # - therefore, stronger influence of (non-linear effects of) variability
  # - therefore, also bigger effects of variability increase
  # - but negligible effects of climate non-linearity once DD becomes over-compensating
  # (dark red lines in same place as lighter red)
  # (median N can still increase with clim sd, but not much, and may decline)
  # - median X increases *exponentially* with zsd, suggesting even greater effects
  # than on population growth rate
  # - however, with overcompensating DD and nonlinear clim, zsd can actually reduce
  # median population size (massive increase -> disproportional decrease)
  # - with overcompensating DD, no clear ordering of DD strength effects on clim sd
  # effects (i.e., stronger DD doesn't mean weaker clim sd effects)

# Extinction risk

xtau <- apply(xmat,c(2,3),function(x,tau=-0.05) sum(x<tau)/length(x))
matplot(zsd,xtau,type="l",col=tim.colors(np),lty=1)
  # - non-linear effects of climate variability on extinction risk
  # - stronger density-dependence -> maximum extinction risk occurs at higher 
  # climate variabilty value
  # - extinction can be highest at intermediate zsd for weak DD 
  # (i.e., at some point, increasing zsd reduces extinction risk), 
  # but at max zsd for strong DD (i.e., zsd always increases extinction risk)

# Distributions
b3eg <- c(1,33,66,100)
b3seq[b3eg]
dplot(xmat[,50,b3eg],nz=1,np=2,bw=0.1,xmin=-5,xmax=7.5,lty=1,
      col=tim.colors(length(b3eg)))

# Example time series
matplot(xmat[9000:9900,50,b3eg],type="l",col=tim.colors(length(b3eg)),lty=1)
plot(xmat[9000:9900,50,min(b3eg)],type="l")
plot(xmat[9000:9900,50,max(b3eg)],type="l")
  # weak DD -> really autocorrelated

acf(xmat[9000:9900,50,b3eg[1]])
acf(xmat[9000:9900,50,b3eg[2]])
acf(xmat[9000:9900,50,b3eg[3]])
acf(xmat[9000:9900,50,b3eg[4]])
  # weak DD -> positive lag-1 ac; overcompensating -> negative lag-1 ac
