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
mycols <- c("red","blue")

zmu <- c(0,0)
zsd <- c(0.25,0.50)
zdem <- zmu # demonstrated relationships
nz <- length(zmu)
nt <- 10^4
zsim <- matrix(NA,nr=nt,nc=nz)
zsim[] <- rnorm(nt*nz,rep(zmu,each=nt),sd=rep(zsd,each=nt))
# zsim[] <- filter(rnorm(nt*nz,rep(zmu,each=nt),sd=zsd), filter=rep(1,3), circular=TRUE)

# 0.5,0.75; -0.25,-1


# Linear climate response -------------------------------------------------

xpars <- rbind(
  data.frame(b0=0,b1=0.5,b2=0,b3=-0.1,b4=0),
  data.frame(b0=0,b1=0.75,b2=0,b3=-0.75,b4=0)
)
# try intercept 0.1

xmat <- xsim(zsim,xpars,nt=nt)

rplot_3eg(zmu,zsd,xpars,xmin=-1,xmax=1)
abline(v=0,lty=3)
dplot(xmat,col=rep(c("blue","red"),each=2),lty=rep(1:2,times=2),bw=0.1)
# shows that same population fluctuations can be made with different
# resistance / resilience combos
# (because stronger climate response can subsitute for more variability)

# Squared climate response ------------------------------------------------

source("Source/gompertz_functions.R")

xpars <- rbind(
  #data.frame(b0=0,b1=0.5,b2=0.75,b3=-0.1,b4=0),
  data.frame(b0=0,b1=0.5,b2=0.75,b3=-0.75,b4=0),
  #data.frame(b0=0,b1=0.5,b2=0.75,b3=-1.5,b4=0),
  data.frame(b0=0,b1=0.5,b2=0.75,b3=-2,b4=0)
)
# try intercept 0.1

xmat <- xsim(zsim,xpars,nt=nt)

rplot_3eg(zmu,zsd,xpars,xmin=-1,xmax=1)
abline(v=0,lty=3)
dplot(xmat,bw=0.1)
  # stronger effects of nonlinearity when DD weaker
  # but DD strength > 1, may help less to have nonlinear effects
  # when DD strength > 2, continually climbs to higher population sizes and never comes back

apply(xmat,2,median)
apply(xmat,2,mean)

Kcalc(z=zmu[1],pars=xpars)

# K1 <- Kcalc(z=zmu,pars=pars1)
# K2 <- Kcalc(z=zmu,pars=pars2)

xsdc <- apply(xmatc,2,sd)
xsdc[2]/xsdc[1]
xsdc[4]/xsdc[3]
xsdc[6]/xsdc[5]
  # ratio between sd change seems identical

# Simulations - squared climate response ----------------------------------

nb3 <- 100
b3seq <- seq(-0.1,-2,length.out=nb3)
xpars <- data.frame(
  b0=rep(0,nb3),b1=rep(0.5,nb3),b2=rep(0.75,nb3),b3=b3seq,b4=rep(0,nb3)
)
xmat <- xsim(zsim,xpars,nt=nt)
xmed <- apply(xmat,2,median)

xmed1 <- xmed[seq(2,200,by=2)]
xmed0 <- xmed[seq(1,199,by=2)]

matplot(b3seq,cbind(xmed0,xmed1),type="l")
  # weaker DD ->
  # - stronger influence of environment
  # - therefore, stronger influence of (non-linear effects of) variability
  # - therefore, also bigger effects of variability increase

# Generalise plotting function (via dimensions) to extra climate combinations?

### Try with other (nonlinear) DD forms



