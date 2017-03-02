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

xpars <- rbind(
  data.frame(b0=0,b1=0.5,b2=0.75,b3=-0.25,b4=0),
  data.frame(b0=0,b1=0.75,b2=0.75,b3=-1,b4=0)
)
# try intercept 0.1

xmat <- xsim(zsim,xpars,nt=nt)

rplot_3eg(zmu,zsd,xpars,xmin=-1,xmax=1)
dplot(xmat,col=rep(c("blue","red"),each=2),lty=rep(1:2,times=2),bw=0.1)
  # shows that same population fluctuations can be made with different
  # resistance / resilience combos

apply(xmat,2,median)
apply(xmat,2,mean)

# K1 <- Kcalc(z=zmu,pars=pars1)
# K2 <- Kcalc(z=zmu,pars=pars2)

xsdc <- apply(xmatc,2,sd)
xsdc[2]/xsdc[1]
xsdc[4]/xsdc[3]
xsdc[6]/xsdc[5]
  # ratio between sd change seems identical

