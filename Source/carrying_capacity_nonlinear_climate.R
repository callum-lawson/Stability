#####################################################################################
# How do nonlinear population growth responses to climate affect carrying capacity? #
#####################################################################################

source("Source/royama_functions.R")

set.seed(5)

zmu <- c(1,2)
zsd <- c(0.25,0.25)
zdem <- zmu # demonstrated relationships
nz <- length(zmu)
nt <- 10^4
zsim <- matrix(NA,nr=nt,nc=nz)
zsim[] <- rnorm(nt*nz,rep(zmu,each=nt),sd=zsd)
# zsim[] <- filter(rnorm(nt*nz,rep(zmu,each=nt),sd=zsd), filter=rep(1,3), circular=TRUE)

myxmin <- 0.5 # 0
myxmax <- 5.5 # 7

cseq <- c("red","blue")
lseq <- 2:1

pars1 <- data.frame(b0=0,b1=1,b2=0,b3=-1/2,b4=0)
pars2 <- data.frame(b0=2,b1=-1/2,b2=0,b3=-5/4,b4=1/2)
pars3 <- data.frame(b0=-2,b1=2.5,b2=0,b3=1/4,b4=-1/2)
# DD WAS -1; interaction WAS -0.5

K1 <- Kcalc(z=zmu,pars=pars1)
K2 <- Kcalc(z=zmu,pars=pars2)
K3 <- Kcalc(z=zmu,pars=pars3)

### SIMS

xmat1 <- xsim(zsim,pars1,nt=nt)
xmat2 <- xsim(zsim,pars2,nt=nt)
xmat3 <- xsim(zsim,pars3,nt=nt)

### PLOTS

myxmin <- 0.5 # 0
myxmax <- 5.5 # 7

rplot(zdem,pars1,ylim=lineylim,sarr=sarr,yval=yval)
rplot(zdem,pars2,ylim=lineylim,sarr=sarr,yval=yval)
rplot(zdem,pars3,ylim=lineylim,sarr=sarr,yval=yval)
