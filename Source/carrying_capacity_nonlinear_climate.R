#####################################################################################
# How do nonlinear population growth responses to climate affect carrying capacity? #
#####################################################################################

source("Source/royama_functions.R")

# b0 = intercept
# b1 = linear climate 	
# b2 = squared climate 	
# b3 = density 
# b4 = climate*density

set.seed(5)

zmu <- c(0,0)
zsd <- c(0.25,0.50)
zdem <- zmu # demonstrated relationships
nz <- length(zmu)
nt <- 10^4
zsim <- matrix(NA,nr=nt,nc=nz)
zsim[] <- rnorm(nt*nz,rep(zmu,each=nt),sd=rep(zsd,each=nt))
# zsim[] <- filter(rnorm(nt*nz,rep(zmu,each=nt),sd=zsd), filter=rep(1,3), circular=TRUE)

myxmin <- 0.5 # 0
myxmax <- 5.5 # 7

cseq <- c("red","blue")
lseq <- 2:1

pars1 <- data.frame(b0=0,b1=0.5,b2=0,b3=-0.25,b4=0)
pars2 <- data.frame(b0=0,b1=0.75,b2=0,b3=-1,b4=0)
pars3 <- data.frame(b0=0,b1=1,b2=0,b3=-1.75,b4=0)
# try intercept 0.1

# K1 <- Kcalc(z=zmu,pars=pars1)
# K2 <- Kcalc(z=zmu,pars=pars2)

### SIMS

xmat1 <- xsim(zsim,pars1,nt=nt)
xmat2 <- xsim(zsim,pars2,nt=nt)
xmat3 <- xsim(zsim,pars3,nt=nt)

xmatc <- cbind(xmat1,xmat2,xmat3)
xdensc <- apply(xmatc,2,density)

source("Source/royama_functions.R")
rplot_3eg(zmu,zsd,rbind(pars1,pars2,pars3),xmin=-1,xmax=1)
abline(v=0,lty=3)

plot(xdensc[[1]],col="blue")
lines(xdensc[[2]],lty=2,col="blue")
lines(xdensc[[3]],col="red")
lines(xdensc[[4]],lty=2,col="red")
lines(xdensc[[5]],col="black")
lines(xdensc[[6]],lty=2,col="black")
  # shows that same population fluctuations can be made with different
  # resistance / resilience combos

xsdc <- apply(xmatc,2,sd)
xsdc[2]/xsdc[1]
xsdc[4]/xsdc[3]
  # ratio between sd change seems identical

### PLOTS

myxmin <- 0.5 # 0
myxmax <- 5.5 # 7

rplot(zdem,pars1,ylim=lineylim,sarr=sarr,yval=yval)
