### Dynamics of three- and four-species food chains ###

require(deSolve)
require(fields)
source("Source/predprey_functions_general.R")

# Input parameters --------------------------------------------------------

tmax <- 60^2 * 24 * 7 * 52 * 10 # maximum length of time in seconds
tf <- 10001
tseq <- seq(0,tmax,length.out=tf)

zmu <- 0 # =20°C 
zsig <- 5 # wave amplitude
zf <- 10^1 # wave frequency over whole time series
zl <- tmax/zf 

e0 <- c(
  m = 10^-6, # very high value relative to consumer
  k = 100,
  a = 6*10^-7, # 3.181989*10^-9, # estimated from data
  h = 60^2, # 0.61, # 1685.586,     # estimated from data
  mu = 2.689*10^-7, # 2.689*10^-6,
  alpha = 0.5,
  phi = 0.1 # relative death rate of eggs
)

e1 <- c(
  m = 0, # 0.639,
  k = 0, # -0.772,
  a = 0.5091663, # -0.03,   # estimated from data
  h = 0, # -1.9, # -0.19, # -0.4660012, # estimated from data
  mu = 0, # 0.639
  alpha = 0,
  phi = 0
)

e2 <- c(
  m = 0,
  k = 0, 
  a = 0, 
  h = 0, 
  mu = -1/4, # metabolic rate per unit mass (could also be -1/3)
  alpha = 0,
  phi = 0
)

M <- 10^c(0,2)

zparms <- list(zmu=zmu,zsig=zsig,zl=zl)
eparms <- list(M=M,e0=e0,e1=e1,e2=e2)

R_0 <- 10^1
C1_0 <- 10^1
C2_0 <- 10^1

y0 <- c(R=R_0,C1=C1_0,C2=C2_0)

# Increasing temperature variability --------------------------------------

### Consumer present

lvar <- ode(y=y0,times=tseq,func=dRCt3,parms=c(eparms,zl=zl,zmu=zmu,zsig=0))
mvar <- ode(y=y0,times=tseq,func=dRCt3,parms=c(eparms,zl=zl,zmu=zmu,zsig=3))
hvar <- ode(y=y0,times=tseq,func=dRCt3,parms=c(eparms,zl=zl,zmu=zmu,zsig=6))

par(mfrow=c(1,1))
matplot(tseq,log(hvar[,-1]),type="l",col="red",bty="n")
matplot(tseq,log(mvar[,-1]),type="l",col="black",add=TRUE)
matplot(tseq,log(lvar[,-1]),type="l",col="blue",add=TRUE)

start <- 1:round(99*tf/100,0)
matplot(tseq[-start],log(hvar[-start,2]),type="l",col="red",bty="n")

### Calculations

hiter <- log(hvar[-(1:round(tf/10,0)),-1])
liter <- log(lvar[-(1:round(tf/10,0)),-1])
# apply(hiter,2,range) - rep(apply(liter,2,median),each=2)
apply(hiter,2,median) - apply(liter,2,median)
  # might be biased at low wave frequencies because slice off part of wave

C1t <- hvar[,3]
C2t <- hvar[,4]

plot(diff(log(C1t[-(1:100)]))~log(C1t[-(1:100)][-1]))
plot(diff(log(C2t[-(1:100)]))~log(C2t[-(1:100)][-1]))
