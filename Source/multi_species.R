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
  m = 10^-5, # very high value relative to consumer
  k = 100,
  a = 6*10^-7, # 3.181989*10^-9, # estimated from data
  h = 60^2, # 0.61, # 1685.586,     # estimated from data
  w = 60^2,
  mu = 2.689*10^-7, # 2.689*10^-6,
  alpha = 0.5,
  phi = 0.1 # relative death rate of eggs
)

e1 <- c(
  m = 0, # 0.639,
  k = 0, # -0.772,
  a = 0.5091663, # -0.03,   # estimated from data
  h = 0, # -1.9, # -0.19, # -0.4660012, # estimated from data
  w = 0,
  mu = 0, # 0.639
  alpha = 0,
  phi = 0
)

e2 <- c(
  m = 0,
  k = 0, 
  a = 0, 
  h = -0.15, # -0.15, # Rall et al 2012
  w = 0, # -0.15,
  mu = -1/4, # metabolic rate per unit mass (could also be -1/3)
  alpha = 0,
  phi = 0
)

M <- 10^c(-2,0)

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

matplot(tseq,log(hvar[,-1]),type="l",col="red",bty="n",ylim=c(3,4))
matplot(tseq,log(mvar[,-1]),type="l",col="black",add=TRUE)
matplot(tseq,log(lvar[,-1]),type="l",col="blue",add=TRUE)

# start <- 1:round(99*tf/100,0)
# matplot(tseq[-start],log(hvar[-start,2]),type="l",col="red",bty="n")

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

### Growth rate curves

library(reshape2)
C2seq <- exp(seq(-5,5,length.out=100))
zseq <- seq(zmu-5,zmu+5,length.out=10)
rd2 <- expand.grid(C2=C2seq,z=zseq)
Rstar2 <- with(rd2, t( Rstarcalc2(C2,z,eparms=eparms) ) )

plot(Rstar2)

parms2 <- eparms
parms2$M <- eparms$M[2]
rd2$r <- with(rd2, dCCdt(R=Rstar2[,"C1"],C=C2,z=z,eparms=parms2))

ra2 <- acast(melt(rd2,id=c("C2","z")),C2~z)
matplot(C2seq,ra2,type="l",lty=1,col=heat.colors(10))
abline(h=0,lty=3,col="grey")

# Three-species - discrete ------------------------------------------------

R_0 <- 10^1
C1_0 <- 10^1
C2_0 <- 10^1

y0 <- c(R=R_0,C1=C1_0,C2=C2_0)

sf <- 5000 + 1 # number of seasons (+1 because new season starts right at end)
sl <- tmax/sf # less than one day delay
sstart <- seq(0,tmax,length.out=sf)
sseq <- sseqgen(tseq,sstart)

discrete <- DRCt_disc3(y0,tseq,sf,parms=c(eparms,zl=zl,zmu=zmu,zsig=1))
matplot(tseq,log(discrete[,-1]),type="l",col="orange")

# Four-species ------------------------------------------------------------

R_0 <- 10^1
C1_0 <- 10^1
C2_0 <- 10^1
C3_0 <- 10^1

y0 <- c(R=R_0,C1=C1_0,C2=C2_0,C3=C3_0)
eparms$M <- 10^c(-2,0,2)

zf <- 10^3 # wave frequency over whole time series
zl <- tmax/zf 

### Dynamics

lvar4 <- ode(y=y0,times=tseq,func=dRCt4,parms=c(eparms,zl=zl,zmu=zmu,zsig=0))
mvar4 <- ode(y=y0,times=tseq,func=dRCt4,parms=c(eparms,zl=zl,zmu=zmu,zsig=3))
hvar4 <- ode(y=y0,times=tseq,func=dRCt4,parms=c(eparms,zl=zl,zmu=zmu,zsig=6))

# par(mfrow=c(1,1))
# matplot(tseq,log(hvar4[,-1]),type="l",col="red",bty="n")
# matplot(tseq,log(mvar4[,-1]),type="l",col="black",add=TRUE)
# matplot(tseq,log(lvar4[,-1]),type="l",col="blue",add=TRUE)

par(mfrow=c(1,1))
matplot(tseq,log(hvar4[,-1]),type="l",col="red",bty="n",ylim=c(4.5,5.5))
matplot(tseq,log(mvar4[,-1]),type="l",col="black",add=TRUE)
matplot(tseq,log(lvar4[,-1]),type="l",col="blue",add=TRUE)

### Growth rate curves

dCCdt <- Vectorize(
  function(R,C,z,eparms){
    parms <- with(eparms, as.list(c( R=R, arrrate(z,M,e0,e1,e2) )) )
    dCt_cons(t=0,y=C,parms=parms)/C
  }, 
  vectorize.args=c("R","C","z")
)

DCC <- Vectorize(
  function(R,C,z,eparms,sl){
    parms <- with(eparms, as.list( arrrate(z,M,e0,e1,e2) ) )
    N <- ode(y=c(R=R,C=C,E=0),
             times=c(0,sl),
             func=dRCt_disc_cons,
             parms=parms
    )
    N[2,"C"] + N[2,"E"]
  },
  vectorize.args=c("R","C","z")
)

eparms$e0["phi"] <- 0

rd2 <- rd
rd2$r <- with(rd2, DCC(Rstar,C,z,eparms,sl*100)/C)

ra2 <- acast(melt(rd2,id=c("C","z")),C~z)
par(mfrow=c(1,1))
matplot(log(Cseq),log(ra2),type="l",lty=1,col=heat.colors(10))
abline(h=0,lty=3,col="grey")

