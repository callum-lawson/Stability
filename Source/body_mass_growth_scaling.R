require("rootSolve")
require("deSolve")
require("phaseR")

source("Source/consumer_resource_functions.R")
bhat <- readRDS("Output/rate_parameters_marginal_brms_06Sep2018.rds")

a <- 10^bhat$mu$bm

bgrow <- function(x,a,eps,u){
  2/(1/a + x) - a * exp(eps*(1-u)) + (u * eps)/x
}

Nseq <- exp(seq(-1,1,length.out=100))

aseq <- c(1,a)
epseq <- c(-0.1,0,0.1)

small1 <- bgrow(Nseq,a=aseq[1],eps=epseq[1],u=0)
small2 <- bgrow(Nseq,a=aseq[1],eps=epseq[2],u=0)
small3 <- bgrow(Nseq,a=aseq[1],eps=epseq[3],u=0)

big1 <- bgrow(Nseq,a=aseq[2],eps=epseq[1],u=0)
big2 <- bgrow(Nseq,a=aseq[2],eps=epseq[2],u=0)
big3 <- bgrow(Nseq,a=aseq[2],eps=epseq[3],u=0)

matplot(log(Nseq),cbind(small1,small2,small3,big1,big2,big3),type="l",
  lty=rep(c(2,1,2),times=2),col=rep(1:2,each=3))
abline(h=0,col="grey",lty=3)

small1b <- bgrow(Nseq,a=aseq[1],eps=epseq[1],u=1)
small2b <- bgrow(Nseq,a=aseq[1],eps=epseq[2],u=1)
small3b <- bgrow(Nseq,a=aseq[1],eps=epseq[3],u=1)

big1b <- bgrow(Nseq,a=aseq[2],eps=epseq[1],u=1)
big2b <- bgrow(Nseq,a=aseq[2],eps=epseq[2],u=1)
big3b <- bgrow(Nseq,a=aseq[2],eps=epseq[3],u=1)

matplot(log(Nseq),cbind(small1b,small2b,small3b,big1b,big2b,big3b),type="l",
        lty=rep(c(2,1,2),times=2),col=rep(1:2,each=3))
abline(h=0,col="grey",lty=3)
abline(h=0.1,col="grey",lty=3)
abline(h=-0.1,col="grey",lty=3)

# Integrate to get distributions ------------------------------------------

Cdt_f <- function(t,y,parms){
  with(parms, {
    eps <- sigma * sin(t*omega)
    dN <- bgrow(x=y,a,eps,u) * y
    list(c(dN = dN))
  })
}

freq <- 1/10

parms <- list(sigma = 0.1, omega = 2*pi*freq) # waves per unit time

parms0 <- c(parms,u=0)
parms1 <- c(parms,u=1)

# Integrate ---------------------------------------------------------------

C0 <- 1
tT <- 10^4
tseq <- seq(0,10^4,length.out=tT)

Ct1a <- ode(y=C0,times=tseq,func=Cdt_f,parms=c(parms0,a=aseq[1]))[,2]
Ct2a <- ode(y=C0,times=tseq,func=Cdt_f,parms=c(parms0,a=aseq[2]))[,2] 
Ct1b <- ode(y=C0,times=tseq,func=Cdt_f,parms=c(parms1,a=aseq[1]))[,2]
Ct2b <- ode(y=C0,times=tseq,func=Cdt_f,parms=c(parms1,a=aseq[2]))[,2] 

plot(Ct1a,Ct2a)
plot(log(Ct1a),log(Ct2a))
matplot(tseq,cbind(Ct1a,Ct2a),type="l")
matplot(tseq,log(cbind(Ct1a,Ct2a)),type="l")

sd(Ct1a[-(1:100)])/sd(Ct2a[-(1:100)])
sd(log(Ct1a[-(1:100)]))/sd(log(Ct2a[-(1:100)]))

plot(Ct1b,Ct2b)
plot(log(Ct1b),log(Ct2b))
matplot(tseq,cbind(Ct1b,Ct2b),type="l")
matplot(tseq,log(cbind(Ct1b,Ct2b)),type="l")

sd(Ct1b[-(1:100)])/sd(Ct2b[-(1:100)])
sd(log(Ct1b[-(1:100)]))/sd(log(Ct2b[-(1:100)]))

# for environmental and immigration, smaller fluctuates less in absolute terms but more in log terms

# Small fluctuations ------------------------------------------------------

### Proportional

par(mfrow=c(1,2))
curve(bgrow(exp(x),a=aseq[1],eps=epseq[2]*0.1,u=0),
  xlim=c(-0.05,0.05),ylim=c(-0.02,0.02))
curve(bgrow(exp(x),a=aseq[1],eps=epseq[1]*0.1,u=0),add=T,lty=2)
curve(bgrow(exp(x),a=aseq[1],eps=epseq[3]*0.1,u=0),add=T,lty=2)
abline(h=0,col="grey",lty=3)

curve(bgrow(exp(x),a=aseq[2],eps=epseq[2]*0.1,u=0),
  xlim=c(0.7,0.8),ylim=c(-0.02,0.02),col="red")
curve(bgrow(exp(x),a=aseq[2],eps=epseq[1]*0.1,u=0),add=T,lty=2,col="red")
curve(bgrow(exp(x),a=aseq[2],eps=epseq[3]*0.1,u=0),add=T,lty=2,col="red")
abline(h=0,col="grey",lty=3)

  # immigration effects probably complex - depend on where k is

