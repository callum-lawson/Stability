### Explore MacArthur R* models in C-space ### 

require("rootSolve")
require("deSolve")
require("phaseR")

# Functions ---------------------------------------------------------------
# i = resource index; j = consumer index
# a = v * a
# mu = V * mu

Rstar_f <- function(C,i,K,A){
  K[i] / (1 + A[i,1]*C[1] + A[i,2]*C[2])
}
  # abiotic

# Rstar_f <- function(C,i,K,A){
#   K[i] * (1 - A[i,1]*C[1] - A[i,2]*C[2])
# }
  # logistic

Cin_f <- function(C,j,K,A){
  A[1,j] * Rstar_f(C,i=1,K,A) + A[2,j] * Rstar_f(C,i=2,K,A) 
}

dC_f <- function(C,j,K,A,M){
  Cin_f(C,j,K,A) - M[j]
}

lndCdt_f <- function(t,y,parms){
  with(parms, {
    list( c( dC1 = dC_f(C=y,1,K,A,M), dC2 = dC_f(C=y,2,K,A,M) ) )
  })
}

# Parameters --------------------------------------------------------------

k1 <- 1
k2 <- 1
a11 <- 1
a12 <- 1
a21 <- 1
a22 <- 1
mu1 <- 1
mu2 <- 1

parms <- list(
  K = c(k1,k2),
  A = matrix(c(a11,a12,a21,a22),nr=2,nc=2),
  M = c(mu1,mu2)
)

# Integrate ---------------------------------------------------------------

lnC0 <- c(lnC1=0,lnC2=0)
tT <- 1000
tseq <- seq(0,100,length.out=tT)

lnCt <- ode(y=lnC0,times=tseq,func=lndCdt_f,parms=parms) 
matplot(lnCt[,1],lnCt[,-1],type="l")

# Eigenvalues -------------------------------------------------------------

equ <- steady(y = lnC0, times = c(0,Inf), func = lndCdt_f, parms = parms, method = "runsteady")$y

jac <- jacobian.full(y=equ,fun=lndCdt_f,parms=parms,time=0)
eig <- eigen(jac)
vec <- eig$vectors
val <- eig$values

inv <- solve(vec)
lnc0 <- lnC0 - equ
lnc0mat <- matrix(rep(inv %*% lnc0, each=tT),nr=tT,nc=2)
mmat <- lnc0mat * exp(outer(tseq,val,"*"))
lnnmat <- t(vec %*% t(mmat))

equdiff <- lnCt[,-1] - rep(equ,each=tT)
ihat <- vec %*% c(1,0)
jhat <- vec %*% c(0,1) 

ihatA <- ihat + equ
jhatA <- jhat + equ
ihatB <- -ihat + equ
jhatB <- -jhat + equ

# Phase space plots -------------------------------------------------------

cxlim <- c(-1,1)
cylim <- c(-1,1)

par(mfrow=c(1,1))
flowField(lndCdt_f,xlim=cxlim,ylim=cylim,parameters=parms,points=30,add=FALSE)
abline(v=equ[1])
abline(h=equ[2])
clines <- list()
clines[[1]] <- nullclines(lndCdt_f, xlim=cxlim, ylim=cylim,
  parameters=parms, points=100,col=rep("blue",2),add.legend=FALSE
)

points(lnC0[1],lnC0[2])
arrows(x0=equ[1],y0=equ[2],x1=ihatA[1],y1=ihatA[2],length=0,col="red",lty=2) # i-hat
arrows(x0=equ[1],y0=equ[2],x1=jhatA[1],y1=jhatA[2],length=0,col="red",lty=3) # j-hat
arrows(x0=equ[1],y0=equ[2],x1=ihatB[1],y1=ihatB[2],length=0,col="red",lty=2,angle=-180) # i-hat 2
arrows(x0=equ[1],y0=equ[2],x1=jhatB[1],y1=jhatB[2],length=0,col="red",lty=3,angle=-180) # j-hat 2

lines(x=lnCt[,"lnC1"],y=lnCt[,"lnC2"],col="purple",lwd=2)
lines(lnnmat+rep(equ,each=tT),col="red",lwd=2)
