### Calculate eigenvectors and eigenvalues for independently-growing species ###

# Trial runs --------------------------------------------------------------

require("rootSolve")
require("deSolve")
require("phaseR")

dNLV <- function(N1,N2,r,a,adash){
  N1 * ( r + a * N1 + adash * N2 )
}

dLV <- function(t,y,parms){
  with(parms, {
    N1 <- y[1]
    N2 <- y[2]
    dN1 <- dNLV(N1,N2,r1,a[1,1],a[2,1])
    dN2 <- dNLV(N2,N1,r2,a[2,2],a[1,2])
    list(c(dN1=dN1,dN2=dN2))
  })
}

N0 <- c(N1=1,N2=1)
tT <- 1000
tseq <- seq(0,100,length.out=tT)

r1 <- 1
r2 <- 1
a11 <- 0 # -1
a22 <- 0 # -1
a12 <- -0.1
a21 <- -0.1

a <- matrix(c(a11,a12,a21,a22),nr=2,nc=2)

parms <- list(r1=r1, r2=r2, a=a)
Nt <- ode(y=N0,times=tseq,func=dLV,parms=parms) 
matplot(Nt[,1],log(Nt[,-1]),type="l")

# Eigenvalues - continuous ------------------------------------------------

equ <- solve(-a,c(r1,r2))

jac <- jacobian.full(y=equ,fun=dLV,parms=parms,time=0)
eig <- eigen(jac)
vec <- eig$vectors
val <- eig$values

inv <- solve(vec)
n0 <- N0 - equ
n0mat <- matrix(rep(inv %*% n0, each=tT),nr=tT,nc=2)
mmat <- n0mat * exp(outer(tseq,val,"*"))
nmat <- t(vec %*% t(mmat))

equdiff <- Nt[,-1] - rep(equ,each=tT)
ihat <- vec %*% c(1,0)
jhat <- vec %*% c(0,1) 

ihatA <- ihat + equ
jhatA <- jhat + equ
ihatB <- -ihat + equ
jhatB <- -jhat + equ

# Phase space plots -------------------------------------------------------

cxlim <- c(0,20)
cylim <- c(0,20)

par(mfrow=c(1,1))
flowField(dLV,xlim=cxlim,ylim=cylim,parameters=parms,points=30,add=FALSE)
abline(v=equ[1])
abline(h=equ[2])
clines <- list()
clines[[1]] <- nullclines(dLV, xlim=cxlim, ylim=cylim,
  parameters=parms, points=100,col=rep("blue",2),add.legend=FALSE
)

points(N0[1],N0[2])
arrows(x0=equ[1],y0=equ[2],x1=ihatA[1],y1=ihatA[2],length=0,col="red",lty=2) # i-hat
arrows(x0=equ[1],y0=equ[2],x1=jhatA[1],y1=jhatA[2],length=0,col="red",lty=3) # j-hat
arrows(x0=equ[1],y0=equ[2],x1=ihatB[1],y1=ihatB[2],length=0,col="red",lty=2,angle=-180) # i-hat 2
arrows(x0=equ[1],y0=equ[2],x1=jhatB[1],y1=jhatB[2],length=0,col="red",lty=3,angle=-180) # j-hat 2

lines(x=Nt[,"N1"],y=Nt[,"N2"],col="purple",lwd=2)
lines(nmat+rep(equ,each=tT),col="red",lwd=2)
