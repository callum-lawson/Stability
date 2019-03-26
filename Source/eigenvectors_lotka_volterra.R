### Calculate eigenvectors and eigenvalues for Lotka-Volterra interaction models ###

require("rootSolve")
require("deSolve")
require("phaseR")

# Setup -------------------------------------------------------------------

dN_Ndt <- function(N1,N2,r,a,adash){
  r + a * N1 + adash * N2
}
  # r + exp(N1 + a_e) + exp(N2 + adash_e)
  # exp makes it look like Ricker

dlnLV <- function(t,y,parms){
  with(parms, {
    N <- exp(y)
    N1 <- N[1]
    N2 <- N[2]
    dlnN1 <- dN_Ndt(N1,N2,r1,A[1,1],A[2,1])
    dlnN2 <- dN_Ndt(N2,N1,r2,A[2,2],A[1,2])
    list(c(dlnN1=dlnN1,dlnN2=dlnN2))
  })
}

lnN0 <- c(lnN1=0,lnN2=0)
tT <- 1000
tseq <- seq(0,100,length.out=tT)

### Competition
v <- 100
r1 <- v*1
r2 <- 1
a11 <- v*-1
a22 <- -1
a12 <- v*-0.1
a21 <- -0.1

### Predation
# r1 <- 1
# r2 <- -0.01
# a11 <- -1
# a22 <- -0.1
# a12 <- -1
# a21 <- 1

A <- matrix(c(a11,a12,a21,a22),nr=2,nc=2)

parms <- list(r1=r1, r2=r2, A=A)
lnNt <- ode(y=lnN0,times=tseq,func=dlnLV,parms=parms) 
matplot(lnNt[,1],lnNt[,-1],type="l")

# Eigenvalues -------------------------------------------------------------

equ <- log(solve(-t(A),c(r1,r2)))

jac <- jacobian.full(y=equ,fun=dlnLV,parms=parms,time=0)
eig <- eigen(jac)
vec <- eig$vectors
val <- eig$values

inv <- solve(vec)
lnn0 <- lnN0 - equ
lnn0mat <- matrix(rep(inv %*% lnn0, each=tT),nr=tT,nc=2)
mmat <- lnn0mat * exp(outer(tseq,val,"*"))
lnnmat <- t(vec %*% t(mmat))

equdiff <- lnNt[,-1] - rep(equ,each=tT)
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
flowField(dlnLV,xlim=cxlim,ylim=cylim,parameters=parms,points=30,add=FALSE)
abline(v=equ[1])
abline(h=equ[2])
clines <- list()
clines[[1]] <- nullclines(dlnLV, xlim=cxlim, ylim=cylim,
  parameters=parms, points=100,col=rep("blue",2),add.legend=FALSE
)

points(lnN0[1],lnN0[2])
arrows(x0=equ[1],y0=equ[2],x1=ihatA[1],y1=ihatA[2],length=0,col="red",lty=2) # i-hat
arrows(x0=equ[1],y0=equ[2],x1=jhatA[1],y1=jhatA[2],length=0,col="red",lty=3) # j-hat
arrows(x0=equ[1],y0=equ[2],x1=ihatB[1],y1=ihatB[2],length=0,col="red",lty=2,angle=-180) # i-hat 2
arrows(x0=equ[1],y0=equ[2],x1=jhatB[1],y1=jhatB[2],length=0,col="red",lty=3,angle=-180) # j-hat 2

lines(x=lnNt[,"lnN1"],y=lnNt[,"lnN2"],col="purple",lwd=2)
lines(lnnmat+rep(equ,each=tT),col="red",lwd=2)
