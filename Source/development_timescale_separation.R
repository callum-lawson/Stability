### Explore consumer development models with timescale separation of resource ### 

require("rootSolve")
require("deSolve")
require("phaseR")

# Functions ---------------------------------------------------------------
# a = v * a

Rstar_f <- function(C,k,a){
  k / (1 + a*C)
}
  # abiotic
 
# Rstar_f <- function(C,k,a){
#   k * (1 - a*C)
# }
  # logistic

dlnJ_f <- function(C,k,a,mu,phi){
  a * Rstar_f(C,k,a) - (phi + mu)
}

dlnC_f <- function(J,C,mu,phi){
  phi * J/C - mu
}

dlnNdt_f <- function(t,y,parms){
  with(parms, {
    N <- exp(y)
    dlnJ <- dlnJ_f(C=N[2],k,a,mu,phi)
    dlnC <- dlnC_f(J=N[1],C=N[2],mu,phi)
    list( c( dlnJ = dlnJ, dlnC = dlnC ) )
  })
}

# Parameters --------------------------------------------------------------

k <- 1
a <- 2
mu <- 0.1
phi <- 1

parms <- list(k=k,a=a,mu=mu,phi=phi)

# Integrate ---------------------------------------------------------------

equ <- steady(y = c(lnJ=0,lnC=0), times = c(0,Inf), func = dlnNdt_f, parms = parms, method = "runsteady")$y

lnN0 <- equ + log(1.5) # 1.5 times above N*
tT <- 1000
tseq <- seq(0,100,length.out=tT)

lnNt <- ode(y=lnN0,times=tseq,func=dlnNdt_f,parms=parms) 
# matplot(lnNt[,1],lnNt[,-1],type="l")

# Eigenvalues -------------------------------------------------------------

jac <- jacobian.full(y=equ,fun=dlnNdt_f,parms=parms,time=0)
eig <- eigen(jac)
vec <- eig$vectors
val <- eig$values

vecdash <- cbind(Re(vec[,1]),Im(vec[,1]))
  # Rotation-scaling theorem

inv <- solve(vecdash)
# lnn0 <- lnN0 - equ
# lnn0mat <- matrix(rep(inv %*% lnn0, each=tT),nr=tT,nc=2)
# mmat <- lnn0mat * exp(outer(tseq,val,"*"))
# lnnmat <- t(vecdash %*% t(mmat))

ihat <- vecdash %*% c(1,0)
jhat <- vecdash %*% c(0,1)

ihatA <- ihat + equ
jhatA <- jhat + equ
ihatB <- -ihat + equ
jhatB <- -jhat + equ

# Phase space plots -------------------------------------------------------

cxlim <- c(-4,-2.25) # c(equ[1]-1,equ[1]+1)
cylim <- c(-1.5,-0.25) # c(equ[2]-1,equ[2]+1)

par(mfrow=c(1,1))
flowField(dlnNdt_f,xlim=cxlim,ylim=cylim,parameters=parms,points=30,add=FALSE)
abline(v=equ[1])
abline(h=equ[2])
clines <- list()
clines[[1]] <- nullclines(dlnNdt_f, xlim=cxlim, ylim=cylim,
  parameters=parms, points=100,col=rep("blue",2),add.legend=FALSE
)

points(lnN0[1],lnN0[2])
arrows(x0=equ[1],y0=equ[2],x1=ihatA[1],y1=ihatA[2],length=0.25,col="red") # i-hat
arrows(x0=equ[1],y0=equ[2],x1=jhatA[1],y1=jhatA[2],length=0.25,col="red") # j-hat
arrows(x0=equ[1],y0=equ[2],x1=ihatB[1],y1=ihatB[2],length=0.25,col="red") # i-hat 2
arrows(x0=equ[1],y0=equ[2],x1=jhatB[1],y1=jhatB[2],length=0.25,col="red") # j-hat 2

points(x=lnNt[,"lnJ"],y=lnNt[,"lnC"],col="purple",type="b",lwd=1)
#points(lnnmat+rep(equ,each=tT),col="red",type="b")

# Perturbation ------------------------------------------------------------

parms2 <- parms
parms2$mu <- parms$mu + 0.1
clines[[2]] <- nullclines(dlnNdt_f, xlim=cxlim, ylim=cylim,
  parameters=parms2, points=100,col=rep("orange",2),add.legend=FALSE
)

### Back

equ2 <- steady(y = c(lnJ=0,lnC=0), times = c(0,Inf), func = dlnNdt_f, parms = parms2, method = "runsteady")$y
lnNt2 <- ode(y=equ2,times=tseq,func=dlnNdt_f,parms=parms) 
points(x=lnNt2[,"lnJ"],y=lnNt2[,"lnC"],col="orange",type="b")

### There

lnNt2b <- ode(y=equ,times=tseq,func=dlnNdt_f,parms=parms2) 
points(x=lnNt2b[,"lnJ"],y=lnNt2b[,"lnC"],col="orange",type="b")

### Equilibria

# equ_all <- rbind(equ,equ2,equ3,equ4)
# points(x=equ_all[,"lnJ"],y=equ_all[,"lnC"],col="black",pch="+",cex=1.1)

# Eigenspace --------------------------------------------------------------

# equvec <- rep(equ,each=tT)
# tran <- t(inv %*% t(lnNt[,-1] - equvec))
# tran2 <- t(inv %*% t(lnNt2[,-1] - equvec))
# tran3 <- t(inv %*% t(lnNt3[,-1] - equvec))
# tran4 <- t(inv %*% t(lnNt4[,-1] - equvec))
# 
# dtran <- rbind(tran2,tran3,tran4)
# plot(dtran,col=rep(c("orange","green","brown"),each=tT),pch=16)

# Simulations with sin waves ----------------------------------------------

# dlnNdt_f2 <- function(t,y,parms){
#   with(parms, {
#     C <- exp(y)
#     eps1 <- sin(t)
#     eps2 <- sin(t - (1 - rho)/2*pi)
#     dJ = dlnC_f(C,j=1,K,A,M,U) + eps1
#     dC = dlnC_f(C,j=2,K,A,M,U) + eps2
#     list( c( dJ=dJ, dC=dC ) )
#   })
# }
# 
# parms_sin <- parms
# parms_sin$rho <- 1
# lnNt_sin <- ode(y=equ,times=tseq,func=dlnNdt_f2,parms=parms_sin) 
# plot(lnNt_sin[,"lnJ"],lnNt_sin[,"lnC"])


