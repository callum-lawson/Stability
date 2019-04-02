### Explore MacArthur R* models in C-space ### 

require("rootSolve")
require("deSolve")
require("phaseR")

# Functions ---------------------------------------------------------------
# i = resource index; j = consumer index
# a = v * a
# mu = V * mu

# Rstar_f <- function(C,i,K,A){
#   K[i] / (1 + A[i,1]*C[1] + A[i,2]*C[2])
# }
#   # abiotic
 
Rstar_f <- function(C,i,K,A){
  K[i] * (1 - A[i,1]*C[1] - A[i,2]*C[2])
}
  # logistic

Cin_f <- function(C,j,K,A){
  A[1,j] * Rstar_f(C,i=1,K,A) + A[2,j] * Rstar_f(C,i=2,K,A) 
}

dC_f <- function(C,j,K,A,M,U){
  Cin_f(C,j,K,A) - M[j] + U[j]/C[j]
}

lndCdt_f <- function(t,y,parms){
  with(parms, {
    C <- exp(y)
    list( c( dC1 = dC_f(C,j=1,K,A,M,U), dC2 = dC_f(C,j=2,K,A,M,U) ) )
  })
}

# Parameters --------------------------------------------------------------

k1 <- 1
k2 <- 1
a11 <- 2
a12 <- 0.5
a21 <- 0.5
a22 <- 2
mu1 <- 1
mu2 <- 0.1
u1 <- 0
u2 <- 0

parms <- list(
  K = c(k1,k2),
  A = matrix(c(a11,a12,a21,a22),nr=2,nc=2),
  M = c(mu1,mu2),
  U = c(u1,u2)
)

# Integrate ---------------------------------------------------------------

equ <- steady(y = c(lnC1=0,lnC2=0), times = c(0,Inf), func = lndCdt_f, parms = parms, method = "runsteady")$y

lnC0 <- equ + log(1.5) # 1.5 times above N*
tT <- 1000
tseq <- seq(0,100,length.out=tT)

lnCt <- ode(y=lnC0,times=tseq,func=lndCdt_f,parms=parms) 
# matplot(lnCt[,1],lnCt[,-1],type="l")

# Eigenvalues -------------------------------------------------------------

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

cxlim <- c(-2.5,-1.5) # c(equ[1]-1,equ[1]+1)
cylim <- c(-0.75,-0.5) # c(equ[2]-1,equ[2]+1)

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

points(x=lnCt[,"lnC1"],y=lnCt[,"lnC2"],col="purple",type="b",lwd=1)
points(lnnmat+rep(equ,each=tT),col="red",type="b")

# Perturbation ------------------------------------------------------------

### Vertical
parms2 <- parms
parms2$M <- parms$M - 0.25
clines[[2]] <- nullclines(lndCdt_f, xlim=cxlim, ylim=cylim,
  parameters=parms2, points=100,col=rep("orange",2),add.legend=FALSE
)

### Horizontal
parms3 <- parms
parms3$M <- parms$M * 0.75
clines[[3]] <- nullclines(lndCdt_f, xlim=cxlim, ylim=cylim,
  parameters=parms3, points=100,col=rep("green",2),add.legend=FALSE
)

### Immigration
parms4 <- parms
parms4$U <- parms$U + 0.05
clines[[4]] <- nullclines(lndCdt_f, xlim=cxlim, ylim=cylim,
  parameters=parms4, points=100,col=rep("brown",2),add.legend=FALSE
)

### Back

equ2 <- steady(y = c(lnC1=0,lnC2=0), times = c(0,Inf), func = lndCdt_f, parms = parms2, method = "runsteady")$y
lnCt2 <- ode(y=equ2,times=tseq,func=lndCdt_f,parms=parms) 
points(x=lnCt2[,"lnC1"],y=lnCt2[,"lnC2"],col="orange",type="b")

equ3 <- steady(y = c(lnC1=0,lnC2=0), times = c(0,Inf), func = lndCdt_f, parms = parms3, method = "runsteady")$y
lnCt3 <- ode(y=equ3,times=tseq,func=lndCdt_f,parms=parms) 
points(x=lnCt3[,"lnC1"],y=lnCt3[,"lnC2"],col="green",type="b")

equ4 <- steady(y = c(lnC1=0,lnC2=0), times = c(0,Inf), func = lndCdt_f, parms = parms4, method = "runsteady")$y
lnCt4 <- ode(y=equ4,times=tseq,func=lndCdt_f,parms=parms) 
points(x=lnCt4[,"lnC1"],y=lnCt4[,"lnC2"],col="brown",type="b")

### There

lnCt2b <- ode(y=equ,times=tseq,func=lndCdt_f,parms=parms2) 
points(x=lnCt2b[,"lnC1"],y=lnCt2b[,"lnC2"],col="orange",type="b")

lnCt3b <- ode(y=equ,times=tseq,func=lndCdt_f,parms=parms3) 
points(x=lnCt3b[,"lnC1"],y=lnCt3b[,"lnC2"],col="green",type="b")

lnCt4b <- ode(y=equ,times=tseq,func=lndCdt_f,parms=parms4) 
points(x=lnCt4b[,"lnC1"],y=lnCt4b[,"lnC2"],col="brown",type="b")

### Equilibria

equ_all <- rbind(equ,equ2,equ3,equ4)
points(x=equ_all[,"lnC1"],y=equ_all[,"lnC2"],col="black",pch="+",cex=1.1)

# Eigenspace --------------------------------------------------------------

equvec <- rep(equ,each=tT)
tran <- t(inv %*% t(lnCt[,-1] - equvec))
tran2 <- t(inv %*% t(lnCt2[,-1] - equvec))
tran3 <- t(inv %*% t(lnCt3[,-1] - equvec))
tran4 <- t(inv %*% t(lnCt4[,-1] - equvec))

dtran <- rbind(tran2,tran3,tran4)
plot(dtran,col=rep(c("orange","green","brown"),each=tT),pch=16)

# Other eigenvalues -------------------------------------------------------

### Other points

jac2 <- jacobian.full(y=equ2,fun=lndCdt_f,parms=parms,time=0)
eig2 <- eigen(jac2)


# Immigration growth curve ------------------------------------------------

im1 <- function(u,N,alpha=1,beta=1){
  u/N + alpha - beta*N
}

im2 <- function(u,N,alpha=1,beta=1){
  u/N + alpha - beta*log(N)
}

implot <- function(FUN=im1,xtype="log"){
  x <- exp(seq(-2,2,length.out=100))
  y <- sapply(seq(-1,1,by=1),FUN,N=x)
  if(xtype=="log") matplot(log(x),y,type="l",col="black",lty=c(2,1,2))
  if(xtype=="abs") matplot(x,y,type="l",col="black",lty=c(2,1,2))
  abline(h=0,lty=3,col="red")
}

par(mfrow=c(2,2))
implot(im1,xtype="log")
implot(im1,xtype="abs")
implot(im2,xtype="log")
implot(im2,xtype="abs")

# Simulations with sin waves ----------------------------------------------

lndCdt_f2 <- function(t,y,parms){
  with(parms, {
    C <- exp(y)
    eps1 <- sin(t)
    eps2 <- sin(t - (1 - rho)/2*pi)
    dC1 = dC_f(C,j=1,K,A,M,U) + eps1
    dC2 = dC_f(C,j=2,K,A,M,U) + eps2
    list( c( dC1=dC1, dC2=dC2 ) )
  })
}

parms_sin <- parms
parms_sin$rho <- 1
lnCt_sin <- ode(y=equ,times=tseq,func=lndCdt_f2,parms=parms_sin) 
plot(lnCt_sin[,"lnC1"],lnCt_sin[,"lnC2"])


