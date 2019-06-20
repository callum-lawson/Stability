### Simulation of structured consumer-resource chains ###

# Trial runs --------------------------------------------------------------

require("rootSolve")
require("phaseR")

### Inputs (parms)
source("Source/consumer_resource_functions_log.R")

sparms = list(
  chainlength = 2,
    # single number or vector of length Ya
    # eggs (storage structure) always start at 0
  nchain = 1,
  store = FALSE,
    # if FALSE, births are allocated directly to feeders
  slevel = c("consumer","resource")[1],
  storetype = c("diffuse","feeding")[2],
    # no "selective" because too complex (ESS work)
  movetype = c("diffuse","feeding","selective")[3],
    # if FALSE, is mixed specialist
  generalist = TRUE,
    # does storage operate through births (TRUE) or diffusion (FALSE)?
  discrete = FALSE,
  mbase = 1, # 1 *GRAM*
  morder = 2,
  tT = 24 * 365 * 50,
  nt = 100,
  sS = 7*52, 
    # number of seasons over time series
  bdt = NULL,   
  nstart = c(0,0)
)

zparms <- list(
  zmu = 20, 
  zsig = 0,
  zl = 24 * 365 
)

bc <- c(
  v = 0.00001,   # max input rate = vk *grams* per m^2 per hour
  k = 0.1,    # grams per m^2
  psi = 0,   # interference:handling time ratio
  omega = 1, # relative feeding rate of y2
  phi_E = 0, # relative death rate of eggs
  phi_m = 1, # relative death rate of y2
  u_E = 1,   # odds ratio of y1:y2 at equilibrium
  m_E = 1,   # migration rate 
  u_m = 0,
  m_m = 0,
  tau_E = 0, # lags in migration functions
  tau_m = 0
)
  # phi and omega could instead by controlled by body masses
  #   (in this case, phi can be fraction of adult body mass)

# bhat <- readRDS("Output/rate_parameters_simulated_03Sep2018.rds")
# bhat <- bdselect(bhat,bpos=rep(1,2))

bhat <- readRDS("Output/rate_parameters_marginal_brms_06Sep2018.rds")
bhat <- bdselect(bhat,bpos=rep(1,2))

# bhat$a$bz <- 0
# bhat$h$bz <- 0
# bhat$alpha$bz <- 0
# bhat$mu$bz[2] <- 0

# bhat$a$b0 <- bhat$a$b0 - 2
bhat$h$b0 <- -100 # 1.25 * 3.2089
# bhat$alpha$b0 <- bhat$alpha$b0 + 10 # bhat$alpha$b0 - 1.5

iparms <- iparmf(bhat,sparms)
parms <- c(sparms,iparms,zparms,bc)

trial <- popint(parms)
# matplot(trial[,-1],type="l")

# Eigenvalues - continuous ------------------------------------------------

# equ <- steady(y=c(R=1,C1=1),fun=d_web,parms=parms,time=0)$y
equ <- with(parms, {
  steady(y=y0,
    parms=parms,
    fun=d_web,
    times=c(0,Inf),
    method="runsteady",
    hold=FALSE
  )$y
})

jac <- jacobian.full(y=equ,fun=d_web,parms=parms,time=0)
eig <- eigen(jac)
vec <- eig$vectors
val <- eig$values
inv <- solve(vec)
n0 <- parms$nstart - equ
n0mat <- matrix(rep(inv %*% n0, each=parms$nt),nr=parms$nt,nc=parms$chainlength)
mmat <- n0mat * exp(outer(parms$tseq,val,"*"))
nmat <- t(vec %*% t(mmat))

# for(i in 1:parms$nt){
#   nmat[i,] <- vec %*% diag(exp(val*parms$tseq[i])) %*% inv %*% n0
# }

equdiff <- trial[,-1] - rep(equ,each=parms$nt)
ihat <- vec %*% c(1,0)
jhat <- vec %*% c(0,1) 

ihatA <- ihat + equ
jhatA <- jhat + equ
ihatB <- -ihat + equ
jhatB <- -jhat + equ

parms2 <- parms3 <- parms
parms2$y0[1] <- parms2$y0[1] * 2
parms3$y0[1] <- parms3$y0[1] * 0.5

trial2 <- popint(parms2)
trial3 <- popint(parms3)

# equdiff2 <- trial2[,-1] - rep(equ,each=parms$nt)
# equdiff3 <- trial3[,-1] - rep(equ,each=parms$nt)

# matplot(cbind(trial[,-1],trial2[,-1],trial3[,-1],nmat+rep(equ,each=parms$nt)),
#   type="l",
#   col=rep(c("black","red","blue","grey"),each=2),
#   lty=rep(1:2,times=4)
#   )

# Phase space plots -------------------------------------------------------

d_chain_express <- function(t,y,parameters){
  with(parameters,{
    dy <- y # convenient way to assign same names as y
    y <- exp(y) # hack to put on log scale
    fy <- f(R=y[-Y],C=y[-1],a,h,alpha,psi)
    xy <- x(y[-1],mu) + c(fy[-1]/alpha[-1],0)
    dy[1] <-  (g(y[1],v,k) - fy[1]/alpha[1]) / y[1] # more hacks
    dy[-1] <- (fy - xy)  / y[-1] # more hacks
    list(dy=dy,fy=fy,xy=xy)
  })
}

bddd <- with(parms, btf(t=0,bd,M,parms))
bddda <- c(bc,bddd,Y=parms$Yc)

cxlim <- c(equ[1]-1,equ[1]+1)
cylim <- c(equ[2]-1,equ[2]+1)

par(mfrow=c(1,1))
flowField(d_chain_express,xlim=cxlim,ylim=cylim,parameters=bddda,points=30,add=FALSE)
clines <- list()
clines[[1]] <- nullclines(d_chain_express, xlim=cxlim, ylim=cylim,
  parameters=bddda, points=100,col=rep("blue",2),add.legend=FALSE
)
# clines[[2]] <- nullclines(romac_phaser,x.lim = cxlim,y.lim = cylim,
#   parameters = pdm[2,], points = 100,
#   colour=rep("red",2)
# )

# tradj <- list()
# tradj[[1]] <- with(parms, 
#   trajectory(d_chain_express, y0=y0, tlim=c(0,tT), tstep=diff(tseq[1:2]), 
#     parameters=bddda)
# )

points(parms$y0[1],parms$y0[2])
abline(v=equ[1])
abline(h=equ[2])
arrows(x0=equ[1],y0=equ[2],x1=ihatA[1],y1=ihatA[2],col="red",lty=2) # i-hat
arrows(x0=equ[1],y0=equ[2],x1=jhatA[1],y1=jhatA[2],col="red",lty=3) # j-hat
arrows(x0=equ[1],y0=equ[2],x1=ihatB[1],y1=ihatB[2],col="red",lty=2,angle=-180) # i-hat 2
arrows(x0=equ[1],y0=equ[2],x1=jhatB[1],y1=jhatB[2],col="red",lty=3,angle=-180) # j-hat 2

lines(x=trial[,"R"],y=trial[,"C1"],col="purple")
# lines(mmat,col="red")
lines(nmat+rep(equ,each=parms$nt),col="red")
# x-squish much faster than y-squish
# overshoot of joint equilibrium occurs because R first goes to its own R*,
# then R* changes as C changes

lines(x=trial2[,"R"],y=trial2[,"C1"],lty=2,col="purple")
lines(x=trial3[,"R"],y=trial3[,"C1"],lty=3,col="purple")

abline(h=parms$y0[2],lty=3) # abline(h=equdiff[1,"C1"],lty=3)
points(equ[1],equ[2],col="red",pch="+")

parms4 <- parms5 <- parms6 <- parms
parms4$y0 <- equ + 0.1 
parms5$y0 <- equ + 0.05 
parms6$y0 <- equ + 0.15 

trial4 <- popint(parms4)
trial5 <- popint(parms5)
trial6 <- popint(parms6)

points(trial4[,-1])
points(trial5[,-1])
points(trial6[,-1])

# Oscillations ------------------------------------------------------------

### Rotation-scale

Bre <- Re(val)[1] # taking first eigenvalue (arbitrary)
Bim <- Im(val)[1]
Bmat <- matrix(c(Bre,-Bim,Bim,Bre),nr=2,nc=2)
atan(Bmat[2,1]/Bmat[1,1]) + 2*pi
# but this is the angle of the reaction vector, not the rotation?

Cmat <- cbind(Re(vec[,1]),Im(vec[,1]))

ihat <- Cmat %*% c(1,0)
jhat <- Cmat %*% c(0,1) 

ihatA <- ihat + equ
jhatA <- jhat + equ
ihatB <- -ihat + equ
jhatB <- -jhat + equ

arrows(x0=equ[1],y0=equ[2],x1=ihatA[1],y1=ihatA[2],col="red",lty=2) # i-hat
arrows(x0=equ[1],y0=equ[2],x1=jhatA[1],y1=jhatA[2],col="red",lty=3) # j-hat
arrows(x0=equ[1],y0=equ[2],x1=ihatB[1],y1=ihatB[2],col="red",lty=2,angle=-180) # i-hat 2
arrows(x0=equ[1],y0=equ[2],x1=jhatB[1],y1=jhatB[2],col="red",lty=3,angle=-180) # j-hat 2

### Cycle speeds
# How does timescale similarity of R and C affect oscillation speed?

# R speed

v_len <- 10
v_seq <- 10^seq(-7,-5,length.out=v_len)
v_parms <- v_jac <- vector("list",length=v_len)
v_equ <- matrix(nr=v_len,nc=2)
v_val <- vector("numeric",length=v_len)

for(i in 1:v_len){
  v_parms[[i]] <- parms
  v_parms[[i]]$v <- v_seq[i]
  v_equ[i,] <- steady(y=parms$y0,
           parms=v_parms[[i]],
           fun=d_web,
           times=c(0,Inf),
           method="runsteady",
           hold=FALSE
           )$y
  v_jac[[i]] <- jacobian.full(y=v_equ[i,],fun=d_web,parms=v_parms[[i]],time=0)
  v_val[i] <- eigen(v_jac[[i]])$values[1] # taking first eigenvalue
}

par(mfrow=c(2,2))
plot(v_seq,Re(v_val),type="b",ylab="shrink speed")
plot(v_seq,Im(v_val),type="b",ylab="spin speed")
plot(log10(v_seq),Re(v_val),type="b",ylab="shrink speed")
plot(log10(v_seq),Im(v_val),type="b",ylab="spin speed")
  # closer timescales -> shrink speed decreases, spin speed increases

# C speed

M_len <- 10
M_seq <- 10^seq(-3,0,length.out=M_len)
M_parms <- M_jac <- vector("list",length=M_len)
M_equ <- matrix(nr=M_len,nc=2)
M_val <- vector("numeric",length=M_len)

for(i in 1:M_len){
  M_parms[[i]] <- parms
  M_parms[[i]]$M <- M_seq[i]
  M_equ[i,] <- steady(y=parms$y0,
                      parms=M_parms[[i]],
                      fun=d_web,
                      times=c(0,Inf),
                      method="runsteady",
                      hold=FALSE
  )$y
  M_jac[[i]] <- jacobian.full(y=M_equ[i,],fun=d_web,parms=M_parms[[i]],time=0)
  M_val[i] <- eigen(M_jac[[i]])$values[1] # taking first eigenvalue
}

par(mfrow=c(2,2))
plot(M_seq,Re(M_val),type="b",ylab="shrink speed")
plot(M_seq,Im(M_val),type="b",ylab="spin speed")
plot(log10(M_seq),Re(M_val),type="b",ylab="shrink speed")
plot(log10(M_seq),Im(M_val),type="b",ylab="spin speed")

# Perturbations -----------------------------------------------------------

### Perturation curves

parms_eps <- parms

n_eps <- 10
eps_seq <- parms$zmu  + seq(-10,10,length.out=n_eps)
eps_eq_seq <- matrix(nr=n_eps,nc=parms$Yc)

for(i in 1:n_eps){
  parms_eps$zmu <- eps_seq[i]
  eps_eq_seq[i,] <- steady(
    y=parms_eps$y0,
    parms=parms_eps,
    fun=d_web,
    times=c(0,Inf),
    method="runsteady",
    hold=FALSE
  )$y
}

points(eps_eq_seq)

### Perturbation isoclines

parmsB <- parmsC <- parms
parmsB$zmu <- parms$zmu + 10
bdddB <- with(parmsB, btf(t=0,bd,M,parmsB))
bdddBa <- c(bc,bdddB,Y=parms$Yc)
equB <- with(parms, {
  steady(y=y0,
         parms=parmsB,
         fun=d_web,
         times=c(0,Inf),
         method="runsteady",
         hold=FALSE
  )$y
})

parmsC$zmu <- parms$zmu - 10
bdddC <- with(parmsC, btf(t=0,bd,M,parmsC))
bdddCa <- c(bc,bdddC,Y=parms$Yc)
equC <- with(parms, {
  steady(y=y0,
         parms=parmsC,
         fun=d_web,
         times=c(0,Inf),
         method="runsteady",
         hold=FALSE
  )$y
})

clines[[2]] <- nullclines(d_chain_express, xlim=cxlim, ylim=cylim,
                          parameters=bdddBa, points=100,col=rep("blue",2),
                          add.legend=FALSE,lty=rep(2,2)
                          )
clines[[3]] <- nullclines(d_chain_express, xlim=cxlim, ylim=cylim,
                          parameters=bdddCa, points=100,col=rep("blue",2),
                          add.legend=FALSE,lty=rep(3,2)
)
points(equB[1],equB[2],col="red",pch="+")
points(equC[1],equC[2],col="red",pch="+")

# plotCircle <- function(x, y, r) {
#   angles <- seq(0,2*pi,length.out=360)
#   lines(r*cos(angles)+x,r*sin(angles)+y,col="darkgrey")
# }
# plotCircle(equ[1],equ[2],r=0.25)
# plotCircle(equ[1],equ[2],r=1)

### Comparing flowfields

# par(mfrow=c(2,2))
# flowField(d_chain_express,xlim=cxlim,ylim=cylim,parameters=bddda,points=30,add=FALSE)
# flowField(d_chain_express,xlim=cxlim,ylim=cylim,parameters=bdddBa,points=30,add=FALSE)
# flowField(d_chain_express,xlim=cxlim,ylim=cylim,parameters=bdddCa,points=30,add=FALSE)

# Timescale separation ----------------------------------------------------

outer(val,val,"/")

with(bddd, parms$k * a / mu) # should be not to small compared to 1
with(bddd, mu / parms$v)     # should be as close to 0 as possible

# Other eigenstuff --------------------------------------------------------

er <- eigen(A)$vectors # right
el <- eigen(t(A))$vectors # left

sum(Conj(t.default(el[,1])) *  v[1] *  t(er)[,1])

e3 <- solve(e1)
e3 <- ginv(e1)

# http://www.wolframalpha.com/input/?i=Eigenvectors%5B%7B%7B-(v%2Ba+C),-a+R%7D,+%7B+a+C,+a+R+-m%7D%7D%5D

# http://www.wolframalpha.com/input/?i=Eigenvectors%5B%7B%7B-a+v+k+%2Fm,v+(k+%2F+m+-+1)%7D,+%7B-m,+0%7D%7D%5D

  # still fluctuates when transitioning from one equilibrium to another
  # [separate comment]

# Growth curves - continuous ----------------------------------------------

Cmin <- -1 
Cmax <- log10(cylim[2])
nC <- 50
Cseq <- 10^seq(Cmin,Cmax,length.out=nC)

sparms1 <- sparms
zparms1 <- zparms
sparms1$discrete = FALSE
sparms1$store = FALSE
zparms1$zsig = 0
iparms1 <- iparmf(bhat,sparms1)
parms1 <- c(sparms1,iparms1,zparms1,bc)

dC1 <- rCfv(Cseq,parms1) / Cseq # *per-capita* growth

parms2 <- parms1
parms2$zmu <- parms1$zmu + 5
dC2 <- rCfv(Cseq,parms2) / Cseq

parms3 <- parms3b <- parms1
parms3b$zsig <- parms1$zsig + 5
parms3$bdt <- with(parms3b, rate_int_l(bd=bd,bn=names(bd),parms=parms3b))
dC3 <- rCfv(Cseq,parms3) / Cseq 
  # parms1 because don't want fluctuating for C* calculation

parms4 <- parms1
parms4$zmu <- parms1$zmu - 5
dC4 <- rCfv(Cseq,parms4) / Cseq

matplot(log10(Cseq), cbind(dC1,dC2,dC3,dC4), 
        type="l", col=c("orange","red","green","blue"))
abline(h=0,col="black",lty=2)
lines(log10(Cseq),apply(cbind(dC2,dC4),1,mean),col="purple",lty=1)

Cstar <- sapply(list(parms1,parms2,parms3,parms4),Cstarf)
# points(log(Cstar),rep(0,length(Cstar)))
points(mean(log10(Cstar[c(2,4)])),0)
  # next up: function for Cstar calculation over vector of different z

spectrum <- function(Rstar,alpha,a,h,mu){
  lrmax <- log( f(R=Rstar,C=1,a,h,alpha) ) # 1 -> per-capita (doesn't include w effects)
  lrmin <- log ( mu )
  lrmax - lrmin
}

lapply(list(parms4,parms1,parms2), function(x) btf(t=0, x$bd, x$M, parms=x))

# R Trajectories ----------------------------------------------------------

parmsR1 <- parmsR2 <- parmsR3 <- parms
parmsR1$y0 <- c(R=0.015,C=1)
parmsR2$y0 <- c(R=0.02,C=1)
parmsR3$y0 <- c(R=0.025,C=1)

R1t <- popint(parmsR1)
R2t <- popint(parmsR2)
R3t <- popint(parmsR3)

R1dt <- diff(R1t)[,"R"]
R2dt <- diff(R2t)[,"R"]
R3dt <- diff(R3t)[,"R"]

plot(R1dt ~ R1t[-1,"R"],type="l")
lines(R2dt ~ R2t[-1,"R"],col="red")
lines(R3dt ~ R3t[-1,"R"],col="blue")


