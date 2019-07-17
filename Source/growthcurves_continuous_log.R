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

### Rotation-scale plots

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

### Cycle eigenvectors and values
  # How does timescale similarity of R and C affect oscillation speed?

oparms <- parms
# oparms$bd$alpha$b0 <- 100
# oparms$bd$alpha$h <- -2.373758 + 2

p_len <- 5
p_tot <- p_len^3
v_seq <- 10^seq(-7,-3,length.out=p_len)
m_seq <- seq(-11,-9,length.out=p_len) # b0 (log scale)
a_seq <- seq(-7,-5,length.out=p_len)  # b0 (log scale)
p_dat <- expand.grid(v=v_seq,m=m_seq,a=a_seq)

p_parms <- p_jac <- vector("list",length=p_tot)
p_equ <- matrix(nr=p_tot,nc=2)
p_vec <- matrix(nr=p_tot,nc=4) # 2 * 2D vectors = 4
p_val <- matrix(nr=p_tot,nc=2)
p_bd <- data.frame(a=rep(NA,p_len),h=rep(NA,p_len),alpha=rep(NA,p_len),mu=rep(NA,p_len))

for(i in 1:p_tot){
  p_parms[[i]] <- oparms
  p_parms[[i]]$v <- p_dat$v[i]
  p_parms[[i]]$bd$mu$b0 <- p_dat$m[i]
  p_parms[[i]]$bd$a$b0 <- p_dat$a[i]
  p_equ[i,] <- steady(y=oparms$y0,
                      parms=p_parms[[i]],
                      fun=d_web,
                      times=c(0,Inf),
                      method="runsteady",
                      hold=FALSE
  )$y
  p_jac[[i]] <- jacobian.full(y=p_equ[i,],fun=d_web,parms=p_parms[[i]],time=0)
  p_eig <- eigen(p_jac[[i]]) # replaced each loop
  p_vec[i,] <- p_eig$vectors # storing both eigenvectors
  p_val[i,] <- p_eig$values # storing both eigenvalues
  p_bd[i,] <- btf(t=0,p_parms[[i]]$bd,oparms$M,p_parms[[i]])
}

p_TS <- log10( p_bd$mu / p_dat$v )
p_R0 <- log10( parms$k * p_bd$a / p_bd$mu)
p_TR <- p_TS - p_R0

p_val_Re <- Re(p_val)
p_val_Im <- Im(p_val)
is_feas <- apply(p_equ > -20,1,function(x) sum(x)==2)
iso <- p_val_Re[,1]==p_val_Re[,2] & is_feas
iss <- p_val_Re[,1]!=p_val_Re[,2] & is_feas

mypchs <- rep(rep(15:19,each=p_len),times=p_len)
mycols <- rep(c("black","purple","blue","orange","red"),each=p_len^2)

### Shrink and spin speeds

p_val_o <- p_val[iso,1] # taking first eigenvalue
p_TS_o <- p_TS[iso]
p_R0_o <- p_R0[iso]
p_v_o <- p_bd$mu[iso] # p_dat$v[iso]
  # adjusts timescale of response (scale-down faster systems more)
  # I.e. making Thieme's timescale adjustment by C lifespan
p_TR_o <- p_TR[iso]

cr <- colorRamp(c("blue", "red"))
scale01 <- function(x) scale(x,center=min(x),scale=max(x)-min(x))
p_R0_o_scaled <- scale01(p_R0_o)
mycols_R0 = rgb(cr(p_R0_o_scaled), max=255)

par(mfrow=c(1,2),mar=c(4.5,4.5,2,2))
plot(p_TR_o,log(-Re(p_val_o)/p_v_o),pch=mypchs[iso],col=mycols_R0,
  xlab=expression(R[0]/p), ylab="shrink speed")
plot(p_TR_o,log(Im(p_val_o)/p_v_o),pch=mypchs[iso],col=mycols_R0,
  xlab=expression(R[0]/p), ylab="spin speed")
# plot(p_R0_o,Re(p_val_o)/p_v_o,pch=mypchs[iso],col=mycols_R0,
#   xlab="Interaction strength", ylab="shrink speed")
# plot(p_R0_o,Im(p_val_o)/p_v_o,pch=mypchs[iso],col=mycols_R0,
#   xlab="Interaction strength", ylab="spin speed")

### Analytical vs eigenvalue TS

val1 <- p_val_Re[iss,1]
val2 <- p_val_Re[iss,2]
p_TS_s <- log10( val1 / val2 )
p_TR_s <- p_TR[iss]
p_R0_s <- p_R0[iss]
p_R0_s_scaled <- scale01(p_R0_s)
mycols_R0_s = rgb(cr(p_R0_s_scaled), max=255)

par(mfrow=c(1,2))
plot(p_TS[iss],-p_TS_s,pch=mypchs[iss],col=mycols[iss],
  xlab="TS (analytic)", ylab="TS (eigenvector)")
abline(0,1,lty=3,col="gray")

plot(p_TR_s,-p_TS_s,pch=mypchs[iss],col=mycols[iss],
  xlab=expression(R[0]/p), ylab="TS (eigenvector)")
abline(0,1,lty=3,col="gray")

par(mfcol=c(2,3))
plot(p_TS_s,val1,pch=mypchs[iss],col=mycols[iss],
     xlab=expression(p), ylab=expression(lambda[1]))
plot(p_TS_s,val2,pch=mypchs[iss],col=mycols[iss],
     xlab=expression(p), ylab=expression(lambda[2]))
plot(p_R0_s,val1,pch=mypchs[iss],col=mycols[iss],
     xlab=expression(R[0]), ylab=expression(lambda[1]))
plot(p_R0_s,val2,pch=mypchs[iss],col=mycols[iss],
     xlab=expression(R[0]), ylab=expression(lambda[2]))
plot(p_TR_s,val1,pch=mypchs[iss],col=mycols[iss],
     xlab=expression(R[0]/p), ylab=expression(lambda[1]))
plot(p_TR_s,val2,pch=mypchs[iss],col=mycols[iss],
     xlab=expression(R[0]/p), ylab=expression(lambda[2]))

### Eigenvectors

p_angle_fast <- ( atan( Re(p_vec[iss,2]) / Re(p_vec[iss,1]) ) + pi ) / pi
  # y always positive, x always negative (so second quadrant)
  
p_angle_slow <- atan( Re(p_vec[iss,4]) / Re(p_vec[iss,3]) ) / pi
  # y always negative, x always positive (so fourth quadrant)

  # using 2*pi to convert to fractions of a circle
  # https://en.wikipedia.org/wiki/Polar_coordinate_system#Converting_between_polar_and_Cartesian_coordinates

par(mfcol=c(2,4))
plot(p_TS[iss],p_angle_fast,pch=mypchs[iss],col=mycols_R0_s,xlab="TS (analytic)")
plot(p_TS[iss],p_angle_slow,pch=mypchs[iss],col=mycols_R0_s,xlab="TS (analytic)")
plot(p_TS_s,p_angle_fast,pch=mypchs[iss],col=mycols_R0_s,xlab="TS (eigenvalue)")
plot(p_TS_s,p_angle_slow,pch=mypchs[iss],col=mycols_R0_s,xlab="TS (eigenvalue)")
plot(p_R0_s,p_angle_fast,pch=mypchs[iss],col=mycols_R0_s,xlab=expression(R[0]))
plot(p_R0_s,p_angle_slow,pch=mypchs[iss],col=mycols_R0_s,xlab=expression(R[0]))
plot(p_TR[iss],p_angle_fast,pch=mypchs[iss],col=mycols_R0_s,xlab=expression(R[0]/p))
plot(p_TR[iss],p_angle_slow,pch=mypchs[iss],col=mycols_R0_s,xlab=expression(R[0]/p))

par(mfrow=c(1,1))
plot(p_angle_fast,p_angle_slow,pch=mypchs[iss],col=mycols[iss])

### Equilibrium check

par(mfrow=c(1,2))
p_TS_feas <- p_TS[is_feas]
plot(p_TS_feas,p_equ[is_feas,1],pch=mypchs[is_feas],col=mycols[is_feas],xlab=expression(log[10](v/mu)))
plot(p_TS_feas,p_equ[is_feas,2],pch=mypchs[is_feas],col=mycols[is_feas],xlab=expression(log[10](v/mu)))

### Anaytical oscillations
# See "Scaling" OneNote file

TS_abiotic_f <- function(R_0,p){
  delta_lambda <- sqrt(R_0^2 - 4*p*(R_0 - 1))
  (-R_0 - delta_lambda) / (-R_0 + delta_lambda) # 1/2's cancel
}

TS_biotic_f <- function(R_0,p){
  delta_lambda <- sqrt(1 - 4*p*R_0*(R_0 - 1))
  (-1 - delta_lambda) / (-1 + delta_lambda) # 1/(2*R_0)'s cancel
}

R_0_seq <- 10^seq(0.01,0.7,length.out=100)
p_seq <- seq(0.01,0.99,length.out=100)
TS_abiotic_mat <- sapply(p_seq, function(p) TS_abiotic_f(R_0=R_0_seq,p=p))
TS_biotic_mat  <- sapply(p_seq, function(p) TS_biotic_f(R_0=R_0_seq,p=p))

par(mfrow=c(2,3),mar=c(4.5,4.5,1.5,1.5))
matplot(R_0_seq,log(TS_abiotic_mat),type="l",col=heat.colors(100),lty=1)
matplot(p_seq,t(log(TS_abiotic_mat)),type="l",col=heat.colors(100),lty=1)
image(R_0_seq,p_seq,log(TS_abiotic_mat))
matplot(R_0_seq,log(TS_biotic_mat),type="l",col=heat.colors(100),lty=1)
matplot(p_seq,t(log(TS_biotic_mat)),type="l",col=heat.colors(100),lty=1)
image(R_0_seq,p_seq,log(TS_biotic_mat))

# Conditions for oscillations:
par(mfrow=c(1,2))
curve(x/(4*(1-1/x)),xlim=c(1,10),
  xlab=expression(R[0]),
  ylab=expression(p[crit]),
  main="abiotic"
  )
curve(1/(4*x*(x-1)),xlim=c(1,10),
  xlab=expression(R[0]),
  ylab=expression(p[crit]),
  main="biotic"
)

# Eigenvector angle difference

dtheta_abiotic_f <-  function(R_0,p){
  delta_lambda <- sqrt(R_0^2 - 4*p*(R_0 - 1))
  vec1 <- 1/(2*(R_0-1)) * (- R_0 - delta_lambda)
  vec2 <- 1/(2*(R_0-1)) * (- R_0 + delta_lambda)
  theta1 <- atan2(vec1,1)
  theta2 <- atan2(vec2,1)
  propdiff <- (theta2 - theta1)/(2*pi)
  return(propdiff)
}

dtheta_biotic_f <-  function(R_0,p){
  delta_lambda <- sqrt(1 - 4*p*R_0*(R_0 - 1))
  vec1 <- 1/(2*R_0*(R_0-1)) * (- 1 - delta_lambda)
  vec2 <- 1/(2*(R_0-1)) * (- 1 + delta_lambda)
  theta1 <- atan2(vec1,1)
  theta2 <- atan2(vec2,1)
  propdiff <- (theta2 - theta1)/(2*pi)
  return(propdiff)
}

dtheta_abiotic_mat <- sapply(p_seq,function(p) dtheta_abiotic_f(R_0_seq,p))
dtheta_biotic_mat <- sapply(p_seq,function(p) dtheta_biotic_f(R_0_seq,p))

par(mfrow=c(2,3),mar=c(4.5,4.5,1.5,1.5))
matplot(R_0_seq,dtheta_abiotic_mat,type="l",col=heat.colors(100),lty=1)
matplot(p_seq,t(dtheta_abiotic_mat),type="l",col=heat.colors(100),lty=1)
image(R_0_seq,p_seq,dtheta_abiotic_mat)
matplot(R_0_seq,dtheta_biotic_mat,type="l",col=heat.colors(100),lty=1)
matplot(p_seq,t(dtheta_biotic_mat),type="l",col=heat.colors(100),lty=1)
image(R_0_seq,p_seq,dtheta_biotic_mat)

par(mfrow=c(1,2))
plot(log(TS_abiotic_mat),dtheta_abiotic_mat,main="abiotic")
plot(log(TS_biotic_mat),dtheta_biotic_mat,main="biotic")

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


