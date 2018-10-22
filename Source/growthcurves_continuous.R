### Simulation of structured consumer-resource chains ###

# Trial runs --------------------------------------------------------------

require("rootSolve")

### Inputs (parms)
source("Source/consumer_resource_functions.R")

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
  mbase = 0.001, # 1mg
  morder = 2,
  tT = 24 * 365,
  nt = 100,
  sS = 7*52, 
    # number of seasons over time series
  bdt = NULL,   
  nstart = c(1,1) # c(1,0.01)
)

zparms <- list(
  zmu = 20, 
  zsig = 0,
  zl = 24 * 365 
)

bc <- c(
  v = 0.01,   # max input rate = vk *grams* per m^2 per hour
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
bhat$h$b0 <- bhat$h$b0 - 100 # 1.25 * 3.2089
# bhat$alpha$b0 <- bhat$alpha$b0 - 1.5
# bhat$alpha$b0 <- bhat$alpha$b0 + 10

iparms <- iparmf(bhat,sparms)
parms <- c(sparms,iparms,zparms,bc)

Cmin <- -1 
Cmax <- -0.5
nC <- 50
Cseq <- 10^seq(Cmin,Cmax,length.out=nC)

trial <- popint(parms)
matplot(trial[,-1],type="l")

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

plot(1,1,type="n",xlim=c(-1,1),ylim=c(-1,2),xlab="R",ylab="C")
points(n0[1],n0[2])
abline(h=0)
abline(v=0)
arrows(x0=0,y0=0,x1=ihat[1],y1=ihat[2],col="red",lty=2) # i-hat
arrows(x0=0,y0=0,x1=jhat[1],y1=jhat[2],col="red",lty=3) # j-hat
arrows(x0=0,y0=0,x1=-ihat[1],y1=-ihat[2],col="red",lty=2,angle=-180) # i-hat 2
arrows(x0=0,y0=0,x1=-jhat[1],y1=-jhat[2],col="red",lty=3,angle=-180) # j-hat 2

lines(x=equdiff[,"R"],y=equdiff[,"C1"],col="black")
# lines(mmat,col="red")
lines(nmat,col="red")
  # x-squish much faster than y-squish
  # overshoot of joint equilibrium occurs because R first goes to its own R*,
  # then R* changes as C changes

# Resource isocline
Ctseq <- seq(0,1,length.out=10)
RCtmat <- t(Rstarfv(Ctseq,parms)) - rep(equ,each=10)
lines(x=RCtmat[,1],y=RCtmat[,2],col="blue")
# plot(x=RCtmat[,1],y=RCtmat[,2],col="blue")

bddd <- with(parms, btf(t=0,bd,M,parms))
outer(val,val,"/")

with(bddd, parms$k * a / mu) # should be not to small compared to 1
with(bddd, mu / parms$v)     # should be as close to 0 as possible

abline(h=equdiff[1,"C1"],lty=3)

parms2 <- parms3 <- parms
parms2$y0[1] <- parms2$y0[1] * 2
parms3$y0[1] <- parms3$y0[1] * 0.5

trial2 <- popint(parms2)
trial3 <- popint(parms3)

equdiff2 <- trial2[,-1] - rep(equ,each=parms$nt)
equdiff3 <- trial3[,-1] - rep(equ,each=parms$nt)

lines(x=equdiff2[,"R"],y=equdiff2[,"C1"],lty=2)
lines(x=equdiff3[,"R"],y=equdiff3[,"C1"],lty=3)

matplot(cbind(trial[,-1],trial2[,-1],trial3[,-1],nmat+rep(equ,each=parms$nt)),
  type="l",
  col=rep(c("black","red","blue","grey"),each=2),
  lty=rep(1:2,times=4)
  )

# Phase space plots -------------------------------------------------------

d_chain_express <- function(t,y,parameters){
  with(parameters,{
    dy <- y # convenient way to assign same names as y
    fy <- f(R=y[-Y],C=y[-1],a,h,alpha,psi)
    xy <- x(y[-1],mu) + c(fy[-1]/alpha[-1],0)
    dy[1] <-  g(y[1],v,k) - fy[1]/alpha[1]
    dy[-1] <- fy - xy
    list(dy=dy,fy=fy,xy=xy)
  })
}

bddda <- c(bc,bddd,Y=parms$Yc)
require(phaseR)

cxlim <- c(0,1)
cylim <- c(0,2)

par(mfrow=c(1,1))
flowField(d_chain_express,xlim=cxlim,ylim=cylim,parameters=bddda,points=30,add=FALSE)
clines <- list()
clines[[1]] <- nullclines(d_chain_express, xlim=cxlim, ylim=cylim,
  parameters=bddda, points=100,colour=rep("blue",2)
)
# clines[[2]] <- nullclines(romac_phaser,x.lim = cxlim,y.lim = cylim,
#   parameters = pdm[2,], points = 100,
#   colour=rep("red",2)
# )

tradj <- list()
tradj[[1]] <- with(parms, 
  trajectory(d_chain_express, y0=y0, tlim=c(0,tT), tstep=diff(tseq[1:2]), 
    parameters=bddda)
)
tradj[[2]] <- with(parms, 
  trajectory(d_chain_express, y0=c(0.01,0.01), tlim=c(0,tT), tstep=diff(tseq[1:2]), 
    parameters=bddda)
)
tradj[[3]] <- with(parms, 
  trajectory(d_chain_express, y0=c(0.01,0.5), tlim=c(0,tT), tstep=diff(tseq[1:2]), 
    parameters=bddda)
)

# next step: trajectories for different climates

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
