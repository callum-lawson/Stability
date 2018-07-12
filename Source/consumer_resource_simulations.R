### Simulation of structured consumer-resource chains ###

# Trial runs --------------------------------------------------------------

### Inputs (parms)
source("Source/Consumer_resource_functions.R")

sparms = list(
  chainlength = 2,
    # single number or vector of length Ya
    # eggs (storage structure) always start at 0
  nchain = 2,
  store = FALSE,
    # if FALSE, births are allocated directly to feeders
  slevel = c("consumer","resource")[1],
  storetype = c("diffuse","feeding")[2],
    # no "selective" because too complex (ESS work)
  movetype = c("diffuse","feeding","selective")[3],
    # if FALSE, is mixed specialist
  generalist = FALSE,
    # does storage operate through births (TRUE) or diffusion (FALSE)?
  discrete = FALSE,
  tT = 24*7*52,
  nt = 24*7*52,
  sS = 7*52, 
    # number of seasons over time series
  bdt = NULL,   
  nstart = 1
    # c(1,1,2, 1,1,1)
)

zparms <- list(
  zmu = 0, 
  zsig = 0,
  zl = 24
)

bc <- c(
  v = 0.1,   # max flow rate = k grams per m^2 per hour
  k = 10,    # 10g per m^2
  psi = 0,   # interference:handling time ratio
  omega = 1, # relative feeding rate of y2
  phi_E = 0, # relative death rate of eggs
  phi_m = 1, # relative death rate of y2
  u_E = 1,   # odds ratio of y1:y2 at equilibrium
  m_E = 1,   # migration rate 
  u_m = 1,
  m_m = 0.01,
  tau_E = 24 * 7, # lags in migration functions
  tau_m = 0
)
  # phi and omega could instead by controlled by body masses
  #   (in this case, phi can be fraction of adult body mass)

bhat <- readRDS("Output/rate_parameters_simulated_21Jun2018.rds")
bhat <- bdselect(bhat,bpos=c(1,1)) # same params for top consumer

iparms <- iparmf(bhat,sparms)
parms <- c(sparms,iparms,zparms,bc)

attach(parms)
y <- y0
t <- t0
trial <- popint(parms)
matplot(log(trial[,-1]),type="l")

# Fluctuation speed -------------------------------------------------------

sparms$discrete <- FALSE # TRUE

zlmin <- 0
zlmax <- 2
nzl <- 10
zlseq <- 24 * 10 ^ seq(zlmin,zlmax,length.out=nzl)
Cpos <- parms$nchain + 2

Narr <- with(parms, array(dim=c(nt,Cpos+1,nzl)))

for(i in 1:nzl){
  zparms <- list(
    zmu = 0,
    zsig = 5,
    zl = zlseq[i]
  )
  parms <- c(sparms,iparms,zparms,bc)
  Narr[,,i] <- popint(parms)
}

### Plot

gmf <- function(x) mean(log(x))
gsf <- function(x) sd(log(x))
require(RColorBrewer)
nshow <- nzl
mypalette <- rev(brewer.pal(nshow,"RdBu"))
nburn <- round(parms$nt/2,0)
rem <- 1:nburn
Cam <- apply(Narr[-rem,Cpos,],2,mean)
Cas <- apply(Narr[-rem,Cpos,],2,sd)
Cgm <- apply(Narr[-rem,Cpos,],2,gmf)
Cgs <- apply(Narr[-rem,Cpos,],2,gsf)
par(mfrow=c(1,1),mar=c(2,2,2,2))
matplot(parms$tseq[-rem],log(Narr[-rem,Cpos,1:nshow]),type="l",lty=1,col=mypalette)
par(mfrow=c(3,3))
for(i in 1:9){
  plot(parms$tseq[-rem],log(Narr[-rem,Cpos,i]),col=mypalette[i],type="l")
  lines(parms$tseq[-rem],log(Narr2[-rem,Cpos,i]),col=mypalette[i],type="l",lty=2)
}
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(Cam~log(zlseq/24),type="b")
plot(Cas~log(zlseq/24),type="b")
plot(Cgm~log(zlseq/24),type="b")
plot(Cgs~log(zlseq/24),type="b")

parmsB <- parms
parmsB$zl <- 24
parmsB$t0 <- -12
Nshift <- popint(parmsB)
lines(log(Nshift[-rem,Cpos])~parms$tseq[-rem],col="green")

zparms <- list(
  zmu = 0,
  zsig = 0,
  zl = 1 # irrelevant
)
parms <- c(sparms,iparms,zparms,bc)
Cstarcons <- Cstarf(parms)

par(mfrow=c(1,1),mar=c(2,2,2,2))
matplot(parms$tseq[-rem],log(Narr[-rem,Cpos,]),type="l",lty=1,col=mypalette)
abline(h=log(Cstarcons),lty=2)
  # NB: only works for continuous

### No-lag comparison

Narr2 <- with(parms, array(dim=dim(Narr)))
bc$tau_E <- 0
for(i in 1:nzl){
  zparms <- list(
    zmu = 0,
    zsig = 5,
    zl = zlseq[i]
  )
  parms <- c(sparms,iparms,zparms,bc)
  Narr2[,,i] <- popint(parms)
}
matplot(parms$tseq[-rem],log(Narr2[-rem,Cpos,1:nshow]),type="l",lty=2,col=mypalette,add=TRUE)

### Cons-env comparison

bc$tau_E <- 24 * 7
zparms <- list(
  zmu = 0,
  zsig = 0,
  zl = 0
)
parms <- c(sparms,iparms,zparms,bc)
Narr3 <- popint(parms)
par(mfrow=c(1,1))
plot(log(Narr3[-rem,3]),type="l")

### Check stability of other lags

sparms$chainlength <- 3
sparms$tT <- 24 * 7 * 52 * 100
sparms$nt <- 7 * 52 * 100
iparms <- iparmf(bhat,sparms)
bc$tau_E <- 24 * 30
parms <- c(sparms,iparms,zparms,bc)
Narr4 <- popint(parms)
par(mfrow=c(1,1))
plot(log(Narr4[((sparms$nt-1000):sparms$nt),4]),type="l")
  # very slow underdamping to equilibrium (not cycles)

# Growth curves - continuous ----------------------------------------------

sparms1 <- sparms
zparms1 <- zparms
sparms1$discrete = FALSE
sparms1$store = FALSE
zparms1$zsig = 0
iparms1 <- iparmf(bhat,sparms1)
parms1 <- c(sparms1,iparms1,zparms1,bc)

Cmin <- -3 # 1.4
Cmax <- 3 # 1.7
nC <- 100
Cseq <- 10^seq(Cmin,Cmax,length.out=nC)

dC1 <- rCfv(Cseq,parms1) / Cseq # *per-capita* growth

parms2 <- parms1
parms2$zmu <- 10
dC2 <- rCfv(Cseq,parms2) / Cseq

parms3 <- parms3b <- parms1
parms3b$zsig <- 5
parms3$bdt <- with(parms3b, rate_int_l(bd=bd,bn=names(bd),parms=parms3b))
dC3 <- rCfv(Cseq,parms3) / Cseq 
  # parms1 because don't want fluctuating for C* calculation

parms4 <- parms1
parms4$zmu <- -10
dC4 <- rCfv(Cseq,parms4) / Cseq

plot(dC1~log(Cseq),type="l",col="orange")
lines(dC2~log(Cseq),col="red")
lines(dC3~log(Cseq),col="green")
lines(dC4~log(Cseq),col="blue")
abline(h=0,col="black",lty=2)

Cstar <- sapply(list(parms1,parms2,parms3,parms4),Cstarf)
points(log(Cstar),rep(0,length(Cstar)))
  # next up: function for Cstar calculation over vector of different z

# Growth curves - discrete ------------------------------------------------

# requires at least 3 species
parmsX <- parms
parmsX$discrete <- TRUE
parmsX$tT <- 24 * 7 * 52 / 12
parmsX$zsig <- 0
parmsX$k <- 10
# parmsX$bd$mu$b0 <- parmsX$bd$mu$b0 - 5
parmsX$phi_E <- 0.1

now1 <- log( RCfv(Cseq,parmsX) / Cseq ) # *per-capita* growth

plot(now1~log(Cseq),type="l")
abline(4,-1,lty=2,col="red")
abline(h=0,lty=2,col="blue")

parmsX$tT <- 24*7*52
parmsX$nt <- 24*7*52
parmsX$sS <- 12
NarrD <- popint(parmsX)
par(mfrow=c(1,1))
plot(log(NarrD[,3]),type="l")
matplot(lo(NarrD[,-1]),type="l")

# bdt <- with(parmsX, btf(t=0, bd, M, parmsX))
# exp(-bdt$mu * parmsX$tT/parmsX$sS)
# function to plot functional responses? d_chain(y=c())


