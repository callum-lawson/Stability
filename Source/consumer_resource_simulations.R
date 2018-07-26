### Simulation of structured consumer-resource chains ###

# Trial runs --------------------------------------------------------------

### Inputs (parms)
source("Source/Consumer_resource_functions.R")

sparms = list(
  chainlength = 3,
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
  tT = 24*7,
  nt = 24*7, # 52,
  sS = 7*52, 
    # number of seasons over time series
  bdt = NULL,   
  nstart = 1
    # c(1,1,2, 1,1,1)
)

zparms <- list(
  zmu = 20, 
  zsig = 0,
  zl = 24*7*52
)

bc <- c(
  v = 1,   # max flow rate = k *grams* per m^2 per hour
  k = 100,    # grams per m^2
  psi = 0,   # interference:handling time ratio
  omega = 1, # relative feeding rate of y2
  phi_E = 0, # relative death rate of eggs
  phi_m = 1, # relative death rate of y2
  u_E = 1,   # odds ratio of y1:y2 at equilibrium
  m_E = 1,   # migration rate 
  u_m = 1,
  m_m = 0.1,
  tau_E = 24 * 7, # lags in migration functions
  tau_m = 0
)
  # phi and omega could instead by controlled by body masses
  #   (in this case, phi can be fraction of adult body mass)

# bhat <- readRDS("Output/rate_parameters_simulated_26Jul2018.rds")
# bhat <- bdselect(bhat,bpos=rep(5,2))
bhat <- readRDS("Output/rate_parameters_marginal_26Jul2018.rds")
bhat <- bdselect(bhat,bpos=rep(1,10))

iparms <- iparmf(bhat,sparms)
parms <- c(sparms,iparms,zparms,bc)

Cmin <- -3 
Cmax <- 3
nC <- 50
Cseq <- 10^seq(Cmin,Cmax,length.out=nC)

trial <- popint(parms)
matplot(log(trial[,-1]),type="l")

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
parms2$zmu <- parms1$zmu + 10
dC2 <- rCfv(Cseq,parms2) / Cseq

parms3 <- parms3b <- parms1
parms3b$zsig <- parms1$zsig + 5
parms3$bdt <- with(parms3b, rate_int_l(bd=bd,bn=names(bd),parms=parms3b))
dC3 <- rCfv(Cseq,parms3) / Cseq 
  # parms1 because don't want fluctuating for C* calculation

parms4 <- parms1
parms4$zmu <- parms1$zmu - 10
dC4 <- rCfv(Cseq,parms4) / Cseq

plot(dC1~log(Cseq),type="l",col="orange")
lines(dC2~log(Cseq),col="red")
lines(dC3~log(Cseq),col="green")
lines(dC4~log(Cseq),col="blue")
abline(h=0,col="black",lty=2)

Cstar <- sapply(list(parms1,parms2,parms3,parms4),Cstarf)
points(log(Cstar),rep(0,length(Cstar)))
  # next up: function for Cstar calculation over vector of different z

spectrum <- function(Rstar,alpha,a,h,mu){
  lrmax <- log( f(R=Rstar,C=1,a,h,alpha) ) # 1 -> per-capita (doesn't include w effects)
  lrmin <- log ( mu )
  lrmax - lrmin
}

with(with(parms, btf(t=0,bd,M,parms1)),spectrum(Rstar=bc["k"],alpha,a,h,mu))
with(with(parms, btf(t=0,bd,M,parms2)),spectrum(Rstar=bc["k"],alpha,a,h,mu))
with(with(parms, btf(t=0,bd,M,parms4)),spectrum(Rstar=bc["k"],alpha,a,h,mu))

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

# Generalism --------------------------------------------------------------

CseqG <- 10^seq(-3,3,length.out=100)

rCseqf <- function(parms){
  Cmin <- -3 
  Cmax <- 3
  nC <- 100
  Cseq <- 10^seq(Cmin,Cmax,length.out=nC)
  rCfv(Cseq,parms) / Cseq
}

rCquick <- function(zmu,parms){
  newparms <- parms
  # newparms$bd$h$b0 <- newparms$bd$h$b0 + 5
  newparms$k <- newparms$k * 5
  newparms$zmu <- zmu
  rCseqf(newparms)
}

parmsG1 <- parmsG2 <- parmsS1 <- parmsS2 <- parms
parmsG2$generalist <- FALSE
  # no need to prevent migration as dC = instantaneous rate

bhatS1 <- bdselect(bhat,bpos=c(1,2))
bhatS2 <- bdselect(bhat,bpos=c(1,3))
sparmsS <- sparms
sparmsS$nchain <- 1
iparmsS1 <- iparmf(bhatS1,sparmsS)
iparmsS2 <- iparmf(bhatS2,sparmsS)
parmsS1 <- c(sparmsS,iparmsS1,zparms,bc)
parmsS2 <- c(sparmsS,iparmsS2,zparms,bc)

zmuseq <- c(-5,0,5)

dCG1 <- sapply(zmuseq,rCquick,parms=parmsG1)
dCG2 <- sapply(zmuseq,rCquick,parms=parmsG2)
dCS1 <- sapply(zmuseq,rCquick,parms=parmsS1)
dCS2 <- sapply(zmuseq,rCquick,parms=parmsS2)

pdf(paste0("Plots/generalist_growthcurves_",format(Sys.Date(),"%d%b%Y"),".pdf"),width=14,height=14)
matplot(log(CseqG),cbind(dCG1,dCG2,dCS1,dCS2),type="l",
        col=rep(c("orange","green","red","blue"),each=3),
        lty=rep(c(2,1,2),times=4)
)
abline(h=0,col="black",lty=2)
dev.off()

# Fluctuation speed -------------------------------------------------------

sparms$discrete <- FALSE # TRUE

zlmin <- 0
zlmax <- 2
nzl <- 10
zlseq <- 24 * 10 ^ seq(zlmin,zlmax,length.out=nzl)
Cpos <- parms$Ya + 1 # +1 for time variable

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

# Cycles ------------------------------------------------------------------

# Rough exploration

bddd <- with(parms, btf(t=0,bd,M,parms))

# bhat <- bdselect(bhat,bpos=rep(5,2))
# bhat$alpha$bm <- 0.1 
# bhat$a$bm <- 0.1 
# bhat$h$bm <- bhat$h$bm + 1
# bhat$mu$b0 <- bhat$mu$b0 + 1
# bhat$mu$bm <- bhat$mu$bm + 1

# attach(parms)
# y <- y0
# t <- t0
trial <- popint(parms)
matplot(log(trial[,-c(1,2)]),type="l")
lines(ifelse(bigger,-10,10),col="blue",lty=2)
strial <- trial[-parms$nt,]

yo <- apply(trial[,-c(1,2)],2,function(x) diff(log(x)))
matplot(yo,type="l")
abline(h=0,col="blue",lty=3)
bigger <- abs(yo[,1]) > abs(yo[,2])

# plot(trial[,4]~trial[,3],type="l")
plot(log(trial[,4])~log(trial[,3]),type="l")

plot(yo[,2]~log(strial[,3]),type="l")
plot(yo[,2]~log(strial[,4]),type="l")

plot(yo[,2]~log(strial[,3]/strial[,4]),type="l")
peep1 <- 1:1000
peep2 <- 1000:2000
peep3 <- 2000:3000
peep4 <- 2000:2500
peep5 <- 2500:3000
plot(yo[peep3,2]~log(strial[peep3,3]/strial[peep3,4]),type="l",col="red",lty=3)
lines(yo[peep1,2]~log(strial[peep1,3]/strial[peep1,4]),lty=1)
lines(yo[peep2,2]~log(strial[peep2,3]/strial[peep2,4]),col="blue",lty=2)
lines(yo[peep4,2]~log(strial[peep4,3]/strial[peep4,4]),col="orange",lty=2)
lines(yo[peep5,2]~log(strial[peep5,3]/strial[peep5,4]),col="green",lty=2)

lRt <- log(trial[-parms$nt,3])
lRmed <- median(lRt)
bigR <- lRt > -3 # lRt > lRmed
plot(yo[bigR,2]~log(strial[bigR,4]),col="blue",type="l")
plot(yo[!bigR,2]~log(strial[!bigR,4]),col="red",type="l")

plot(yo[,2]~log(strial[,4]),type="n")
lines(yo[bigR,2]~log(strial[bigR,4]),col="blue")
lines(yo[!bigR,2]~log(strial[!bigR,4]),col="red")

parmsR <- parms
parmsR$tT = 24*7
parmsR$nT = 2
parmsR1 <- parmsR2 <- parmsR3 <- parmsR
parmsR1$y0 <- c(1,0.1,1)
parmsR2$y0 <- c(1,1,1)
parmsR3$y0 <- c(1,10,1)

hi1 <- popint(parmsR1)
hi2 <- popint(parmsR2)
hi3 <- popint(parmsR3)

trial <- function(y1,y2){
  ds <- unlist(romac_phaser(t = 0, y = c(y1,y2), parameters = pdm[1,]))
  diff(abs(ds))
}
  
uniroot.all(trial, interval=c(0,10^10), y1=0.1)

hi <- ode(y = c(1,1), times = seq(0,10^7,length.out=1000), func = romac_phaser, parms = pdm[1,])
matplot(log(hi[,-1]),type="l")

pstarsearch <- Vectorize(
  function(N,C1,C2){
    require(rootSolve)
    try1 <- try(
      root1 <- uniroot.all(hi5, interval=c(0,1), N=N, C1=C1, C2=C2)
    )
    if(class(try1)=="try-error" | length(try1)==0){
      return(NA)
    }
    else{
      return(root1) # returns p
    } 
  },
  vectorize.args=c("N")
)

RCof <- function(C,parms){
  require(rootSolve)
  with(parms, {
    y0[Ycc] <- C / nchain # split C over chains (p=0.5)
    Rstar <- steady(y=y0,
                    parms=parms,
                    fun=d_web,
                    times=c(0,Inf),
                    method="runsteady",
                    hold=TRUE
    )$y
    sum( unlist(d_web(t=0, y=Rstar, parms))[Ycc] )
  })
}


