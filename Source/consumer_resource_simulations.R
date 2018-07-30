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
  tT = 24 * 7,
  nt = 24 * 7 * 10,
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
  v = 10,   # max input rate = vk *grams* per m^2 per hour
  k = 10,    # grams per m^2
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

# bhat <- readRDS("Output/rate_parameters_simulated_27Jul2018.rds")
# bhat <- bdselect(bhat,bpos=c(1,2,1,3))

bhat <- readRDS("Output/rate_parameters_marginal_27Jul2018.rds")
bhat <- bdselect(bhat,bpos=rep(1,2))

bhat$a$bz <- 0
bhat$h$bz <- 0
bhat$alpha$bz <- 0
bhat$mu$bz[1] <- 0

bhat$h$b0 <- bhat$h$b0 - 1.25 * 3.2089

iparms <- iparmf(bhat,sparms)
parms <- c(sparms,iparms,zparms,bc)

Cmin <- -3 
Cmax <- 3
nC <- 100
Cseq <- 10^seq(Cmin,Cmax,length.out=nC)

trial <- popint(parms)
matplot(log(trial[,-1]),type="l")

# parms2 <- parms
# newstart <- runsteady(y = parms$y0, time = c(0,Inf), func = d_web, parms = parms)
# parms2$y0 <- newstart$y
# parms2$zmu <- 25
# trial2 <- popint(parms2)
# matplot(log(trial2[,-1]),type="l")
  # still fluctuates when transitioning from one equilibrium to another

bddd <- with(parms, btf(t=0,bd,M,parms))

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

# Growth rates of basal consumer ------------------------------------------

sat <- function(N,k,a,alpha,mu){
  alpha * a * ( k - a*N ) - mu 
}
conk <- function(k,a,alpha,mu){
  1/a * (k-mu/(alpha*a))
}
conmax <- function(k,a,alpha,mu){
  alpha*a*k - mu
}
relcurve <- function(N,k,a,alpha,mu){
  newN <- N * conk(k,a,alpha,mu)
  sat(newN,k,a,alpha,mu)/conmax(k,a,alpha,mu)
}

v <- 50
k <- 50
a <- bddd[1,]$a
alpha <- bddd[1,]$alpha
mu <- bddd[1,]$mu

curve(sat(N=x,k,a,alpha,mu)*x,xlim=c(0,200))
curve(sat(N=x,0.5*k,a,alpha,mu)*x,add=T,col="red")
curve(sat(N=x,k,a,alpha,mu)*x,add=T,col="blue")

r <- 1
K <- 1000
curve(x*r*(1-x/K),add=TRUE,col="grey")

curve(relcurve(N=x,k,a,alpha,mu),xlim=c(0,1),ylim=c(0,1))
curve(relcurve(N=x,0.5*k,a,alpha,mu),add=T,col="red")    
curve(relcurve(N=x,k,0.5*a,alpha,mu),add=T,col="blue")

# Cycles ------------------------------------------------------------------

# Rough exploration

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


