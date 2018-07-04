### Simulation of structured consumer-resource chains ###

# Trial runs --------------------------------------------------------------

### Inputs (parms)
source("Source/Consumer_resource_functions.R")

sparms = list(
  chainlength = 2,
  # single number or vector of length Ya
  # eggs (storage structure) always start at 0
  nchain = 1,
  store = FALSE,
  slevel = c("consumer","resource")[1],
  storetype = c("diffuse","feeding")[2],
  # if FALSE, births are allocated directly to feeders
  movetype = c("diffuse","selective","feeding")[3],
  # if FALSE, is mixed specialist
  generalist = FALSE,
  # does storage operate through births (TRUE) or diffusion (FALSE)?
  discrete = FALSE, # taus become season lengths
  tT = 24*7*52,
  nt = 1000,
  sS = 1, # number of seasons over time series
  bdt=NULL,   #  bdt can be supplied here
  nstart = 1 # c(1,1,2, 1,1,1),
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
  phi_E = 0.1, # relative death rate of eggs
  phi_m = 1, # relative death rate of y2
  u_E = 1,   # rates in migration functions
  m_E = 1,   # u = odds ratio of y1:y2 at equilibrium
  u_m = 1,
  m_m = 0.01,
  tau_E = 0, # lags in migration functions
  tau_m = 0
)
  # phi and omega could instead by controlled by body masses
  #   (in this case, phi can be fraction of adult body mass)

bhat <- readRDS("Output/rate_parameters_simulated_21Jun2018.rds")
# bhat <- bdselect(bhat,bpos=rep(1:2,2))

iparms <- iparmf(bhat,sparms)
parms <- c(sparms,iparms,zparms,bc)

# Population growth curves ------------------------------------------------

# attach(parms)

# Cmin <- -5 # 1.4
# Cmax <- 5 # 1.7
# nC <- 100
# Cseq <- 10^seq(Cmin,Cmax,length.out=nC)
# 
# wow1 <- rCfv(Cseq,parms) / Cseq # *per-capita* growth 
# 
# parms2 <- parms
# parms2$zmu <- 5
# wow2 <- rCfv(Cseq,parms2) / Cseq 
# 
# parms3 <- parms
# parms3$zsig <- 5
# parms3$bdt <- rate_int_l(bd=bd,bn=names(bd),parms=parms3)
# wow3 <- rCfv(Cseq,parms3) / Cseq 
# 
# parms4 <- parms
# parms4$zmu <- -5
# wow4 <- rCfv(Cseq,parms4) / Cseq 
# 
# plot(wow1~log10(Cseq),type="l",col="orange")
# lines(wow2~log10(Cseq),col="red")
# lines(wow3~log10(Cseq),col="green")
# lines(wow4~log10(Cseq),col="blue")
# abline(h=0,col="black",lty=2)
# 
# Cstar <- sapply(list(parms,parms2,parms3,parms4),Cstarf)
# points(log10(Cstar),rep(0,length(Cstar)))
#   # next up: function for Cstar calculation over vector of different z

# Fluctuation speed -------------------------------------------------------

zlmin <- -1
zlmax <- 2
nzl <- 10
zlseq <- 24 * 10 ^ seq(zlmin,zlmax,length.out=nzl)
Cpos <- parms$nchain + 2

Narr <- with(parms, array(dim=c(nt,Cpos,nzl)))

for(i in 1:nzl){
  zparms <- list(
    zmu = 0, 
    zsig = 5,
    zl = zlseq[i]
  )
  parms <- c(sparms,iparms,zparms,bc)
  Narr[,,i] <- popint(parms)
}

gmf <- function(x) mean(log(x))
require(RColorBrewer)
mypalette <- rev(brewer.pal(nzl,"RdBu"))
nburn <- 500
rem <- 1:nburn
Cgm <- apply(Narr[-rem,Cpos,],2,gmf)
Csd <- apply(Narr[-rem,Cpos,],2,sd)
par(mfrow=c(1,1),mar=c(2,2,2,2))
matplot(parms$tseq[-rem],log10(Narr[-rem,Cpos,]),type="l",lty=1,col=mypalette)
par(mfrow=c(3,3))
for(i in 1:9){
  plot(parms$tseq[-rem],log10(Narr[-rem,Cpos,i]),col=mypalette[i],type="l")
}
par(mfrow=c(1,2),mar=c(2,2,2,2))
plot(Cgm~log10(zlseq/24),type="b")
plot(Csd~log10(zlseq/24),type="b")

zparms <- list(
  zmu = 0, 
  zsig = 0,
  zl = 1 # irrelevant
)
parms <- c(sparms,iparms,zparms,bc)
Cstarcons <- Cstarf(parms)

par(mfrow=c(1,1),mar=c(2,2,2,2))
matplot(parms$tseq[-rem],log10(Narr[-rem,Cpos,]),type="l",lty=1,col=mypalette)
abline(h=log10(Cstarcons),lty=2)

# Other -------------------------------------------------------------------

# popint <- function(y0,tseq,parms){
#   # if(nrow(bhat[])!=length(M)) stop("wrong masses or params")
#   # ! generalist + "resource"
# }
# 
# tseq <- seq(0,24*60,length.out=100)
# require(deSolve)
# if(parms$tau==0) lvar <- ode(y=y0,times=tseq,func=d_web,parms=parms)
# if(parms$tau>0) lvar <- dede(y=y0,times=tseq,func=d_web,parms=parms)