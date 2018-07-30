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


