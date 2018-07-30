### Simulation of structured consumer-resource chains ###

# Trial runs --------------------------------------------------------------

### Inputs (parms)
source("Source/Consumer_resource_functions.R")

sparms = list(
  chainlength = 3,
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

bhat <- readRDS("Output/rate_parameters_simulated_27Jul2018.rds")
bhat <- bdselect(bhat,bpos=c(1,2,1,3))

# bhat <- readRDS("Output/rate_parameters_marginal_27Jul2018.rds")
# bhat <- bdselect(bhat,bpos=rep(1,4))

bhat$a$bz <- 0
bhat$h$bz <- 0
bhat$alpha$bz <- 0
bhat$mu$bz[1] <- 0

# bhat$h$b0 <- bhat$h$b0 - 1.25 * 3.2089

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
  # newparms$k <- newparms$k * 5
  newparms$zmu <- zmu
  rCseqf(newparms)
}

parmsG1 <- parmsG2 <- parmsS1 <- parmsS2 <- parms
parmsG2$generalist <- FALSE
  # no need to prevent migration as dC = instantaneous rate

bhatS1 <- bdselect(bhat,bpos=1:2)
bhatS2 <- bdselect(bhat,bpos=3:4)
sparmsS <- sparms
sparmsS$nchain <- 1
iparmsS1 <- iparmf(bhatS1,sparmsS)
iparmsS2 <- iparmf(bhatS2,sparmsS)
parmsS1 <- c(sparmsS,iparmsS1,zparms,bc)
parmsS2 <- c(sparmsS,iparmsS2,zparms,bc)

zmuseq <- c(15,20,25)

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



