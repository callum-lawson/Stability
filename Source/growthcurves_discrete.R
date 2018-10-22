### Simulation of structured consumer-resource chains ###

# Trial runs --------------------------------------------------------------

### Inputs (parms)
source("Source/consumer_resource_functions.R")

sparms = list(
  chainlength = 3,
    # single number or vector of length Ya
    # eggs (storage structure) always start at 0
  nchain = 1,
  store = TRUE,
    # if FALSE, births are allocated directly to feeders
  slevel = c("consumer","resource")[1],
  storetype = c("diffuse","feeding")[2],
    # no "selective" because too complex (ESS work)
  movetype = c("diffuse","feeding","selective")[3],
    # if FALSE, is mixed specialist
  generalist = FALSE,
    # does storage operate through births (TRUE) or diffusion (FALSE)?
  discrete = TRUE,
  mbase = 0.001, # 1mg
  morder = 2,
  tT = 24 * 365 * 10 * 10,
  nt = 24 * 365,
  sS = 10 * 10, 
    # *number* of seasons over time series
  bdt = NULL,   
  nstart = 0.1
    # c(1,1,2, 1,1,1)
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
  phi_E = 0.001, # relative death rate of eggs
  phi_m = 1, # relative death rate of y2
  u_E = 1,   # odds ratio of y1:y2 at equilibrium
  m_E = 0,   # migration rate 
  u_m = 1,
  m_m = 0.1,
  tau_E = 0, # lags in migration functions - 
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

# bhat$a$b0 <- bhat$mu$b0 + 2
# bhat$h$b0 <- bhat$h$b0 - 10 # 1.25 * 3.2089

iparms <- iparmf(bhat,sparms)
parms <- c(sparms,iparms,zparms,bc)

Cmin <- -0.5 
Cmax <- -0.35
nC <- 100
Cseq <- 10^seq(Cmin,Cmax,length.out=nC)

trial <- popint(parms)
matplot(log10(trial[,-1]),type="l")

bddd <- with(parms, btf(t=0,bd,M,parms))

# Growth curves - discrete ------------------------------------------------

# requires at least 3 species
parmsh <- parmsm <- parmsl <- parms
parmsh$zmu <- parms$zmu + 5
parmsl$zmu <- parms$zmu - 5

# *per-capita* growth
nowh <- log( RCfv(Cseq,parmsh) / Cseq ) 
nowm <- log( RCfv(Cseq,parmsm) / Cseq )
nowl <- log( RCfv(Cseq,parmsl) / Cseq ) 

nowjj <- cbind(nowh,nowm,nowl)

matplot(log(Cseq), nowjj, type="l",col=c("red","orange","blue"))
abline(-1,-1,lty=2,col="green")
abline(h=0,lty=2,col="grey")
lines(log(Cseq), apply(nowjj[,-2],1,mean),col="purple")

# bdt <- with(parmsX, btf(t=0, bd, M, parmsX))
# exp(-bdt$mu * parmsX$tT/parmsX$sS)
# function to plot functional responses? d_chain(y=c())

nh <- popint(parmsh)
nm <- popint(parmsm)
nl <- popint(parmsl)
njj <- data.frame(t=parms$tseq,h=log10(nh[,"C2"]),m=log10(nm[,"C2"]),l=log10(nl[,"C2"]))

require(tidyr)
require(ggplot2)

ndd <- gather(njj, "Tr", "n", 2:4) 

matplot(parmsh$tseq,njj[,-1],type="l",col=c("red","orange","blue"))


