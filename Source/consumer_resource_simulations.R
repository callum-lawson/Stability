### Simulation of structured consumer-resource chains ###

# Trial runs --------------------------------------------------------------

### Inputs (parms)
source("Source/Consumer_resource_functions.R")

zparms <- list(
  zmu = 0, 
  zsig = 0,
  zl = 24*7
  )

sparms = list(
  chainlength = 3,
  # single number or vector of length Ya
  # eggs (storage structure) always start at 0
  nchain = 2,
  store = FALSE,
  slevel = c("consumer","resource")[1],
  storetype = c("diffuse","feeding")[1],
  # if FALSE, births are allocated directly to feeders
  movetype = c("diffuse","selective")[2],
  # if FALSE, is mixed specialist
  generalist = TRUE,
  # does storage operate through births (TRUE) or diffusion (FALSE)?
  nstart = c(1,1,2, 1,1,1)
)

# later parms
# lagtype <- c("none","random","continuous","discrete")

bc <- c(
  v = 0.1,     # max flow rate = k grams per m^2 per hour
  k = 10,    # 10g per m^2
  psi = 0,   # interference:handling time ratio
  phi = 1,   # relative death rate of eggs / y2
  omega = 1, # relative feeding rate of eggs
  u_E = 1,   # rates in migration functions
  m_E = 1,   # u = odds ratio of y1:y2 at equilibrium
  u_m = 1,
  m_m = 1,
  tau_E = 0, # lags in migration functions
  tau_m = 0
)
# phi and omega could instead by controlled by body masses
#   (in this case, phi can be fraction of adult body mass)

# if(nchain==2){
#   omega_new <- rep(1,iparms$Yb)
#   omega_new[iparms$Yr] <- bc$omega # Yr because 
#   bc$omega <- omega_new
# }
bhat <- readRDS("Output/rate_parameters_simulated_21Jun2018.rds")
# bhat <- bdselect(bhat,bpos=rep(1:2,2))

tparms <- list(
  t0 = 0,
  tT = 24*7*12,
  nt = 1000
)

iparms <- iparmf(bhat,sparms)
parms <- c(bc,zparms,sparms,iparms)
attach(parms)
t <- 0
y <- y0

tseq <- with(tparms, seq(t0,tT,length.out=nt))
require(deSolve)
hi <- ode(y=y0,times=tseq,func=d_web,parms=parms)
matplot(log(hi[,-1]),type="l")

# popint <- function(y0,tseq,parms){
#   # if(nrow(bhat[])!=length(M)) stop("wrong masses or params")
#   # ! generalist + "resource"
# }
# 
# tseq <- seq(0,24*60,length.out=100)
# require(deSolve)
# if(parms$tau==0) lvar <- ode(y=y0,times=tseq,func=d_web,parms=parms)
# if(parms$tau>0) lvar <- dede(y=y0,times=tseq,func=d_web,parms=parms)