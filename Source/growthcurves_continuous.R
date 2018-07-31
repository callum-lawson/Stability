### Simulation of structured consumer-resource chains ###

# Trial runs --------------------------------------------------------------

### Inputs (parms)
source("Source/consumer_resource_functions.R")

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
  zl = 24
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

# bhat$a$bz <- 0
bhat$h$bz <- 0
bhat$alpha$bz <- 0
# bhat$mu$bz[2] <- 0
bhat$a$bz <- bhat$mu$bz
bhat$mu$bz[2] <- 0

# bhat$h$b0 <- bhat$h$b0 - 1.25 * 3.2089

iparms <- iparmf(bhat,sparms)
parms <- c(sparms,iparms,zparms,bc)

Cmin <- -3 
Cmax <- 3
nC <- 100
Cseq <- 10^seq(Cmin,Cmax,length.out=nC)

trial <- popint(parms)
newparms <- parms
newparms$zsig <- 5
newtrial <- popint(newparms)
matplot(log(newtrial[,-1]),type="l")
matplot(log(trial[,-1]),type="l",add=T)

# parms2 <- parms
# newstart <- runsteady(y = parms$y0, time = c(0,Inf), func = d_web, parms = parms)
# parms2$y0 <- newstart$y
# parms2$zmu <- 25
# trial2 <- popint(parms2)
# matplot(log(trial2[,-1]),type="l")
  # still fluctuates when transitioning from one equilibrium to another

bddd <- with(parms, btf(t=0,bd,M,parms))

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

matplot(log(Cseq), cbind(dC1,dC2,dC3,dC4), type="l", col=c("orange","red","green","blue"))
abline(h=0,col="black",lty=2)
lines(log(Cseq),apply(cbind(dC2,dC4),1,mean),col="purple",lty=1)

Cstar <- sapply(list(parms1,parms2,parms3,parms4),Cstarf)
# points(log(Cstar),rep(0,length(Cstar)))
points(mean(log(Cstar[c(2,4)])),0)
  # next up: function for Cstar calculation over vector of different z

spectrum <- function(Rstar,alpha,a,h,mu){
  lrmax <- log( f(R=Rstar,C=1,a,h,alpha) ) # 1 -> per-capita (doesn't include w effects)
  lrmin <- log ( mu )
  lrmax - lrmin
}

with(with(parms, btf(t=0,bd,M,parms1)),spectrum(Rstar=bc["k"],alpha,a,h,mu))
with(with(parms, btf(t=0,bd,M,parms2)),spectrum(Rstar=bc["k"],alpha,a,h,mu))
with(with(parms, btf(t=0,bd,M,parms4)),spectrum(Rstar=bc["k"],alpha,a,h,mu))


