### Explore functional response parameter distributions from Rall et al. 2010 ###

require(plyr)
require(lme4)
source("Source/consumer_resource_functions.R")

# Load feeding rates ------------------------------------------------------

# frf <- read.csv("Data/DatabaseXX2_corr.csv",header=T)
frf <- read.csv("Data/DatabaseXX2_Apr2016_rechecked.csv",header=T)

fr <- subset(frf,select=c("publication.short",
                         "ecosystem.type",
                         "predator.met.group",
                         "predator.mass.g",
                         "prey.mass.g",
                         "attack.rate",
                         "handling.time",
                         "half.saturation.density",
                         "maximum.ingestion.rate",
                         "hill.exponent",
                         "temperature.degree.celcius"
                         ))

fr <- rename(fr,replace=c(
             "publication.short" = "pub",
             "ecosystem.type" = "system",
             "predator.met.group" = "metgroup",
             "predator.mass.g"="cmass",
             "prey.mass.g"="rmass",
             "attack.rate"="a",
             "handling.time"="h",
             "half.saturation.density"="Rhalf",
             "maximum.ingestion.rate"="fmax",
             "hill.exponent"="hill",
             "temperature.degree.celcius"="Tc"
             ))


# Load metabolic rates ----------------------------------------------------

mrf <- read.csv("Data/Lang_et_al_2017_data.csv",skip=1) 
  # assimilation efficiency and respiration data
mr <- subset(mrf,select=c(
  "taxonomic.group.consumer",
  "consumer.type",
  "body.size.gram",
  "temperature.degree.C",
  "assimilation.efficiency",
  "metabolism.J.h",
  "reference.short"
  )) 
mr <- rename(mr,replace=c(
  "taxonomic.group.consumer" = "group",
  "consumer.type" = "guild",
  "body.size.gram" = "cmass",
  "temperature.degree.C" = "Tc",
  "assimilation.efficiency" = "alpha",
  "metabolism.J.h" = "mu",
  "reference.short" = "pub"
))

# Analyse feeding rates ---------------------------------------------------

fr <- subset(fr,!metgroup %in% c("endovert","unicell"))
  # temporarily remove small-sample-size groups
fr$group <- with(fr, droplevels(system:metgroup))
fr$Tr <- with(fr, arrtemp(Tc) / log(10))
  # using log10 to improve model convergence
  # log_a(x) = log_b(x)/log_b(a)
  # changes intercept but not other parameters
fr$al <- with(fr, log10(a * 60^2) )
  # attack rate in Yuanheng's database in m^2 per second
fr$hl <- with(fr, log10(h / 60^2) )
  # handling times per hour, instead of per second as in Yuanheng's database
fr$cml <- with(fr, log10(cmass) ) 
  # mass in mg
sigma <- 0.01
fr$rml <- with(fr, log10(sigma*cmass/rmass) )
  # sets intercept with consumer 100 times larger than resource

ma <- lmer(al ~ offset(cml) + cml + rml + Tr + (1|pub) + (1+cml+Tr|group), data=fr)
  # mass as offset -> attack and max feeding rates per gram of consumer
  # so attack rate in m^2 per hour per consumer gram
  # each gram of consumer covers 0.03 cm^2 per hour!
mh <- lmer(hl ~ offset(-(rml + log10(sigma))) 
           + cml + rml + Tr + (1|pub) + (1+cml+Tr|group), 
           data=fr)
  # handling time per resource gram per consumer gram, i.e. 
  #   time for 1 feeding consumer gram to process 1 resource gram
  # each resource individual yields *M[R] grams*, so handling time needs to be 
  #   scaled down by amount of grams of in each resource individual
  # fmax = R/C * M[R]/M[C]
  #   so handling time scaled by inverse (negative log)
  
# Analyse metabolic rates -------------------------------------------------

mr[mr==-999.9] <- NA
mr$Tr <- arrtemp(mr$Tc) / log(10)
mr$cml <- log10(mr$cmass)
  # original mass in grams
qlogis10 <- function(x) log10(x/(1-x))
mr$el <- qlogis10(mr$alpha)
mr$ml <- log10(mr$mu)

me <- lmer(el ~ cml + Tr + (1+Tr+cml|guild) + (1|group) + (1|pub), data=mr)
  # no offset because percentage, not rate
mm <- lmer(ml ~ offset(cml) + cml + Tr + (1+Tr+cml|guild) + (1|group) + (1|pub), data=mr)
  # J used by each consumer gram in each hour

# Extract and save parameters ---------------------------------------------

fixadj <- function(m){
  fixcoef <- fixef(m)
  fixcoef["(Intercept)"] <- 10 ^ fixcoef["(Intercept)"] 
    # * 1 ^ fixcoef["cml"]
    # check whether mass power is correct here (was 1000)
  if(! "rml" %in% names(fixcoef) ){
    fixcoef <- c(fixcoef,0)
    names(fixcoef)[length(fixcoef)] <- "rml"
      # add on zero body mass ratio effect if effect is absent
  }
  fixcoef <- rename(fixcoef,replace=c(
    "(Intercept)" = "b0",
    "cml" = "bm",
    "Tr" = "bz",
    "rml" = "br"
  ))
  return(as.list(fixcoef))
}
  # changes intercept to absolute (un-logged) scale
  # shifts intercept to consumer mass of 1g instead of 1mg
  # (needs to be done for any intercept-related parameters)

mlist <- list(a=ma,h=mh,alpha=me,mu=mm)
( fixlist <- lapply(mlist,fixadj) )

fixlist$mu$b0 <- fixlist$mu$b0 / 16
  # protein = 16 kJ/g = 16 J/mg
  # so lose 1 mg of mass for every 16 J of energy
  # (original mortality rates in J)
  # https://en.wikipedia.org/wiki/Specific_energy#Energy_density_of_food
  
### Mean parameters

saveRDS(fixlist,
        paste0("Output/rate_parameters_marginal_",
               format(Sys.Date(),"%d%b%Y"),
               ".rds")
        )

### Randomly-drawn parameters

sdadj <- function(m){
  sdgroup <- attr(VarCorr(mlist$a)$group, "stddev")
  fixcoef["(Intercept)"] <- 10 ^ fixcoef["(Intercept)"] 
  # * 1 ^ fixcoef["cml"]
  # check whether mass power is correct here (was 1000)
  fixcoef <- rename(fixcoef,replace=c(
    "(Intercept)" = "b0",
    "cml" = "bm",
    "Tr" = "bz",
    "rml" = "br"
  ))
  return(as.list(fixcoef))
}
 # incomplete - finish later

# Plot residuals ----------------------------------------------------------

### Functions

resd <- function(m){
  dat <- m@frame
  fix <- fixef(m)
  varnames <- names(fix)[-1] # delete intercept
  res <- residuals(m)
  adjd <- lapply(varnames,function(n){
    x = dat[,n]
    y = res + fix[n] * x 
    # total difference = residuals + fixed effect
    return(data.frame(x,y))
  })
  return(adjd)
}

quickplot <- function(x,y){
  plot(y~x,pch="+")
  lines(supsmu(x,y,bass=5),col="red")
}

resplots <- function(l){
  lapply(l,function(d) with(d, quickplot(x,y)))
} 

### Plots

reslist <- lapply(mlist,resd)

par(mfrow=c(2,3),mar=c(3,3,3,3))
resplots(reslist$a)
resplots(reslist$h)
par(mfrow=c(2,2))
resplots(reslist$alpha)
resplots(reslist$mu)

### Covariance

with(fr, quickplot(log10(a),log10(h)))
with(mlist, quickplot(resid(a),resid(h)))
  # in data, higher attack rates associated with smaller handling times
  # but this association largely accounted for by mass and temperature

# Judge the size of temperature effects -----------------------------------

tempchange <- function(m){
  exp(diff(fixef(m)["Tr"]*arrtemp(c(0,5),z0=z0)))
}

tm <- sapply(mlist,tempchange)

baseparrs <- function(m){
  exp(fixef(m)["(Intercept)"])
}

bp <- sapply(mlist,baseparrs)

# curve(10^-1*bp[3]*bp[1]*exp(x)/(1+bp[1]*(1/bp[2])*exp(x))-bp[4],xlim=c(-2,2))
# curve(10^-1*bp[3]*tm[3]*bp[1]*tm[1]*exp(x)/(1+bp[1]*tm[1]*(1/(bp[2]*tm[2]))*exp(x))-bp[4]*tm[4],add=T,col="red")
#   # 10^-1 is aribrary scaler to match intake and mortality rates
# 
# tm[4] # mortality
# tm[2] * tm[3] # high food
# tm[1] * tm[3] # low food
#   # when little food, increasing efficiency / attack makes little difference 
#   #   relative to mortality
#   # but when lots of food, handling time can make a big difference 
#   #   (big food multiplier, whereas mortality still has same base rate)

# Predictions for simulations ---------------------------------------------

E0a <- as.numeric(
  exp(predict(ma,data.frame(cmass=10^-4,rmass=10^-6,Tr=arrtemp(z0)),re.form=~0))
  )
E0h <- as.numeric(
  exp(predict(mh,data.frame(cmass=10^-4,rmass=10^-6,Tr=arrtemp(z0)),re.form=~0))
  )
  # 100 mg consumer, 1 mg prey
E1a <- as.numeric(fixef(ma)[names(fixef(ma))=="Tr"])
E1h <- as.numeric(fixef(mh)[names(fixef(mh))=="Tr"])

