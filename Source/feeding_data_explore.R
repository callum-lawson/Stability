### Explore functional response parameter distributions from Rall et al. 2010 ###

require(plyr)
require(lme4)
source("Source/consumer_resource_functions.R")

# Load feeding rates ------------------------------------------------------

# frf <- read.csv("Data/DatabaseXX2_corr.csv",header=T)
frf <- read.csv("Data/DatabaseXX2_Apr2016_rechecked.csv",header=T)

T0 <- 273.15 # 0°C = 293.15 K 
z0 <- 293.15 # intercept for models (K)
frf$temperature.kelvin <- frf$temperature.degree.celcius + T0
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
                         "temperature.kelvin"
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
             "temperature.kelvin"="Tk"
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
fr$Tr <- with(fr, arrtemp(Tk,z0=z0) / log(10))
  # using log10 to improve model convergence
  # log_a(x) = log_b(x)/log_b(a)
  # changes intercept but not other parameters
fr$al <- with(fr, log10(a * 60^2) )
fr$hl <- with(fr, log10(1/(h * 60^2)) )
  # rates per hour instead of per second
fr$cml <- with(fr, log10(cmass * 1000) ) 
# mass in mg
fr$rml <- with(fr, log10(cmass/rmass) )

ma <- lmer(al ~ cml + rml + Tr + (1|pub) + (1+cml+Tr|group), data=fr)
mh <- lmer(hl ~ cml + rml + Tr + (1|pub) + (1+cml+Tr|group), data=fr)

# Analyse metabolic rates -------------------------------------------------

mr[mr==-999.9] <- NA
mr$Tk <- mr$Tc + T0
mr$Tr <- arrtemp(mr$Tk, z0=z0) / log(10)
mr$cml <- log10(mr$cmass * 1000)
  # mass in mg
qlogis10 <- function(x) log10(x/(1-x))
mr$el <- qlogis10(mr$alpha)
mr$ml <- log10(mr$mu)

me <- lmer(el ~ cml + Tr + (1+Tr+cml|guild) + (1|group) + (1|pub), data=mr)
mm <- lmer(ml ~ cml + Tr + (1+Tr+cml|guild) + (1|group) + (1|pub), data=mr)
  # mass as offset?

# Plot residuals ----------------------------------------------------------

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

reslist <- lapply(mlist,resd)

par(mfrow=c(2,3))
resplots(reslist$ma)
resplots(reslist$mh)
par(mfrow=c(2,2))
resplots(reslist$me)
resplots(reslist$mm)

# Extract and save parameters ---------------------------------------------

fixadj <- function(m){
  fixcoef <- fixef(m)
  fixcoef["(Intercept)"] <- fixcoef["(Intercept)"] * log(10)
  return(fixcoef)
}
# convert intercept from base 10 to base e
# (needs to be done for any intercept-related parameters)

mlist <- list(ma=ma,mh=mh,me=me,mm=mm)
( fixlist <- lapply(mlist,fixadj) )

saveRDS(fixlist,
        paste0("Output/rate_parameters_marginal_",
               format(Sys.Date(),"%d%b%Y"),
               ".rds")
        )


# Judge the size of temperature effects -----------------------------------

tempchange <- function(m){
  exp(diff(fixef(m)["Tr"]*arrtemp(c(0,5),z0=z0)))
}

tm <- sapply(mlist,tempchange)

baseparrs <- function(m){
  exp(fixef(m)["(Intercept)"])
}

bp <- sapply(mlist,baseparrs)

curve(10^-1*bp[3]*bp[1]*exp(x)/(1+bp[1]*(1/bp[2])*exp(x))-bp[4],xlim=c(-2,2))
curve(10^-1*bp[3]*tm[3]*bp[1]*tm[1]*exp(x)/(1+bp[1]*tm[1]*(1/(bp[2]*tm[2]))*exp(x))-bp[4]*tm[4],add=T,col="red")
  # 10^-1 is aribrary scaler to match intake and mortality rates

tm[4] # mortality
tm[2] * tm[3] # high food
tm[1] * tm[3] # low food
  # when little food, increasing efficiency / attack makes little difference 
  #   relative to mortality
  # but when lots of food, handling time can make a big difference 
  #   (big food multiplier, whereas mortality still has same base rate)


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

