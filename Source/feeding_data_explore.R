### Explore functional response parameter distributions from Rall et al. 2010 ###

require(plyr)
require(lme4)
require(brms)
library(ggplot2)
source("Source/consumer_resource_functions.R")

# Load data ---------------------------------------------------------------

### Feeding rates

# frf <- read.csv("Data/DatabaseXX2_corr.csv",header=T)
frf <- read.csv("Data/DatabaseXX2_Apr2016_rechecked.csv",header=T)

fr <- subset(frf,select=c("publication.short",
                         "ecosystem.type",
                         "search.space",
                         "predator.met.group",
                         "predator.Class",
                         "consumption.type",
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
             "search.space" = "dimension",
             "predator.met.group" = "metgroup",
             "predator.Class" = "group",
             "consumption.type"="guild",
             "predator.mass.g"="cmass",
             "prey.mass.g"="rmass",
             "attack.rate"="a",
             "handling.time"="h",
             "half.saturation.density"="Rhalf",
             "maximum.ingestion.rate"="fmax",
             "hill.exponent"="hill",
             "temperature.degree.celcius"="Tc"
             ))

### Metabolic rates

mrf <- read.csv("Data/Lang_et_al_2017_data.csv",skip=1) 
  # assimilation efficiency and respiration data
mra <- subset(mrf,select=c(
  "taxonomic.group.consumer",
  "consumer.type",
  "body.size.gram",
  "temperature.degree.C",
  "assimilation.efficiency",
  "metabolism.J.h",
  "reference.short"
  )) 
mra <- rename(mra,replace=c(
  "taxonomic.group.consumer" = "group",
  "consumer.type" = "guild",
  "body.size.gram" = "cmass",
  "temperature.degree.C" = "Tc",
  "assimilation.efficiency" = "alpha",
  "metabolism.J.h" = "mu",
  "reference.short" = "pub"
))

# Collate variables -------------------------------------------------------

### Feeding data

# fr <- subset(fr,!metgroup %in% c("endovert","unicell"))
#   # temporarily remove small-sample-size groups
# fr$metsyst <- with(fr, droplevels(system:metgroup))
# frf$predator.met.group[frf$prey.Genus=="Paramecium"] <- "unicell"

fr$guild <- as.character(fr$guild)
fr$guild[fr$guild %in% c("bacterivore","fungivore")] <- "detritivore"
fr$guild[fr$guild=="parasitoid"] <- "carnivore"
fr$guild <- as.factor(fr$guild)
fr$gg <- with(fr, droplevels(group:guild))
fr$Tr <- with(fr, arrtemp(Tc))
fr$cml <- with(fr, log(cmass) ) 
  # mass in *grams*
fr$rml <- with(fr, cml - log(rmass) )
fr$al <- with(fr, log(a * 60^2) - cml )
  # attack rate in m^2 per hour instead of per second as in Yuanheng's database
  # -cml -> attack rate per consumer gram instead of per consumer individual
  # no need to adjust for resource weight because this accounted for in density units
  #   (i.e. R in g/m^2 instead of N/m^2)
fr$hl <- with(fr, log(h / 60^2) + rml )
  # handling times per hour instead of per second as in Yuanheng's database
  # fmax = R/C * M[R]/M[C] / h, so handling time (inverse) needs to be multiplied by *M[C]/M[R]*
  # each R individual yields *M[R] grams*, so if e.g. one R yields 2g, h per g = h/2
  # if each C individual weighs 2g, 1g of C would take twice as long to handle it (h'=2*h)
  # or: a*h*R needs to be dimensionless, a*R is in g[R]/(h*g[C]), so h needs to be in h*g[C]/g[R]

### Metabolic data

mra[mra==-999.9] <- NA
mra$gg <- with(mra, group:guild)
mra$Tr <- arrtemp(mra$Tc) 
mra$cml <- log(mra$cmass)
  # mass in *grams*

ef <- droplevels(subset(mra,!is.na(alpha)))
mr <- droplevels(subset(mra,!is.na(mu)))

ef$el <- qlogis(ef$alpha)
  # no offset because percentage, not rate
mr$ml <- log(mr$mu) - log(16000) - mr$cml 
  # original data on mortality rates (mu) in J/h
  # using 16 kJ of energy burns 1g of mass
  # so mortality in g/h = (rates in J/h) / 16000
  # protein = 16 kJ/g https://en.wikipedia.org/wiki/Specific_energy#Energy_density_of_food
  # glucose = 16 kJ/g https://en.wikipedia.org/wiki/Glucose#Energy_source 
  # (fat provides more as has higher energy density)
  # different estimate: 7 kJ/g (Yodzis & Innes 1992; de Roos, Box 3.4)
  # metabolic rates in database are given per *individual* - checked by comparison with:
  #   Anderson & Prestwich (1982) "Respiratory gas exchange in spiders", Table 3
  #   Nephila clavipes: 157 / 1000 [per ml instead of micro l] * 20.1 [ml 02 to J] 
  #     = 3.15 J/h per individual (848 mg); c.f. 3.03 J/h in Lang et al.

# Units:
# a: m^2 / (h g[C]) (area searched per g[C] per hour)
# h: h g[C] / g[R] (time taken by each g[C] to process each g[R])
# alpha: % (fraction of killed g[R] that is converted to new g[C])
# mu: /h (g[C] used per g[C] per hour)

# Feeding rate models -----------------------------------------------------

### Uncorrelated (lme4)

ma <- lmer(al ~ Tr + cml + rml + (1|gg) + (1|pub), data=fr)
  # 10 AIC better than (1|group) + (1|pub) 
  # 20 AIC better than (Tr + cml + rml|group) + (1|pub) 
  # not enough observations per group to estimate intercept-slope correlations
  # bigger body masses:
  #   - increase attack rates per consumer individual (Rall et al.)
  #   - decrease attack rates per consumer gram (splitting up more effective)

mh <- lmer(hl ~ Tr + cml + rml + (1|gg) + (1|pub), data=fr)
  # 20 AIC better than (1|guild) + (1|group) + (1|pub)
  # 60 AIC worse than (Tr + cml + rml|group) + (1|pub), and intercept / cml very different:
  # (Intercept)          Tr         cml         rml 
  # 1.25624752 -0.27798215  0.04896353  0.47241528 
  # large among-group variance; unchanged by fitting log(rmass) instead of log(cmass/rmass)

# mh2 <- lmer(hl ~ Tr + I(Tr^2) + cml + rml + (1 + Tr + I(Tr^2) + cml + rml | group) + (1 | pub), data=fr)
# AIC(mh,mh2)
#  # no evidence for U-shaped temperature response (Englund et al. 2011)
# Fitting rmass instead of (cmass/rmass) makes very little difference

### Uncorrelated (brms)

mah <- brm(
  cbind(al, hl) ~ Tr + cml + rml + (1 | p | gg) + (1 | q | pub),
  # prior = set_prior("normal(0,8)"),
  warmup = 1000, iter = 2000,
  data = fr, chains = 3, cores = 3
)
  # how to fit crossed random effects here?

# Assimilation models -----------------------------------------------------

me <- lmer(el ~ Tr + cml + (1 | gg) + (1|pub), data=ef, REML=F)
  # 20 AIC worse than (1 + Tr + cml | guild) + (1 + Tr + cml | group) + (1|pub) about 

# Metabolic rate models ---------------------------------------------------

mm <- lmer(ml ~ Tr + cml + (1 | gg) + (1|pub), data=mr)
  # 3 AIC better than (1 | group) + (1 | guild) + (1|pub)
  # 170 AIC worse than (1 + Tr + cml | gg) + (1|pub) but marginal params very simliar:
  # (Intercept)           Tr          cml  
  # -9.3173       0.6312      -0.3003  

# use brms instead?

# Extract and save parameters ---------------------------------------------

### Mean parameters

fixadj <- function(m){
  fixcoef <- fixef(m)
  if(! "rml" %in% names(fixcoef) ){
    fixcoef <- c(fixcoef,0)
    names(fixcoef)[length(fixcoef)] <- "rml"
      # add on zero body mass ratio effect if effect is absent
  }
  fixcoef <- rename(fixcoef,replace=c(
    "(Intercept)" = "b0",
    "Tr" = "bz",
    "cml" = "bm",
    "rml" = "br"
  ))
  return(as.data.frame(t(fixcoef)))
}
  # changes intercept to absolute (un-logged) scale
  # shifts intercept to consumer mass of 1g instead of 1mg
  # (needs to be done for any intercept-related parameters)

mlist <- list(a=ma,h=mh,alpha=me,mu=mm)
( fixlist <- lapply(mlist,fixadj) )
  
# Save parameters

saveRDS(mah,
        paste0("Output/mcmc_brms_",
               format(Sys.Date(),"%d%b%Y"),
               ".rds")
)

saveRDS(fixlist,
        paste0("Output/rate_parameters_marginal_",
               format(Sys.Date(),"%d%b%Y"),
               ".rds")
        )

# brms parameters

fixlist_brms <- fixlist
brms_coefs <- fixef(mah)[,"Estimate"]
fixlist_brms$a[1:4] <- brms_coefs[c(1,3:5)]
fixlist_brms$h[1:4] <- brms_coefs[c(2,6:8)]

saveRDS(fixlist_brms,
        paste0("Output/rate_parameters_marginal_brms_",
               format(Sys.Date(),"%d%b%Y"),
               ".rds")
)


cbind(t(fixlist$a),t(fixlist_brms$a))
cbind(t(fixlist$h),t(fixlist_brms$h))

### Randomly-drawn parameters

sdadj <- function(m,n){
  require(MASS)

  fixcoef <- unlist(fixadj(m))
  
  vcgroup <- VarCorr(m)$group
  sdgroup <- attr(vcgroup, "stddev")
  corrgroup <- attr(vcgroup, "correlation")
  
  sd1 <- sdgroup["(Intercept)"]
  sd2 <- sdgroup["Tr"]
  rho <- corrgroup["(Intercept)","Tr"]
  if(rho < -0.99 | rho > 0.99){
    rho <- 0
    # if can't estimate correlation, assume uncorrelated
  }
  
  sig <- matrix(c(sd1^2,rep(rho*sd1*sd2,2),sd2^2),nr=2,nc=2)
  simd1 <- mvrnorm(n=n, mu=c(fixcoef["b0"],fixcoef["bz"]), Sigma=sig)
  simd <- data.frame(simd1,bm=rep(fixcoef["bm"],n),br=rep(fixcoef["br"],n))

  return(simd)
}
 # incomplete - finish later
 
set.seed(1)
simdlist <- lapply(mlist,sdadj,n=10) 

# Save parameters

saveRDS(simdlist,
        paste0("Output/rate_parameters_simulated_",
               format(Sys.Date(),"%d%b%Y"),
               ".rds")
)

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
with(mlist, cor.test(resid(a),resid(h)))
  # in data, higher attack rates associated with smaller handling times
  # this association partly accounted for by mass and temperature
  # but even after accounting for this, higer a associated with lower h

# Exploratory tables ------------------------------------------------------

frtab <- with(fr, table(pub,gg))
mrtab <- with(mr, table(pub,gg))

table(rowSums(frtab>0))
table(rowSums(mrtab>0))
# almost all studies only have one taxononmic group
colSums(frtab)
colSums(mrtab)
# but each group contains many separate studies
frtab[frtab>0]
mrtab[mrtab>0]
# 25% of groups only have one study to them, but some have many studies
with(fr,table(group))

ucount <- function(x) length(unique(x))
utemp <- with(fr, tapply(Tr,pub,ucount))
umass <- with(fr, tapply(cml,pub,ucount))
uboth <- with(fr, tapply(paste(cml,Tr),pub,ucount))
table(umass>1)
table(utemp>1)
# among all studies, 50% have one pred, 60% have one temp, and 33% have one pred at one temp

# Exploratory graphs ------------------------------------------------------

frg <- subset(fr,group %in% c("Actinopterygii","Arachnida","Insecta"))
frg$group <- droplevels(frg$group)

fr$ares <- resid(ma)
fr$hres <- resid(mh)
with(fr, plot(ares~hres))
ef$eres <- resid(me)
mr$mres <- resid(mm)

ggplot(fr, aes(x=cml, y=ares)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method = "lm",alpha=0.1) + 
  theme_bw()
ggplot(fr, aes(x=cml, y=hres)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method = "lm",alpha=0.1) + 
  theme_bw()
ggplot(ef, aes(x=cml, y=eres)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method = "lm",alpha=0.1) + 
  theme_bw()
ggplot(mr, aes(x=cml, y=mres)) + 
  geom_point(aes(color=group,alpha=0.1,stroke=0)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method = "lm",alpha=0.1) + 
  coord_cartesian(ylim=c(-3.75,3.75)) +
  theme_bw()

ggplot(frg, aes(x=Tr, y=al)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method = "lm",alpha=0.1) + 
  theme_bw()
ggplot(frg, aes(x=cml, y=al)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method = "lm",alpha=0.1) + 
  theme_bw()
ggplot(frg, aes(x=rml, y=al)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method = "lm",alpha=0.1) + 
  theme_bw()
# facet_wrap(~ metsyst) +

ggplot(fr, aes(x=Tr, y=hl)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method = "lm",alpha=0.1) + 
  theme_bw()
ggplot(fr, aes(x=cml, y=hl)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method = "lm",alpha=0.1) + 
  theme_bw()
ggplot(fr, aes(x=rml, y=hl)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method = "lm",alpha=0.1) + 
  theme_bw()
# high handling times of Philippova.1988.IntRevueGesHydrobiol checked in paper, approximately correct

ggplot(mr, aes(x=Tr, y=ml)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method="lm", alpha=0.1) + 
  theme_bw()
ggplot(mr, aes(x=cml, y=ml)) + 
  geom_point(aes(color=group)) + 
  geom_smooth(aes(group=group,color=group,fill=group), method="lm", alpha=0.1) + 
  theme_bw()

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

# E0a <- as.numeric(
#   exp(predict(ma,data.frame(cmass=10^-4,rmass=10^-6,Tr=arrtemp(z0)),re.form=~0))
#   )
# E0h <- as.numeric(
#   exp(predict(mh,data.frame(cmass=10^-4,rmass=10^-6,Tr=arrtemp(z0)),re.form=~0))
#   )
#   # 100 mg consumer, 1 mg prey
# E1a <- as.numeric(fixef(ma)[names(fixef(ma))=="Tr"])
# E1h <- as.numeric(fixef(mh)[names(fixef(mh))=="Tr"])

