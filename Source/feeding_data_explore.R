### Explore functional response parameter distributions from Rall et al. 2010 ###

require(plyr)
require(lme4)

# frf <- read.csv("Data/DatabaseXX2_corr.csv",header=T)
frf <- read.csv("Data/DatabaseXX2_Apr2016_rechecked.csv",header=T)

k <- 8.6173303 * 10^-5
T0 <- 293.15
frf$temperature.kelvin <- frf$temperature.degree.celcius + T0
fr <- subset(frf,select=c("publication.short",
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
             "predator.mass.g"="cmass",
             "prey.mass.g"="rmass",
             "attack.rate"="a",
             "handling.time"="h",
             "half.saturation.density"="Rhalf",
             "maximum.ingestion.rate"="fmax",
             "hill.exponent"="hill",
             "temperature.kelvin"="Tk"
             ))

fr$Tr <- with(fr, (Tk-T0)/(k*Tk*T0))
ma <- lmer(log(a) ~  log(cmass) + log(rmass) + Tr + (1|pub), data=fr)
mh <- lmer(log(h) ~  log(cmass) + log(rmass) + Tr + (1|pub), data=fr)

coef_a <- as.matrix(coef(ma)[[1]][as.numeric(fr$pub),])
coef_h <- as.matrix(coef(mh)[[1]][as.numeric(fr$pub),])

modmat_a <- model.matrix(ma)
modmat_h <- model.matrix(mh)

varnames <- colnames(coef_a)
notT <- varnames!="Tr"
notc <- varnames!="log(cmass)"
notr <- varnames!="log(rmass)"
notcr <- !varnames %in% c("log(cmass)","log(rmass)")

ares_T <- log(fr$a) - rowSums(coef_a[,notT] * modmat_a[,notT] )
ares_c <- log(fr$a) - rowSums(coef_a[,notc] * modmat_a[,notc] )
ares_r <- log(fr$a) - rowSums(coef_a[,notr] * modmat_a[,notr] )
ares_cr <- log(fr$a) - rowSums(coef_a[,notcr] * modmat_a[,notcr] )

hres_T <- log(fr$h) - rowSums(coef_h[,notT] * modmat_h[,notT] )
hres_c <- log(fr$h) - rowSums(coef_h[,notc] * modmat_h[,notc] )
hres_r <- log(fr$h) - rowSums(coef_h[,notr] * modmat_h[,notr] )
hres_cr <- log(fr$h) - rowSums(coef_h[,notcr] * modmat_h[,notcr] )

quickplot <- function(x,y){
  plot(y~x,pch="+")
  lines(supsmu(x,y),col="red")
}

par(mfrow=c(1,4))
with(fr, quickplot(Tk,ares_T))
with(fr, quickplot(log(cmass),ares_c))
with(fr, quickplot(log(rmass),ares_r))
with(fr, quickplot(Tk,hres_T))
with(fr, quickplot(log(cmass),hres_c))
with(fr, quickplot(log(rmass),hres_r))

par(mfrow=c(1,2))
with(fr, quickplot(log(cmass/rmass),ares_r))
with(fr, quickplot(log(cmass/rmass),hres_r))
  # the bigger the consumer relative to its prey,
  # the more slowly it attacks it but the faster it consumes it?

