##########################################################
# Royama plots for different forms of density-dependence #
##########################################################

r_expsurv <- function(X,z,Smax=0.25,d=0.05,B=10){
  log(Smax) - d*exp(X) + log(1+B*exp(z))
}
  # Barraquand & Yoccoz 2013, equation (6)
  # S decreases exponentially with population size (rate d)
  # B increases exponentially [my choice] with climate

curve(r_expsurv(x,z=0),xlim=c(-0.5,5))
curve(r_expsurv(x,z=-1),add=T,col="blue")
curve(r_expsurv(x,z=1),add=T,col="red")
abline(h=0,lty=3)
