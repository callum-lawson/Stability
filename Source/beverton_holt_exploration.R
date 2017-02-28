####################################################################
# Explore how random variation in Beverton-Holt mortality function #
# affects the distribution of population sizes                     #
####################################################################

# Functions ---------------------------------------------------------------

BH <- function(n_in,m0,m1){
  n_out <- n_in * exp(-m0) / ( 1 + (m1/m0)*(1-exp(-m0))*n_in )
  return(n_out)
}
# assuming T = 1

# Exploration -------------------------------------------------------------

ln0_min <- -5
ln0_max <- 0
nseq <- 10^3
neps <- 10^3
alpha_m <- 0.1
beta_m <- 1

# n0_seq <- exp(seq(ln0_min,ln0_max,length.out=nseq))
# eps_m <- seq(-1,1,length.out=neps)
n0_seq <- exp(rnorm(nseq,0,1.5))
eps_m <- rnorm(neps,0,2)
  
m0 <- m1 <- n0 <- n1 <- matrix(NA,nr=nseq,nc=neps)

m0[] <- rep(exp(outer(alpha_m,eps_m,"+")),each=nseq)
m1[] <- rep(exp(outer(beta_m,eps_m,"+")),each=nseq)
n0[] <- rep(n0_seq,times=neps)
n1 <- BH(n0,m0,m1)

par(mfrow=c(2,1))
hist(n0,breaks=100,prob=T)
hist(n1,breaks=100,prob=T)

hist(n1[1,],breaks=100)
hist(n1[2,],breaks=100)
hist(n1[3,],breaks=100)

matplot(n0[,1],n1,type="l")
abline(0,1,col="purple",lty=3)

(K <- BH(10^10,m0[1,],m1[1,]))
