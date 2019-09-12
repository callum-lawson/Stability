### Do random interaction matrices produce interaction patterns from Barbier et al. 2019? ###

require("rootSolve")
require("phaseR")

dLV <- function(t,y,parms){
  with(parms, {
    dy <- y * r * (1 - alpha %*% y/k) # top row = ii,ij,ik.. (i.e. i effects)
    list(dy)
  })
}

n_sim <- 1000
n_species <- 20
mu_alpha <- -3
sd_alpha <- 0.1

r <- 1 # sets the timescale; no effect on equilibria
mu_k <- 0
sd_k <- 0.1
set.seed(1)
seq_k <- exp(rnorm(n_species,mu_k,sd_k))
mat_krat <- t(outer(seq_k,seq_k,"/")) # alpha weighted by j/i, not i/j, so transpose
y0 <- rep(1,n_species)

rank_delta_alpha <- rank_delta_beta <- array(dim=c(n_species,n_species,n_sim))

for(i in 1:n_sim){
  
  set.seed(i)
  seq_alpha <- plogis(rnorm(n_species^2,mu_alpha,sd_alpha))
  alpha <- mat_alpha <- matrix(seq_alpha,nr=n_species,nc=n_species)
  mat_beta <- mat_krat * mat_alpha
  diag(alpha) <- 1
  diag(mat_alpha) <- NA
  diag(mat_beta) <- NA
  mat_delta_alpha <- log(mat_alpha) - log(mean(mat_alpha,na.rm=TRUE))
  mat_delta_beta  <- log(mat_beta)  - log(mean(mat_beta,na.rm=TRUE))
  
  parms <- list(r=r, k=seq_k, alpha=alpha)
  
  ystar <- runsteady(y=y0,parms=parms,fun=dLV,times=c(0,Inf))$y

  log_ystar <- log(ystar)
  ord_ystar <- order(log_ystar)
  log_eta <- log(ystar) - log(seq_k)
  ord_eta <- order(log_eta)
  
  rank_delta_alpha[,,i] <- mat_delta_alpha[ord_ystar,ord_ystar]
  rank_delta_beta[,,i]  <- mat_delta_beta[ord_eta,ord_eta]
  
}

par(mfrow=c(1,2))
mean_delta_alpha <- apply(rank_delta_alpha,c(1,2),mean)
image(mean_delta_alpha)

mean_delta_beta  <- apply(rank_delta_beta,c(1,2),mean)
image(mean_delta_beta)

cor_delta_beta <- array(dim=c(n_species,n_species,n_sim))
for (i in 1:n_sim){
  cor_delta_beta[,,i] <- cor(rank_delta_beta[,,i],use="pairwise.complete.obs")
  cor_delta_beta[,,i][upper.tri(cor_delta_beta[,,i], diag=TRUE)] <- NA
}
cor_delta_beta[cor_delta_beta==1] <- NA
mean_cor_delta_beta <- apply(cor_delta_beta,c(1,2),mean)
image(mean_cor_delta_beta)

# extinct <- vector("logical",length=n_sim)
# tau <- 10^-16
# extinct[i] <- TRUE %in% (ystar < tau)
# ord_delta_beta_coexist <- rank_delta_beta
# ord_delta_beta_coexist[,,extinct] <- NA
# mean_delta_beta_coexist <- apply(rank_delta_beta_coexist,c(1,2),mean,na.rm=TRUE)
# image(mean_delta_beta_coexist)