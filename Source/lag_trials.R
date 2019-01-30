### Exploring one-species system with a lag, in 2D ###

require("deSolve")
require("rootSolve")
require("phaseR")

# Abiotic -----------------------------------------------------------------

derivs <- function(t, y, parms) {
  if (t < 1)
    dy <- 1
  else
    dy <- 1 - lagvalue(t - 1)
  list(c(dy))
}

yinit <- 1
times <- seq(0, 30, 0.1)
times_start <- times<=(30-1)
times_end <- times>=1

yout <- dede(y = yinit, times = times, func = derivs, parms = NULL)[,2]
yout_lag <- c(rep(NA,sum(!times_end)),yout[times_end])

plot(yout, type = "l", lwd = 2, main = "dy/dt = -y(t-1)")
plot(yout,yout_lag)
  # plotting flowfield not worthwhile becuase y is fully determined by x?

# plot(yout[times_end],yout_lag[times_start])

# Biotic ------------------------------------------------------------------

derivs2 <- function(t, y, parms) {
  if (t < 1)
    dy <- - y
  else
    dy <- lagvalue(t - 1) - y
  list(c(dy))
}

yout2 <- dede(y = yinit, times = times, func = derivs2, parms = NULL)[,2]
yout2_lag <- c(rep(NA,sum(!times_end)),yout2[times_end])

points(yout2,yout2_lag,col="red")


