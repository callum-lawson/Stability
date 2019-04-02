# ratio of distances travelled (exponential of integral towards zero)

decay <- function(t,a,b){
  a * (1 - exp(-b*t))
}

curve(decay(t=x,a=1,b=1),xlim=c(0,5))
curve(decay(t=x,a=0.5,b=2),col="red",add=T)

# => similar under blue noise, but shallower DD fluctuates more under red noise
# no simple function for difference / ratio between these

curve(decay(t=x,a=1,b=1)-decay(t=x,a=0.5,b=2),xlim=c(0,5))

# Horizontal perturbations:
curve(decay(t=x,a=1,b=1),xlim=c(0,5))
curve(decay(t=x,a=0.5*2,b=2),col="red",add=T)
  # order reversed: more stable species fluctuates more