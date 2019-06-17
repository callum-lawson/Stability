### Comparing self-regulation around k between biotic and abiotic models ###

abiotic <- function(x,v,k){
  v * k * exp(-x) - v
}

biotic <-  function(x,v,k){
  v - v/k * exp(x)
}

abcurve <- function(v,k,scale){
  curve(abiotic(x,v,k),xlim=c(log(k)-scale,log(k)+scale))
  curve(biotic(x,v,k),col="red",add=T)
  abline(h=0,col="gray",lty=3)
}

par(mfcol=c(2,2),mar=c(3,3,3,3))
abcurve(v=1,k=2,scale=1)
abcurve(v=1,k=2,scale=0.001)
abcurve(v=1,k=0.5,scale=1)
abcurve(v=1,k=0.5,scale=0.001)
