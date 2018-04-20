model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
  
  R = 1/(1+A^n)               # Repressor
  dA = M*L - delta*A - v*M*A  # Allolactose
  dM = c0 + c*(1-R) - d*M     # mRNA
  
  return(list(c(dA, dM)))  
  }) 
}

p <- c(L=1,c=1,c0=0.05,d=1,delta=0.2,n=5,v=0.25)
s <- c(A=0,M=0)
plane(xmax=4,ymax=1.2)
low <- newton(s,plot=T)
mid <- newton(c(A=0.8,M=0.2),plot=T)
hig <- newton(c(A=2,M=1),plot=T)
continue(mid,x="L",y="A",xmax=2,ymax=4)

