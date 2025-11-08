#MCMC with a transformation
library(coda); library(MCMCpack); library(lattice)

#MCMC generalisations leading up to RJ MCMC

###
#Detailed balance with a Jacobian term
x=1; 
N=100000; X=rep(NA,N)
for (k in 1:N) {
  u=runif(1,min=0.5,max=2)
  xp=u*x
  if (log(runif(1))<x-xp-log(u)) {
    x=xp
  }
  X[k]=x
}
hist(X,freq=F,100); lines(v<-seq(0,5,0.2),dexp(v,1),lwd=2,col=2)
effectiveSize(X)

###
#Detailed balance with a Jacobian and matched pairs of moves
x=1; 
N=100000; X=rep(NA,N)
for (k in 1:N) {
  move=sample(1:2,1);
  if (move==1) {
    u=runif(1,min=0.5,max=1)
    xp=u*x
    if (log(runif(1))<x-xp-log(u)-log(2)) {
      x=xp
    }
  } else {
    u=runif(1,min=1,max=2)
    xp=u*x
    if (log(runif(1))<x-xp-log(u)+log(2)) {
      x=xp
    }
  }
  X[k]=x
}
hist(X,freq=F,100); lines(v<-seq(0,5,0.2),dexp(v,1),lwd=2,col=2)
effectiveSize(X) #slightly more efficient?

