#ABC

#some standard mcmc and plotting libraries 
library(coda);library(lattice)

######################################################
#Basic ABC - example estimate Poisson mean
#from a small sample set

#Use a small synthetic data set to illustrate
lambda.true=2; 
n.dat=5
#Y=sort(rpois(n.dat,lambda=lambda.true)) #sorted data is sufficient
Y=c(1,1,2,2,3) #example used had these values
alpha=1; beta=1;
#5 poisson numbers, estimate Poisson mean; Prior Gamma(alpha,beta)

#this is just rejection not ABC
K=100000; LR=c();
for (k in 1:K) {
  lambda=rgamma(1,shape=alpha,rate=beta)
  y.sim=sort(rpois(n.dat,lambda)) #sorting increases our chance of getting a hit
  if (all(y.sim==Y)) {LR=c(LR,lambda)}
}

#pdf('ABCpoisson.pdf')
ho=hist(LR,20,plot=F)
plot(ho$mids,ho$density,ylim=c(0,1.1),type='l',col=2,xlab='LAMBDA',ylab='posterior and approx posterior',lwd=2)
lines(ho$mids,dgamma(ho$mids,shape=alpha+sum(Y),rate=beta+n.dat),lwd=2)
#above shows this sort of rejection works as expected

#ABC - choose mean(Y) as statistic - propose lambda.sim using prior, 
#simulate y.sim~Poisson(lambda.sim) accept if mean(y.sim) is
#within delta of mean(Y); first delta is 1
K=10000; LA=SA=c(); my=mean(Y); delta=1
for (k in 1:K) {
  lambda=rgamma(1,shape=alpha,rate=beta)
  y.sim=sort(rpois(n.dat,lambda)) #sorting increases our chance of getting a hit
  if (abs(mean(y.sim)-my)<delta) {LA=c(LA,lambda); SA=c(SA,mean(y.sim))}
}
ho=hist(LA,50,plot=F)
lines(ho$mids,ho$density,type='l',col=3,lwd=2)

#Regression correction
b=coef(lm(LA~SA))[2]
LA.adj=LA+b*(mean(Y)-SA)
ho=hist(LA.adj,50,plot=F)
lines(ho$mids,ho$density,type='l',col=3,lwd=2,lty=2)

#repeat with smaller tolerance delta
K=100000; LB=SB=c(); my=mean(Y); delta=0.5
for (k in 1:K) {
  lambda=rgamma(1,shape=alpha,rate=beta)
  y.sim=sort(rpois(n.dat,lambda)) #sorting increases our chance of getting a hit
  if (abs(mean(y.sim)-my)<delta) {LB=c(LB,lambda); SB=c(SB,mean(y.sim))}
}
ho=hist(LB,50,plot=F)
lines(ho$mids,ho$density,type='l',col=4,lwd=2)

#Regression correction
b=coef(lm(LB~SB))[2]
LB.adj=LB+b*(mean(Y)-SB)
ho=hist(LB.adj,50,plot=F)
lines(ho$mids,ho$density,type='l',col=4,lwd=2,lty=2)

legend('topright',
       c('True posterior','Rejection','delta=1','delta 1, adjusted','delta=0.5','delta=0.5, adjusted'),
       col=c(1,2,3,3,4,4),lwd=c(2,2,2,2,2,2),lty=c(1,1,1,2,1,2))
#dev.off()
######################################################################
#Ising model example

n.disagree<-function(X) {
  #count number disagreeing nbrs
  sum(abs(diff(X))+abs(diff(t(X))))
}

mcmc.ising<-function(theta,T,n) {
  #MCMC simulation of Y~Ising(theta); nxn lattice, T-steps
  Y<-matrix(runif(n^2)>1/2,n,n) #random start state
  hashY<-n.disagree(Y) 
  for (j in 1:T) {
    i<-1+floor(runif(1)*n^2)
    Yp<-Y
    Yp[i]<-1-Y[i]
    hashYp<-n.disagree(Yp)
    if (runif(1)<exp(theta*(hashY-hashYp))) {
      Y<-Yp
      hashY<-hashYp
    }
  }
  return(Y)
}
#example call
Yp8=mcmc.ising(theta=0.8,T=200000,n=32) #T chosen by checking convergence

#ran above with different thetas and saved
Yp4=mcmc.ising(theta=0.4,T=200000,n=32)
Yp1p2=mcmc.ising(theta=1.2,T=200000,n=32)

#pdf('IsingExamples.pdf',10,4)
par(mfrow=c(1,3))
image(Yp4,col=gray(0:255/255),axes=F,main='Ising model, theta=0.4'); box()
image(Yp8,col=gray(0:255/255),axes=F,main='Ising model, theta=0.8'); box()
image(Yp1p2,col=gray(0:255/255),axes=F,main='Ising model, theta=1.2'); box()
#dev.off()

#Ising ABC
#Use synthetic data for known theta to illustrate
theta.true=0.8
Y=mcmc.ising(theta=theta.true,T=100000,n=8)
image(Y,col=gray(0:255/255),axes=F); box()

#ABC part
n=8; max.dis=(4*2+4*(n-2)*3+(n-2)^2*4)/2 #the maximum number of disagreeing nbrs
s=n.disagree(Y)/max.dis #s=S(Y) scaled so 0<s<1 convenient for setting delta
theta=S=c(); K=1000

delta=0.1; #run this twice for delta=0.1/0.05
for (k in 1:K) {
  th=rexp(1,2)                #Exp(2) prior
  Yp=mcmc.ising(th,T=1000,n=8)#simulate Y'
  sp=n.disagree(Yp)/max.dis   #Statistic S(y)
  if (abs(sp-s)<delta) {      #if simulated data is close to Y accept
    theta=c(theta,th)
    S=c(S,sp)
  }
}
thetap1=theta; Sp1=S     #save results for different delta values
#thetap05=theta; Sp08=S

#pdf('IsingABC.pdf',6,6)
plot(
  density(thetap05),
  main='',xlab='THETA',ylab='Approx Posterior density',xlim=c(0,2),ylim=c(0,3))
lines(density(thetap1),col=2)

#regression adjustment
adjust<-function(theta,S) {
  b=coef(lm(theta~S))[2]
  theta.adj=theta+(s-S)*b
  return(theta.adj)
}
lines(density(adjust(thetap1,Sp1)),col=2,lty=2)
lines(density(adjust(thetap05,Sp05)),lty=2)
abline(v=theta.true,lwd=2,col=3)
legend('topright',c('delta=0.05','delta=0.1','0.05 adjusted','0.1 adjusted','true'),
col=c(1,2,1,2,3),lty=c(1,1,2,2,1),lwd=c(1,1,1,1,2))
#dev.off()

#Remark - there is a library(abc) package with many useful tools
library(abc)
help(package=abc)

