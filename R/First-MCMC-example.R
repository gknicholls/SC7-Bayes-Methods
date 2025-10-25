
#SC7 Bayes Methods First examples for basic MCMC up to Gibbs and Probit

######################################################
# First MCMC example
# MCMC simulating X_t according to p=HyperGeom(K=10,N=20,n=10).
K<-10; N<-20; n<-10; Bm<-max(0,n+K-N); Bp<-min(n,K);
T<-100000; X<-rep(NA,T); X[1]<-Bm #start at lower limit X[1]=0
for (t in 1:(T-1)) {
   i<-X[t]
   j<-sample(c(i-1,i+1),1)
   if (j<Bm | j>Bp) {
     X[t+1]<-i         #must reject if outside SP
   } else {
     a<-min(1,(choose(K,j)*choose(N-K,n-j))/(choose(K,i)*choose(N-K,n-i)))
     U<-runif(1)
     if (U<=a) {
        X[t+1]<-j
     } else {
        X[t+1]<-i
     }
  }
}

#pdf(file="mcmcHG.pdf",3.5,3.5)
plot(X[1:200],type='l',xlab='MCMC step',ylab='X_t')
#dev.off()

(phat=table(X)/T)

#pdf(file="mcmcHistHG.pdf",3.5,3.5)
plot(Bm:Bp,round(dhyper(Bm:Bp,m=K,n=N-K,k=n),6),col=1,lwd=2,lty=2,type='l',
     xlab="Omega={0,1,...,10}",ylab="HyperGeom(i;K=10,N=20,n=10)")
lines(as.numeric(names(phat)),phat,type='l',col=2,lwd=2)
legend('topright',c('True PMF','MCMC estimate'),col=c(1,2),lwd=c(2,2),lty=c(2,1))
#dev.off()

######################################################
#MCMC simulating X_t according to a mixture of normals
f<-function(x,mu1,mu2,S1i,S2i,p1=0.5) {
  #mixture of normals, density up to constant factor
  c1<-exp(-t(x-mu1)%*%S1i%*%(x-mu1))
  c2<-exp(-t(x-mu2)%*%S2i%*%(x-mu2))
  return(p1*c1+(1-p1)*c2)
}

mcmc.binorm<-function(n,a,x0,mu1,mu2,S1i,S2i) {
  X=matrix(NA,n,2); X[1,]=x=x0;
  for (t in 1:(n-1)) {
    y<-x+(2*runif(2)-1)*a
    MHR<-f(y,mu1,mu2,S1i,S2i)/f(x,mu1,mu2,S1i,S2i)
    if (runif(1)<MHR) x<-y
    X[t+1,]<-x
  }
  return(X)
}

mu1=c(1,1); mu2=c(4,4); S=diag(2); S1i=S2i=solve(S);
X<-mcmc.binorm(a=4,n=10000,x0=mu1,mu1,mu2,S1i,S2i)

plot(X[1:2000,1],pch=16,type='l',main='',xlab='MCMC step t=1,2,3,...',ylab='MCMC state X[1]_t')
nshow=500; plot(X[1:nshow,],pch=16,ann=F)
lines(X[1:nshow,1],X[1:nshow,2],col=3)
plot(X,pch=16,ann=F)
m=30; u=seq(-2,7,length.out=m); v=u; z=matrix(NA,m,m); 
for (i in 1:m) { 
  for (j in 1:m) { 
    z[i,j]=f(c(u[i],v[j]),mu1,mu2,S1i,S2i) 
  } 
}
contour(u,v,z,col=2,nlevels=7,add=T,lwd=2)

## One coordinate at a time
mcmc.oct.binorm<-function(n,a,x0,mu1,mu2,S1i,S2i,xi=0.5) {
  X=matrix(NA,n,2); X[1,]=x=x0;
  for (t in 1:(n-1)) {
    k<-sample(x=1:2,size=1,prob=c(xi,1-xi)) #pick a component k to update
    y=x; y[k]<-x[k]+(2*runif(1)-1)*a #overwrite entry k
    MHR<-f(y,mu1,mu2,S1i,S2i)/f(x,mu1,mu2,S1i,S2i) #otherwise the same
    if (runif(1)<MHR) x<-y
    X[t+1,]<-x
  }
  return(X)
}
mu1=c(1,1); mu2=c(4,4); S=diag(2); S1i=S2i=solve(S);
set.seed(4)
X<-mcmc.oct.binorm(a=4,n=10000,x0=mu1,mu1,mu2,S1i,S2i)

plot(X[1:2000,1],pch=16,type='l',main='',xlab='MCMC step t=1,2,3,...',ylab='MCMC state X[1]_t')
nshow=1000; plot(X[1:nshow,],pch=16,ann=F,xlim=c(-1,6),ylim=c(-1,6))
lines(X[1:nshow,1],X[1:nshow,2],col=3,lwd=1.5)
plot(X,pch=16,ann=F)
m=30; u=seq(-2,7,length.out=m); v=u; z=matrix(NA,m,m); 
for (i in 1:m) { 
  for (j in 1:m) { 
    z[i,j]=f(c(u[i],v[j]),mu1,mu2,S1i,S2i) 
  } 
}
contour(u,v,z,col=2,nlevels=7,add=T,lwd=2)

######################################################
#Convergence and mixing

#small jumps a=2 inefficient
X<-mcmc.binorm(n=3000,a=2,x0=c(-6,0),mu1,mu2,S1i,S2i)
#pdf(file="mcmcMX2.pdf",3.5,3.5)
plot(X[,1],pch=16,type='l',main='a=2',xlab='MCMC step t=1,2,3,...',ylab='MCMC state X[1]_t')
T<-100
abline(v=T,lwd=2,lty=2)
#dev.off()

#larger jumps a=4
X<-mcmc.binorm(n=3000,a=4,x0=c(-6,0),mu1,mu2,S1i,S2i)
#pdf(file="mcmcMX4.pdf",3.5,3.5)
plot(X[,1],pch=16,type='l',main='a=4',xlab='MCMC step t=1,2,3,...',ylab='MCMC state X[1]_t')
T<-100
abline(v=T,lwd=2,lty=2)
#dev.off()

#many useful MCMC analysis tools are in the 'coda' package
library(coda)

#the autocorrelation function and ESS
X<-mcmc.binorm(n=(n<-10000),a=2,x0=c(-6,0),mu1,mu2,S1i,S2i)
#pdf(file="acfMX2.pdf",3.5,3.5)
acf(X[T:n,1],main='',lag.max=200)
effectiveSize(X)
#dev.off()

X<-mcmc.binorm(n=10000,a=4,x0=c(-6,0),mu1,mu2,S1i,S2i)
#pdf(file="acfMX4.pdf",3.5,3.5)
acf(X[T:n,1],main='',lag.max=200)
effectiveSize(X)
#dev.off()

#summary has a method for MCMC output - the std(fbar) is 'time series SE'
summary(as.mcmc(X[T:n,]))

#We can roughly work out "tau" (the IACT) from a glance at the ACF
#Then n/tau is the ESS. Once we have tau_f we can compute var(fbar)
#using var(fbar)=var(f)*tau_f/n. The summary command gives us std(fbar)
#directly

#the 'brute force' way to estimate var(fbar) repeats multiple runs 
K<-100; m<-matrix(0,K,2); for (i in 1:K) m[i,]<-apply(mcmc.binorm(n=1000,a=4,x0=mu1,mu1,mu2,S1i,S2i),2,mean)
apply(m,2,sd)

#compare
X<-mcmc.binorm(n=1000,a=4,x0=mu1,mu1,mu2,S1i,S2i)
summary(as.mcmc(X))

#One straightforward approach to convergence checking is to make multiple independent runs
#and check the the distribution of each statistic of interest is stable across runs. 
#Here I look at the distribution of X_1

X<-list(); K=3; xinit=seq(from=-5,to=10,length.out=K)
for (k in 1:K) {
  X[[k]]<-mcmc.binorm(n=3000,a=2,x0=xinit[k],mu1,mu2,S1i,S2i)
}
#pdf(file="mcmcMR.pdf",height=3.5,width=8)
par(mfrow=c(1,2));
plot(X[[1]][,1],type='l',ylim=c(-5,10),xlab='MCMC step t=1,2,3,...',ylab='MCMC state X[1]_t'); for (k in 2:K) {lines(X[[k]][,1],col=k)}
plot(density(X[[1]][,1]),type='l',main='',xlab='X_1',xlim=c(-2,7),ylim=c(0,0.5)); for (k in 2:K) {lines(density(X[[k]][,1]),col=k)}
#dev.off()


######################################################
# MCMC for Probit using DA as in lect 3
# code based on http://www4.stat.ncsu.edu/~reich/ST740/code/probit.R

#######  Code to fit the model    ########:
#y[i] =  Bern[Phi(eta[i])] where eta[i]=x[i,]%*%theta
#x[i,] = covariates for sub i
#theta-prior is theta~N(0,I_p*sd.theta^2) for theta[i], i=1:p

#see at end for example

#Draw samples from a truncated normal
#useful for z_i~p(z_i|theta,y_i)

rtnorm<-function(n,mu,sigma,lower,upper){ 
  lp<-pnorm(lower,mu,sigma) 
  up<-pnorm(upper,mu,sigma)  
  qnorm(runif(n,lp,up),mu,sigma) 
}

#main probit DA MCMC

probit<-function(y,x,sd.theta=100,iters=10000,z.save=10){
  
  n<-length(y)
  p<-ncol(x)
  # low < z|y < high - if y=0 z<0, if y=1, z>0
  low<-ifelse(y==1,0,-Inf)
  high<-ifelse(y==1,Inf,0)
  
  #Initial values
  z<-y-.5          #satisfy low < z|y < high
  theta<-rep(0,p)  #start at prior mean
  
  #store samples here
  keep.theta<-matrix(0,iters,p)
  keep.z<-matrix(0,iters,z.save)
  
  #Do some matrix manipulations offline
  txx<-t(x)%*%x
  V<-solve(txx+diag(p)/sd.theta^2) #V from L3 theta|z~N(mu,V)
  P1<-V%*%t(x)                     #P1=VX^T so mu=P1%*%z from L3
  P2<-t(chol(V))                   #if Z~N(0,I_p) then P2%*%Z~N(0,V)
  
  #Let's go!
  for(i in 1:iters){
    
    #update the latent probit variables, z:
    eta<-x%*%theta
    z<-rtnorm(n,eta,1,low,high)
    
    #update theta:
    theta<-P1%*%z+P2%*%rnorm(p) #p(theta|z)=N(theta; VX^Tz,V)
    
    keep.theta[i,]<-theta
    keep.z[i,]<-z[1:z.save]
    
  }
  
  
  list(theta=keep.theta,z=keep.z)}


#Simulated binary data
#y~b_1+b_2*x_2+...+b_10*x_10
#with effect sizes b_1=0, b_2=0.11, b_3=0.22, ..., b_10=1
#and std normal covariates x_2,...,x_10

set.seed(0820)

n<-500
p<-10
true.theta<-seq(0,1,length=p)
x<-matrix(rnorm(n*p),n,p)
x[,1]<-1
xi<-pnorm(x%*%true.theta)
y<-rbinom(n,1,xi)


#Fit the model
fit<-probit(y,x)

#Plot results
boxplot(fit$theta,main="Posterior of theta")
lines(true.theta)
boxplot(fit$z,main="Posterior of z for 10 observations")

######################################################

