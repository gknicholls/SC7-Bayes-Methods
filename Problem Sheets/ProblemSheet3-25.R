#PS3 Bayes Methods MT 2025
library(coda); library(MCMCpack); library(lattice)

#Code illustrating some of the methods from PS3 - this is all optional but may help with understanding

#Q3 ABC-MCMC for Poisson
#Use a small synthetic data set to illustrate
Y=c(1,1,2,2,3) #ABC example in lecture had these values
n.dat=length(Y)
alpha=1; beta=1; #prior hyperparameters
#5 poisson numbers, estimate Poisson mean; Prior Gamma(alpha,beta)

#ABC settings
my=mean(Y); delta=0.5

#ABC-MCMC
K=10000; SS=10; 
Lambda=rep(NA,K/SS); lambda=rgamma(1,shape=alpha,rate=beta); w=1
for (k in 1:K) {
  lp=lambda+w*rnorm(1);
  if (lp>0) {
    y.sim=rpois(n.dat,lp)
    if (abs(mean(y.sim)-my)<delta) {
      MHR=dgamma(lp,alpha,beta,log=T)-dgamma(lambda,alpha,beta,log=T)
      if (log(runif(1))<MHR) {
        lambda=lp
      }
    }
  }
  
  if (k%%SS==0) {
    Lambda[k/SS]=lambda;
  }
}
ho=hist(Lambda,sqrt(K/SS),plot=F)
plot(ho$mids,ho$density,ylim=c(0,1.1),type='s',col=2,xlab='LAMBDA',ylab='posterior and approx posterior',lwd=2)
lines(ho$mids,dgamma(ho$mids,shape=alpha+sum(Y),rate=beta+n.dat),lwd=2)
lines(ho$mids,dgamma(ho$mids,shape=alpha,rate=beta),lwd=2,lty=2)
legend('topright',c('prior','posterior','MCMC-ABC'),col=c(1,1,2),lwd=2,lty=c(2,1,1))

#
#Q4 and 8 on model averaging with outliers

#load data
library(MASS)
data(hills)
a=hills
str(a)

#transform data
a$y=sqrt(a$time)
a$climb=scale(a$climb)
a$dist=scale(a$dist)

X=model.matrix(y~climb+dist,data=a)
n=dim(a)[1]

llk<-function(z,b,s) {
  #z outlier index, b regression beta, s regression sigma
  #fixed scale for outliers 
  rho=3 
  sum(log(z*dnorm(a$y,mean=X%*%b,sd=rho*s)+(1-z)*dnorm(a$y,mean=X%*%b,sd=s)))
}

lprior<-function(p,b,z,s) {
  #p is probability a point is an outlier
  tot=dbeta(p,shape1=1,shape2=9,log=TRUE)
  tot=tot+sum(z*log(p)+(1-z)*log(1-p))
  tot=tot+sum(dt(b/2.5,df=1,log=TRUE))
  tot=tot-log(s)
  return(tot)
}

#inits
p=0.1
b=coef(lm(y~climb+dist,data=a))
z=rep(0,n)
s=1
olp=lprior(p,b,z,s); oll=llk(z,b,s)

#Run length, subsampling - these choices give a run time around 3 minutes on my laptop
T=100000; SS=100; 
Theta=matrix(NA,T/SS,length(c(p,b,s))); Z=matrix(NA,T/SS,n);
#Jump sizes
u=0.05; w=0.2; v=0.1;

for (i in 1:T) {
  #update p
  pp=p+u*rnorm(1);
  if (pp>0 & pp<1) {
    nlp=lprior(pp,b,z,s);
    if (log(runif(1))<nlp-olp) {
      p=pp; olp=nlp;
    }
  }
  #update b (ie beta)
  bp=b+w*rnorm(3);
  nll=llk(z,bp,s); nlp=lprior(p,bp,z,s);
  if (log(runif(1))<nlp+nll-olp-oll) {
    b=bp; olp=nlp; oll=nll
  }
  
  #update s (ie sigma)   
  sp=s+v*rnorm(1);
  if (sp>0) {
    nll=llk(z,b,sp); nlp=lprior(p,b,z,sp);
    if (log(runif(1))<nlp+nll-olp-oll) {
      s=sp; olp=nlp; oll=nll
    }
  }
  
  #update z
  for (j in 1:n) {
    zp=z; zp[j]=1-z[j];
    nll=llk(zp,b,s); nlp=lprior(p,b,zp,s);
    if (log(runif(1))<nlp+nll-olp-oll) {
      z=zp; olp=nlp; oll=nll
    }
  }
  
  #save the sample state
  if (i%%SS==0) {
    Theta[i/SS,]=c(p,b,s)
    Z[i/SS,]=z
  }
}
plot(as.mcmc(Theta))
effectiveSize(Theta)
image(Z)
#####
# Use the output Theta and Z to identify outliers etc

op=apply(Z,2,mean);
#pdf('PS2outlier.pdf',5,5)
plot(1:n,op,pch=16,ylim=c(0,1.5),xlab='data index',ylab='outlier probability'); 
text(1:n,op,row.names(a),srt=90,pos=4,cex=0.8)
abline(h=1,lty=2)
#dev.off()
summary(as.mcmc(Theta))
HPDinterval(as.mcmc(Theta))

########################################################################
#Q5 and 9 illustrating MCMC with a Jacobian

lf<-function(x) {
  -log(x)/2-log(1+x^2)
}

m<-function(T=10000,nu=8,x0=2) {
 
  X=rep(NA,T);
  X[1]=x=x0;
  for (k in 2:T) {
    u=rt(1,nu); xp=x^u
    MHR=lf(xp)-lf(x)+dt(1/u,nu,log=TRUE)-dt(u,nu,log=TRUE)+log(abs(x^(u-1)/u))
    if (is.na(MHR)) MHR=-Inf
    if (log(runif(1))<MHR) {x=xp}
    X[k]=x
  }
  return(X)
}
X<-m(T<-1e4)

plot(log(X))
effectiveSize(as.mcmc(X))

#check density correct in [L,U] any 0<L<U<infinity
#care needed as singular at 0 and heavy tail

dx=0.01; L=0.5; U=4
z=seq(L,U,dx)
Z=sum(exp(lf(z)))*dx; 
Z 		#this is just the approx norm const for x|x in [L,U] - for ease of estimating Z and plotting
hist(X[L<X & X<U],round(sqrt(T)),freq=FALSE); 
lines(z,exp(lf(z))/Z,col=2)

#choice of nu?
nu=c(1,2,4,8,16,32); es=0*nu
for (k in 1:length(nu)) {
  X<-m(T=1e6,nu=nu[k]); es[k]=effectiveSize(as.mcmc(X));  
}
plot(nu,es) 
#alot of variation but nu=1 tends to be lower than nu ~ 8

#remark - actually the CDF and Z for this distribution can be computed exactly. For example
#$\int_0^\infty x^(1/a)/(1+x^(1+b/a))dx=a\pi/((a+b)\sin(\pi(b-1)/(a+b)))$ for a>1 and b>1


########################################################################
#Q10 illustrating RJ-MCMC targeting the Inhomogeneous Poisson Process

dmix<-function(x,mu,sig,w) {
  #calculate lambda(x) at locations in x
  nx=length(x)
  dout=apply(sapply(1:nx,function(i) w*dnorm(x[i],mu,sig)),2,sum)
  return(dout)
}

pmix<-function(x,mu,sig,w) {
  #calculate cumulative integral from -inf to each location in x
  nx=length(x)
  pout=apply(sapply(1:nx,function(i) w*pnorm(x[i],mu,sig)),2,sum)
  return(pout)
}

rIPP<-function(T,SS,mu,sig,w,L,U,y0=c()) {
  #MCMC targeting y~PPP(lambda,[L,U])
  Y=vector('list',T/SS) #T MCMC steps save every SS steps
  y=y0                  #initialise - y will be the point set
  N=length(y)           
  for (t in 1:T) {
    if (runif(1)<0.5) { #wp 1/2 propose to add a point
      np=runif(1,min=L,max=U) #propose a new pt loc
      #p(y')/p(y)=lambda(np) as everything else cancels
      #q(y'|y)=0.5 x 1/(U-L)
      #q(y|y')=0.5 x 1/(N+1) if currently N pts 
      MHR=log(dmix(np,mu,sig,w))+log(U-L)-log(N+1)
      if (log(runif(1))<MHR) {
        N=N+1
        y=c(y,np) #tack the point on the end
      }
    } else {            #wp 1/2 propose to delete a point
      if (N>0) {        #can only delete if N>0 - if N=0 reject proposal
        i=sample(1:N,1) #choose a oint to delete
        #p(y')/p(y)=1/lambda(y[i]) 
        #q(y'|y)=0.5 x 1/N
        #q(y|y')=0.5 x 1/(U-L)  
        MHR=-log(dmix(y[i],mu,sig,w))-log(U-L)+log(N)
        if (log(runif(1))<MHR) {
          N=N-1
          y=y[-i]
        }
      }
    }
    if (!(t%%SS)) {Y[[t/SS]]=y} #save every SS steps
  }
  return(Y)
}

#initialise
set.seed(21)

#interval over which the IPPP runs
L=0; U=100

#define the true rate lambda(t) as a mixture
M=5;                        #number of components
muT=sort(L+(U-L)*runif(M))  #locations of means
sigT=rgamma(M,(U-L)/10,1)   #sd.dev at each mean about 10% of support
Wav=6
wT=Wav*rgamma(M,1,1)        #component weights - gives E(N) around MxW

#plot lambda(x)
x=seq(L,U,length.out=1000)
lamT=dmix(x,muT,sigT,wT)    #evaluate lambda(x) at points in x
clamT=pmix(x,muT,sigT,wT)-pmix(L,muT,sigT,wT)   #evaluate cumulative integral of lambda(x)
interval=c(L,U)             #take the whole range
(mn=diff(pmix(interval,muT,sigT,wT))) #int(lam,-inf,U)-int(lam,-inf,L)

#ignore the black lines they are checking my pmix function is right
plot(x,clamT,type='l',ylim=c(0,1.1*max(clamT))); abline(h=mn,lty=2)
lam.sc=mn/max(lamT)
lines(x,lamT*lam.sc,col=2) #this is lambda(t)

#run MCMC to sample y~PPP(lambda,[L,U])
#T steps, subsample every SS steps
T=100000; SS=10; NS=T/SS
Y<-rIPP(T,SS,muT,sigT,wT,L,U,y0=c()) #returns samples as a list

rp<-Y[[NS]];  #last sampled state
points(rp,0*rp,pch=16,col=3,cex=1) #more points where lambda is larger

N=sapply(Y,length) #number of points in each of the NS samples
plot(N) #this and next quick convergence check
effectiveSize(N) #ideally > 1000

#Check we got that all correct
#number of events in any interval [a,b] should be Poisson(int(lam,a,b))
interval=c(25,75) #c(L,U) #can take any interval we like
(mn=diff(pmix(interval,muT,sigT,wT))) #expected number from PPP
Nint=sapply(Y,function(y) sum(y>interval[1] & y<interval[2])) #MCMC sampled number of events in interval for each sample
mean(Nint) #MCMC mean should match mn above
#more demandingly the dbn of Nint should be Poisson(mn)
brk=-0.5:max(Nint+0.5); hs=hist(Nint,breaks=brk,freq=FALSE); 
lines(hs$mids,dpois(hs$mids,mn),col=2)


