
#RJ-MCMC
library(coda); library(MCMCpack); library(lattice)

#################################################################################
#fitting mixture models with an unknown number of components
#In the following example we keep the model as simple as possible
#(maybe even simpler than possible) to highlight the RJ aspects.

#This example is inspired by (but is simpler than)
#"On Bayesian analysis of mixtures with an unknown number of components" 
#Richardson & Green, JRSSB, 59, 731-792 (1997). 

#Load the galaxy velocity data
y=scan("http://www.stats.ox.ac.uk/~nicholls/BayesMethods/galaxies.txt")

#a couple of prior hyperparameters - 
lambda=10  #clusters ~ Poisson(lambda)
alpha=1    #cluster weights Dirichlet(alpha,...,alpha)

DEBUG=FALSE #set to TRUE to check get Poisson>0 number of clusters
log.lkd<-function(y,w,mu,sigma,debug.check=DEBUG) {
  if (debug.check==TRUE) return(0)
  tot=0
  for (i in 1:length(y)) {
    tot=tot+log(sum(w*dnorm(y[i],mean=mu,sd=sigma)))
  }
  return(tot) #debug check get prior
}


#emphasis here on illustrating RJ so _very_ simple prior
#Suppose we have information informing the following prior
log.prior<-function(w,mu,sigma) {                    
  mdp=dpois(nc<-length(w),lambda,log=TRUE)           #80 data pts, lambda=10 => 1-20 clusters
  sp=sum(dgamma(sigma,shape=1.5,rate=0.5,log=TRUE))  #mean 3, rule out very dense clusters, sd about 2.5
  wp=log(ddirichlet(w,rep(alpha,length(w))))         #cluster weight sum to one else uniform
  mup=sum(dnorm(mu,mean=20,sd=10,log=TRUE))          #mu is 0-40 at 2 sigma
  return(sum(mup+sp+wp+mdp+lfactorial(nc)))          #remember to add nc! term as theta was a set not a vector
}

#initialise MCMC with one component
mu=mean(y)
sigma=sd(y)
w=1
oll=log.lkd(y,w,mu,sigma)
olp=log.prior(w,mu,sigma)

K=100000; SS=100; Nsamp=K/SS;  #lecture used 10 x this
PL=matrix(NA,Nsamp,3); colnames(PL)<-c("d","olp","oll"); TH=list();
for (k in 1:K) {
  OK=TRUE
  move=sample(1:5,1) #choose a MCMC move (2 birth death, 3 fixed dimension)
  d=length(w)
  if (move==1) {
    #add a cluster
    i=sample(1:d,1) #pick a component to split - the weight w[i] is shared with the new component
    wn=runif(1,min=0,max=w[i]); wp=c(w,wn); wp[i]=w[i]-wn
    mun=rnorm(1,mean=20,sd=10); mup=c(mu,mun)
    sn=rgamma(1,shape=1.5,rate=0.5); sigmap=c(sigma,sn)
    qd=dnorm(mun,mean=20,sd=10,log=TRUE)+dgamma(sn,shape=1.5,rate=0.5,log=TRUE)-log(w[i])
    qn=0
    rhod=-log(d)
    rhon=-log((d+1)*d)
  }
  if (move==2) {
    #kill a cluster
    if (d>1) {
      ij=sample(1:d,2,replace=FALSE); i=ij[1]; j=ij[2]
      #pick component to kill and a component to increment with the weight from the killed cluster
      wp=w; wp[j]=w[j]+w[i]; wp=wp[-i]
      mup=mu[-i]; sigmap=sigma[-i]
      qd=0
      qn=dnorm(mu[i],mean=20,sd=10,log=TRUE)+dgamma(sigma[i],shape=1.5,rate=0.5,log=TRUE)-log(w[i]+w[j])
      rhod=-log(d*(d-1))
      rhon=-log(d-1)
    } else {
      OK=FALSE
    }
  }
  if (move==3) {
    #fixed dimension mu - simple RW MH
    i=sample(1:d,1); 
    mup=mu; mup[i]=rnorm(1,mean=mu[i],sd=1);
    sigmap=sigma; wp=w;
    qd=qn=rhod=rhon=0      
  }
  if (move==4) {
    #fixed dimension sigma - scaling move
    i=sample(1:d,1); 
    sigmap=sigma; u=runif(1,0.5,2); sigmap[i]=sigma[i]*u
    mup=mu; wp=w;
    qn=-log(u); qd=rhod=rhon=0      
  }
  if (move==5) {
    #fixed dimension w - redistribute the weights between two clusters chosen UAR 
    if (d>1) {
      ij=sample(1:d,2,replace=FALSE); i=ij[1]; j=ij[2]
      wp=w; wp[i]=runif(1,min=0,max=w[i]+w[j]); wp[j]=w[i]+w[j]-wp[i]
      mup=mu; sigmap=sigma;
      qn=qd=rhod=rhon=0
    } else {
      OK=FALSE
    }
  }

  if (OK) {
    nll=log.lkd(y,wp,mup,sigmap)
    nlp=log.prior(wp,mup,sigmap)  
    MHR=nll+nlp-oll-olp-rhod+rhon-qd+qn
    if (log(runif(1))<MHR) {
      w=wp; mu=mup; sigma=sigmap
      oll=nll; olp=nlp
    }
  }
  
  if (k%%SS==0) {
    TH[[k/SS]]=list(w=w,mu=mu,sigma=sigma)
    PL[k/SS,]=c(length(w),olp,oll)
  }
}

#Convergence checks - repeat runs omitted 
#pdf('RJmixtureGmcmc.pdf',6,6)
xyplot(out.mcmc<-as.mcmc(PL)) #reasonbly flat, no obvious trend
#dev.off()
effectiveSize(out.mcmc)       #adequate
par(mfrow=c(1,3),oma=c(1,1,1,1)); #acf plots VG
for (i in 1:dim(out.mcmc)[2]) {
  par(mai=0.2*c(1,1,1,1)); mxl=100
  plot(acf(out.mcmc[,i],lag.max=mxl,plot=F),type='l',ann=F,xaxp=c(0,mxl,2),yaxp=c(0,1,1)); 
  text(mxl/10,0.8,colnames(out.mcmc)[i])
}

#posterior dbn over number of components
#pdf('RJmixtureGhistM.pdf',8,4)
par(mfrow=c(1,1))
hist(d<-PL[,1],breaks=0.5:(0.5+max(PL[,1])),main='',freq=FALSE,cex.lab=1.3,cex.axis=1.3,
     xlab='Number of components m in mixture',ylab='density',xlim=c(0,1+max(PL[,1])))
if (DEBUG) { 
  #if we are debugging check we reproduce the Poisson>0 dbn for the number of clusters
  mn=min(d); mx=max(d); 
  lines(mn:mx,dpois(mn:mx,lambda)/(1-dpois(0,lambda)))
}
#dev.off()
table(d)/Nsamp

#compute posterior predictive distributions - weights (w) vs. means (mu) scater plot
#pdf('RJmixtureGscatterMUW.pdf',8,4)
J=200
x=seq(0,40,length.out=J)
den=matrix(0,max(d)-min(d)+1,J)
plot(c(),c(),xlim=c(0,50),ylim=c(0,1),
     xlab='MU, colored by number of clusters in state',ylab='W component weight',
     cex.lab=1.3,cex.axis=1.3)
for (k in 1:Nsamp) { 
  w=TH[[k]]$w; mu=TH[[k]]$mu; sigma=TH[[k]]$sigma
  ind=d[k]-min(d)+1
  points(mu,w,pch='.',col=ind)
  den[ind,]=den[ind,]+exp(apply(t(x),2,log.lkd,w,mu,sigma))
}
#next line can throw an error if dims in sample are not contigious
cden=den/as.vector(table(d)) #den/c(as.vector(table(d)),0)
mden=apply(den,2,sum)/Nsamp
#dev.off()


#plot posterior predictive distributions
#pdf('RJmixtureGppd.pdf',8,4)
#data
hist(y,breaks=x,freq=FALSE,main='',xlab='velocity',ylab='density')
#posterior predicitve mean
for (i in c(4,6,8)) {
  #posterior predictive dbn conditioned on "i" components
  lines(x,cden[i,],lwd=2,col=1+(i-2)/2)
}
lines(x,mden,col=1,lwd=3)
#dev.off()

#pdf('RJmixtureGdata.pdf',8,4)
hist(y,breaks=x,freq=FALSE,main='',xlab='velocity',ylab='density')
#dev.off()

#shows the label switching if we dont sort mu at plot
#pdf('RJmixtureLS.pdf',8,4)
plot(c(1,Nsamp),c(0,0),type='n',ylim=c(0,50),ylab='mu',xlab='MCMC step (x 100)')
for (i in 1:Nsamp) {
  points(rep(i,PL[i,1]),(TH[[i]]$mu),pch='.',cex=2,col=1:PL[i,1])
}
#dev.off()

#shows the tendency for clusters with larger means to have larger variance
plot(0,0,type='n',xlab="means mu",ylab="st.devs sigma",xlim=c(0,50),ylim=c(0,15))
for (i in 1:Nsamp) {
  o=order(TH[[i]]$mu)
  points(TH[[i]]$mu[o],TH[[i]]$sigma[o],pch=1,cex=4*TH[[i]]$w[o],col=1:PL[i,1],
         cex.lab=1.3,cex.axis=1.3)
}

