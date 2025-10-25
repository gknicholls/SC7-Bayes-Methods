
#Radiocarbon dating example for Lecture 2

rm(list=ls())
#some standard mcmc and plotting libraries 
library(coda);library(lattice)

#load the data
rcc.dat<-read.table("http://www.stats.ox.ac.uk/~nicholls/BayesMethods/SHRCC2013.txt",sep=",",header=TRUE)
str(rcc.dat)
attach(rcc.dat)

Earliest=1000; Latest=1;
y.BP=Latest:Earliest
mu=spline(CAL.BP, X14C.age, xmin = 0.9*Latest, xmax = 1.1*Earliest, xout=y.BP)$y
err=spline(CAL.BP,Error, xmin = 0.9*Latest, xmax = 1.1*Earliest, xout=y.BP)$y

#Data - come from a Moari rivermouth settlement in NZ
#NZ 7758 (1, 1) 580 47
#NZ 7761 (2, 1) 600 50
#NZ 7757 (3, 1) 537 44
#NZ 7756 (4, 1) 670 47
#NZ 7755 (5, 1) 646 47
#WK 2589 (5, 2) 630 35
#NZ 7771 (6, 1) 660 46

#radiocarbon dates
y=c(580,600,537,670,646,630,660); nd=length(y)
#measurement errors
d=c(47,50,44,47,47,35,46)
#layers
s=c(1,2,3,4,5,5,6)
#the laboratory date identifiers
nm=c("NZ 7758","NZ 7761","NZ 7757","NZ 7756","NZ 7755","WK 2589","NZ 7771")

#these are very conservative bounds on the settlement dates
#current best estimates suggest NZ was only settled around 700 BP 
#so a lower bound of 1000 is quite conservative. The upper bound
#comes from information independent of the dates
L=500; U=1000; 


#the log likelihood for calendar dates theta
llk<-function(theta) {
  sig=sqrt(d^2+err[theta]^2)
  llk=sum(dnorm(y,mean=mu[theta],sd=sig,log=T))
  return(llk)
}

#the log prior for the shrinkage prior
lpr<-function(psi) {
  S=psi[2]-psi[1]
  lpr=-nd*log(S)-log(U-L-S)
  return(lpr)
}

#illustrate the calibration function and likelihood
#pdf('rcdcallkd.pdf',8,8)
x=L:U;
plot(x,mu[x],type='l',xlab='Calendar Years Before the Present (Theta)',
     ylab='Uncalibrated Radiocarbon Age (Y)'); 
lines(x,mu[x]+err[x],lty=2); lines(x,mu[x]-err[x],lty=2)
i=6; 
abline(h=y[i],col=3); abline(h=y[i]-d[i],lty=2,col=3); abline(h=y[i]+d[i],lty=2,col=3)
z=dnorm(y[i],mean=mu[x],sd=sqrt(d[i]^2+err[x]^2)); z=(z-min(z))/abs(max(z)); 
z=par()$usr[3]+100*z
lines(x,z,col=2,lwd=2)
legend('topleft',c('calibration curve mu','mu +/- err',
       'uncalibrated carbon date y','y +/- sigma','(scaled) likelihood for theta'),
       lty=c(1,2,1,2,1),col=c(1,1,3,3,2))
#dev.off()

#checking the min max and span for the uniform prior
K=10000
theta.sim=matrix(runif(K*nd,min=L,max=U),K,nd)
lower=apply(theta.sim,1,min)
upper=apply(theta.sim,1,max)
span=upper-lower
#pdf('unifpriorstats.pdf',12,4)
par(mfrow=c(1,3))
hist(lower,30,freq=FALSE); hist(upper,30,freq=FALSE); bins=hist(span,30,freq=FALSE); 
span.uprior<-function(x) {nd*(nd-1)*x^(nd-2)*(U-L-x)/(U-L)^nd}
lines(bins$mids,span.uprior(bins$mids))
#dev.off()
#far from non-informative

#checking the min max and span for the shrinkage prior
span=runif(K,min=0,max=U-L); 
for (k in 1:K) {lower[k]=runif(1,min=L,max=U-span[k]); upper[k]=lower[k]+span[k];}
#pdf('shrinkpriorstats.pdf',12,4)
par(mfrow=c(1,3))
hist(lower,30,freq=FALSE); hist(upper,30,freq=FALSE); hist(span,30,freq=FALSE); 
#dev.off()

#MCMC sampler very simple
rcd.mcmc<-function(nd,L,U,T,Subs,shrink=TRUE,w,scale.move=TRUE) {
  #initialise the start state for the MCMC
  theta=sort(runif(nd,min=L,max=U))
  psi=c((L+theta[1])/2,(theta[nd]+U)/2)
  llo=llk(theta)
  lpo=lpr(psi)
  OP=matrix(NA,T/SubS,1+2+nd+2)
  colnames(OP)<-c('span','psi1','psi2','thet1','thet2','thet3','thet4','thet5','thet6','thet7','llkd','lpri')

  #run the MCMC T steps
  for (t in 1:T) {
    #go through each of the parameters in turn
    for (i in 1:nd) {
      thetap=theta; 
      #the proposal is uniform in allowed range
      if (shrink) {
        thetap[i]=runif(1,min=psi[1],max=psi[2]);
      } else {
        thetap[i]=runif(1,min=L,max=U)
      }
      lln=llk(thetap); #no change to prior
      logMHR=lln-llo
      if (log(runif(1))<logMHR) {
        theta=thetap
        llo=lln; 
      }
    }

    #if we are using the shrinkage prior (psi etc) update the extra parameters
    if (shrink) {
      psip=psi; psip[1]=runif(1,min=L,max=min(theta));
      lpn=lpr(psip)
      logMHR=lpn-lpo
      if (log(runif(1))<logMHR) {
        psi=psip
        lpo=lpn
      }
      psip=psi; psip[2]=runif(1,min=max(theta),max=U);
      lpn=lpr(psip)
      logMHR=lpn-lpo
      if (log(runif(1))<logMHR) {
        psi=psip
        lpo=lpn
      }
      
      if (scale.move) {
        
        delta=runif(1,w,1/w)
        m=mean(c(psi,theta))
        psip=m+(psi-m)*delta
        thetap=m+(theta-m)*delta
        if ( (psip[2]<U) & (psip[1]>L) ) {
          lpn=lpr(psip)
          lln=llk(thetap)
          lqq=(nd-1)*log(delta)
          logMHR=lln-llo+lpn-lpo+lqq
          if (log(runif(1))<logMHR) {
            psi=psip
            theta=thetap
            lpo=lpn
            llo=lln
          }
        }
      }
    }
    #subsample - take a sample every SubS steps - store it in OutPut
    if (!t%%SubS) {
      span=diff(psi)
      OP[t/SubS,]=c(span,psi,theta,llo,lpo)
    }
  }
  return(OP)
}

#set up the MCMC run parameters - Number of samples returned is T/SubS
T=1000000; SubS=100; Nsamp=T/SubS

#run the MCMC for the two models
#About 10 minutes for about 1e7 samples, ACT about is about 100 hence SubS.

#run the MCMC for the two models
OP=rcd.mcmc(nd,L,U,T,SubS,shrink=TRUE,w=0.25,scale.move=TRUE)
OP.shrink=as.mcmc(OP)

OP=rcd.mcmc(nd,L,U,T,SubS,shrink=FALSE)
OP[,1]=apply(OP[,4:10],1,max)-apply(OP[,4:10],1,min);
OP.unif=as.mcmc(OP[,c(1,4:11)])

nvar=length(colnames(OP))

#standard output analysis for both runs - convergence looks OK
OP=OP.shrink
xyplot(OP)
plot(as.matrix(OP)[,2],as.matrix(OP)[,3],pch='.',main='scatterplot of psi1 v. psi2',xlab='psi1',ylab='psi2')
levelplot(OP,scales = list(x = list(rot = 45))) 
effectiveSize(OP)
par(mfrow=c(3,4),oma=c(1,1,1,1));                   #acf plots
for (i in 1:nvar) {
  par(mai=0.2*c(1,1,1,1)); 
  plot(acf(OP[,i],lag.max=200,plot=F),type='l',ann=F,xaxp=c(0,200,2),yaxp=c(0,1,1)); 
  text(50,0.8,colnames(OP)[i])
}

OP=OP.unif
xyplot(OP)
levelplot(OP,scales = list(x = list(rot = 45))) 
effectiveSize(OP)
par(mfrow=c(3,3),oma=c(1,1,1,1));                   #acf plots
for (i in 1:(nvar-2)) {
  par(mai=0.2*c(1,1,1,1)); 
  plot(acf(OP[,i],lag.max=200,plot=F),type='l',ann=F,xaxp=c(0,200,2),yaxp=c(0,1,1)); 
  text(50,0.8,colnames(OP)[i])
}

#the span is of particular interest - do the two priors lead to different conclusions?
#pdf('unifshrinkspanhist.pdf',8,4)
par(mfrow=c(1,2))
hist(OP.shrink[,1],breaks=seq(0,U-L,length.out=50),xlab='span',main='Posterior span, shringage prior')
hist(OP.unif[,1],breaks=seq(0,U-L,length.out=50),xlab='span',main='Posterior span, uniform prior')
#dev.off()

#HPD interval for span
HPDinterval(OP.shrink[,1])
HPDinterval(OP.unif[,1])

#harmonic mean estimators for marginal likelihood
#p(y|unif prior)
a=1/(mean(1/exp(OP.unif[,9])))      
#p(y|shrinkage prior)
b=1/(mean(1/exp(OP.shrink[,11])))
b/a
#the shrinkage prior is preferred but we shouldnt trust 
#this estimator it has high variance

#Bridge estimate Marginal lkd of shrinkage prior = unif Span
#sim prior
n.samp=dim(OP.shrink)[1]; n.samp.p=100*n.samp;
TH.pri=matrix(NA,n.samp.p,nd)
for (k in 1:n.samp.p) {
  span=runif(1,min=0,max=U-L); 
  lower=runif(1,min=L,max=U-span); 
  upper=lower+span;
  TH.pri[k,]=runif(nd,min=lower,max=upper); 
}
TH.pos=OP.shrink[,4:10]
L.pri=rep(NA,n.samp.p); L.pos=rep(NA,n.samp)
for (k in 1:n.samp) {L.pos[k]=llk(TH.pos[k,]);}
for (k in 1:n.samp.p) {L.pri[k]=llk(TH.pri[k,]);}
ml.shrink=mean(exp(L.pri/2))/mean(exp(-L.pos/2))

#Bridge estimate Marginal lkd of unif prior
#sim prior
n.samp=dim(OP.unif)[1]; n.samp.p=100*n.samp;
TH.pri=matrix(NA,n.samp.p,nd)
for (k in 1:n.samp.p) {
  TH.pri[k,]=runif(nd,min=L,max=U)
}
TH.pos=OP.unif[,2:8]
L.pri=rep(NA,n.samp.p); L.pos=rep(NA,n.samp)
for (k in 1:n.samp) {L.pos[k]=llk(TH.pos[k,]);}
for (k in 1:n.samp.p) {L.pri[k]=llk(TH.pri[k,]);}
ml.unif=mean(exp(L.pri/2))/mean(exp(-L.pos/2))

(Bayes.Factor=ml.shrink/ml.unif)

#Favours shrinkage fairly clearly.




