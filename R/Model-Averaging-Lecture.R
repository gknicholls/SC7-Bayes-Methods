#Model averaging

library(MCMCpack)

############
#example model averaging over choice of link functions
#Challenger O-ring failures (fail=1, OK=0)
fail=c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0)
temp=scale(c(53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81))
plot(temp,fail)

#this is classical so MCMC code targeting posterior is available in R library "MCMCpack"
#results in lectures use run length N x 10.

N=1000000
prior.std.dev=3
v=1/prior.std.dev^2 #precision v=1/9 is sigma^2=3 in prior
b.p=MCMCprobit(fail~temp,b0=c(0,0),B0=c(v,v),seed=4,mcmc=N,marginal.likelihood='Laplace')
b.l=MCMClogit(fail~temp,b0=c(0,0),B0=c(v,v),seed=5,mcmc=N,marginal.likelihood='Laplace')

logit<-function(beta,x=temp) {eta=beta[1]+beta[2]*x; exp(eta)/(1+exp(eta))}
probit<-function(beta,x=temp) {eta=beta[1]+beta[2]*x; pnorm(eta)}

lkd<-function(beta,logit=T) {
  if (logit) {pb=logit(beta)} else {pb=probit(beta)}
  return(sum(dbinom(fail,1,pb,log=T)))
}

#First compute Bayes factors - gives weights for the two models
l.pr=l.lo=pr.lo=pr.pr=l.lo.pr=l.pr.lo=rep(NA,N)
for (k in 1:N) { #takes about a minute for N=1e6
  l.lo[k]=lkd(b.l[k,],logit=T); 
  l.pr[k]=lkd(b.p[k,],logit=F);
  l.lo.pr[k]=lkd(b.p[k,],logit=T); 
  l.pr.lo[k]=lkd(b.l[k,],logit=F);
  pr.lo[k]=sum(dnorm(b.l[k,],0,9,log=T))
  pr.pr[k]=sum(dnorm(b.p[k,],0,9,log=T))
}

#and Bridge using h=1/sqrt(p1*p2) - recommended by MW as good place to start
BF=mean(exp(l.lo.pr/2-l.pr/2))/mean(exp(l.pr.lo/2-l.lo/2))

#now compute pi(theta|y) averged over the two models
#pdf('linkfnMA.pdf',8,8)
i=1
bp=seq(min( c(b.p[,i],b.l[,i]) ), max(c(b.p[,i],b.l[,i])) ,length.out=100)
d.p=hist(b.p[,i],breaks=bp,plot=F)
d.l=hist(b.l[,i],breaks=bp,plot=F)
plot(d.p$mids,d.p$density,type='l',xlim=c(-2,0.5),cex.lab=1.4,
     xlab='beta_1 values',ylab='marginal posterior density',lwd=2); 
lines(d.l$mids,d.l$density,col=2,lwd=2)
legend('topright',c('probit posterior','logit','model average'),col=1:3, lty=c(1,1,2), lwd=c(2,2,2))

w1=BF/(1+BF); w2=1/(1+BF);
lines(d.l$mids,w1*d.l$density+w2*d.p$density,col=3,lty=2,lwd=2)
#dev.off()


##########################################################################
##########################################################################
# A more substantial example of model averaging

library(lattice)

###

#Data: swiss data set relating family fertility of people

library(MASS)
sw<-swiss
#head(sw)
sw[,-1]<-log((swiss[,-1]+1)/(101-swiss[,-1]))
sw[,1]<-scale(sw[,1])
n<-dim(sw)[1]; p<-dim(sw)[2] #counting intercept

#remove some outliers - see PS3 for Bayesian treatment of outliers
sw.lm<-lm(Fertility~Infant.Mortality+Examination+Education+Catholic+Agriculture,
           data=sw)
ol<-cooks.distance(sw.lm)>(8/(n-2*p))
sw<-sw[-which(ol),]

sw.lm<-lm(Fertility~Infant.Mortality+Examination+Education+Catholic+Agriculture,
          data=sw)
summary(sw.lm)
y=sw$Fertility
X=model.matrix(lm(Fertility~Infant.Mortality+Examination+Education+Catholic+Agriculture,
                         data=sw)) #include intercept
#  cbind(rep(1,n),sw[,c("Agriculture","Examination","Education","Catholic","Infant.Mortality")])
n=dim(X)[1]; p=dim(X)[2];

#set known prior and lkd parameters
xi=3/p        #prior probability covariate is in model - 1/2 is uniform

#given the model and regression pars just standard normal linear model
log.lkd<-function(th,z,sig=sigma,yd=y) {
  return(sum(dnorm(yd,X%*%(z*th),sd=sig,log=TRUE)))
}

#here I take simple priors - the focus is on the MCMC
log.prior.theta <- function(th) {
  #return(sum(dnorm(th,mean=0,sd=3,log=TRUE))) #covariates and response are O(1)
  return(sum(dt(th,df=3,log=TRUE))) #covariates and response are O(1)
}
log.prior.sigma <- function(sigma) {
  return(-log(sigma)) #Jeffreys
}

#prior over models - binomial(p,xi) (say, puts each covariate in model wp xi)
log.prior.model<-function(z,pr=xi,nc=p) {
  ne=sum(z); 
  return(ne*log(xi)+(nc-ne)*log(1-xi))
}

#setup MCMC
K=10000000; SS=1000; Nsamp=K/SS; #run for lect was 100 times longer
TH=Z=THZ=matrix(NA,Nsamp,p); 
colnames(TH)<-colnames(THZ)<-colnames(Z)<-colnames(X);
SG=matrix(NA,Nsamp,1); colnames(SG)<-'sigma'
colnames(SG)<-c('sigma')
#fixed dimension rw MH jump size
w=1;

#initialise at the MLE with all pars in the model
z=rep(1,p);                      #parameter i=1,...,p is in z[i]=1 or out z[i]=0
theta=coef(sw.lm)           #init all in and MLE
sigma=summary(sw.lm)$sigma  #init sigma using residual standard error

#run MCMC
for (k in 1:K) {
  #vary the z's
  i=sample(1:p,1);    #choose a parameter to toggle in/out
  zp=z; zp[i]=1-z[i]; #propose the new model
  #if (all(diff(zp)!=1)) {
    MHR_num=log.lkd(theta,zp,sigma,y)+log.prior.model(zp);
    MHR_den=log.lkd(theta,z,sigma,y)+log.prior.model(z);
    if (log(runif(1))<(MHR_num-MHR_den)) {
      z=zp
    }
  #}
    
  #vary theta
  i=sample(1:p,1);
  thetap=theta; thetap[i]=theta[i]+w*(2*runif(1)-1) #propose new theta 
  MHR_num=log.prior.theta(thetap)+log.lkd(thetap,z,sigma,y);
  MHR_den=log.prior.theta(theta)+log.lkd(theta,z,sigma,y);
  if (log(runif(1))<(MHR_num-MHR_den)) {
    theta=thetap
  }
  
  u=runif(1,min=0.5,max=2); sigmap=sigma*u; #propose new sigma by scaling 
  MHR_num=log.lkd(theta,z,sigmap,y)+log.prior.sigma(sigmap);
  MHR_den=log.lkd(theta,z,sigma,y)+log.prior.sigma(sigma);
  if (log(runif(1))<(MHR_num-MHR_den-log(u))) {
    sigma=sigmap
  }
  
  if (k%%SS==0) {
    TH[k/SS,]=theta;
    Z[k/SS,]=z;
    THZ[k/SS,]=theta*z;
    SG[k/SS]=sigma;
  }
}

#superficial convergence check based on inspection of MCMC traces and ESS
Nsamp=dim(Z)[1]
s=seq(10,Nsamp,by=5)

xyplot(as.mcmc(cbind(TH[s,],SG[s,,drop=FALSE])))
#pdf('IVswissRmcmc24.pdf',6,6)
xyplot(as.mcmc(cbind(THZ[s,],SG[s,,drop=FALSE])))
#dev.off()
(es.ss<-effectiveSize(as.mcmc(cbind(THZ,SG))))
#posterior probabilities for the different models
(Zpt.ss<-sort(round(100*table(Zs<-apply(Z,1,paste,collapse=''))/Nsamp,1),decreasing=TRUE))
#marginal posterior probability for each effect
#apply(Z,2,mean)
dZ=apply(Z,1,sum); (dZt.ss=100*table(dZ)/Nsamp)
hz=hist(dZ,xlab='dim=sum_i z_i',main='dim distribution S&S',
        breaks=seq(min(dZ)-0.5,max(dZ)+0.5,1))

THZs<-THZ[s,]; Zs<-Z[s,]; SGs<-SG[s,,drop=FALSE]

#how does the joint posterior dbn for the top two effects change between model av and fixing z^*111000
#pdf('VswissScat24.pdf',6,6)
is=which(apply(Zs,1,function(x){all(x==c(1,1,1,0,0,0))}))
plot(jitter(THZs[,2],100),jitter(THZs[,3],100),pch=3,cex=0.5,col=1,xlab="Infant.Mortality",ylab="Examination"); 
points(THZs[is,2],THZs[is,3],pch=16,cex=0.5,col=2)
legend('topright',pch=c(3,16),col=c(1,2),legend=c('theta_z|y','theta_z|y,z^*'))
#dev.off()

#how do the marginal posterior dbns change between model av and fixing z^*111000
boxplot(cbind(THZ,SG),at=(0:p)-0.15,boxwex=0.2)
i=which(apply(Z,1,function(x){all(x==c(1,1,1,0,0,0))}))
boxplot(cbind(THZ[i,],SG[i]),add=TRUE,at=(0:p)+0.15,col=2,boxwex=0.2,xaxt='n')

#the linear predictor is eta=X*theta - get the HPD for eta 
#using model averaging and compare if select model first
th.ma=apply(THZ,2,mean) #posterior mean with model averaging
eta.ma=X%*%th.ma
i=which(apply(Z,1,function(x){all(x==c(1,1,1,0,0,0))}))
th.zs=apply(THZ[i,],2,mean) #posterior mean for the MAP model
eta.zs=X%*%th.zs

(hpd.ma=HPDinterval(as.mcmc(t(X%*%t(THZ)))))
(hpd.zs=HPDinterval(as.mcmc(t(X%*%t(THZ[i,])))))

#they look similar but a few differ a bit
library(EnvStats)
errorBar(x=4*(1:n)-2.25,y=eta.zs,lower=hpd.zs[,1],upper=hpd.zs[,2],
         col=2,bar.ends.size=0.1,gap=FALSE,pch='.',
         xlab='observation index i=1:45',ylab='HPD for linear predictor eta_i')
#points(4*(1:n)-2,eta.zs,pch=16,col=2)
errorBar(x=4*(1:n)-1.75,y=eta.ma,lower=hpd.ma[,1],upper=hpd.ma[,2],bar.ends.size=0.1,gap=FALSE,pch='.',add=TRUE)
#points(4*(1:n)-2,eta.ma,pch=16,col=1)
