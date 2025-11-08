##########################################################################
##########################################################################
# RJ-MCMC for spike and slab prior, repeats the model averaging analysis 
# of the swiss dataset

library(lattice)
library(MCMCpack)

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
y=sw$Fertility
X=model.matrix(lm(Fertility~Infant.Mortality+Examination+Education+Catholic+Agriculture,
                  data=sw)) #include intercept
n=dim(X)[1]; p=dim(X)[2];

#set known prior and lkd parameters
xi=3/p        #prior probability covariate is in model - 1/2 is uniform

#given the model and regression pars just standard normal linear model
#mo is the subset of covariates in the model
log.lkd<-function(th,mo,sig=sigma,yd=y) {
  #return(0) #for debig check
  if (length(mo)==0) {
    llk=sum(dnorm(yd,mean=0,sd=sig,log=TRUE)) #can have empty model - lin pred=0
  } else {
    if (length(mo)==1) {
      llk=sum(dnorm(yd,X[,mo]*th,sd=sig,log=TRUE))
    } else {
      llk=sum(dnorm(yd,X[,mo]%*%th,sd=sig,log=TRUE))
    }
  }
  return(llk)
}

#here I take simple priors - the focus is on the MCMC
log.prior.theta <- function(th) {
  #return(sum(dnorm(th,mean=0,sd=3,log=TRUE))) #covariates and response are O(1)
  return(sum(dt(th,df=3,log=TRUE))) #covariates and response are O(1)
}

#simulate a single draw from the theta prior
rprior.theta <- function() {
  #return(rnorm(1,mean=0,sd=3))
  return(rt(1,df=3))
}

log.prior.sigma <- function(sigma) {
  return(-log(sigma))
}

#prior over models - binomial(p,xi) (say, puts each covariate in model wp xi)
log.prior.model<-function(m,pr=xi,nc=p) {
  ne=length(m); 
  return(ne*log(xi)+(nc-ne)*log(1-xi))
}

#handle interesting sample() behavior when length(x)=1
resample <- function(x, ...) x[sample.int(length(x), ...)]

#setup MCMC
K=100000; SS=10; Nsamp=K/SS; #run for lect was 100 times longer
TH=M=vector('list',Nsamp)      #save theta and m - variable dim so go in list
LLK=SG=matrix(NA,Nsamp,1)      #save the log-likelihood and sample sigma-values
colnames(LLK)<-c('log-lkd')
colnames(SG)<-c('sigma')

#fixed dimension rw MH jump size
w=1;

#initialise at the MLE with all pars in the model
m=1:p;                           #parameter i=1,...,p  
theta=coef(sw.lm)            #init all in and MLE
sigma=summary(sw.lm)$sigma  #init sigma using residual standard error

#run MCMC
for (k in 1:K) {
  #vary the model
  nm=length(m);
  
  if (runif(1)<0.5) { #delete an effect
    if (nm==0) {
      OK=FALSE
    } else {
      OK=TRUE
      i=resample(1:nm,1)
      mp=m[-i]
      thetap=theta[-i]
      MHR_num=log.lkd(thetap,mp,sigma,y)+log.prior.model(mp)+
        log(nm);#+log.prior.theta(thetap)+log.prior.theta(theta[i]); #prior factors cancel num/denom
      MHR_den=log.lkd(theta,m,sigma,y)+log.prior.model(m)+
        log(p-nm+1);#+log.prior.theta(theta);
    }
  } else { #add an effect
    if (nm==p) {
      OK=FALSE
    } else {
      OK=TRUE
      i=resample(setdiff(1:p,m),1)
      thip=rprior.theta()
      mp=sort(c(i,m)) #put it in the right place so it matches its theta
      loc=which(mp==i); 
      if (loc==1) thetap=c(thip,theta)
      if (loc==(nm+1)) thetap=c(theta,thip)
      if (loc>1 & loc<(nm+1)) thetap=c(theta[1:(loc-1)],thip,theta[loc:nm])
      MHR_num=log.lkd(thetap,mp,sigma,y)+log.prior.model(mp)+
        log(p-nm);#+log.prior.theta(thetap);
      MHR_den=log.lkd(theta,m,sigma,y)+log.prior.model(m)+
        log(nm+1);#+log.prior.theta(theta)+log.prior.theta(thip);
    }
  }
  
  if (OK && log(runif(1))<(MHR_num-MHR_den)) {
    m=mp; nm=length(m); theta=thetap;
  }
  
  
  #vary theta by itself - help mixing
  if (nm>0) {
    i=resample(1:nm,1);
    thetap=theta; thetap[i]=theta[i]+w*(2*runif(1)-1) #propose new theta 
    MHR_num=log.prior.theta(thetap)+log.lkd(thetap,m,sigma,y);
    MHR_den=log.prior.theta(theta)+log.lkd(theta,m,sigma,y);
    if (log(runif(1))<(MHR_num-MHR_den)) {
      theta=thetap
    }
  }
  
  u=runif(1,min=0.5,max=2); sigmap=sigma*u; #propose new sigma by scaling 
  MHR_num=log.lkd(theta,m,sigmap,y)+log.prior.sigma(sigmap);
  MHR_den=log.lkd(theta,m,sigma,y)+log.prior.sigma(sigma);
  if (log(runif(1))<(MHR_num-MHR_den-log(u))) {
    sigma=sigmap
  }
  
  if (k%%SS==0) {
    TH[[k/SS]]=theta;
    M[[k/SS]]=m;
    SG[k/SS]=sigma;
    LLK[k/SS]=log.lkd(theta,m,sigma,y)
  }
}

#below superficial convergence check based on inspection of MCMC traces
#for ease of plotting map back to the theta_z notation used in spike and slab
#can recycle output analysis code from that example

Zr=t(sapply(M,function(x) {z=rep(0,p); z[x]=1; return(z)}))
Zpt.rj<-sort(round(100*table(Zs<-apply(Zr,1,paste,collapse=''))/Nsamp,1),
           decreasing=TRUE)

#posterior probabilities for the different models
sort(round(100*table(Zs<-apply(Zr,1,paste,collapse=''))/Nsamp,1),decreasing=TRUE)
#marginal posterior probability for each effect
apply(Zr,2,mean)
dZr=apply(Zr,1,sum); (dZt.rj=100*table(dZr)/Nsamp)

THZ=t( sapply( 1:Nsamp, 
               function(x) {
                 theta=rep(0,p); 
                 theta[M[[x]]]=TH[[x]]; 
                 return(theta)
               }
)
)
colnames(THZ)<-colnames(X);
DM=sapply(M,length)

#rport ESS values for MCMC
(es.rj<-effectiveSize(as.mcmc(cbind(THZ,SG))))

#now some graphs

#burn-in 10 and possible subsampling for plotting when doing long runs
s=seq(10,Nsamp,by=1) #make by=5 or 10 if doing a very long run
THZs<-THZ[s,]; Zs<-Zr[s,]; SGs<-SG[s,,drop=FALSE]; DMs<-DM[s]

#classic MCMC parameter run-traces
xyplot(as.mcmc(cbind(THZs,SGs)))

#pdf('RJswissMdim24.pdf',8,4)
par(mfrow=c(1,2))
plot(DMs,type='l',xlab='MCMC step (thinned)',ylab='Dimension')#,main='RJ-MCMC/Spike and Slab prior swiss data')
hist(DM,breaks=seq(min(DM)-0.5,max(DM)+0.5,by=1),main='',xlab='dim(theta)',freq=FALSE)
#dev.off()

#how does the joint posterior dbn for the top two effects change between model av and fixing z^*111000
#pdf('RJswissScat24.pdf',6,6)
is=which(apply(Zs,1,function(x){all(x==c(1,1,1,0,0,0))}))
plot(jitter(THZs[,2],100),jitter(THZs[,3],100),pch=3,cex=0.5,col=1,xlab="Infant.Mortality",ylab="Examination"); 
points(THZs[is,2],THZs[is,3],pch=16,cex=0.5,col=2)
legend('topright',pch=c(3,16),col=c(1,2),legend=c('theta_z|y','theta_z|y,z^*'))
#dev.off()

#how do the marginal posterior dbns change between model av and fixing z^*111000
boxplot(cbind(THZ,SG),at=(0:p)-0.15,boxwex=0.2)
i=which(apply(Zr,1,function(x){all(x==c(1,1,1,0,0,0))}))
boxplot(cbind(THZ[i,],SG[i]),add=TRUE,at=(0:p)+0.15,col=2,boxwex=0.2,xaxt='n')

#the linear predictor is eta=X*theta - get the HPD for eta 
#using model averaging and compare if select model first
th.ma=apply(THZ,2,mean) #posterior mean with model avergaing
eta.ma=X%*%th.ma
i=which(apply(Zr,1,function(x){all(x==c(1,1,1,0,0,0))}))
th.zs=apply(THZ[i,],2,mean) #posterior mean for the MAP model
eta.zs=X%*%th.zs

(hpd.ma=HPDinterval(as.mcmc(t(X%*%t(THZ)))))
(hpd.zs=HPDinterval(as.mcmc(t(X%*%t(THZ[i,])))))

#they look similar but a few differ a bit
library(EnvStats)
errorBar(x=4*(1:n)-2.25,y=eta.zs,lower=hpd.zs[,1],upper=hpd.zs[,2],col=2,bar.ends.size=0.1,gap=FALSE,pch='.')
#points(4*(1:n)-2,eta.zs,pch=16,col=2)
errorBar(x=4*(1:n)-1.75,y=eta.ma,lower=hpd.ma[,1],upper=hpd.ma[,2],bar.ends.size=0.1,gap=FALSE,pch='.',add=TRUE)
#points(4*(1:n)-2,eta.ma,pch=16,col=1)


