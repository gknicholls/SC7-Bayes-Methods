
#PS4 - sample answers to the RJ-MCMC questions

library(coda) #MCMC output analysis stuff

#######################################################################################
#PS4 Q1

#pi(x,m)
pr<-function(x,M) {
  if (M==1) {
    return(1/3)
  } else {
    return(4*x/3)
  }
}

#MCMC setup
K=100000; SS=10; Nsamp=K/SS; X=M=rep(0,Nsamp)
x=1/2; m=1;

for (k in 1:K) {
  if (m==1) {
    #we are adding a dimension
    xp=rbeta(1,1/2,1/2); mp=2;
    MHR=pr(xp,mp)/(pr(x,m)*dbeta(xp,1/2,1/2))
  } else {
    #we are deleting a dimension
    xp=1/2; mp=1;
    MHR=pr(xp,mp)*dbeta(x,1/2,1/2)/pr(x,m)
  }
  #MHR's from notes; proposals use Beta(1/2,1/2)
  if (runif(1)<MHR) {
    x=xp
    m=mp
  }
  if (k%%SS==0) {X[k/SS]=x; M[k/SS]=m} 
}

#display the output
#pdf('RJMCMCsimple.pdf',10,4)
par(mfrow=c(1,3)); 
plot(X[1:30],type='l',main='MCMC output'); points(X[1:30],pch=16); abline(h=0.5,lty=2)
hist(X,breaks=seq(0,1,length.out=20),main='');
plot(ecdf(X),main='CDF, estimated/exact',xlab='X')
xa=seq(0,1,0.02); lines(xa,(xa<0.5)*2*xa^2/3+(xa>=0.5)*(1/3+2*xa^2/3),col=2)
#dev.off()

########################################################################


########################################################################
#PS4 Q3

###
#RJ MCMC - is it skew?
#
library(fGarch) #skew normal
#take a look at SN dbn
par(mfrow=c(1,1))
hist(rsnorm(10000,mean=0,sd=1,xi=3),100) #try different values to see effect of xi

#About the data: https://gksmyth.github.io/ozdasl/general/shoshoni.html
y=unlist(read.table("https://gksmyth.github.io/ozdasl/general/shoshoni.txt",header=TRUE))
hist(y)

#in case that website moves again
#y <- c(0.693,0.670,0.654,0.749,0.606,0.553,0.601,0.609,0.672,0.662,
#       0.606,0.615,0.844,0.570,0.933,0.576,0.668,0.628,0.690,0.611)

#muT=0; sgT=1; xiT=2; isT=TRUE; n=30; y=rsnorm(n,mean=muT,sd=sgT,xi=xiT) 
#check it all works using synth data if you like

#given the model and regression pars just standard normal linear model
log.lkd<-function(mu,sg,xi,is,yd=y) {
  #here and below is=1/0 means skew/not (same as m=2/1 in Q)
  if (is) {
    llk=sum(dsnorm(yd,mu,sg,xi,log=TRUE)) #you need library(fGarch) for this
  } else {
    llk=sum(dnorm(yd,mu,sg,log=TRUE))
  }
  return(llk)
}

#here I take simple priors - for hyperparameters see below
log.prior.param <- function(mu,sg,xi,is) {
  mp=dnorm(mu,mean=PM,sd=PSD,log=T)				#prior on MU
  sgp=dgamma(sg,shape=PS,rate=PR,log=T)               #prior on SIGMA
  if (is) {
    xp=dgamma(xi,shape=SHx,rate=RTx,log=T)		#prior on XI
  } else {
    xp=0
  }
  return(mp+sgp+xp)
}

#prior over models - with probability PrS=1/2 there is skew
log.prior.model<-function(is) {
  return(log(PrS*is+(1-PrS)*(1-is))) 
}

#the density for the skew/xi proposal distribution q(theta'[i]|m')
lq<-function(xi) {
  return(dgamma(xi,shape=SHx,rate=RTx,log=T))
}

#hyper pars
PM=0; PSD=3; PS=1; PR=1; SHx=1.5; RTx=1/2; PrS=1/2;

#setup MCMC
K=10000; SS=10; Nsamp=K/SS; 
TH=matrix(NA,Nsamp,4); colnames(TH)<-c('mu','sg','xi','is');
L=matrix(NA,Nsamp,2); colnames(L)<-c('lpr','llk')
#fixed dimension rw MH jump size
w=1;

#initialise standard normal
mu=0; sg=1; xi=1; is=FALSE;

#run MCMC
for (k in 1:K) {
  
  #two variable dim updates - notice we always propose the 'other' model
  if (!is) {                          #is=FALSE means normal/m=1
    #we are adding a dimension
    xip=rgamma(1,shape=SHx,rate=RTx); isp=TRUE #note this has density lq() above
    MHR_num=log.lkd(mu,sg,xip,isp,y)+log.prior.param(mu,sg,xip,isp)+log.prior.model(isp);
    MHR_den=log.lkd(mu,sg,xi,is,y)+log.prior.param(mu,sg,xi,is)+log.prior.model(is)+lq(xip);
    if (log(runif(1))<(MHR_num-MHR_den)) {
      xi=xip
      is=isp
    }
  } else {                            #is=TRUE so currently in m=2/skew-normal
    #we are killing a dimension
    xip=1; isp=FALSE;   
    MHR_num=log.lkd(mu,sg,xip,isp,y)+log.prior.param(mu,sg,xip,isp)+log.prior.model(isp)+lq(xi);
    MHR_den=log.lkd(mu,sg,xi,is,y)+log.prior.param(mu,sg,xi,is)+log.prior.model(is);
    if (log(runif(1))<(MHR_num-MHR_den)) {
      xi=xip
      is=isp
    }
  }
  
  #fixed dimension updates
  if (is) {
    xip=rgamma(1,shape=SHx,rate=RTx);
    MHR_num=log.lkd(mu,sg,xip,is,y)+log.prior.param(mu,sg,xip,is)+dgamma(xi,shape=SHx,rate=RTx,log=T);
    MHR_den=log.lkd(mu,sg,xi,is,y)+log.prior.param(mu,sg,xi,is)+dgamma(xi,shape=SHx,rate=RTx,log=T);
    if (log(runif(1))<(MHR_num-MHR_den)) {
      xi=xip
    }
  }
  
  u=runif(1,min=0.5,max=2); sgp=sg*u; #propose new sigma by scaling 
  MHR_num=log.lkd(mu,sgp,xi,is,y)+log.prior.param(mu,sgp,xi,is)
  MHR_den=log.lkd(mu,sg,xi,is,y)+log.prior.param(mu,sg,xi,is)
  if (log(runif(1))<(MHR_num-MHR_den-log(u))) {
    sg=sgp
  }
  
  mup=mu+rnorm(1,mean=0,sd=w)
  MHR_num=log.lkd(mup,sg,xi,is,y)+log.prior.param(mup,sg,xi,is)
  MHR_den=log.lkd(mu,sg,xi,is,y)+log.prior.param(mu,sg,xi,is);
  if (log(runif(1))<(MHR_num-MHR_den)) {
    mu=mup
  }
  
  if (k%%SS==0) {
    TH[k/SS,]=c(mu,sg,xi,is);
    L[k/SS,]=c(log.prior.param(mu,sg,xi,is),log.lkd(mu,sg,xi,is,y))
  }
}

#here is=TRUE is skew so if model 1 is m=1, we want !is to be true
p1.hat=mean(!TH[,4]); p1.hat #p(m=1|y)
B12=p1.hat/(1-p1.hat); B12; 1/B12 #Bayes factor (2nd number is evidence in favor of skew)

#FOLLOWING NOT NEEDED FOR PS4 QUESTION
#superficial convergence check based on inspection of MCMC traces
plot(as.mcmc(cbind(TH,L)))

#show variation in fitted distribution
hist(y,16,freq=F,xlab='Aspect ratio',ylab='',main=''); 
vx=seq(0.5,1,0.005); 
for (k in sample(1:Nsamp,20)) {
  vy=dsnorm(vx,mean=TH[k,1],sd=TH[k,2],xi=TH[k,3])
  lines(vx,vy,col=grey(0.8))
}

#overlay posterior predictive distribution
thp=apply(TH,2,mean)
vy=dsnorm(vx,mean=thp[1],sd=thp[2],xi=thp[3])
lines(vx,vy,col=2)

########################################################################


########################################################################
#PS4 Q8 - Bayesian non-parametrics - The coal-mining disaster data

library(boot) #happens to have coal-mining data

#coal mining disasters - each event records an explosion 
#in a coal mine leading to significant loss of life

min(coal); 
max(coal); 
L<-1851
U<-1951
coal.dat<-coal[coal<U] #focus on a century 1850-1950
par(mfrow=c(1,2))

#pdf('cmf.pdf')
hist(coal.dat,breaks=seq(L,U,1),main='',cex.lab=1.4) #number of events per year
#dev.off()


#interesting to look at the fitted mean rate for events over the century.
rf<-function(L,U,lam,tm) {
  #turn set of CP params into a rate fn - give rate in each year
  yr=seq(L,U-1,1)
  v=c(tm,U)
  r=numeric(length(yr))
  #for each year what was the rate in that year
  for (i in 1:length(yr)) {
    j=which.max(v>yr[i]) #which v-interval does the year fall in
    r[i]=lam[j]          #use the lambda for that interval
  }
  return(r) #a vector giving the rate in each year
}


#####################################################################################################
#Answer to the RJ-MCMC exercise - just the birth death move as it is irreducible by itself

#model CP events as PP(mu.cp)
mu.cp=5
#model yearly rate - constant between CP's with mean mu.lam events/yr
mu.lam=1

log.post<-function(n,lam,tm,prior=FALSE){
  t1=dpois(n,mu.cp,log=TRUE)
  t2=lfactorial(n)-n*log(U-L)
  t3=sum(dexp(lam,rate=1/mu.lam,log=TRUE))
  if (prior) {
    t4=0
  } else {
    counts=hist(coal.dat,breaks=tm,plot=FALSE)$counts
    t4=sum(lfactorial(counts)-counts*log(diff(tm))+dpois(counts,lambda=lam*diff(tm),log=TRUE))
  } 
  lp=t1+t2+t3+t4
  return(lp)
}
resample <- function(x, ...) x[sample.int(length(x), ...)]

#MCMC setup
T=100000; SS=100;
TM<-LM<-rep(list(),T/SS); N<-rep(NA,T/SS); LPOST<-rep(NA,T/SS); 
n.accept=0
DEBUG=FALSE #check code by checking we recover the prior when we switch off the LKD.

#init state
n<-rpois(1,lambda=mu.cp);                   #number of CP events
lam<-rexp(n+1,rate=1/mu.lam);               #rates in each interval between CP's
tm<-c(L,sort(runif(n,min=L,max=U)),U);      #times of the CP's including end times
olp=log.post(n,lam,tm,prior=DEBUG)

for (t in 1:T) {
  p.up=(n==0)+0.5*(n>0); p.down=0.5; #cant delete if no CP's - we only use p.down from n>0
  if (runif(1)<p.up) {
    #add a CP
    move=1 #'add'
    
    n.new=n+1
    i=resample(1:(n+1),1) #choose an interval to split - the first interval is 1, the last is n+1
    tm.p=runif(1,min=tm[i],max=tm[i+1])    #pick a new time uniformly in the split interval
    
    lam.p=rexp(2,rate=1/mu.lam)            #generate 2 new lambda rates in this interval
    tm.new=c(tm[1:i],tm.p,tm[(i+1):(n+2)]) #tm's are 1,...,n+2 as 1 is L and n+2 is U
    if (n==0) {lam.new=lam.p} else {       #put the new times and lambda in the right places to form "new"
      if (i==1) {lam.new=c(lam.p,lam[2:(n+1)])} 
      if (i==(n+1)) {lam.new=c(lam[1:n],lam.p)}
      if (i>=2 & i<=n) {lam.new=c(lam[1:(i-1)],lam.p,lam[(i+1):(n+1)])}
    }
    nlp=log.post(n.new,lam.new,tm.new,prior=DEBUG)
    
    #work out the Q/Q factors
    lam.p.ld=sum(dexp(lam.p,rate=1/mu.lam,log=TRUE))
    q.up=-log(tm[i+1]-tm[i])-log(n+1)+log(p.up)+lam.p.ld
    
    #when we delete we choose a tm to delete - there are n+1 of these in the new state
    lam.ld=dexp(lam[i],rate=1/mu.lam,log=TRUE)
    q.down=-log(n+1)+log(p.down)+lam.ld
    
    qq=q.down-q.up
    
  } else {
    
    #delete a CP
    move=2 #'delete'
    
    n.new=n-1                            #should never be -1 as we cant propose it - see p.down formula
    i=resample(2:(n+1),1)                #choose a CP to delete - dont choose an endpoint L/U at i=1/n+2
    lam.p=rexp(1,rate=1/mu.lam)          #when we delete a CP we pick a new lambda for the merged interval
    tm.new=tm[-i]                        #tm's are 1,...,n+2 as 1 is L and n+2 is U
    if (n.new==0) {lam.new=lam.p} else { #form the new state putting things in the right places
      if (i==2) {lam.new=c(lam.p,lam[3:(n+1)])} 
      if (i==(n+1)) {lam.new=c(lam[1:(n-1)],lam.p)}
      if (i>=3 & i<=n) {lam.new=c(lam[1:(i-2)],lam.p,lam[(i+1):(n+1)])}
    }
    nlp=log.post(n.new,lam.new,tm.new,prior=DEBUG)
    
    #when we add to get back there are n-1 CP's in the new state so n intervals to choose from
    lam.ld=sum(dexp(lam[c(i-1,i)],rate=1/mu.lam,log=TRUE))  #this is like g_{m',m}(u')
    q.up=-log(tm[i+1]-tm[i-1])-log(n)+log(p.up)+lam.ld
    
    #when we delete we choose a tm to delete - there are n of these in the current state
    lam.p.ld=dexp(lam.p,rate=1/mu.lam,log=TRUE)             #this is like g_{m,m'}(u)
    q.down=-log(n)+log(p.down)+lam.p.ld
    
    qq=q.up-q.down
    
  }
  
  MHR=nlp-olp+qq
  if (log(runif(1))<MHR) {
    olp=nlp
    lam=lam.new
    tm=tm.new
    n=n.new
    n.accept=n.accept+1
  }
  
  if (t%%SS==0) {
    LM[[t/SS]]=lam #save into list as dim varies
    TM[[t/SS]]=tm[c(-1,-(n+2))] 
    N[t/SS]=n
    LPOST[t/SS]=olp #save the evolving log posterior values as a summary stat for convergence check
  }
  
}

#Output analysis
N.samples=length(N)
plot(LPOST)
(ESS<-effectiveSize(N))
#very rudimentary convergence diagnostics

#compare marginal posteriors with priors
#N - number of CP's
par(mfrow=c(1,2))
cp.n=sort(unique(N))
plot(cp.n,table(N)/length(N),xlab='number of CPs',ylab='prob',type='l')
x=0:max(cp.n)
lines(x,dpois(x,mu.cp),col=2)
legend('topright',lty=c(1,1),col=c(1,2),legend=c('in samples','in prior'))
#the prior and posterior are quite similar, sometimes a cause for concern
#as it can suggest the prior is dominating the data. Here I think we just 
#guessed right

#lambda[1] - the rate in the first interval
lam1=sapply(LM,function(x) x[1])
x=hist(lam1,freq=FALSE,floor(sqrt(ESS)), #posterior
       xlab='lambda in first interval',ylab='density',main='lambda samples and prior')$mids
lines(x,dexp(x,rate=1/mu.lam)) #the prior

#rf() function needed here - defined above

#display the variation in the estimated rate function
par(mfrow=c(1,1))
hist(coal.dat,breaks=seq(L,U,1))

#average the sampled rates - plot each rate function
yr=seq(L,U-1,1) #the vector of years
r.rj=rep(0,length(yr))
for (m in 1:N.samples) {
  lambda.rate.fn.rj=rf(L,U,LM[[m]],TM[[m]])
  lines(yr,lambda.rate.fn.rj,col=grey(0.5))
  r.rj=r.rj+lambda.rate.fn.rj
}
#plot the posterior mean rate
mean.rate.rj=r.rj/N.samples
lines(yr,mean.rate.rj,col=2,lwd=2)

#estimate prob for CP in each year
tot.rj=rep(0,length(yr))
for (i in 1:N.samples) {
  #TMo[[i]]) is the lists of CP-years in the i'th sample
  index=floor(TM[[i]])-L+1      #get the index of each year in the yr vector
  tot.rj[index]=tot.rj[index]+1 #add one to those year
}
cp.prob.rj=tot.rj/N.samples
plot(yr,smooth(cp.prob.rj),pch=16,type = 'l',col=2,lwd=2,ylab='prob for a CP per year')
points(coal.dat,runif(length(coal.dat),0,0.025),pch='+')



