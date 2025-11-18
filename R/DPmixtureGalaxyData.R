#Bayes Methods - Normal-Dirichlet-Process mixture
#fit to galaxy data

library(MCMCpack) #dinvgamma for transparency

###############################################
#MCMC fitting the galaxy data

y=scan("http://www.stats.ox.ac.uk/~nicholls/BayesMethods/galaxies.txt") #82 galaxy transverse velocities
nd=length(y)           #number of data points

#the priors for the means are not quite the same as those I use in the RJ-MCMC example
#because I want the priors to be conjugate I use inverse gamma for sigma^2 rather than
#gamma on sigma as I used before mu~N(20,10) as before, but sigma^2~IGamma(alpha0,beta0)
#the choice of alpha0 and beta0 is intended to give a similar prior for sigma

log.prior.param<-function(mu,sg,mu0=20,sigma0=10,alpha0=2,beta0=1/9) {
  a=sum(dnorm(mu,mean=mu0,sd=sigma0,log=TRUE))
  b=sum(log(dinvgamma(sg^2,shape=alpha0,scale=beta0)))
  return(a+b)
}

sigma0=10; mu0=20
alpha0=2; beta0=1/9 #IGamma sg^2 mean=beta/(alpha-1)=9 so sg~3

#choice of alpha based on simulations above
alpha=1

#initialise MCMC state - S is known given mu but it is convenient to maintain both
S=S.init=rep(1,nd)
mu=mu.init=mean(y)
sg=sg.init=sd(y)

#run parameters and output storage
J=10000; SS=10;  #pretty minimal but enough for plausible estimates
CL=matrix(NA,J/SS,nd); PL=matrix(NA,J/SS,3); TH=list();

#debug check
db=1 #set to 0 to get CRP prior over dbns

#MCMC J steps each step sweeps mu, sigma and the cluster membership
set.seed(1)
for (j in 1:J) {

  nS=table(S) #the number of pts in each cluster
  K=max(S)    #the number of clusters

  #go through each cluster and update mu[k] and sg[k] given the cluster membership S 
  for (k in 1:K) {
    i.in.k=which(S==k); yk=y[i.in.k] #which data points in this cluster
    nk=nS[k]                         #how many pts in this cluster
    if (nk>0) {
      #sample mu[k] given mu[-k], sg and S
      sgk=sg[k]
      bm=sqrt(1/(nk/sgk^2+1/sigma0^2))
      am=bm^2*(sum(yk)/sgk^2+mu0/sigma0^2)
      mu[k]=rnorm(1,mean=am,sd=bm)
      
      #sample sg[k] given sg[-k], mu and S
      as=alpha0+nk/2
      bs=beta0+(1/2)*sum((yk-mu[k])^2)
      sg[k]=1/sqrt(rgamma(1,shape=as,rate=bs))
    }
  } 
  
  #go through each point y[i] and update its cluster membership S[i]
  for (i in 1:nd) {
    #the number of pts per cluster and the number of clusters may have changed 
    nS=table(S)
    K=max(S)

    #we need the weights with i removed 
    nSp=nS; 
    k=S[i]; 
    if (nS[k]==1) {
      nSp[k]=alpha
      pold=dnorm(y[i],mean=mu,sd=sg)^db*nSp
      pnew=0
    } else {
      nSp[k]=nSp[k]-1
      pold=dnorm(y[i],mean=mu,sd=sg)^db*nSp
      
      #if we move i to a new cluster we need mu and sg values for the new cluster
      #only possible if nS[i]>1
      mkn=rnorm(1,mean=mu0,sd=sigma0)
      sgn=1/sqrt(rgamma(1,shape=alpha0,rate=beta0))
      
      #alpha*p(y[i]|mu.new,sig.new)
      pnew=alpha*dnorm(y[i],mean=mkn,sd=sgn)^db
    } 
    
    #the new cluster is either one of the old ones or a new one
    #if nS[i]=1 then pr[1]=0 (cant move to a new K+1-cluster)
    pr=c(pnew,pold);
    
    #Normalise
    pr=pr/sum(pr)
    
    #pick a new cluster
    kn=sample(0:K,1,prob=pr)
    
    if (kn==0) { #new cluster, append mu.new and sig.new to mu, sg 
      S[i]=K+1
      mu=c(mu,mkn)
      sg=c(sg,sgn)
    } else {     #put the pt in cluster kn
      S[i]=kn
    }

    #keeping the cluster labels packed (ie 1:K) and sorted (ie min(S_k)<min(S_k') -> k<k')
    #probably a bad choice of data structure it cant be this hard! Bugzone.
    if (nS[k]==1 & kn!=k) {
      #pack the labels
      ib=which(S>k); 
      S[ib]=S[ib]-1; 
      mu=mu[-k]; 
      sg=sg[-k]
      # #standardise labeling
      # b=rep(NA,max(S)); b[S[1]]=1
      # d=1; So=S;
      # for (a in 1:length(S)) {
      #   x=S[a];
      #   if (is.na(b[x])) {
      #     d=d+1
      #     b[x]=d;
      #   }
      #   So[a]=b[x]
      # }
      # o=order(b)
      # mu=mu[o]; sg=sg[o]
    }
  }
  
  #collect samples from the MCMC every SS steps
  if (j%%SS==0) {
    #collect a certain amount of redundant stuff
    TH[[j/SS]]=list(w=table(S),mu=mu,sigma=sg) #was something like list(w=rdirichlet(1,alpha/max(S)+table(S)),mu=mu,sigma=sg)
    CL[j/SS,]=S
    PL[j/SS,]=c(length(mu),min(S),max(S))  #check we are keeping the cluster labels 'packed'
  }

}

##
if (db==0) {
  summary(as.mcmc(PL[,1]))$statistics[c(1,4)]
  n=length(y)
  EK=sum(alpha/(alpha+0:(n-1))); EK #should match in 2 x sigma
}

###
#Superficial Output Analysis
effectiveSize(as.mcmc(CL)) #cluster membership for each data point
library(lattice); 
plot(as.mcmc(PL[,1])) #number of components K in fitted DP mixture

###

###
#visualise co-occurence matrix p(i,j) = prob i,j in same cluster 
Nsamp=J/SS;
com=matrix(0,nd,nd)
for (i in 1:nd) {
  for (j in 1:nd) {
    for (m in 1:Nsamp) {
      com[i,j]=com[i,j]+(CL[m,i]==CL[m,j])
    }
  }
}
com=com/(Nsamp)
image(com)
###


###
#posterior dbn over number of components - exactly the same as for RJ-MCMC
#pdf('DPmixtureGhistM.pdf',8,4)
hist(d<-PL[,1],breaks=0.5:(0.5+max(PL[,1])),main='',freq=FALSE,
     xlab='Number of components m in mixture',ylab='density',xlim=c(0,1+max(PL[,1])))
#dev.off()
table(d)/Nsamp

#loglkd including the weights (not used in MCMC)
log.lkd<-function(y,w,mu,sigma) {
  tot=0
  for (i in 1:length(y)) {
    tot=tot+log(sum(w*dnorm(y[i],mean=mu,sd=sigma)))
  }
  return(tot)
}



###
#compute posterior predictive distributions - weights (w) vs. means (mu) scatter plot
#pdf('DPmixtureGscatterMUW.pdf',8,4)
L=200
x=seq(0,40,length.out=L)
den=matrix(0,max(d)-min(d)+1,L)
plot(c(),c(),xlim=c(0,50),ylim=c(0,1),xlab='MU, colored by number of clusters in state',ylab='W component weight')
for (k in 1:Nsamp) { 
  #hat.p(y'|y) = E_{theta^*,S|y}[ ap(y')/(a+n)+sum_k p(y'|theta^*_k)n_k/(a+n) ]
  #now E_{theta^*,S} is estimated using posterior samples
  #p(y')=E_{theta'}[ p(y'|theta') ] is estimated using prior samples theta' ~ h()
  #hat.p(y')=p(y'|theta')
  #simulate theta'
  m.prior=rnorm(1,mean=mu0,sd=sigma0)
  sg.prior=1/sqrt(rgamma(1,shape=alpha0,rate=beta0))
  #put theta' in the "alpha" position of the w vector
  mu=c(TH[[k]]$mu,m.prior); sigma=c(TH[[k]]$sigma,sg.prior)
  w=c(TH[[k]]$w,alpha); 
  w<-w/sum(w);
  #now w=(n_1,...,n_K,alpha)/(n+alpha) and the prior sample theta' picks up the alpha weight.
  ind=d[k]-min(d)+1
  points(mu,w,pch=16,cex=0.5,col=ind) #not very helpful to me - locations of means colored by number of clusters
  den[ind,]=den[ind,]+exp(apply(t(x),2,log.lkd,w,mu,sigma))
}
cden=den/as.vector(table(d)) #useful if you want to look at conditional predictive given K clusters
mden=apply(den,2,sum)/Nsamp  #this is the estimate of the posterior predictive.
#dev.off()



#plot posterior predictive distributions
#pdf('DPmixtureGppd.pdf',8,4)
#data
hist(y,breaks=x,freq=FALSE,main='',xlab='velocity',ylab='density')
#posterior preditive dbn conditioned on "i" components
for (i in c(1,3,5)) {
  lines(x,cden[i,],lwd=2,col=i+1)
}
#posterior predictive mean
lines(x,mden,lwd=2,col=1)
#dev.off()

#shows the label switching if we dont sort mu at plot
#pdf('DPmixtureLS.pdf',8,4)
plot(c(1,Nsamp),c(0,0),type='n',ylim=c(0,50),ylab='mu',xlab='MCMC step (x 100)')
for (i in 1:Nsamp) {
  points(rep(i,PL[i,1]),(TH[[i]]$mu),pch='.',cex=2,col=1:PL[i,1])
}
#dev.off()

###
#debugging relics

# Nsamp=max(which(!is.na(PL[,1])))
# THs=TH; CLs=CL; PLs=PL
# TH<-TH[1:Nsamp]
# CL<-CL[1:Nsamp,]
# PL<-PL[1:Nsamp,]
i=42
mv=sapply(100:Nsamp,function(x) TH[[x]]$mu[CL[x,i]])
sv=sapply(100:Nsamp,function(x) TH[[x]]$sigma[CL[x,i]])
plot(sv,mv,main=paste('mu, sigma for the cluster y[i] with i=',as.character(i),' is in',sep=""),
                      xlab='mu',ylab='sigma')
abline(h=y[i],col=2,lwd=2)
