
#Estimating Bayes Factors/prior simulation example

############################################################################
#bayes factor example - see also code for radiocarbon dating case study

library(MCMCpack)

##################################################################################
#example - Bayesian model selection of link function
#Challenger O-ring failures (fail=1, OK=0)
fail=c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0)
temp=scale(c(53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81))
plot(temp,fail)

#this is classical so MCMC code targeting posterior is available in R library "MCMCpack"
#results in lectures use run length N x 10.

N=1000
prior.std.dev=3
v=1/prior.std.dev^2 #precision v=1/9 is sigma^2=3 in prior
b.p=MCMCprobit(fail~temp,b0=c(0,0),B0=c(v,v),seed=4,mcmc=N,marginal.likelihood='Laplace')
b.l=MCMClogit(fail~temp,b0=c(0,0),B0=c(v,v),seed=5,mcmc=N,marginal.likelihood='Laplace')

plot(b.p)

logit<-function(beta,x=temp) {eta=beta[1]+beta[2]*x; exp(eta)/(1+exp(eta))}
probit<-function(beta,x=temp) {eta=beta[1]+beta[2]*x; pnorm(eta)}

lkd<-function(beta,logit=T) {
  if (logit) {pb=logit(beta)} else {pb=probit(beta)}
  return(sum(dbinom(fail,1,pb,log=T)))
}

l.pr=l.lo=pr.lo=pr.pr=l.lo.pr=l.pr.lo=rep(NA,N)
for (k in 1:N) {
  l.lo[k]=lkd(b.l[k,],logit=T); 
  l.pr[k]=lkd(b.p[k,],logit=F);
  l.lo.pr[k]=lkd(b.p[k,],logit=T); 
  l.pr.lo[k]=lkd(b.l[k,],logit=F);
  pr.lo[k]=sum(dnorm(b.l[k,],0,9,log=T))
  pr.pr[k]=sum(dnorm(b.p[k,],0,9,log=T))
}

#harmonic
ml.lo=1/mean(1/exp(l.lo))
ml.pr=1/mean(1/exp(l.pr))
ml.lo/ml.pr

#Bridge estimator just using h=1 - simple, performs well
mean(exp(l.lo.pr)*pr.pr)/mean(exp(l.pr.lo)*pr.lo)
#and BE using h=1/sqrt(p1*p2) - recommended by Meng & Wong as good place to start
mean(exp(l.lo.pr/2-l.pr/2))/mean(exp(l.pr.lo/2-l.lo/2))

#builtin function uses a Laplace estimator 
exp(attr(b.l,"logmarglike"))/exp(attr(b.p,"logmarglike"))

#the naive estimator is usually hopeless
b.prior=matrix(rnorm(2*N,0,3),N,2)
l.pr.n=l.lo.n=rep(NA,N)
for (k in 1:N) {l.lo.n[k]=lkd(b.prior[k,],logit=T); l.pr.n[k]=lkd(b.prior[k,],logit=F)}
mean(exp(l.lo.n))/mean(exp(l.pr.n))

#I did a run with 1E+6 samples got Bridge1&2/Laplace/Harmonic/Naive = 2.76/2.67/1.92/2.76

########################################################
#Using simulation in prior elicitation Challenger/logistic example
#Now discussed in Chapter 1

logit<-function(beta,x) {eta=beta[1]+beta[2]*x; exp(eta)/(1+exp(eta))}

#pdf('logitpriorsample.pdf',12,6)
x=seq(min(temp),max(temp),length.out=100)
N=100
plot.std.dev=c(prior.std.dev,10,1000)
par(mfrow=c(1,length(plot.std.dev)),mai=c(1,0.75,1,0.1))
for (k in 1:length(plot.std.dev)) {
  b.prior=matrix(rnorm(2*N,0,plot.std.dev[k]),N,2) #simulate prior
  y=rep(NA,N); for (j in 1:N) {y[j]=logit(b.prior[j,2],x[j])}
  plot(x,y,type='l',xlim=c(min(temp),max(temp)),ylim=c(0,1),
       xlab='scale(temp)',ylab='failure probability',
     main=paste('Prior std. dev. ',plot.std.dev[k]),cex.main=1.5,cex.lab=2)
  for (j in 1:N) {if (b.prior[j,2]>0) {y=logit(b.prior[j,],x); lines(x,y);}} 
  #plot failure probability functions
}
#dev.off()
