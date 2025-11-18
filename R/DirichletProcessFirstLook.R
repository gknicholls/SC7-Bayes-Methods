#Bayes Methods - CRP/Dirichlet Process simple stuff

#example simulate DP(alpha,H=normal(mu)*gamma(sigma))

#DP parameter
alpha=1
#(theta_1,...,theta_n)
n=1e2

DP.sim<-function(alpha,n) {
  #one person at table 1
  S=1
  K=1
  mu=rnorm(1,mean=20,sd=10)
  sg=rgamma(1,shape=1.5,rate=0.5)
  #bring in the customers 1 by 1
  for (i in 2:n) {
    #n_k=nS[k], table() counts the number at each table k=1..K
    nS=table(S)
    #the probability vector is propto [n_1,...,n_K, alpha]
    pr=c(nS,alpha); pr=pr/sum(pr)
    #pick a table
    k=sample(1:(K+1),1,prob=pr)
    if (k==(K+1)) { #new cluster
      K=K+1
      S=c(S,K)
      mu=c(mu,rnorm(1,mean=20,sd=10))
      sg=c(sg,rgamma(1,shape=1.5,rate=0.5))
    } else {
      S=c(S,k)
    }
  }
  return(list(S=S,mu=mu,sg=sg))
}

#examples with different alpha showing the discrete DP dbn on mu
#pdf('DPalpha.pdf',6,6)
ro=2; co=2; ex=ro*co
alpha=10^seq(-1,2,length.out=ex)
par(mfrow=c(ro,co))

for (i in 1:ex) {
  #simulate DP
  O=DP.sim(alpha[i],n)

  #display discrete character of DP
  lf=0.1+0.9*(i%%co==1); up=0.1+0.9*(i<=ro); rg=0.1+0.9*(i%%co==0); dn=0.1+0.9*(i>(ex-co));
  xt=ifelse(i>(ex-co),"s","n")
  yt=ifelse(i%%co==1,"s","n")
    par(mai=0.75*c(dn,lf,up,rg))
    plot(O$mu,O$mu^0,xlab='mu',ylab='w',ylim=c(0,1),xlim=c(0,40),
       type='n',xaxt=xt,yaxt=yt,cex.lab=1.4)
  text(30,0.6,paste("a=",as.character(signif(alpha[i],2))))
  w=table(O$S); w=w/sum(w)
  for (k in 1:length(O$mu)) { lines( c(1,1)*O$mu[k],c(0,w[k]),lwd=2)}
}
#dev.off()

#check we agree on the mean number of clusters - n=82, alpha=1
J=1000; K=rep(NA,J); alpha=1; n=82 #this will be n for the main example application
for (j in 1:J) {
  O=DP.sim(alpha,n)
  K[j]=max(O$S)
}
mean(K)
EK=sum(alpha/(alpha+0:(n-1))); EK
#reasonable agreement
