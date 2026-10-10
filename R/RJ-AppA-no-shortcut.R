#experiment exploring RJ giving example of case where shortcut doesn't work
#See Section on RJ-MCMC appendix A

rm(list=ls())
library(rgl)

#map for adding a dimension
#if the current state is m=1, x=(x_1,x_2) and we draw u=(u_1,u_2)~N(0,I_2) 
#then set (x',u')=(x'_1,x'_2,x'_3,u') using the linear map (x',u')=psi*(x,u)
psi=matrix(c(1,0,1/2,0,0,1,1/2,0,1,1,0,-1,1,0,1,0),4,4)

#map for deleting a dimension
#if the current state is m'=2, x'=(x'_1,x'_2,x'_3) then draw u'~N(0,1) 
#and set (x,u)=psi_inv*(x',u')
psi_inv=solve(psi)

d=log(det(psi)) #J_\psi

gq<-function(x){
  #make a proposal and pass back the proposal and the values of g_{1,2} and g_{2,1}
  if (length(x)==2) {           # if x in R^2 then propose x' in R^3
    u=rnorm(2)                  # u in R^2
    q=sum(dnorm(u,log=TRUE))    # g_{1,2}(u)
    dt=d                        # Jacobian for adding
    xu=psi%*%matrix(c(x,u),4,1) # (x',u')
    xp=xu[1:3]                  # x'=(x'_1,x'_2,x'_3)
    up=xu[4]                    # u' in R^1
    qp=dnorm(up,log=TRUE)       # g_{2,1}(u')
  } else {                      # if x in R^3 then propose x' in R^2
    if (length(x)!=3) stop('wierdness')
    u=rnorm(1)                  # u in R^1
    q=dnorm(u,log=TRUE)         # g_{2,1}(u)
    dt=-d                       # Jacobian for deleting
    xu=psi_inv%*%matrix(c(x,u),4,1) 
    xp=xu[1:2]                  # x'=(x'_1,x'_2)
    up=xu[3:4]                  # u'=(u'_1,u'_2)
    qp=sum(dnorm(up,log=TRUE))  # g_{1,2}(u')
  }
  return(list(xp=xp,q=q,qp=qp,dt=dt)) #all the stuff we need to evaluate the acc. prob.
}

p<-function(x,w=1/2){
  #the target dbn is p(x,m)=p(x|m)p(m) with p(x|m)=N(x;0,I_{m+1}) and p(m)=1/2
  wm=ifelse(length(x)==2,w,1-w)
  pr=log(wm)+sum(dnorm(x,log=TRUE))
  return(pr)
}

K=100000 #run length
SS=10    #collect a sample every SS steps
M=rep(NA,K/SS) #save m (can be 1 or 2, but use TRUE/FALSE for m=1)
X=matrix(0,K/SS,3) #save x can be 2 or 3 components
x=c(0,0,0) #start state

for (k in 1:K) {
  xq=gq(x) #propose
  MHR=xq$dt+p(xq$xp)+xq$qp-p(x)-xq$q #Hasting ratio 
  if (log(runif(1))<MHR) x=xq$xp #if accept then overwrite x with x'
  if (k%%SS==0) { #every SS steps
    M[k/SS]=(m_is_1<-length(x)==2) #M[k]=TRUE or FALSE for dim(x)=2 or 3 so m=1 or 2
    if (m_is_1) {
      X[k/SS,1:2]=x
    } else {
      X[k/SS,]=x
    }
  }
  #ensure irreducibility - trivial target prod_{i=1}^m N(x_1;0,1)
  x[1]=rnorm(1)
}

#show the 3D structure with m=1 samples in a plane (in cyan)
#and the m=2 samples in 3D (in pink) obvious normal structure
mycolors <- c('hotpink3', 'darkcyan')
plot3d( 
  x=X[,1], y=X[,2], z=X[,3], 
  col = mycolors[M+1], 
  type = 's', 
  radius = .05,
  xlab="X1", ylab="X2", zlab="X3")

#show the RJ update is not irreducible by itself
K=100
X=matrix(0,K,3)
x=c(0,0)
for (k in 1:K) {
  #make two steps, m=1->m=2->m=1 (but using different rv each time)
  u=rnorm(2); v=c(x,u); vp=psi%*%v; xp=vp[1:3]
  #X[k,]=xp
  up=rnorm(1); vp=c(xp,up); v=psi_inv%*%vp; x=v[1:2]
  X[k,1:2]=x
}
X[,1]/X[,2]
#so after two steps m=1->m=2->m=1 from x to y,
#x_1-3x_2=y_1-3y_2 (so this is invariant => not irreducible)
#here x_1=x_2=0 so y_1/y_2=3 as seen in X[,1]/X[,2]

#ahat=summary(lm(X[,1]~X[,2]))$coefficients[2,1]; a=3; ahat-a
