set.seed(1993,kind=NULL)
library(VIM)

## model setting ##
nsim=1000 #number of simulation runs
N=10000 #Population size
n=200   #sample size
x=1+rexp(N)
e=rnorm(N)
y=4+x+e ##
muy=mean(y)
mux=mean(x)

## missing data ##
eta_i<-cbind(rep(1,N),x) %*% t(t(c(2,0.1)))
expit<-function(x){
  return(exp(x)/(1+exp(x)))
}
R<-unlist(lapply(eta_i,FUN = function(x){
  return(rbinom(1,1,expit(x)))
}))
Y<-R*y
Y[Y==0]<-NA
muY<-mean(Y,na.rm = T)
muY
muy
z=cbind(x,Y)
miss.vals = which(apply(z,1,function(z){return(sum(is.na(z)))})!=0) 
c(nrow(z), nrow(z) - length(miss.vals), length(miss.vals)/nrow(z) )

## take a sample by PPS ##
z=x/sum(x)
syspps=function(x,n){
  N=length(x)
  U=sample(N,N)
  xx=x[U]
  z=rep(0,N)
  for(i in 1:N) z[i]=n*sum(xx[1:i])/sum(x)
  r=runif(1)
  s=numeric()
  for(i in 1:N){
    if(z[i]>=r){
      s=c(s,U[i])
      r=r+1
    }
  }
  return(s[order(s)])
}
s=syspps(z,n)
zs=z[s]
Ys=Y[s]
xs=x[s]
pis=n*zs
zz=cbind(xs,Ys)
miss.vals = which(apply(zz,1,function(zz){return(sum(is.na(zz)))})!=0) 
c(nrow(zz), nrow(zz) - length(miss.vals), length(miss.vals)/nrow(zz) )
si.nni=kNN(zz,k=1)
si.nni=subset(si.nni,select=xs:Ys)
plot(zz,pch=20,main="Figure 1")
points(si.nni[miss.vals,],pch=4,col=c(4))
ys=si.nni$Ys
muyhat_ht=sum(ys/pis)/N
muxhat_ht=sum(xs/pis)/N
muyhat_gr=(muyhat_ht/muxhat_ht)*mux
s2x=var(xs)
sxy=cov(xs,ys)
betahat=sxy/s2x
muyhat_greg = muyhat_ht+betahat*(mux-muxhat_ht)
muyhat_gr
muyhat_greg

##variance estimation for generalized ratio estimator ##
r=muyhat_ht/muxhat_ht
li=ys-r*xs
Ri=li/zs
s2=sum((Ri-mean(Ri))^2)/(n-1)
VL_gr=s2/n/N^2

VBoot_gr=function(xs,ys,zs,mux,n,B)
{
  V=rep(0,B)
  for(i in 1:B)
  {
    bsam=sample(n,n,replace=T)
    bx=xs[bsam]
    by=ys[bsam]
    bz=zs[bsam]
    bpi=n*bz
    muyhat_ht=sum(by/bpi)/N
    muxhat_ht=sum(bx/bpi)/N
    V[i]=(muyhat_ht/muxhat_ht)*mux
  }
  VB_gr=(B-1)*var(V)/B
  return(VB_gr)
}
VB_gr=VBoot_gr(xs,ys,zs,mux,n,1000)

VJack_gr=function(xs,ys,pis,mux,n)
{
  V=rep(0,n)
  for (i in 1:n)
  {
    muyhat_ht=sum(ys[-i]/pis[-i])/N
    muxhat_ht=sum(xs[-i]/pis[-i])/N
    V[i]=(muyhat_ht/muxhat_ht)*mux
  }
  VJ_gr=((n-1)^2/n)*var(V)
  return(VJ_gr)
}
VJ_gr=VJack_gr(xs,ys,pis,mux,n)

VL_gr
VB_gr
VJ_gr

## CIs for generalized ratio estimator ##
mu1=muyhat_gr-1.96*sqrt(VL_gr)
mu2=muyhat_gr+1.96*sqrt(VL_gr)
mu3=muyhat_gr-1.96*sqrt(VB_gr)
mu4=muyhat_gr+1.96*sqrt(VB_gr)
mu5=muyhat_gr-1.96*sqrt(VJ_gr)
mu6=muyhat_gr+1.96*sqrt(VJ_gr)
CI_gr1=c(mu1,mu2)
CI_gr2=c(mu3,mu4)
CI_gr3=c(mu5,mu6)
CI_gr1
CI_gr2
CI_gr3

## variance estimation for generalized regression estimator ##
ei=ys-betahat*xs
Ri=ei/zs
s2=sum((Ri-mean(Ri))^2)/(n-1)
VL_greg=s2/n/N^2

VBoot_greg=function(xs,ys,zs,mux,n,B)
{
  V=rep(0,B)
  for(i in 1:B)
  {
    bsam=sample(n,n,replace=T)
    bx=xs[bsam]
    by=ys[bsam]
    bz=zs[bsam]
    bpi=n*bz
    muyhat_ht=sum(by/bpi)/N
    muxhat_ht=sum(bx/bpi)/N
    sxy=cov(bx,by)
    s2x=var(bx)
    betahat=sxy/s2x
    V[i]=muyhat_ht+betahat*(mux-muxhat_ht)
  }
  VB_greg=(B-1)*var(V)/B
  return(VB_greg)
}
VB_greg=VBoot_greg(xs,ys,zs,mux,n,1000)

VJack_greg=function(xs,ys,pis,mux,n)
{
  V=rep(0,n)
  for (i in 1:n)
  {
    sxyj=cov(xs[-i],ys[-i])
    s2xj=var(xs[-i])
    betahatj=sxyj/s2xj
    muyhat_ht=sum(ys[-i]/pis[-i])/N
    muxhat_ht=sum(xs[-i]/pis[-i])/N
    V[i]=muyhat_ht+betahatj*(mux-muxhat_ht)
  }
  VJ_greg=((n-1)^2/n)*var(V)
  return(VJ_greg)
}
VJ_greg=VJack_greg(xs,ys,pis,mux,n)

VL_greg
VB_greg
VJ_greg

## CIs for generalized regression estimator ##
mu1=muyhat_greg-1.96*sqrt(VL_greg)
mu2=muyhat_greg+1.96*sqrt(VL_greg)
mu3=muyhat_greg-1.96*sqrt(VB_greg)
mu4=muyhat_greg+1.96*sqrt(VB_greg)
mu5=muyhat_greg-1.96*sqrt(VJ_greg)
mu6=muyhat_greg+1.96*sqrt(VJ_greg)
CI_greg1=c(mu1,mu2)
CI_greg2=c(mu3,mu4)
CI_greg3=c(mu5,mu6)
CI_greg1
CI_greg2
CI_greg3

## bias, MSE, length, cp for generalized ratio ##
bias=c(0)
mse=c(0)
cp=matrix(0,3,3)
len=c(0,0,0)
for(m in 1:nsim){
  s=syspps(z,n)
  zs=z[s]
  Ys=Y[s]
  xs=x[s]
  pis=n*zs
  zz=cbind(xs,Ys)
  si.nni=kNN(zz,k=1)
  si.nni=subset(si.nni,select=xs:Ys)
  ys=si.nni$Ys
  muyhat_ht=sum(ys/pis)/N
  muxhat_ht=sum(xs/pis)/N
  muyhat_gr=(muyhat_ht/muxhat_ht)*mux
  r=muyhat_ht/muxhat_ht
  li=ys-r*xs
  Ri=li/zs
  s2=sum((Ri-mean(Ri))^2)/(n-1)
  VL_gr=s2/n/N^2
  VB_gr=VBoot_gr(xs,ys,zs,mux,n,1000)
  VJ_gr=VJack_gr(xs,ys,pis,mux,n)
  mu1=muyhat_gr-1.96*sqrt(VL_gr)
  mu2=muyhat_gr+1.96*sqrt(VL_gr)
  mu3=muyhat_gr-1.96*sqrt(VB_gr)
  mu4=muyhat_gr+1.96*sqrt(VB_gr)
  mu5=muyhat_gr-1.96*sqrt(VJ_gr)
  mu6=muyhat_gr+1.96*sqrt(VJ_gr)
  bias[1]=bias[1]+muyhat_gr-muy
  mse[1]=mse[1]+(muyhat_gr-muy)^2
  len[1]=len[1]+mu2-mu1
  len[2]=len[2]+mu4-mu3
  len[3]=len[3]+mu6-mu5
  cp[1,1]=cp[1,1]+(muy<=mu1)
  cp[1,2]=cp[1,2]+(muy>mu1)*(muy<mu2)
  cp[1,3]=cp[1,3]+(muy>=mu2)
  cp[2,1]=cp[2,1]+(muy<=mu3)
  cp[2,2]=cp[2,2]+(muy>mu3)*(muy<mu4)
  cp[2,3]=cp[2,3]+(muy>=mu4)
  cp[3,1]=cp[3,1]+(muy<=mu5)
  cp[3,2]=cp[3,2]+(muy>mu5)*(muy<mu6)
  cp[3,3]=cp[3,3]+(muy>=mu6)
}
bias_gr=bias/(muy*nsim)
mse_gr=mse/nsim
len_gr=len/nsim
cp_gr=cp/nsim
bias_gr
mse_gr
len_gr
cp_gr

## bias, MSE, length, cp for generalized regression ##
bias=c(0)
mse=c(0)
cp=matrix(0,3,3)
len=c(0,0,0)
for(m in 1:nsim){
  s=syspps(z,n)
  zs=z[s]
  Ys=Y[s]
  xs=x[s]
  pis=n*zs
  zz=cbind(xs,Ys)
  si.nni=kNN(zz,k=1)
  si.nni=subset(si.nni,select=xs:Ys)
  ys=si.nni$Ys
  muyhat_ht=sum(ys/pis)/N
  muxhat_ht=sum(xs/pis)/N
  s2x=var(xs)
  sxy=cov(xs,ys)
  betahat=sxy/s2x
  muyhat_greg=muyhat_ht+betahat*(mux-muxhat_ht)
  ei=ys-betahat*xs
  Ri=ei/zs
  s2=sum((Ri-mean(Ri))^2)/(n-1)
  VL_greg=s2/n/N^2
  VB_greg=VBoot_greg(xs,ys,zs,mux,n,1000)
  VJ_greg=VJack_greg(xs,ys,pis,mux,n)
  mu1=muyhat_greg-1.96*sqrt(VL_greg)
  mu2=muyhat_greg+1.96*sqrt(VL_greg)
  mu3=muyhat_greg-1.96*sqrt(VB_greg)
  mu4=muyhat_greg+1.96*sqrt(VB_greg)
  mu5=muyhat_greg-1.96*sqrt(VJ_greg)
  mu6=muyhat_greg+1.96*sqrt(VJ_greg)
  bias[1]=bias[1]+muyhat_greg-muy
  mse[1]=mse[1]+(muyhat_greg-muy)^2
  len[1]=len[1]+mu2-mu1
  len[2]=len[2]+mu4-mu3
  len[3]=len[3]+mu6-mu5
  cp[1,1]=cp[1,1]+(muy<=mu1)
  cp[1,2]=cp[1,2]+(muy>mu1)*(muy<mu2)
  cp[1,3]=cp[1,3]+(muy>=mu2)
  cp[2,1]=cp[2,1]+(muy<=mu3)
  cp[2,2]=cp[2,2]+(muy>mu3)*(muy<mu4)
  cp[2,3]=cp[2,3]+(muy>=mu4)
  cp[3,1]=cp[3,1]+(muy<=mu5)
  cp[3,2]=cp[3,2]+(muy>mu5)*(muy<mu6)
  cp[3,3]=cp[3,3]+(muy>=mu6)
}
bias_greg=bias/(muy*nsim)
mse_greg=mse/nsim
len_greg=len/nsim
cp_greg=cp/nsim
bias_greg
mse_greg
len_greg
cp_greg