library(expm)
library(mapfit)
library(MCMCpack)
library(pracma)
library(FAdist)
library(Rsolnp)
library(pracma)
automobile<-read.table(file = "C:/Users/e0265296/Desktop/automobile.txt",header = TRUE)
##change the time scale from month to year
lags<-automobile[,2]
lags<-lags/12
ts<-automobile[,3]
ts<-ts/12
##warr is the warranty length
warr<-18/12
##censor is the end-of-study date
censor<-54/12
##ss0 is the total number of products under study
ss0<-589
##ncensor is the number of censored products
ncensor<-ss0-length(ts)


mixture_Erlang<-function(mypara){
  pips<-rep(1/m1,m1)
  pi1ps<-rep(1/m2,m2)
  if(m1>1){
    p0ps<--1*diag(rep(mypara[1],m1))
    for(i in 1:(m1-1)){
      p0ps[i,i+1]<-p0ps[i,i]*(-1)
    }
  }else{
    p0ps<--1*matrix(rep(mypara[1],m1))
  }
  nups<-apply(p0ps,MARGIN = 1,sum)*(-1)
  
  if(m2>1){
    p1ps<--1*diag(rep(mypara[2],m2))
    for(i in 1:(m2-1)){
      p1ps[i,i+1]<-p1ps[i,i]*(-1)
    }
  }else{
    p1ps<--1*matrix(rep(mypara[2],m2))
  }
  nu1ps<-apply(p1ps,MARGIN = 1,sum)*(-1)
  y1<-ts
  y2<-lags
  es=rep(1,m1+m2)
  D0=cbind(rbind(p0ps,matrix(0,m2,m1)),rbind(nups%*%t(pi1ps),p1ps))
  mm1<-m1+m2
  prob<-sum(c(pips,rep(0,m2))*expAtv(D0,es,warr)$eAtv)
  tempmat<-diag(mm1)
  for(k in (m1+1):mm1){
    prob<-prob-sum(c(pips,rep(0,m2))*expAtv(D0,tempmat[k,],warr)$eAtv)*(1-sum(tempmat[k,]*expAtv(D0,es,censor-warr)$eAtv))
  }
  res<-sum(dph(y1,ph=ph(alpha=pips,Q=p0ps,xi=nups),log=T))+sum(dph(y2,ph=ph(alpha=pi1ps,Q=p1ps,xi=nu1ps),log=T))+ncensor*log(prob)
  return(-res)
}

myfunction1<-function(x,alpha=numeric(0),Q=numeric(0),xi=numeric(0)){
  if(length(x)>0){
    res<-rep(0,length(x))
    for(i in 1:length(x)){
      res[i]<-sum(alpha*expAtv(Q,xi,x[i])$eAtv)
    }
  }else{
    res<-sum(alpha*expAtv(Q,xi,x)$eAtv)
  }
  return(res)
}

intfunction1<-function(x,alpha=numeric(0),Q=numeric(0),xi=numeric(0)){
  if(length(x)>0){
    res<-rep(0,length(x))
    for(i in 1:length(x)){
      res[i]<-integrate(f=myfunction1,lower = x[i],upper = Inf,alpha=alpha,Q=Q,xi=xi,rel.tol = .Machine$double.eps^0.5, abs.tol = .Machine$double.eps^0.5,subdivisions = 200L)$value
    }
  }else{
    res<-integrate(f=myfunction1,lower = x,upper = Inf,alpha=alpha,Q=Q,xi=xi,rel.tol = .Machine$double.eps^0.5, abs.tol = .Machine$double.eps^0.5,subdivisions = 200L)$value
  }
  return(res)
}

emcore<-function(y,pips,nups,p0ps,mups,Rmax,
                 y1,pi1ps,nu1ps,p1ps,mu1ps,censor,warr,ncensor,CRmax){
  ss=length(y);m=length(pips)
  if(m==1){
    p0ps<-matrix(p0ps)
  }
  y=sort(y);dy=c(y[1],diff(y));re_p0=p0ps/mups+diag(m)
  fks=matrix(expAtv(t(p0ps),pips,dy[1])$eAtv,m,ss)
  bks=matrix(expAtv(p0ps,nups,dy[1])$eAtv,m,ss)
  for(i in 2:ss){
    fks[,i]=expAtv(t(p0ps),fks[,i-1],dy[i])$eAtv
    bks[,i]=expAtv(p0ps,bks[,i-1],dy[i])$eAtv
  }
  cks=matrix(pips/sum(pips*bks[,ss]),ss,m,byrow=T)
  for(i in (ss-1):1){
    cks[i,]=expAtv(t(p0ps),cks[i+1,],dy[i+1])$eAtv+pips/sum(pips*bks[,i])
  }
  bks=cbind(nups,bks)
  Hs=matrix(0,m,m)
  for(i in 1:ss){
    Hks=matrix(0,m,m)
    alps=matrix(bks[,i],m,Rmax+1)
    for(j in 2:(Rmax+1)) alps[,j]=re_p0%*%alps[,j-1]
    bets=matrix(dpois(Rmax+1,mups*dy[i])*cks[i,],Rmax+1,m,byrow=T)
    for(j in Rmax:1){
      bets[j,]=bets[j+1,]%*%re_p0+dpois(j,mups*dy[i])*cks[i,]
    }
    for(j in 1:(Rmax+1)){
      Hks=Hks+alps[,j]%*%t(bets[j,])
    }
    Hs=Hs+Hks/mups
  }
  Hs=t(Hs)
  Bs=rep(0,m);As=rep(0,m)
  for(i in 1:ss){
    As=As+fks[,i]*nups/sum(pips*bks[,i+1])
    Bs=Bs+pips*bks[,i+1]/sum(pips*bks[,i+1])
  }
  m1=length(pi1ps)
  if(m1==1){
    p1ps<-matrix(p1ps)
  }
  y1=sort(y1);dy1=c(y1[1],diff(y1));re_p01=p1ps/mu1ps+diag(m1)
  fks1=matrix(expAtv(t(p1ps),pi1ps,dy1[1])$eAtv,m1,ss)
  bks1=matrix(expAtv(p1ps,nu1ps,dy1[1])$eAtv,m1,ss)
  for(i in 2:ss){
    fks1[,i]=expAtv(t(p1ps),fks1[,i-1],dy1[i])$eAtv
    bks1[,i]=expAtv(p1ps,bks1[,i-1],dy1[i])$eAtv
  }
  cks1=matrix(pi1ps/sum(pi1ps*bks1[,ss]),ss,m1,byrow=T)
  for(i in (ss-1):1){
    cks1[i,]=expAtv(t(p1ps),cks1[i+1,],dy1[i+1])$eAtv+pi1ps/sum(pi1ps*bks1[,i])
  }
  bks1=cbind(nu1ps,bks1)
  Hs1=matrix(0,m1,m1)
  for(i in 1:ss){
    Hks1=matrix(0,m1,m1)
    alps1=matrix(bks1[,i],m1,Rmax+1)
    for(j in 2:(Rmax+1)) alps1[,j]=re_p01%*%alps1[,j-1]
    bets1=matrix(dpois(Rmax+1,mu1ps*dy1[i])*cks1[i,],Rmax+1,m1,byrow=T)
    for(j in Rmax:1){
      bets1[j,]=bets1[j+1,]%*%re_p01+dpois(j,mu1ps*dy1[i])*cks1[i,]
    }
    for(j in 1:(Rmax+1)){
      Hks1=Hks1+alps1[,j]%*%t(bets1[j,])
    }
    Hs1=Hs1+Hks1/mu1ps
  }
  Hs1=t(Hs1)
  Bs1=rep(0,m1);As1=rep(0,m1)
  for(i in 1:ss){
    As1=As1+fks1[,i]*nu1ps/sum(pi1ps*bks1[,i+1])
    Bs1=Bs1+pi1ps*bks1[,i+1]/sum(pi1ps*bks1[,i+1])
  }
  ##########################################################################################for missing data 
  hc=pips*expAtv(p0ps,rep(1,m),warr)$eAtv;sumhc=sum(hc[1:m])
  CHks=matrix(0,m,m)
  calps=matrix(rep(1,m),m,CRmax+1)
  for(j in 2:(CRmax+1)) calps[,j]=re_p0%*%calps[,j-1]
  cbets=matrix(dpois(CRmax+1,mups*warr)*pips,CRmax+1,m,byrow=T)
  for(j in CRmax:1){
    cbets[j,]=cbets[j+1,]%*%re_p0+dpois(j,mups*warr)*pips
  }
  for(j in 1:(CRmax+1)){
    CHks=CHks+calps[,j]%*%t(cbets[j,])
  }
  CHks=t(CHks)/sumhc/mups
  
  mm1=m+m1;es=rep(1,mm1)
  mmups=max(mups,mu1ps)
  D0=cbind(rbind(p0ps,matrix(0,m1,m)),rbind(nups%*%t(pi1ps),p1ps))
  
  prob<-0
  tempmat<-diag(mm1)
  tempmat1<-diag(m)
  for(k in 1:m){
    prob<-prob+sum(c(pips,rep(0,m1))*expAtv(D0,tempmat[k,],warr)$eAtv)*sum(tempmat[k,]*expAtv(D0,rep(1,mm1),censor-warr)$eAtv)
  }
  p1<-prob
  p2<-sum(c(pips,rep(0,m1))*expAtv(D0,es,censor)$eAtv)
  p3<-sumhc
  q1<-1-p2-p3+p1
  q2<-1-q1
  q3<-(p3-p1)/(1-p2)
  q4<-q1/(1-p2)
  q5<-p2
  q6<-1-p2
  q7<-p1/p3
  q8<-(p3-p1)/p3
  
  re_D0=D0/mmups+diag(mm1)
  
  CHks2=matrix(0,mm1,mm1)
  calps=matrix(es,mm1,CRmax+1)
  for(j in 2:(CRmax+1)) calps[,j]=re_D0%*%calps[,j-1]
  cbets=matrix(dpois(CRmax+1,mmups*censor)*c(pips,rep(0,m1)),CRmax+1,mm1,byrow=T)
  for(j in CRmax:1){
    cbets[j,]=cbets[j+1,]%*%re_D0+dpois(j,mmups*censor)*c(pips,rep(0,m1))
  }
  for(j in 1:(CRmax+1)){
    CHks2=CHks2+calps[,j]%*%t(cbets[j,])
  }
  CHks2=t(CHks2)/p2/mmups
  C2Inf<-rep(0,mm1)
  for(i in 1:mm1){
    C2Inf[i]<-intfunction1(x=censor,alpha = c(pips,rep(0,m1)),Q=D0,xi=tempmat[i,])
  }
  C2Inf<-C2Inf/p2
  CHks2<-CHks2+C2Inf%*%t(es)
  myi0_XplusT<-C2Inf[(m+1):mm1]*nu1ps
  ################till now we have computed E[*|X+T>C], which is denoted as CHks2. (not multiplied by ncensor)
  
  CHks31=matrix(0,mm1,mm1)
  for(k in 1:m){
    temp_CHks=matrix(0,mm1,mm1)
    calps=matrix(tempmat[k,],mm1,CRmax+1)
    for(j in 2:(CRmax+1)) calps[,j]=re_D0%*%calps[,j-1]
    cbets=matrix(dpois(CRmax+1,mmups*(warr))*c(pips,rep(0,m1)),CRmax+1,mm1,byrow=T)
    for(j in CRmax:1){
      cbets[j,]=cbets[j+1,]%*%re_D0+dpois(j,mmups*(warr))*c(pips,rep(0,m1))
    }
    for(j in 1:(CRmax+1)){
      temp_CHks=temp_CHks+calps[,j]%*%t(cbets[j,])
    }
    temp_CHks=t(temp_CHks)/mmups
    CHks31<-CHks31+temp_CHks*sum(tempmat[k,]*expAtv(D0,es,censor-warr)$eAtv)
  }
  CHks31<-CHks31/prob
  CHks3=matrix(0,mm1,mm1)
  for(k in 1:m){
    temp_CHks=matrix(0,mm1,mm1)
    calps=matrix(es,mm1,CRmax+1)
    for(j in 2:(CRmax+1)) calps[,j]=re_D0%*%calps[,j-1]
    cbets=matrix(dpois(CRmax+1,mmups*(censor-warr))*tempmat[k,],CRmax+1,mm1,byrow=T)
    for(j in CRmax:1){
      cbets[j,]=cbets[j+1,]%*%re_D0+dpois(j,mmups*(censor-warr))*tempmat[k,]
    }
    for(j in 1:(CRmax+1)){
      temp_CHks=temp_CHks+calps[,j]%*%t(cbets[j,])
    }
    temp_CHks=t(temp_CHks)/mmups
    CHks3<-CHks3+temp_CHks*sum(pips*expAtv(p0ps,tempmat1[k,],warr)$eAtv)
  }
  CHks3<-CHks3/prob
  CHks3[1:m,1:m]<-CHks31[1:m,1:m]+CHks3[1:m,1:m]
  newC2Inf<-rep(0,mm1)
  for(i in 1:mm1){
    for(j in 1:mm1){
      int_ij<-intfunction1(x=0,alpha = tempmat[j,],Q=D0,xi=tempmat[i,])
      multiplier<-0
      for(k in 1:m){
        multiplier<-multiplier+sum(c(pips,rep(0,m1))*expAtv(D0,tempmat[k,],warr)$eAtv)*sum(tempmat[k,]*expAtv(D0,tempmat[j,],censor-warr)$eAtv)
      }
      newC2Inf[i]<-newC2Inf[i]+multiplier*int_ij
    }
  }
  newC2Inf<-newC2Inf/prob
  CHks3<-CHks3+newC2Inf%*%t(es)
  myi0_both<-newC2Inf[(m+1):mm1]*nu1ps
  ################till now we have computed E[*|X+T>C,T>tau], which is denoted as CHks3. (not multiplied by ncensor)
  
  tau2Inf<-rep(0,mm1)
  for(i in 1:mm1){
    for(k in 1:m){
      res<-intfunction1(x=0,alpha = tempmat[k,],Q=D0,xi=tempmat[i,])
      tau2Inf[i]<-tau2Inf[i]+res*sum(c(pips,rep(0,m1))*expAtv(D0,tempmat[k,],warr)$eAtv)
    }
  }
  tau2Inf<-tau2Inf/p3
  CHks4<-tau2Inf%*%t(es)
  CHks4[1:m,1:m]<-CHks+CHks4[1:m,1:m]
  myi0_T<-tau2Inf[(m+1):mm1]*nu1ps
  ################till now we have computed E[*|T>tau], which is denoted as CHks4. (not multiplied by ncensor)
  
  zero2Inf<-rep(0,mm1)
  for(i in 1:mm1){
    zero2Inf[i]<-intfunction1(x=0,alpha = c(pips,rep(0,m1)),Q=D0,xi=tempmat[i,])
  }
  CHks5<-zero2Inf%*%t(es)
  myi0_none<-zero2Inf[(m+1):mm1]*nu1ps
  ################till now we have computed E[*], which is denoted as CHks5. (not multiplied by ncensor)
  
  w1<-(CHks4-CHks3*q7)/q8
  w2<-(CHks5-CHks2*q5)/q6
  w3<-(w2-w1*q3)/q4
  w4<-(CHks5-w3*q1)/q2
  CHks<-w4*ncensor
  
  w1<-(myi0_T-myi0_both*q7)/q8
  w2<-(myi0_none-myi0_XplusT*q5)/q6
  w3<-(w2-w1*q3)/q4
  w4<-(myi0_none-w3*q1)/q2
  myi0<-w4*ncensor
  
  simple2<-(c(pips,rep(0,m1))*expAtv(D0,es,censor)$eAtv)/sum(c(pips,rep(0,m1))*expAtv(D0,es,censor)$eAtv)
  simple2<-simple2[1:m]
  simple3<-rep(0,mm1)
  for(k in 1:m){
    temp_simple<-c(pips,rep(0,m1))*expAtv(D0,tempmat[k,],warr)$eAtv
    simple3<-simple3+temp_simple*sum(tempmat[k,]*expAtv(D0,es,censor-warr)$eAtv)
  }
  simple3<-simple3/prob
  simple3<-simple3[1:m]
  simple4<-(pips*expAtv(p0ps,rep(1,m),warr)$eAtv)/sum(pips*expAtv(p0ps,rep(1,m),warr)$eAtv)
  simple5<-pips
  
  w1<-(simple4-simple3*q7)/q8
  w2<-(simple5-simple2*q5)/q6
  w3<-(w2-w1*q3)/q4
  w4<-(simple5-w3*q1)/q2
  simple<-w4
  
  ####################################M Step
  Ns=Hs+CHks[1:m,1:m]; Zs=Diag(Ns)
  Ns1=Hs1+CHks[(1+m):mm1,(1+m):mm1]; Zs1=Diag(Ns1)
  subchks=nups%*%t(pi1ps)*CHks[1:m,(1+m):mm1]
  sumsub=apply(subchks,1,sum);sumsub1=apply(subchks,2,sum)
  ##estimate lifetime parameters
  newpips=(Bs+ncensor*simple)/(ss+ncensor)
  newnups=(As+sumsub)/Zs
  Ns=Ns*p0ps
  newp0ps=matrix(0,m,m)
  for(j in 1:m){
    newp0ps[j,]=Ns[j,]/Zs[j]
    newp0ps[j,j]=-(sum(newp0ps[j,-j])+newnups[j])
  }
  newmups=max(abs(Diag(newp0ps)))+0.1
  #estimate delay paramters
  newpi1ps=(Bs1+sumsub1)/sum(Bs1+sumsub1)
  newnu1ps=(As1+myi0)/Zs1
  Ns1=Ns1*p1ps
  newp1ps=matrix(0,m1,m1)
  for(j in 1:m1){
    newp1ps[j,]=Ns1[j,]/Zs1[j]
    newp1ps[j,j]=-(sum(newp1ps[j,-j])+newnu1ps[j])
  }
  newmu1ps=max(abs(Diag(newp1ps)))+0.1
  ####results
  list(newpips,newnups,newp0ps,newmups,newpi1ps,newnu1ps,newp1ps,newmu1ps)
}

canonical<-function(newpips,newnups,newp0ps,newpi1ps,newnu1ps,newp1ps){
  m<-length(newpips);m1<-length(newpi1ps)
  if(m>1){
    for(j in 1:(m-1)){
      for(i in 1:(m-j)){
        if(-newp0ps[i,i]>-newp0ps[i+1,i+1]){
          newpips[i]<-newpips[i]+newpips[i+1]*(1-newp0ps[i+1,i+1]/newp0ps[i,i])
          newpips[i+1]<-newpips[i+1]*newp0ps[i+1,i+1]/newp0ps[i,i]
          temp<-newp0ps[i,i]
          newp0ps[i,i]<-newp0ps[i+1,i+1]
          newp0ps[i+1,i+1]<-temp
        }
      }
    }
    for(i in 1:(m-1)){
      newp0ps[i,i+1]<-(-newp0ps[i,i])
    }
    newnups[m]<-(-newp0ps[m,m])
  }
  if(m1>1){
    for(j in 1:(m1-1)){
      for(i in 1:(m1-j)){
        if(-newp1ps[i,i]>-newp1ps[i+1,i+1]){
          newpi1ps[i]<-newpi1ps[i]+newpi1ps[i+1]*(1-newp1ps[i+1,i+1]/newp1ps[i,i])
          newpi1ps[i+1]<-newpi1ps[i+1]*newp1ps[i+1,i+1]/newp1ps[i,i]
          temp<-newp1ps[i,i]
          newp1ps[i,i]<-newp1ps[i+1,i+1]
          newp1ps[i+1,i+1]<-temp
        }
      }
    }
    for(i in 1:(m1-1)){
      newp1ps[i,i+1]<-(-newp1ps[i,i])
    }
    newnu1ps[m1]<-(-newp1ps[m1,m1])
  }
  
  list(newpips,newnups,newp0ps,newpi1ps,newnu1ps,newp1ps)
}

loglike<-function(pips,pi1ps,p0ps,p1ps,nups,nu1ps,y1,y2,warr,censor,ncensor){
  m=length(pips);m1=length(pi1ps);es=rep(1,m+m1)
  if(m==1){
    p0ps<-matrix(p0ps)
  }
  if(m1==1){
    p1ps<-matrix(p1ps)
  }
  D0=cbind(rbind(p0ps,matrix(0,m1,m)),rbind(nups%*%t(pi1ps),p1ps))
  mm1<-m+m1
  prob<-0
  tempmat<-diag(mm1)
  tempmat1<-diag(m)
  for(k in 1:m){
    prob<-prob+sum(c(pips,rep(0,m1))*expAtv(D0,tempmat[k,],warr)$eAtv)*sum(tempmat[k,]*expAtv(D0,rep(1,mm1),censor-warr)$eAtv)
  }
  p1<-prob
  p2<-sum(c(pips,rep(0,m1))*expAtv(D0,es,censor)$eAtv)
  p3<-sum(pips*expAtv(p0ps,rep(1,m),warr)$eAtv)
  q1<-1-p2-p3+p1
  q2<-1-q1
  return(sum(dph(y1,ph=ph(alpha=pips,Q=p0ps,xi=nups),log=T))+
           sum(dph(y2,ph=ph(alpha=pi1ps,Q=p1ps,xi=nu1ps),log=T))+ncensor*log(q2))
}

m1<-4
m2<-5

ptm<-proc.time()

pips<-rep(1/m1,m1)
pi1ps<-rep(1/m2,m2)
###################use the SEM algorithm to select initial value
aaa<-solnp(pars = c(1,1),fun = mixture_Erlang,LB=rep(1e-2,2))$pars
s1<-aaa[1]
s2<-aaa[2]

mixture_Erlang_time_cost<-as.numeric((proc.time()-ptm)[3])

initial_pips<-rep(1/m1,m1)
initial_pi1ps<-rep(1/m2,m2)

if(m1>1){
  initial_p0ps<-diag(rep(s1,m1))*(-1)
  for(i in 1:(m1-1)){
    initial_p0ps[i,i+1]<-initial_p0ps[i,i]*(-1)
  }
}else{
  initial_p0ps<-matrix(s1*(-1))
}
initial_nups<-apply(initial_p0ps,MARGIN = 1,sum)*(-1)

if(m2>1){
  initial_p1ps<-diag(rep(s2,m2))*(-1)
  for(i in 1:(m2-1)){
    initial_p1ps[i,i+1]<-initial_p1ps[i,i]*(-1)
  }
}else{
  initial_p1ps<-matrix(s2*(-1))
}
initial_nu1ps<-apply(initial_p1ps,MARGIN = 1,sum)*(-1)


ptm<-proc.time()
ites=0
p0ps=array(0,dim=c(m1,m1,2))#intensity matrix
p1ps=array(0,dim=c(m2,m2,2))
nups=matrix(0,2,m1)#intensity to absorbing state
nu1ps=matrix(0,2,m2)
pips=matrix(0,2,m1)#initial probability
pi1ps=matrix(0,2,m2)
mups=rep(0,2)#randomization parameter
mu1ps=rep(0,2)

p0ps[,,1]<-initial_p0ps
p1ps[,,1]<-initial_p1ps
nups[1,]<-initial_nups
nu1ps[1,]<-initial_nu1ps
pips[1,]<-initial_pips
pi1ps[1,]<-initial_pi1ps
if(m1>1){
  mups[1]=max(abs(diag(p0ps[,,1])))+0.1
}else{
  mups[1]=abs(p0ps[,,1])+0.1
}
if(m2>1){
  mu1ps[1]=max(abs(diag(p1ps[,,1])))+0.1
}else{
  mu1ps[1]=abs(p1ps[,,1])+0.1
}

############
Rmax=100
CRmax=100
increment<-Inf
###EM
while(((increment>1e-3)|(ites<10))&(increment>0)){
  ########save the parameter estimation of the last iteration
  pips[2,]=pips[1,]
  nups[2,]=nups[1,]
  p0ps[,,2]=p0ps[,,1]
  mups[2]=mups[1]
  pi1ps[2,]=pi1ps[1,]
  nu1ps[2,]=nu1ps[1,]
  p1ps[,,2]=p1ps[,,1]
  mu1ps[2]=mu1ps[1]
  updates=emcore(ts,pips[1,],nups[1,],p0ps[,,1],mups[1],Rmax,
                 lags,pi1ps[1,],nu1ps[1,],p1ps[,,1],mu1ps[1],censor,warr,ncensor,CRmax)
  pips[1,]=updates[[1]]
  nups[1,]=updates[[2]]
  p0ps[,,1]=updates[[3]]
  mups[1]=updates[[4]]
  pi1ps[1,]=updates[[5]]
  nu1ps[1,]=updates[[6]]
  p1ps[,,1]=updates[[7]]
  mu1ps[1]=updates[[8]]
  increment<-loglike(pips[1,],pi1ps[1,],p0ps[,,1],p1ps[,,1],nups[1,],nu1ps[1,],ts,lags,warr,censor,ncensor)-loglike(pips[2,],pi1ps[2,],p0ps[,,2],p1ps[,,2],nups[2,],nu1ps[2,],ts,lags,warr,censor,ncensor)
  ites=ites+1
  updates<-canonical(pips[1,],nups[1,],p0ps[,,1],pi1ps[1,],nu1ps[1,],p1ps[,,1])
  pips[1,]<-updates[[1]]
  nups[1,]<-updates[[2]]
  p0ps[,,1]<-updates[[3]]
  pi1ps[1,]<-updates[[4]]
  nu1ps[1,]<-updates[[5]]
  p1ps[,,1]<-updates[[6]]
  print(c(ites,increment))
}
Performance1_MLE2<-as.numeric((proc.time()-ptm)[3])
Performance2_MLE2<-ites
##########estimate
pi_est1=pips[1,];pi_est2=pi1ps[1,]
xi_est1=nups[1,];xi_est2=nu1ps[1,]
T_est1=p0ps[,,1];T_est2=p1ps[,,1]
EM_loglikelihood<-loglike(pips[1,],pi1ps[1,],p0ps[,,1],p1ps[,,1],nups[1,],nu1ps[1,],ts,lags,warr,censor,ncensor)
EM_increment<-increment
