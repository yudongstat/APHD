library(expm)
library(mapfit)
library(MCMCpack)
library(pracma)
library(FAdist)
library(Rsolnp)
library(pracma)
furnace<-read.table(file = "C:/Users/e0265296/Desktop/furnace.txt",header = TRUE)
##Following Ghosh(2010), we standardize the data by dividing its std deviation
furnace[,2]<-furnace[,2]/sqrt(var(furnace[,2]))
furnace[,3]<-furnace[,3]/sqrt(var(furnace[,3]))
##warr is the censoring time which acts on the sum of the sales lag and product lifetime
warr<-max(furnace[,2]+furnace[,3])
##ss0 is the total number of products under study
ss0<-400
lags<-furnace[,2]
ts<-furnace[,3]
pre_ts<-ts
pre_lags<-lags
pre_warr<-warr
dm<-length(ts)
########################ncensor is the number of censored products
ncensor<-400-133
censor<-rep(warr,ncensor)

emcore<-function(y,pips,nups,p0ps,mups,Rmax,
                 y1,pi1ps,nu1ps,p1ps,mu1ps,censor,ncensor,CRmax){
  ###for lags date
  ss=length(y);m=length(pips)
  if(m==1){
    p0ps<-matrix(p0ps)
  }
  y=sort(y);dy=c(y[1],diff(y));re_p0=p0ps/mups+diag(m)
  ######compute fk,bk,ck, k=1,..ss
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
  #######compute H_k for each dy[i]
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
  newHs=matrix(0,m,m)
  diag(newHs)<-diag(Hs)
  if(m>1){
    for(i in 1:(m-1)){
      newHs[i,i+1]<-Hs[i,i+1]
    }
  }
  Hs<-newHs
  #########################observed part
  Bs=rep(0,m);As=rep(0,m)
  for(i in 1:ss){
    As=As+fks[,i]*nups/sum(pips*bks[,i+1])
    Bs=Bs+pips*bks[,i+1]/sum(pips*bks[,i+1])
  }
  As<-c(rep(0,m-1),As[m])
  #########for product lifetime
  m1=length(pi1ps)
  if(m1==1){
    p1ps<-matrix(p1ps)
  }
  y1=sort(y1);dy1=c(y1[1],diff(y1));re_p01=p1ps/mu1ps+diag(m1)
  ######compute fk,bk,ck, k=1,..ss
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
  #######compute H_k for each dy[i]
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
  newHs1=matrix(0,m1,m1)
  diag(newHs1)<-diag(Hs1)
  if(m1>1){
    for(i in 1:(m1-1)){
      newHs1[i,i+1]<-Hs1[i,i+1]
    }
  }
  Hs1<-newHs1
  #########################observed part
  Bs1=rep(0,m1);As1=rep(0,m1)
  for(i in 1:ss){
    As1=As1+fks1[,i]*nu1ps/sum(pi1ps*bks1[,i+1])
    Bs1=Bs1+pi1ps*bks1[,i+1]/sum(pi1ps*bks1[,i+1])
  }
  As1<-c(rep(0,m1-1),As1[m1])
  #####for missing data 
  mm1=m+m1;es=rep(1,mm1)
  D0=cbind(rbind(p0ps,matrix(0,m1,m)),rbind(nups%*%t(pi1ps),p1ps))
  myhc<-rep(0,mm1)
  for(i in 1:ncensor){
    hc=c(pips,rep(0,m1))*expAtv(D0,es,censor[i])$eAtv
    myhc<-myhc+hc/sum(hc)
  }
  Bs=Bs+myhc[1:m]
  
  #for the computation of censored part, unnecessary to use bakward-forward technique, directly use (35) in Okamura's paper
  mmups=max(mups,mu1ps)#randomization parameter
  Alpha<-c(pips,rep(0,m1))
  censor=sort(censor);dy=c(censor[1],diff(censor));re_D0=D0/mmups+diag(mm1)
  ######compute fk,bk,ck, k=1,..ncensor
  fks=matrix(expAtv(t(D0),Alpha,dy[1])$eAtv,mm1,ncensor)
  bks=matrix(expAtv(D0,es,dy[1])$eAtv,mm1,ncensor)
  for(i in 2:ncensor){
    fks[,i]=expAtv(t(D0),fks[,i-1],dy[i])$eAtv
    bks[,i]=expAtv(D0,bks[,i-1],dy[i])$eAtv
  }
  cks=matrix(Alpha/sum(Alpha*bks[,ncensor]),ncensor,mm1,byrow=T)
  for(i in (ncensor-1):1){
    cks[i,]=expAtv(t(D0),cks[i+1,],dy[i+1])$eAtv+Alpha/sum(Alpha*bks[,i])
  }
  bks=cbind(es,bks)
  #######compute H_k for each dy[i]
  CHks=matrix(0,mm1,mm1)
  for(i in 1:ncensor){
    Hks=matrix(0,mm1,mm1)
    alps=matrix(bks[,i],mm1,Rmax+1)
    for(j in 2:(Rmax+1)) alps[,j]=re_D0%*%alps[,j-1]
    bets=matrix(dpois(Rmax+1,mmups*dy[i])*cks[i,],Rmax+1,mm1,byrow=T)
    for(j in Rmax:1){
      bets[j,]=bets[j+1,]%*%re_D0+dpois(j,mmups*dy[i])*cks[i,]
    }
    for(j in 1:(Rmax+1)){
      Hks=Hks+alps[,j]%*%t(bets[j,])
    }
    CHks=CHks+Hks/mmups
  }
  CHks=t(CHks)
  #############
  Ns=Hs+CHks[1:m,1:m]; Zs=diag(Ns)
  Ns1=Hs1+CHks[(1+m):mm1,(1+m):mm1]; Zs1=diag(Ns1)
  subchks=nups%*%t(pi1ps)*CHks[1:m,(1+m):mm1]
  sumsub=apply(subchks,1,sum);sumsub1=apply(subchks,2,sum)
  ##estimate lags parameters: mu, pi, p0 and nu
  newpips=Bs/(ss+ncensor)
  newnups=(As+sumsub)/Zs
  Ns=Ns*p0ps
  newp0ps=matrix(0,m,m)
  for(j in 1:m){
    newp0ps[j,]=Ns[j,]/Zs[j]
    newp0ps[j,j]=-(sum(newp0ps[j,-j])+newnups[j])
  }
  newmups=max(abs(diag(newp0ps)))
  
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
  
  #estimate lifetime paramters
  ##sample mu, pi, p0 and nu
  newpi1ps=(Bs1+sumsub1)/sum(Bs1+sumsub1)
  newnu1ps=As1/Zs1
  Ns1=Ns1*p1ps
  newp1ps=matrix(0,m1,m1)
  for(j in 1:m1){
    newp1ps[j,]=Ns1[j,]/Zs1[j]
    newp1ps[j,j]=-(sum(newp1ps[j,-j])+newnu1ps[j])
  }
  newmu1ps=max(abs(diag(newp1ps)))
  
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
  ####results
  list(newpips,newnups,newp0ps,newmups,newpi1ps,newnu1ps,newp1ps,newmu1ps)
}


loglike<-function(apl1,apl2,p0ps,p1ps,nups,nu1ps,y1,y2,censor,ncensor){
  m1=length(apl1);m2=length(apl2);es=rep(1,m1+m2)
  D0=cbind(rbind(p0ps,matrix(0,m2,m1)),rbind(nups%*%t(apl2),p1ps))
  mysum<-0
  for(i in 1:ncensor){
    mysum<-mysum+log(sum(c(apl1,rep(0,m2))*expAtv(D0,es,censor[i])$eAtv))
  }
  return(sum(dph(y1,ph=ph(alpha=apl1,Q=p0ps,xi=nups),log=T))+
           sum(dph(y2,ph=ph(alpha=apl2,Q=p1ps,xi=nu1ps),log=T))+mysum)
}

observed_loglike<-function(mypara){
  pips<-mypara[1:m1]
  if(m1>1){
    p0ps<-diag(mypara[(m1+1):(2*m1)])*(-1)
    for(i in 1:(m1-1)){
      p0ps[i,i+1]<-p0ps[i,i]*(-1)
    }
  }else{
    p0ps<-matrix((mypara[(m1+1):(2*m1)])*(-1)) 
  }
  
  nups<-apply(p0ps,MARGIN = 1,sum)*(-1)
  pi1ps<-mypara[(2*m1+1):(2*m1+m2)]
  
  if(m2>1){
    p1ps<-diag(mypara[(2*m1+m2+1):(2*m1+2*m2)])*(-1)
    for(i in 1:(m2-1)){
      p1ps[i,i+1]<-p1ps[i,i]*(-1)
    }
  }else{
    p1ps<-matrix((mypara[(2*m1+m2+1):(2*m1+2*m2)])*(-1))
  }
  
  nu1ps<-apply(p1ps,MARGIN = 1,sum)*(-1)
  
  apl1<-pips
  apl2<-pi1ps
  y1<-lags
  y2<-ts
  es=rep(1,m1+m2)
  D0=cbind(rbind(p0ps,matrix(0,m2,m1)),rbind(nups%*%t(apl2),p1ps))
  mysum<-0
  for(i in 1:ncensor){
    mysum<-mysum+log(sum(c(apl1,rep(0,m2))*expAtv(D0,es,censor[i])$eAtv))
  }
  my_res<-sum(dph(y1,ph=ph(alpha=apl1,Q=p0ps,xi=nups),log=T))+sum(dph(y2,ph=ph(alpha=apl2,Q=p1ps,xi=nu1ps),log=T))+mysum
  return(-my_res)
}

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
  
  apl1<-pips
  apl2<-pi1ps
  y1<-lags
  y2<-ts
  es=rep(1,m1+m2)
  D0=cbind(rbind(p0ps,matrix(0,m2,m1)),rbind(nups%*%t(apl2),p1ps))
  mysum<-0
  for(i in 1:ncensor){
    mysum<-mysum+log(sum(c(apl1,rep(0,m2))*expAtv(D0,es,censor[i])$eAtv))
  }
  my_res<-sum(dph(y1,ph=ph(alpha=apl1,Q=p0ps,xi=nups),log=T))+sum(dph(y2,ph=ph(alpha=apl2,Q=p1ps,xi=nu1ps),log=T))+mysum
  return(-my_res)
}



m1<-4
m2<-7
ptm<-proc.time()
pips<-rep(1/m1,m1)
pi1ps<-rep(1/m2,m2)
###################use the SEM algorithm to select initial value
aaa<-solnp(pars = c(1,1),fun = mixture_Erlang,LB=rep(1e-2,2),control = list(trace=0))$pars
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
loglv=c(-Inf)
increment<-Inf
###EM
while(((increment>1e-3)|(ites<10))&(increment>0)){
  ########save the parameter estimation of the last iteration
  p0ps[,,2]=p0ps[,,1];nups[2,]=nups[1,];pips[2,]=pips[1,]
  p1ps[,,2]=p1ps[,,1];nu1ps[2,]=nu1ps[1,];pi1ps[2,]=pi1ps[1,]
  updates=emcore(lags,pips[1,],nups[1,],p0ps[,,1],mups[1],Rmax,
                 ts,pi1ps[1,],nu1ps[1,],p1ps[,,1],mu1ps[1],censor,ncensor,CRmax)
  pips[1,]=updates[[1]]
  nups[1,]=updates[[2]]
  p0ps[,,1]=updates[[3]]
  mups[1]=updates[[4]]
  pi1ps[1,]=updates[[5]]
  nu1ps[1,]=updates[[6]]
  p1ps[,,1]=updates[[7]]
  mu1ps[1]=updates[[8]]
  loglv=c(loglv,
          loglike(pips[1,],pi1ps[1,],p0ps[,,1],p1ps[,,1],nups[1,],
                  nu1ps[1,],lags,ts,censor,ncensor))
  ites=ites+1
  increment<-loglv[ites+1]-loglv[ites]
  print(c(ites,increment))
}
Performance1_MLE<-as.numeric((proc.time()-ptm)[3])
Performance2_MLE<-ites
##########estimate
pi_est1=pips[1,];pi_est2=pi1ps[1,]
xi_est1=nups[1,];xi_est2=nu1ps[1,]
T_est1=p0ps[,,1];T_est2=p1ps[,,1]
EM_loglikelihood<-loglv[length(loglv)]
EM_increment<-increment
