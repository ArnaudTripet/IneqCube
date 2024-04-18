
fast.flight.cube.ineq<-function(X,pik,B,r,deepness=1,EPS=0.00000001){
  
  require(sampling)
  require(hitandrun)
  
  Kernel <- function(A) residuals(lm(rnorm(nrow(A))~A-1))
  
  A=as.matrix(X/pik)
  if(nrow(A)==0) A=matrix(0,c(length(pik),1))
  
  ######
  B=as.matrix(B/rep(1,length(pik)))
  TT=!(nrow(B)==0)
  if(TT) c=c(t(B)%*%pik)
  
  
  
  ######
  
  if(TT){
    TE=abs(c(t(B)%*%pik)-r)<=EPS
    A=cbind(A,B[,TE])
    r=r[!TE]
    B=B[,!TE]
    c=c[!TE]
  }
  
  B=as.matrix(B/rep(1,length(pik)))
  TT=!(nrow(B)==0)
  
  if(is.null(A)) prof=deepness else{if(is.matrix(A))prof=ncol(A)+deepness else prof=1+deepness }
  TEST=(EPS<pik) & (pik<1-EPS)
  prof2=min(sum(TEST),prof)
  if(TT & prof2!=0) {BR=matrix(B[TEST,],c(sum(TEST),length(B[TEST,])/sum(TEST)))[1:prof2,];
  if(ncol(B)!=0) BR=matrix(BR,c(length(BR)/ncol(B),ncol(B)));}
  if(prof2==0) a=0 else {pikR=pik[TEST][1:prof2];
  AR=matrix(A[TEST,],c(sum(TEST),length(A[TEST,])/sum(TEST)))[1:prof2,];
  AR=matrix(AR,c(length(AR)/ncol(A),ncol(A)));
  if(nrow(AR)==1 | sum(abs(AR))==0) NN=matrix(1,c(1,1)) else NN=Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))))}
  if(ncol(NN)==0){a <- 0}else{
    a <- 1
    u <- NN[,1]
  }
  
  # print(NN)
  while(a>0.000001){
    
    #####
    #####
    
    
    l1=min(pmax((1-pikR)/u,-pikR/u))
    l2=min(pmax((pikR-1)/u,pikR/u))
    if(TT ){
      c=c(t(B)%*%pik)
      cR=c(t(BR)%*%pikR)
      pet=(r-c)/c(t(BR)%*%u)
      if(sum(pet>0)>0) nu1=min(pet[pet>0])  else nu1=100000000000
      if(sum(pet<0)>0) nu2=min(-pet[pet<0]) else nu2=100000000000
      l1=min(l1,nu1)
      l2=min(l2,nu2)
    }  
    pik[TEST][1:prof2]=if (runif(1)<l2/(l1+l2)) pikR+l1*u else pikR-l2*u
    TEST=(EPS<pik) & (pik<1-EPS)
    
    ##### 
    
    if(TT){
      TE=abs(c(t(B)%*%pik)-r)<=EPS
      A=cbind(A,B[,TE])
      r=r[!TE]
      B=B[,!TE]
      c=c[!TE]
    }
    B=as.matrix(B/rep(1,length(pik)))
    TT=!(nrow(B)==0)
    
    ###
    
    prof=ncol(A)+deepness
    prof2=min(sum(TEST),prof)
    if(TT & prof2!=0) {BR=matrix(B[TEST,],c(sum(TEST),length(B[TEST,])/sum(TEST)))[1:prof2,];
    if(ncol(B)!=0) BR=matrix(BR,c(length(BR)/ncol(B),ncol(B)));}
    if(prof2==0) a=0 else {pikR=pik[TEST][1:prof2];
    AR=matrix(A[TEST,],c(sum(TEST),sum(A[TEST,])/sum(TEST)))[1:prof2,];
    AR=matrix(AR,c(length(AR)/ncol(A),ncol(A)));
    if(sum(abs(AR))==0) {NN=matrix(1,c(nrow(AR),1))} else {u <- Kernel(matrix(AR,c(length(AR)/ncol(A),ncol(A))))}
    a=sum(abs(u))}
    
    # if(sum(abs(AR))==0) {NN=matrix(1,c(nrow(AR),1))} else {NN=Null(matrix(as.numeric(as.factor(AR)),c(length(AR)/ncol(A),ncol(A)))-1)}
    # a=ncol(NN)} 
    
    ####
  }
  
  pik
  
}