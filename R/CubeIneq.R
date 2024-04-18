cube.ineq <- function(X,pik,B,r,index=0.01,deepness=1,EPS=0.0000000001){
  require(StratifiedSampling)
  
  nc <- ncol(B)/length(pik)
  nnc <- length(length(pik)*1:nc)
  rowB <- c()
  
  for(i in 1:nnc){
    if(i==1){
      rowB <- rbind(rowB,B[,1:(length(pik)*(1:nc))[i]])
    }else{
      rowB <- rbind(rowB,B[,(nrow(rowB)+1):(length(pik)*(1:nc))[i]]) 
    }
  }
  
  # s <- fastflightcube(cbind(X,t(B)*pik),pik,comment=F)
  s <- ffphase(cbind(pik,X,t(rowB)*pik),pik)
  
  add.i <- 1
  
  
  while(all(round(s,5)==0 | round(s,5)==1)==F){
    
    # print(round(s,5))
    s <- fast.flight.cube.ineq(s,s,B,sign(r)*abs(r)*add.i,deepness,EPS)
    # print(sum(s))
    add.i <- add.i+index
    # if(length(round(s,5)[which(round(s,5)!=0 & round(s,5)!=1)])==1){s <- round(s)}
    
  }
  return(round(s,5))
  
} 
