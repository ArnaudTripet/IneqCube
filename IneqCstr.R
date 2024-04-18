ineq.cstr <- function(X,pik,cstr="lower",bound=1){
  if(sum(tolower(cstr)==c("greater","less","two.sided"))!=1){
    stop("'cstr' should be one of “two.sided”, “less”, “greater”")
  }
  N=length(pik)
  DD=as.matrix(dist(X))
  
  size <- 1
  
  if(tolower(cstr)=="two.sided"){
    size <- 2
  }
  
  
  IND=array(0,c(N,size*N*length(bound)))
  
  M <- 1
  
  if(cstr=="less"){
    for(bound_i in bound){
      for(i in ( (1:N) + (N*(1:(size*length(bound))-1)[M])) ){
        o <- order(DD[,i==( (1:N) + (N*(1:(size*length(bound))-1)[M])) ])
        cond <- cumsum(pik[o])<=bound_i
        IND[o,i][cond]=1
      }
      M <- M+1
    }
  }
  if(cstr=="greater"){
    for(bound_i in bound){
      for(i in ( (1:N) + (N*(1:(size*length(bound))-1)[M])) ){
        o=order(DD[,i==( (1:N) + (N*(1:(size*length(bound))-1)[M]))])
        cond=cumsum(pik[o])<=bound_i
        z=max(which(cond,TRUE))
        if(z<N & bound_i != sum(pik[cond])){cond[z+1]=TRUE}
        IND[o,i][cond]=-1
      }
      M <- M+1
    }
  }
  if(cstr=="two.sided"){
    
    
    
    for(bound_i in bound){
      
      for(i in ( (1:N) + (N*(1:(size*length(bound))-1)[M])) ){
        o <- order(DD[,i==( (1:N) + (N*(1:(size*length(bound))-1)[M]))])
        cond <- cumsum(pik[o])<=bound_i
        IND[o,i][cond]=1
      }
      M <- M+1
    }
    for(bound_i in bound){
      for(k in ( (1:N) + (N*(1:(size*length(bound))-1)[M])) ){
        o=order(DD[,k==( (1:N) + (N*(1:(size*length(bound))-1)[M]))])
        cond=cumsum(pik[o])<=bound_i
        z=max(which(cond,TRUE))
        if(z<N & bound_i != sum(pik[cond]) ){cond[z+1]=TRUE}
        IND[o,k][cond]=-1
      }
      M <- M+1
    }
  }
  return(IND)
}