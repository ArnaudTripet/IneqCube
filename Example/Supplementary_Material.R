###################################################
# Balanced Sampling With Inequalities
# RCode and Example
###################################################
# Library
require(sampling) # needed to compute inclusion probabilities
require(MASS) # Null fct
require(hitandrun) # inequality
require(StratifiedSampling) # Cube method
###################################################
# Load the function
###################################################
#######################################################
# Computations if several matrices of constraints
# Compute the inequality constraint
# X : the dataset
# pik : your inclusion probabilities
# cstr : the inequality bound the you want to obtain (Either <=, >= or both)
# Bound : bound at 1, 2 ... >, 1 by default
#######################################################
ineq.cstr <- function(X,pik,cstr="less",bound=1){
  if(sum(tolower(cstr)==c("greater","less","two.sided"))!=1){
    stop("’cstr’ should be one of \two.sided", "less", "greater")
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
###################################################################
# Fast cube ineq
# The inequality constraint is t(A)%*%S>=t(A)%*%pik whwrw A = X/pik
# The equality constraint is t(B)%*%S<=r
###################################################################
fast.flight.cube.ineq<-function(X,pik,B,r,deepness=1,EPS=1e-6){
A=as.matrix(X/pik)
if(nrow(A)==0) A=matrix(0,c(length(pik),1))
######
B=as.matrix(B/rep(1,length(pik)))
TT=!(nrow(B)==0)
if(TT) c=c(t(B)%*%pik)
######
if(TT)
{
TE=abs(c(t(B)%*%pik)-r)<=EPS
A=cbind(A,B[,TE])
r=r[!TE]
B=B[,!TE]
c=c[!TE]
}
B=as.matrix(B/rep(1,length(pik)))
TT=!(nrow(B)==0)
if(is.null(A)){
prof=deepness}else{
if(is.matrix(A))prof=ncol(A)+deepness else prof=1+deepness }
TEST=(EPS<pik) & (pik<1-EPS)
prof2=min(sum(TEST),prof)
if(TT & prof2!=0){
BR=matrix(B[TEST,],c(sum(TEST),length(B[TEST,])/sum(TEST)))[1:prof2,];
if(ncol(B)!=0) BR=matrix(BR,c(length(BR)/ncol(B),ncol(B)));}
if(prof2==0) a=0 else {pikR=pik[TEST][1:prof2];
AR=matrix(A[TEST,],c(sum(TEST),length(A[TEST,])/sum(TEST)))[1:prof2,];
AR=matrix(AR,c(length(AR)/ncol(A),ncol(A)));
if(nrow(AR)==1 & sum(AR)==0)
NN=matrix(1,c(1,1)) else NN=Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))));
a=ncol(NN)}
while(a>0.5)
{
#####
#####
u=NN[,1]
l1=min(pmax((1-pikR)/u,-pikR/u))
l2=min(pmax((pikR-1)/u,pikR/u))
if(TT ){
c=c(t(B)%*%pik)
cR=c(t(BR)%*%pikR)
pet=(r-c)/c(t(BR)%*%u)
if(sum(pet>0)>0) nu1=min(pet[pet>0]) else nu1=100000000000
if(sum(pet<0)>0) nu2=min(-pet[pet<0]) else nu2=100000000000
l1=min(l1,nu1)
l2=min(l2,nu2)
}
pik[TEST][1:prof2]=if (runif(1)<l2/(l1+l2)) pikR+l1*u else pikR-l2*u
TEST=(EPS<pik) & (pik<1-EPS)
#####
if(TT)
{
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
if(TT & prof2!=0)
{BR=matrix(B[TEST,],c(sum(TEST),length(B[TEST,])/sum(TEST)))[1:prof2,];
if(ncol(B)!=0) BR=matrix(BR,c(length(BR)/ncol(B),ncol(B)));}
if(prof2==0) a=0 else {pikR=pik[TEST][1:prof2];
AR=matrix(A[TEST,],c(sum(TEST),sum(A[TEST,])/sum(TEST)))[1:prof2,];
AR=matrix(AR,c(length(AR)/ncol(A),ncol(A)));
if(sum(abs(AR))==0)
NN=matrix(1,c(nrow(AR),1)) else
NN=Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))));
a=ncol(NN)}
####
######
}
pik
}
###################################################
# Sampling function
# B : inequality the constraints
# r: The equality constraint is t(B)%*%S<=r
###################################################
cube.ineq <- function(X,pik,B,r,index=0.01,deepness=1,EPS=1e-6){
fast.flight.cube.ineq<-function(X,pik,B,r,deepness=1,EPS=0.00000001){
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
if(is.null(A))
prof=deepness else
{if(is.matrix(A))prof=ncol(A)+deepness else prof=1+deepness }
TEST=(EPS<pik) & (pik<1-EPS)
prof2=min(sum(TEST),prof)
if(TT & prof2!=0)
{BR=matrix(B[TEST,],c(sum(TEST),length(B[TEST,])/sum(TEST)))[1:prof2,];
if(ncol(B)!=0) BR=matrix(BR,c(length(BR)/ncol(B),ncol(B)));}
if(prof2==0) a=0 else {pikR=pik[TEST][1:prof2];
AR=matrix(A[TEST,],c(sum(TEST),length(A[TEST,])/sum(TEST)))[1:prof2,];
AR=matrix(AR,c(length(AR)/ncol(A),ncol(A)));
if(nrow(AR)==1 | sum(abs(AR))==0) NN=matrix(1,c(1,1)) else
NN=Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))))}
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
if(sum(pet>0)>0) nu1=min(pet[pet>0]) else nu1=100000000000
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
if(TT & prof2!=0)
{BR=matrix(B[TEST,],c(sum(TEST),length(B[TEST,])/sum(TEST)))[1:prof2,];
if(ncol(B)!=0) BR=matrix(BR,c(length(BR)/ncol(B),ncol(B)));}
if(prof2==0) a=0 else {pikR=pik[TEST][1:prof2];
AR=matrix(A[TEST,],c(sum(TEST),sum(A[TEST,])/sum(TEST)))[1:prof2,];
AR=matrix(AR,c(length(AR)/ncol(A),ncol(A)));
if(sum(abs(AR))==0)
{NN=matrix(1,c(nrow(AR),1))} else
{u <- Kernel(matrix(AR,c(length(AR)/ncol(A),ncol(A))))}
a=sum(abs(u))}
####
}
pik
}
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
s <- ffphase(cbind(pik,X,t(rowB)*pik),pik)
add.i <- 1
while(all(round(s,5)==0 | round(s,5)==1)==F){
s <- fast.flight.cube.ineq(s,s,B,sign(r)*abs(r)*add.i,deepness,EPS)
add.i <- add.i+index
}
return(round(s,5))
}
###################################################
# Matrix rounding example
###################################################
###################################################
# Set parameters here
###################################################
values <- c(15,21,17,09,
10,08,13,07,
06,09,05,08,
04,03,06,06,
03,02,05,08) # Value orignal matrix
M1 <- matrix(values,5,4,byrow=T) # Changes to create a rounding problem
M <- M1/16.5
M[which(M1/16.5 >1 )] <- c((M1/16.5)[which(M1/16.5 >1 )] -1) # Final Matrix
############################################## Application
X=cbind(disjunctive(rep(1:5,4)),disjunctive(sort(rep(1:4,5)))) # X
pik <- c(M) # Inclusion probabilities
B=cbind(X,-X) # Constraints Upper and less
r <- c(ceiling(t(B)%*%pik))
s <- fast.flight.cube.ineq(pik,pik,B,r) # Application of the algorithm
############################################## Results
(Z1=addmargins(round(matrix(s,dim(M)),4))) # Rouned Matrix with margins
(Z2=addmargins(M)) # Original Matrix with margins
Z1-Z2 # Must only contain values between -1 and 1
###################################################
# Spatial sampling example
###################################################
###################################################
# Set population parameters here
###################################################
### Creater dataset
N <- 100 # Size of the dataset
n <- 20 # Size of the sample
X <- cbind(runif(N),runif(N)) # dataset
pik <- inclusionprobabilities(runif(N),n) # Inclusion probabilities
sum(pik) # should be equal to n
### Inequality Constraints
LG1 <- ineq.cstr(X,pik,"two.sided",1) # inequality less and greater than 1
r <- ceiling(t(LG1) %*% pik)
###################################################
# Draw Sample
###################################################
s <- cube.ineq(X,pik,LG1,r) # Draw sample
###################################################
# A look at the sample
###################################################
plot(X) # plot the population
points(X[s==1,],pch=19) # show selected sample