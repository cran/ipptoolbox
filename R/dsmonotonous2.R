dsmonotonous2 <-
function(fhandle,xlo,xhi,mass,...,samples=NULL){
# Evaluates f(min(a,b)) and f(max(a,b)) boundary function values.
#=========================================================================   
# dsbound<- function(fhandle,xlo,xhi,mass) 
# 
# Propagates focal elements through a function. Used by Monte-Carlo
# propagation algorithm. Most simple and fastest "optimization routine".
# Only for monotonous increasing functions.
# Evaluates f(min(a,b)) and f(max(a,b)) boundary function values. Can
# process one focal element or a vector of focal elements at once.
#
# Input:
# fhandle: Function handle or function name y=f(x). 
# xlo:  Single focal element: lower bound of parameter vector.
#       Focal element vector: matrix containing the lower bounds of the
#       parameter vector.
# xhi:  Single focal element: upper bound of parameter vector.
#       Focal element vector: matrix containing the upper bounds of the
#       parameter vector.
# mass: Single focal element: Mass of parameter vector.
#       Focal element vector: Vector of masses.
#
# Output:
# y: Focal element(s) propagated through f
#
# Usage: Propagating through the function f
# f1=(function(x){ z <- 1-(sum(x));});
# lambda=dsstruct(c(2,3,1));
# x1=dsadf('qexp',100,lambda);
# x2=dsadf('qexp',100,lambda);
# xlo=cbind(as.matrix(x1[,1]),as.matrix(x2[,1]));
# xhi=cbind(as.matrix(x1[,2]),as.matrix(x2[,2]));
# mass=cbind(as.matrix(x1[,3]),as.matrix(x2[,3]));
# y=dsbound(f1,xlo,xhi,mass);
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
y=matrix(0,dim(xlo)[1],3);
twovec=2^(0:(dim(xlo)[2]-1));
iter=2^dim(xlo)[2]-1;
ytemp=matrix(0,dim(xlo)[1],iter+1);
for (i in 0:iter){
xvec=matrix(0,dim(xlo)[1],dim(xlo)[2]);
selection=floor(i/twovec)%%2;
xvec[,selection==0]=xlo[,selection==0];
xvec[,selection==1]=xhi[,selection==1];
ytemp[,i+1]=fhandle(xvec,...);
} 
    y[,1]=apply(ytemp,1,min);
    y[,2]=apply(ytemp,1,max);    
y[,3]=apply(mass,1,prod);
samples=NULL
dsbound=list(y,samples);
}

