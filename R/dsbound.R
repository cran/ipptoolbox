dsbound <-
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
    a=fhandle(xlo,...);
b=fhandle(xhi,...);
c=mass[,1];
if (dim(mass)[2]>1){
for (i in 2:dim(mass)[2]){ c=c*mass[,i]};
}

	temp=cbind(a,b);
if(is.null(samples)){
samples=matrix(NA,0,dim(xlo)[2]+1)
}
samples=rbind(samples,cbind(xlo,a),cbind(xhi,b))
    y[,1]=apply(temp,1,min);
    y[,2]=apply(temp,1,max);    
y[,3]=c;
dsbound=list(y,samples);
}

