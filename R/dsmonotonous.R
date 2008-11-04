dsmonotonous <-
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
xtestlo=matrix(0,1,dim(xlo)[2]);
direction=matrix(NA,1,dim(xlo)[2]);
xveclo=matrix(0,dim(xlo)[1],dim(xlo)[2]);
xvechi=matrix(0,dim(xlo)[1],dim(xlo)[2]);
xtestlo[1,]=xlo[1,];
#if(exists("optdirection")==FALSE || length(optdirection)!=length(direction)){
direction=TRUE;
for (i in 1:dim(xlo)[2]){
xtesthi=xtestlo;
xtesthi[i]=xhi[1,i];
ytestlo=fhandle(xtestlo,...);
ytesthi=fhandle(xtesthi,...);
direction[i]=ytesthi>ytestlo;
}
#optdirection<<-direction
#print("Set optdirection to:")
#print(direction)
#}
#direction=optdirection;
xveclo[,direction]=as.matrix(xlo[,direction]);
xveclo[,!direction]=as.matrix(xhi[,!direction]);
xvechi[,direction]=as.matrix(xhi[,direction]);
xvechi[,!direction]=as.matrix(xlo[,!direction]);
y[,1]=fhandle(xveclo,...);
y[,2]=fhandle(xvechi,...);
y[,3]=mass[,1];
if (dim(mass)[2]>1){
for (i in 2:dim(mass)[2]){ y[,3]=y[,3]*mass[,i]};
}
if(is.null(samples)){
samples=list();
samples$s=matrix(NA,0,dim(xlo)[2]+1)
samples$lo=matrix(NA,0,dim(xlo)[2]+1)
samples$hi=matrix(NA,0,dim(xlo)[2]+1)
samples$m=matrix(NA,0,dim(xlo)[2]+1)
}
samples$s=rbind(samples$s,cbind(xveclo,y[,1]),cbind(xvechi,y[,2]))
samples$lo=rbind(samples$lo,cbind(xlo,y[,1]))
samples$hi=rbind(samples$hi,cbind(xhi,y[,2]))
samples$m=rbind(samples$m,cbind(mass,y[,3]))
dsbound=list(y,samples);
}

