`dsopt` <-
function(fr,xlo,xhi,mass,...,samples=NULL){
# Optimization routine for complex nonlinear functions using Matlabs Optimization Toolbox
#=========================================================================   
# dsopt<- function(fr,xlo,xhi,mass) 
# 
# Propagates focal elements through a function. Used by Monte-Carlo
# propagation algorithm. Time-consuming and
# accurate "optimization routine" for complex nonlinear functions.
# Requires Matlabs Optimization toolbox routine fmincon.
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
# vargin: Parameters of function fhandle
#
# Output:
# y: Focal element(s) propagated through f
#
# Usage: Propagating through Sinus function
# f1=(function(x){ z <- sin(sum(x));});
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

myfminmax=function(f,mylb,myub,startxmin,startxmax,...){
indX=which(mylb!=myub);
indD=which(mylb==myub);
lb=mylb[indX];
ub=myub[indX];
d=mylb[indD];
size=length(mylb);
erg=list();
erg$value=c(0,0);
if(size!=length(indX)){
}
y=optim(par = startxmin[indX], fn = myfun,gr=NULL,f,indX,d,indD,size,ispositive=TRUE,..., method = "L-BFGS-B", lower=lb, upper = ub);
erg$value[1]=y$value;
erg$parmin=startxmin;
erg$parmin[indX]=y$par;
erg$parmin[indD]=d;
y=optim(par = startxmax[indX], fn = myfun,gr=NULL,f,indX,d,indD,size,ispositive=FALSE,..., method = "L-BFGS-B", lower=lb, upper = ub);
erg$value[2]=-y$value;
erg$parmax=startxmax;
erg$parmax[indX]=y$par;
erg$parmax[indD]=d;
erg
}

myfun=function(x,f,indX,d,indD,size,ispositive,...){
if(size!=length(indX)){
print(size)
print(indX)
print(indD)
print(x)
}
myx=matrix(0,1,size);
myx[indX]=x;
myx[indD]=d;
y=f(myx,...);
if(ispositive==FALSE){
y=-y
}
y
}

myconstrfminmax=function(f,mylb,myub,startxmin,startxmax,...){
indX=which(mylb!=myub);
indD=which(mylb==myub);
lb=mylb[indX];
ub=myub[indX];
d=mylb[indD];
size=length(mylb);
sizereal=length(lb);
erg=list();
erg$value=c(0,0);
if(size!=length(indX)){
print(lb)
print(ub)
}

m1=matrix(0,sizereal,sizereal);
m1[cbind(1:(sizereal),1:(sizereal))]=1;
m2=m1;
m2[cbind(1:(sizereal),1:(sizereal))]=-1;
ci=c(lb,-ub)
ui=rbind(m1,m2);
y=constrOptim(startxmin[indX], myfun,grad=NULL,ui, ci,mu=1e-4,control=list(maxit=70),method="Nelder-Mead",outer.iterations=50,outer.eps=1e-5, f,indX,d,indD,size,...);
erg$value[1]=y$value;
erg$parmin=startxmin;
erg$parmin[indX]=y$par;
erg$parmin[indD]=d;
y=constrOptim(startxmin[indX], myfun,grad=NULL,ui, ci,mu=1e-4,control=list(maxit=70,fnscale=-1),method="Nelder-Mead",outer.iterations=50,outer.eps=1e-5, f,indX,d,indD,size,...);
erg$value[2]=y$value;
erg$parmax=startxmax;
erg$parmax[indX]=y$par;
erg$parmax[indD]=d;
erg
}


siz=dim(xlo)[1];
dimx=dim(xlo)[2];
archivexmin=matrix(NA,siz,dimx);
archivexmax=matrix(NA,siz,dimx);
archiveymin=matrix(NA,siz,1);
archiveymax=matrix(NA,siz,1);
amin=matrix(FALSE,siz,1);
amax=matrix(FALSE,siz,1);
erg=matrix(0,siz,3);
for (i in 1:siz){
for (j in 1:dimx){
amin[amin]=archivexmin[amin,j]>=xlo[i,j]&archivexmin[amin,j]<=xhi[i,j];
amax[amax]=archivexmax[amax,j]>=xlo[i,j]&archivexmax[amax,j]<=xhi[i,j];
}
a=which.min(archiveymin[amin,])
if(length(a)>0){
startxmin=archivexmin[amin,,drop=0][a,,drop=0];
} else {
startxmin=(xlo[i,]+xhi[i,])/2;
}
a=which.max(archiveymax[amax,])
if(length(a)>0){
startxmax=archivexmax[amax,,drop=0][a,,drop=0];
} else {
startxmax=(xlo[i,]+xhi[i,])/2;
}
temp=myfminmax(fr,xlo[i,],xhi[i,],startxmin,startxmax,...);
erg[i,1:2]=temp$value;
#erg[i,2]=myfmax(fr,xlo[i,],xhi[i,]);
erg[i,3]=prod(mass[i,]);
archivexmin[i,]=temp$parmin;
archivexmax[i,]=temp$parmax;
archiveymin[i,1]=temp$value[1];
archiveymax[i,1]=temp$value[2];
amin[1:i]=TRUE;
amax[1:i]=TRUE;
}
if(is.null(samples)){
samples=list();
samples$s=matrix(NA,0,dim(archivexmin)[2]+1)
samples$lo=matrix(NA,0,dim(archivexmin)[2]+1)
samples$hi=matrix(NA,0,dim(archivexmin)[2]+1)
samples$m=matrix(NA,0,dim(archivexmin)[2]+1)
}
samples$s=rbind(samples$s,cbind(archivexmin,archiveymin),cbind(archivexmax,archiveymax))
samples$lo=rbind(samples$lo,cbind(archivexmin,archiveymin))
samples$hi=rbind(samples$hi,cbind(archivexmax,archiveymax))
samples$m=rbind(samples$m,cbind(mass,erg[,3]))

dsopt=list(erg,NULL);

}

