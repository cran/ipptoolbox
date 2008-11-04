`dsvariance` <-
function(x){
# Calculates variance of a BPA
# This is an implementation of the fast bisection algorithms of Kreinovich et al.
# The algorithm is described in detail in the following paper:
# Kreinovich, V., G. Xiang, et al. (2006).
# "Computing mean and variance under Dempster–Shafer uncertainty:
# Towards faster algorithms."
# International Journal of Approximate Reasoning 42(3): 212-227.
#=========================================================================   
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

#Function for calculating lower bound of the variance
calcvardown=function(x){
xhilo=rbind(cbind(x[,1],x[,3],0),cbind(x[,2],x[,3],1),Inf);
xhilo=xhilo[order(xhilo[,1]),]
k=bisect(testdown,xhilo);
indhi=xhilo[,3]&(xhilo[,1]<=xhilo[k,1]);
indlo=!xhilo[,3]&(xhilo[,1]>=xhilo[k+1,1]);
S=sum(xhilo[indhi,1]*xhilo[indhi,2])+sum(xhilo[indlo,1]*xhilo[indlo,2])
Sig=sum(xhilo[indhi,2])+sum(xhilo[indlo,2])
if(Sig==0){
Vlo=0;
}else{
rk=S/Sig;
Vlo=sum((xhilo[indhi,1]-rk)^2*xhilo[indhi,2])+sum((xhilo[indlo,1]-rk)^2*xhilo[indlo,2])
}
Vlo
}

#Function for calculating upper bound of the variance
calcvarup=function(x){
xmid=cbind((x[,1]+x[,2])/2,x)
xmid=xmid[order(xmid[,1]),]
n=dim(xmid)[1];
temp=rle(xmid[,1]);
temp$values=cumsum(temp$lengths)
t1=inverse.rle(temp)
temp$values=temp$values-temp$lengths+1
t2=inverse.rle(temp)
xmid=cbind(xmid,t1,t2)
k=bisect(testup,xmid);
if(is.null(k)){
Vup1=-Inf;
}else{
E=0;
Vup1=0;
if(k>1){
E=E+sum(xmid[1:(k-1),2]*xmid[1:(k-1),4])
}
E=E+sum(xmid[k:n,3]*xmid[k:n,4])
if(k>1){
Vup1=sum((xmid[1:(k-1),2]-E)^2*xmid[1:(k-1),4])
}
Vup1=Vup1+sum((xmid[k:n,3]-E)^2*xmid[k:n,4])
}
k=bisect(testup2,xmid);
if(is.null(k)){
Vup2=-Inf;
}else{
Vup2=0;
if(k>1){
Vup2=sum((xmid[1:(k-1),2]-xmid[k,1])^2*xmid[1:(k-1),4])
}
Vup2=Vup2+sum((xmid[k:n,3]-xmid[k,1])^2*xmid[k:n,4])
}
erg=max(Vup1,Vup2)
}

#Test if bisection for lower bound algorithm is successfully terminated
testdown=function(x,k,lb,ub){
if (lb>=ub){
return(c(1,k,lb-1,ub-1))
}
xk=x[k,1];
n=dim(x)[1];
indhi=x[,3]&(x[,1]<=xk);
xk2=x[k+1,1];
indlo=!x[,3]&(x[,1]>=xk2);
erg1=sum((xk-x[indhi,1])*x[indhi,2])
erg2=sum((x[indlo,1]-xk)*x[indlo,2])
erg3=sum((xk2-x[indhi,1])*x[indhi,2])
erg4=sum((x[indlo,1]-xk2)*x[indlo,2])
if(erg1>erg2){
ub=k-1
return(c(0,floor((lb+ub)/2),lb,ub))
}else {
if(erg3<=erg4){
lb=k+1;
return(c(0,floor((lb+ub)/2),lb,ub))
} else {
return(c(1,floor((lb+ub)/2),lb,ub))
}
}
}

#Test if bisection for upper bound algorithm is successfully terminated
testup=function(xmid,k,lb,ub){
xk=xmid[k,1];
xk2=xmid[k-1,1];
n=dim(xmid)[1];
erg1=sum((xmid[k:n,3]-xk)*xmid[k:n,4])
erg2=sum((xk-xmid[1:(k-1),2])*xmid[1:(k-1),4])
erg3=sum((xmid[k:n,3]-xk2)*xmid[k:n,4])
if(k>1){
erg4=sum((xk2-xmid[1:(k-1),2])*xmid[1:(k-1),4])
} else{
erg4=0;
}
if(erg1>=erg2){
lb=k+1;
return(c(0,floor((lb+ub)/2),lb,ub))
}else {
if(erg3<erg4){
ub=k-1
return(c(0,floor((lb+ub)/2),lb,ub))
} else {
return(c(1,floor((lb+ub)/2),lb,ub))
}
}
}
testup2=function(xmid,k,lb,ub){
xk=xmid[k,1];
n=dim(xmid)[1];
lk=xmid[k,5];
if(k<n){
erg1=sum((xmid[(lk+1):n,3]-xk)*xmid[(lk+1):n,4])
} else {
erg1=0
}
erg2=sum((xk-xmid[1:lk,2])*xmid[1:lk,4])
erg3=sum((xmid[k:n,3]-xk)*xmid[k:n,4])
if(k>1){
erg4=sum((xk-xmid[1:(k-1),2])*xmid[1:(k-1),4])
} else{
erg4=0;
}
if(erg1>erg2){
lb=k+1;
k=floor((lb+ub)/2)
return(c(0,xmid[k,6],lb,ub))
}else {
if(erg3<erg4){
ub=k-1
k=floor((lb+ub)/2)
return(c(0,xmid[k,6],lb,ub))
} else {
k=floor((lb+ub)/2)
return(c(1,xmid[k,6],lb,ub))
}
}
}
#Bisection algorithm for function fun and starting position xmid
bisect=function(fun,xmid){
n=dim(xmid)[1];
temp=c(0,floor(n/2+1),1,n+1);
tempold=c(0,0,0,Inf)
nofound=FALSE;
while(temp[1]==0&&!nofound) {
nofound=!(temp[3]>tempold[3]||temp[4]<tempold[4])
tempold=temp;
temp=fun(xmid,temp[2],temp[3],temp[4]);
} 
if(nofound){
return (NULL)}else{
return(temp[2])
}
}
#Calculate variance using lower and upper bound algorithms
varup=calcvarup(x);
vardown=calcvardown(x);
erg=c(vardown,varup);
return(erg)
}

