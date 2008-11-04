dsevalmc <-
function(fhandle,x, mcIT, optimizer=dsbound,corr=NULL,samples=NULL,fnoptions=NULL){
# Propagates Dempster-Shafer structures through arbitrary functions y=f(x)
#=========================================================================
# y=dsevalmc(fhandle,x, mcIT, optimizer,corr, varargin)
# Propagates Dempster-Shafer structures through arbitrary functions y=f(x).
# Samples from the cross product of all focals
#
# Input:
# fhandle: Function handle or function name.
# x: Array of function parameters.
# mcIT: Number of Monte-Carlo samples
# optimizer (optional): Optimizer to use (default: dsbound)
# corr (optional): Correlation coefficients for Gaussian copula
#
# Output:
# y: Result of propagating x through f(x)
#
# Usage:
# f1=(function(x){ z <- sqrt(abs(sum(x)));});
# f2=(function(x){ z <- sin(sum(x));});
# mu1=dsstruct(c(2,3,1));
# mu2=dsstruct(c(4,4,1));
# dss1=dsodf('qnorm',200,mu);
# dss2=dsodf('qnorm',500,mu);
# correlations=0.5;
# y_f1_independent=dsevalmc(f1,list(dss1,dss2),1000)
# y_f1_dependent=dsevalmc(f1,list(dss1,dss2),1000,'dsbound',correlations)
# y_f2=dsevalmc(f2,list(dss1,dss2),1000,'dsopt')
# dscdf(y_f1_independent)
# dscdf(y_f1_dependent)
# dscdf(y_f2)
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
require(copula);
require(AlgDesign);
if (is.list(x[[1]])){
    mult<-length(x);}
else{
    mult<-1;
    x<-list(x);
}
len=vector("numeric",mult);
for (i in 1:mult){
    len[i]=length(x[[i]]);
}
usedependency=1;
#if nargin<4 | isempty(optimizer)
#    optimizer='dsbound';
#end
cop=numeric();
if (is.null(corr)){
    usedependency=0;
} else {
cop <- normalCopula(corr,max(len),dispstr='un');
}
#if nargin<6
    isvectorfun=1;
#end
sample=matrix(0,len[1],mult);
for (j in 1:mult){
    for (i in 1:len[j]){	
        sample[i,j]=dim(x[[j]][[i]])[1];
}
}
for (k in 1:mult){
    for (i in 1:len[k]){
        temp2<-x[[k]][[i]];
ptemp=matrix(0,dim(temp2)[1],5);
ptemp[,1:3]=temp2;
       a=sort(ptemp[,1]+ptemp[,2]);
b=as.matrix(order(as.matrix(ptemp[,1])+as.matrix(ptemp[,2])));
        ptemp[,]=as.matrix(ptemp[b,]);
        ptemp[1,4]=ptemp[1,3];
        ptemp[1,5]=0;
if (dim(ptemp)[1]>1){
        for (j in 2:(dim(ptemp)[1])){
            ptemp[j,4]=ptemp[j-1,4]+ptemp[j,3];
            ptemp[j,5]=ptemp[j-1,4];
        }
}
        x[[k]][[i]]=ptemp;
}

if (usedependency==0){
randnos=matrix(runif(mcIT*max(len)),nrow=mcIT);
}
else
{
    randnos=matrix(0,mcIT,max(len));
    stepsize=min(mcIT, 10000);
    for (i in seq(1,mcIT,stepsize))
            randnos[i:min(mcIT, i+stepsize-1),]=rcopula(cop,stepsize);
}
}
y=list();
toEvalVec=list();
for (nP in 1:mult){
    tempEvalVec=matrix(1,mcIT,len[nP]);
    for (i  in 1:len[nP]){
#if(is.null(samples)==TRUE){
#samples=matrix(NA,0,len[nP]+1)
#}
rvec=matrix(NA,length(randnos[,i])+1,2);
        rvec[,1]=as.matrix(c(randnos[,i],2000));
	  #rvec=cbind(rvec,as.matrix(1:length(rvec)));
rvec[,2]=1:dim(rvec)[1];	  	  
rvec<-rvec[order(rvec[,1]),];
        lb=x[[nP]][[i]][,5];
        ub=x[[nP]][[i]][,4];
        count=1;
        for (j in 1:length(lb)){
            for (k in count:(mcIT+1)){
                if (rvec[k,1]>=ub[j]){
                    count<-k;
                    break;
                }
     tempEvalVec[rvec[k,2],i]<-j;
            }
		}
    }
        toEvalVec[[nP]]=tempEvalVec;
}

for (nP in 1:mult){
    lb=matrix(0,mcIT,len[nP]);
    ub=matrix(0,mcIT,len[nP]);
    mass=matrix(0,mcIT,len[nP]);
    for (i in 1:len[nP]){
temp=x[[nP]][[i]];
vec=toEvalVec[[nP]];
        lb[,i]=x[[nP]][[i]][toEvalVec[[nP]][,i],1];
        ub[,i]=x[[nP]][[i]][toEvalVec[[nP]][,i],2];
        mass[,i]=x[[nP]][[i]][toEvalVec[[nP]][,i],3];
    }
    if (isvectorfun==1){
if(is.null(fnoptions)){
    temp=optimizer(fhandle, lb,ub,mass,samples=samples);
} else {
    temp=optimizer(fhandle, lb,ub,mass,fnoptions,samples=samples);
}
            ytemp=temp[[1]];
samples=temp[[2]];
}
    else{
        ytemp=matrix(0,mcIT,3);
        for (i in 1:mcIT){
if(is.null(fnoptions)){
temp=optimizer(fhandle, lb[i,],ub[i,],mass[i,],samples=samples)

} else {
temp=optimizer(fhandle, lb[i,],ub[i,],mass[i,],fnoptions,samples=samples)
}
            ytemp[i,]=temp[[1]];
samples=temp[[2]];
}
}
    ytemp[,3]=1/mcIT;
temp=dsnorm(ytemp);
    y[[nP]]=temp$ds;
}
if (length(y)==1)
{
    y=y[[1]];
}
dsevalmc<-list(y,samples);
}

