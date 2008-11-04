`dskstest` <-
function(val,ds,...){
# Expected value interval
#=========================================================================   
# dsexpect<-function(x) 
# 
# Expected value of the Dempster-Shafer structure x
#
# Input:
# x: Dempster-Shafer structure
#
# Output:
# y: Expected value (lower & upper bound)
#
# Usage:
# lambda1=dsstruct(c(2,3,1))
# dss=dsodf('qexp',1000,lambda1);
# y=dsexpect(dss)
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

   pkstwo <- function(x, tol = 1e-06) {
        if (is.numeric(x)) 
            x <- as.vector(x)
        else stop("argument 'x' must be numeric")
        p <- rep(0, length(x))
        p[is.na(x)] <- NA
        IND <- which(!is.na(x) & (x > 0))
        if (length(IND) > 0) {
            p[IND] <- .C("pkstwo", as.integer(length(x[IND])), 
                p = as.double(x[IND]), as.double(tol), PACKAGE = "stats")$p
        }
        return(p)
    }

if (is.matrix(val)==FALSE){
val=matrix(val);
}
TIES <- FALSE
kss=matrix(NA,dim(val)[1],dim(val)[2]);
for (i in 1:dim(val)[2]){
xunsorted=val[,i];
x=sort(xunsorted);
n=length(x);
if(is.numeric(ds)){
	fx=dsbelpltests(ds,x);
} else {
pars=(...);
xlo=matrix(NA,length(x),length(pars)+1);
xhi=xlo;
xlo[,1]=x;
xhi[,1]=x;
for (j in 1:length(pars)) {
if (length(pars[[j]])<2) {
pars[[j]]=c(pars[[j]],pars[[j]]);
}
		xlo[,j+1]=pars[[j]][1]
		xhi[,j+1]=pars[[j]][2]
}
	if (length(pars)==1){
		mycdf=function(x){erg=ds(x[,1],x[,2])}
}

	else if (length(pars)==2){
		mycdf=function(x){erg=ds(x[,1],x[,2],x[,3])}
} else {
		mycdf=function(x){erg=ds(x[,1],x[,2],x[,3],x[,4])}
}
mass=xlo;
fx=dsmonotonous2(mycdf,xlo,xhi,mass)[[1]]
}
#ks=pmin(abs(c(1:n)/n-fx[,1]),abs(c(1:n)/n-fx[,2]))
#ks2=pmin(abs(c(0:(n-1))/n-fx[,1]),abs(c(0:(n-1))/n-fx[,2]))
ks=-pmin(c(0:(n-1))/n-fx[,1],fx[,2]-c(1:n)/n,0)
#ks=pmax(ks,ks2);
#ks[fx[,1]<(c(1:n)/n) & fx[,2]>(c(0:n-1)/n)]=0
kss[,i]=ks[order(order(xunsorted))];
}
if(dim(val)[2]>1){
kss=pmax(kss[,1],kss[,2]);
}
ks=max(kss)

if(is.numeric(ds)){
         PVAL <- 1 - pkstwo(sqrt(n * dim(ds)[1]/(n + dim(ds)[1])) * ks)
} else {
PVAL <- 1 - pkstwo(sqrt(n) * ks)
}
PVAL
}

