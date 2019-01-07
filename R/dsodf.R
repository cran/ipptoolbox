#' @import AlgDesign
#' @export
`dsodf` <-
function(fhandle,intervalnumber,...){
# Outer discretization of an inverse probability cdf
#=========================================================================
# Reference : Tonon, F. (2004). "Using random set theory to propagate
# epistemic uncertainty through a mechanical system." Reliability
# Engineering and System Safety 85(1-3): 169-181.
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
#require(AlgDesign)
autoload("factorial","AlgDesign");
mythres=(1-1E-15)/floor(intervalnumber);
count=1;
params=list(...)
l=length(params);
for (i in 1:l){
#If params is no BPA, create a BPA
if (length(params[[i]])==1){
params[[i]]=dsstruct(c(params[[i]],params[[i]],1))
}
}
doestr=numeric();
ind=numeric();

#Create number of all parameter combinations
for (i in 1:l){
    siz=dim(params[[i]])[1];
    if (siz>1){
        doestr=cbind(doestr,siz);
        ind=cbind(ind,i);
}
}
tempdoe=numeric();
if (length(doestr)>0){
tempdoe=gen.factorial(doestr,length(doestr),factors="all");
tempdoe=as.matrix(tempdoe);
 tempdoe=apply(tempdoe,2,as.numeric);
}else{
    tempdoe=1;
}
if (max(tempdoe)>1){
tempdoel=dim(tempdoe)[1];
}else{
tempdoel=1;
}
doe=matrix(1,dim(as.matrix(tempdoe))[1],l);
doe[,ind]=as.matrix(tempdoe)[,1:length(ind)];
erg=matrix(0,(intervalnumber)*tempdoel,3);
for (i in 1:dim(doe)[1]){
    prob=1;
    for (j in 1:l){
        prob=prob*params[[j]][doe[i,j],3];
    }
temp=numeric();
    temp[1:l]=2;
    doe2=gen.factorial(temp,length(temp),factors="all");
doe2=as.matrix(doe2)
doe2=apply(doe2,2,as.numeric)
        t1=seq(1E-8,(1-1E-8),mythres);
	t2=append(seq(1E-8+mythres,(1-1E-8),mythres),1-1E-8);
if (is.matrix(doe2)&min(dim(doe2))>1){ 
tempdoe2=dim(as.matrix(doe2))[1];
}else{
tempdoe2=2;
}
borders=matrix(0,length(t1),tempdoe2*2);
    for (k in 1:tempdoe2){
peval=list();
        for (k2 in 1:l){
            peval[[k2]]=params[[k2]][doe[i,k2],doe2[k,k2]];
}
      borders[,2*k-1]=do.call(fhandle,append(list(t1),peval));
borders[,2*k]=do.call(fhandle,append(list(t2),peval));
}
    erg[count:(count+dim(as.matrix(borders))[1]-1),1]=apply(borders,1,min);
    erg[count:(count+dim(as.matrix(borders))[1]-1),2]=apply(borders,1,max);
    erg[count:(count+dim(as.matrix(borders))[1]-1),3]=mythres*prob;
    count=count+dim(as.matrix(borders))[1];
}
dsodf=dsnorm(erg)$ds;
}

