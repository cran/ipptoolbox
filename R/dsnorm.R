`dsnorm` <-
function(y){
# Normalizes a Dempster-Shafer structure
#=========================================================================   
# dsnorm <- function(y) 
# Normalizes a Dempster-Shafer structure x. Normalizes masses to 1 and
# interchanges lower_bound & upper_bound, if lower_bound>upper_bound.
#
# Input: 
# x: Dempster-Shafer structure to be normalized
#
# Output:
# y: Normalized Dempster-Shafer structure
#
# Usage:
# foo=dsstruct(c(3,2,1))
# bar=dsnorm(foo)
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
vec<-y[,1]>y[,2];
a<-0;
b<-0;
if (sum(vec)>0){
    a<-1;
}
temp<-y[which(vec),1];
y[which(vec),1]<-y[which(vec),2];
y[which(vec),2]<-temp;
s<-sum(y[,3]);
if (s!=1)
    b<-1;
y[,3]<-y[,3]/s;
dsnorm=list(ds=y,a=a,b=b);
}

