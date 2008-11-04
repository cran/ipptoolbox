`dsbelpltests` <-
function(x,a){
# Belief and Plausibility of a set of points a respective to a BPA x.
#=========================================================================
# Reference : Ferson, S., V. Kreinovich, et al. (2003). Constructing
# Probability Boxes and Dempster-Shafer Structures. Albuquerque, Sandia
# National Laboratories.
# Link      : http://citeseer.ist.psu.edu/660030.html
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

xbel=matrix(NA,dim(x)[1],2);
xpl=xbel;
xbel[,c(1,2)]=x[order(x[,1]),c(1,3)];
xpl[,c(1,2)]=x[order(x[,2]),c(2,3)];
asorted=sort(a);

sum=0;
count=1;
n=length(asorted);
nds=dim(xbel)[1];
erg=matrix(1,n,2);
for (i in 1:nds){
while (count<=n && asorted[count]<xbel[i,1]){
erg[count,2]=sum;
count=count+1;
}
sum=sum+xbel[i,2];
}
count=1;
sum=0;
for (i in 1:nds){
while (count<=n && asorted[count]<xpl[i,1]){
erg[count,1]=sum;
count=count+1;
}
sum=sum+xpl[i,2];
}
erg=erg[order(order(a)),];
erg
}

