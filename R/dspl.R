dspl <-
function(x){
# Plausibility function
#=========================================================================   
# y=dspl(x) 
#
# Input:
# x: Dempster-Shafer structure. 
#
# Output:
# y: Plottable plausibility function [x,Bel(-inf,x)]
#
# Example:
#
# lambda1=dsstruct(c(2,3,1));
# dss=dsodf('qexp',1000,lambda1);
# y=dspl(dss);
# plot(y[,1],y[,2],type='l')
#=========================================================================
# Reference : Ferson, S., V. Kreinovich, et al. (2003). Constructing
# Probability Boxes and Dempster-Shafer Structures. Albuquerque, Sandia
# National Laboratories.
# Link      : http://citeseer.ist.psu.edu/660030.html
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
y<-x[,c(1,3),drop=FALSE];
y<-y[order(y[,1]),,drop=FALSE];
if(dim(x)[1]>1){
for (i in 2:dim(x)[1]){
    y[i,2]=y[i-1,2]+y[i,2];
}
}
dspl=y;
}

