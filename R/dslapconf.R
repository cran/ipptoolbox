`dslapconf` <-
function(x,lims=c(-Inf,Inf)){
# Generates BPA from Kolmogorov-Smirnov bounds on data set x (points, intervals)
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

if (is.matrix(x)==FALSE){
x=matrix(x);
}

# If x is consisting of point values, convert to zero-width intervals
if (dim(x)[2]==1){
x=cbind(x,x);
}

# Create n+1 focal elements that are divided by the data
erg=matrix(NA,dim(x)[1]+1,3);
erg[,1]=c(lims[1],sort(x[,1]))
erg[,2]=c(sort(x[,2]),lims[2])
erg[,3]=1/(dim(x)[1]+1);
erg
}

