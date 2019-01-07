#' @export
`dsecdffit` <-
function(x){
# Generates BPA from empirical CDF bounds on data set x (points, intervals)
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#========================================================================

if (is.matrix(x)==FALSE){
x=matrix(x);
}
if (dim(x)[2]==1){
x=cbind(x,x);
}
erg=matrix(NA,dim(x)[1],3);
erg[,c(1,2)]=x
erg[,3]=1/(dim(x)[1]);
erg
}

