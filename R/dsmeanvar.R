#' @export
`dsmeanvar` <-
function(intervalnumber, mean, var){
# Generates a BPA enclosing all CDFs with given mean and variance
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

#Generate number of probability values
p=seq(0.5/intervalnumber,1,1/(intervalnumber));

#Minimum and maximum deviation from mean of the CDFs for probability p
d1=sqrt(var/(p+p^2/(1-p)))
#Add/substractg this to/from the mean
lo=mean-d1;
hi=mean+d1*p/(1-p);

#Form BPA
erg=dsstruct(cbind(lo,hi,1/intervalnumber));
}

