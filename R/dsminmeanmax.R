dsminmeanmax <-
function(intervalnumber,min,mean,max){
# Generates a BPA enclosing all CDFs with given min, mean and max
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
p=seq(0,1,1/(intervalnumber-1));
lo=max-(max-mean)/p;lo[lo<min]=min;
hi=(p*min-mean)/(p-1);hi[hi>max]=max;hi[hi==-Inf]=max
erg=dsstruct(cbind(lo,hi,1/intervalnumber))
}

