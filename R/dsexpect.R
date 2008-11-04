dsexpect <-
function(x){
# Expected value interval
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
dsexpect=c(0,0);
dsexpect[1]=sum(cbind(x[,1]*x[,3]),na.rm = TRUE);
dsexpect[2]=sum(cbind(x[,2]*x[,3]),na.rm = TRUE);
dsexpect
}

