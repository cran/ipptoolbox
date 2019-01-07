#' @export
`dsecdf` <-
function(x){
# Returns empirical CDF of a sample x or set of empirical CDFs of an interval sample x
#=========================================================================   
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

#If x is a set of intervals, calculate empirical CDF bounds
if (is.matrix(x) &&dim(x)[2]==2){
xlo=sort(x[,1]);
xhi=sort(x[,2]);
y=c((1:dim(x)[1])/dim(x)[1]);
erg=cbind(xlo,xhi,y);

}else{
#Calculate empirical CDF
x=sort(x);
y=c((1:length(x))/length(x));
erg=cbind(x,y);
}
}

