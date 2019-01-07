#' @export
`dsstruct` <-
function(x){
# geometric mean calculation
# Returns a new Dempster-Shafer structure
#=========================================================================   
# y=dsstruct(x) 
#
# Input: 
# x (optional): A N-3 array containing [lower_bound,upper_bound,mass]
# 
# Output:
# y: a new Dempster-Shafer structure
# Example: 
# 
# dss=dsstruct(matrix(c(2,3,0.5,3,4,0.5),ncol=3,byrow=TRUE))
# 
# If the sum of masses is not 1, masses will be normalized. If
# lower_bound>upper_bound, they will be interchanged.
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
#if (nargin>0)
if (is.matrix(x)==FALSE){
x=matrix(x,ncol=3,byrow=TRUE);
}
ds<-as.matrix(x);
y=dsnorm(ds);
ds=y$ds
if (y$a>0)
print('Warning, lower bound > upper bound for some elements: BPA normalized')
if (y$b>0)
print('Warning, mass <>1: BPA normalized')
dsstruct<-ds
}

