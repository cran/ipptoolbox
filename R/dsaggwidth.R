#' @export
`dsaggwidth` <-
function(x){
# Aggregated width of all focal elements of a BPA x (area between Bel(x) and Pl(x))
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
dsaggwidth=sum((x[,2]-x[,1])*x[,3])
}

