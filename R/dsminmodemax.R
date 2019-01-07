#' @export
`dsminmodemax` <-
function(intervalnumber,min,mode,max){
# Generates a BPA enclosing all CDFs with given min, mode and max
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
lo=seq(min,mode,(mode-min)/(intervalnumber-1));
hi=seq(mode,max,(max-mode)/(intervalnumber-1))
erg=dsstruct(cbind(lo,hi,1/intervalnumber))
}

