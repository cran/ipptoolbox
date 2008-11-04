`dswmix` <-
function(...,w=NULL) {
# Weighted mixture of BPAs "..." with optional weights w
#=========================================================================
# Reference : Ferson, S., V. Kreinovich, et al. (2003). Constructing
# Probability Boxes and Dempster-Shafer Structures. Albuquerque, Sandia
# National Laboratories.
# Link      : http://citeseer.ist.psu.edu/660030.html
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
erg=numeric()
x=list(...)
if (is.null(w)){
    w=matrix(1,length(x));
}
for (i in 1:length(x)){
    x[[i]][,3]=x[[i]][,3]*w[i];
    erg=rbind(erg,x[[i]]);
}
dswmix=dsnorm(erg)$ds;
}

