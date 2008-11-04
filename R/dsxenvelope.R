dsxenvelope <-
function(...,maxfocals=10000000){
# XEnveloping of BPAs "..." with optional weights w
#=========================================================================
# Reference : Ferson, S., V. Kreinovich, et al. (2003). Constructing
# Probability Boxes and Dempster-Shafer Structures. Albuquerque, Sandia
# National Laboratories.
# Link      : http://citeseer.ist.psu.edu/660030.html
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

# List of BPAs to aggregate
x=list(...)
erg=x[[1]];
# Calculate maximal input BPA length from maximal output BPA length
acc=maxfocals^(1/length(x))
lbs=list();
ubs=list();
ms=list();

# Reduce BPA in case of exceeding length
for (i in 1:length(x)){
if(length(x[[i]])[1]>acc){
        x[[i]]=dsred(x[[i]],1/acc);
}
lbs=c(lbs,list(x[[i]][,1]))
ubs=c(ubs,list(x[[i]][,2]))
ms=c(ms,list(x[[i]][,3]))
}

# Expand grid for creating cross product
tmplo=expand.grid(lbs)
tmphi=expand.grid(ubs)
tmpms=expand.grid(ms)
erg=matrix(NA,dim(tmplo)[1],3)
erg[,1]=Inf
erg[,2]=-Inf
erg[,3]=1

# Calculate envelope for each cross product
# element (using pmax instead of apply to optimize speed)
for (i in 1:length(x)){
erg[,1]=pmin(erg[,1],tmplo[,i])
erg[,2]=pmax(erg[,2],tmphi[,i])
erg[,3]=erg[,3]*tmpms[,i]
}

#Normalize BPA
erg=dsnorm(erg)$ds;
}

