dsbelpl <-
function(x,a){
# Belief and Plausibility of an interval
#=========================================================================   
# [bel,pl]=dsbelpl(x,interval) 
# 
# Belief and Plausibility of an interval
#
# Input:
# x: Dempster-Shafer structure. 
# a: Interval a of which to obtain Bel(x in a) and Pl(x in a).
#
# Hints:
# For obtaining Bel/Pl (x < z) call dsbelpl(x,[-inf,z])
# For obtaining Bel/Pl (x > z) call dsbelpl(x,[z,inf])
#
# Output:
# bel: Bel(x in a)
# pl: Pl(x in a)
#
# Example:
#
# lambda1=dsstruct(c(1/3,1/2,1));
# x=dsodf('qexp',1000,lambda1);
# bel_pl=dsbelpl(x,c(5,10)) 
# bel_pl_2=dsbelpl(x,c(-Inf,5))
# bel_pl_3=dsbelpl(x,c(10,Inf))
#=========================================================================
# Reference : Ferson, S., V. Kreinovich, et al. (2003). Constructing
# Probability Boxes and Dempster-Shafer Structures. Albuquerque, Sandia
# National Laboratories.
# Link      : http://citeseer.ist.psu.edu/660030.html
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
erg=c(0,0);
erg[1]=sum(x[(x[,1]>=a[1]&x[,2]<=a[2]),3]);
erg[2]=sum(x[!(x[,1]>a[2]|x[,2]<a[1]),3]);
dsbelpl=erg;
}

