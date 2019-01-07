#' @export
`dsenvelope` <-
function(...){
# Enveloping of BPAs "..."
#=========================================================================
# Reference : Ferson, S., V. Kreinovich, et al. (2003). Constructing
# Probability Boxes and Dempster-Shafer Structures. Albuquerque, Sandia
# National Laboratories.
# Link      : http://citeseer.ist.psu.edu/660030.html
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

x=list(...)
xlo=list();
xhi=list();
n=dim(x[[1]])[1];
m=length(x);
tochecklo=numeric()
tocheckhi=numeric()

# Convert BPA in Pbox
# Calculate masses for new Pbox
for(i in 1:m){
x[[i]]=x[[i]][order(x[[i]][,1]),,drop=FALSE]
xlo[[i]]=cbind(x[[i]],cumsum(x[[i]][,3]));
tochecklo=c(tochecklo,xlo[[i]][,4]);
x[[i]]=x[[i]][order(x[[i]][,2]),,drop=FALSE]
xhi[[i]]=cbind(x[[i]],cumsum(x[[i]][,3]));
tocheckhi=c(tocheckhi,xhi[[i]][,4]);
}
old=0
tocheck=c(tochecklo, tocheckhi)
tocheck=sort(unique(tocheck))
erg=matrix(NA,length(tocheck),3);

#Calculate bounds for input Pboxes
for(i in 1:length(tocheck)){
lo=numeric(m);
hi=numeric(m);
for(j in 1:m){
indlo=min(which(xlo[[j]][,4]+1E-10>=tocheck[i]))
lo[j]=xlo[[j]][indlo,1]
indhi=min(which(xhi[[j]][,4]+1E-10>=tocheck[i]))
hi[j]=xhi[[j]][indhi,2]
}

# Create new focal element by enveloping
erg[i,]=c(min(lo),max(hi),tocheck[i]-old)
old=tocheck[i];
}
#print(erg)
erg=dsstruct(erg)
}

