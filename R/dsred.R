#' @export
`dsred` <-
function(y,thres=0.001){
# Reduces the size of a Dempster-Shafer structure
#=========================================================================   
# y=dsred(x,thres)
# Reduces the size of a Dempster-Shafer structure x by merging focal
# elements.
#
# Input: 
# x: Dempster-Shafer structure to be reduced
# thres (optional): 0<thres<1, minimal mass of focals. Default: 0.001
#
# Output:
# y: Reduced Dempster-Shafer structure
#
# Usage:
# lambda=dsstruct(c(2,3,1))
# dss=dsodf('qexp',10000,lambda)
# y=dsred(dss,0.05)
# dscdf(dss);
# dscdf(y);
# =========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
y[,]=y[order(y[,1]),];
down=0;
up=0;
siz=nrow(y);
last=0;
step=siz*thres;
if (step < 1){
    return(y)
}
i=1;
num=0;
imin=0;
imax=0;
erg=matrix(0,0,3);    
while (i<=siz){
    imin=i;
    mass=0;
	while (i<=siz && mass<thres)    {
        mass=mass+y[i,3];
        imax=i;
        i=i+1;
}
        down=min(y[imin:imax,1]);
        up=max(y[imin:imax,2]);
erg<-rbind(erg,t(c(down,up,mass)));
}
dsred=erg;
}

