#' @importFrom grDevices dev.new
#' @importFrom graphics barplot legend lines title
#' @export
`dscdf` <-
function(x,xrange=NULL,col=c(3,4),...,newplot=TRUE){
# Plots Bel(-inf,x) and Pl(-inf,x), the analogon to a cdf of x.
#=========================================================================   
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

if(is.matrix(x)==FALSE){
x=matrix(x,ncol=3)
}

#Calculate Bel(x) and Pl(x)
down=dspl(x)
up=dsbel(x)

#Prepare stair plot
xa<-rbind(as.matrix(down[1,1]),as.matrix(down[,1]),as.matrix(up[dim(up)[1],1]));
xb<-rbind(as.matrix(0),as.matrix(down[,2]),as.matrix(1));
xc<-rbind(as.matrix(down[1,1]),as.matrix(up[,1]));
xd<-rbind(as.matrix(0),as.matrix(up[,2]));

#Prepare plot bounds on x-axis
if (is.null(xrange)){
xrange=c(min(xa[xa>-Inf]),max(xc[xc<Inf]));
xrange[1]=xrange[1]-1/20*(xrange[2]-xrange[1])
xrange[2]=xrange[2]+1/20*(xrange[2]-xrange[1])
}
xa[xa==-Inf]=xrange[1]-1/20*(xrange[2]-xrange[1])
xc[xc==-Inf]=xrange[1]-1/20*(xrange[2]-xrange[1])
xa[xa==Inf]=xrange[2]+1/20*(xrange[2]-xrange[1])
xc[xc==Inf]=xrange[2]+1/20*(xrange[2]-xrange[1])

#Plot plausibility and Belief
if(newplot==TRUE){
plot(xa,xb,type="s",col=col[2],lty=1,xlim=xrange,lwd=3,...);
} else {
lines(xa,xb,type="s",col=col[2],lty=1,xlim=xrange,lwd=3,...);
}
lines(xc,xd,type="s",col=col[1],lwd=2,lty=1);
legend("bottomright",c("Bel","Pl"),col = col,lwd=3,lty = c(1, 1),merge = TRUE, bg = 'gray90')
erg=list();
erg$Bel=cbind(xa,xb)
erg$Pl=cbind(xc,xd)
erg=erg
}

