`dsqqplot` <-
function(ds,sample,points=FALSE){
# Plots quantile-quantile plots for a BPA ds and a set of points sample.
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

#Calculate the empirical CDF of the sample
p=dsecdf(sample)

#Calculate Bel(x) and Pl(x) of the sample
p2=dsbelpltests(ds,sort(sample))

x=p[,1];
y=p[,2];

#Plot Bel(x) and Pl(x) respective to quantiles
plot(y,p2[,1],'l',col=4,lwd=2,xlab="Experimental data quantiles",ylab="Theoretical quantiles (Bel/Pl)")
lines(y,p2[,2],'l',col=4,lwd=2)

#Plot empirical quantiles and points
if(points==TRUE){
points(y,y,col=3,lwd=1,pch=20)
legend("bottomright",c("Theor. quantiles (Bel/Pl)"),col = c(4),lty = c(1),merge = TRUE, bg = 'gray90')
} else {
lines(y,y,'l',col=3,lwd=2,lty=2)
legend("bottomright",c("Exp. data quantiles","Theor. quantiles (Bel/Pl)"),col = c(3,4),lty = c(1, 1),merge = TRUE, bg = 'gray90')
}
}

