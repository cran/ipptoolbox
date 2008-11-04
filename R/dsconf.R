dsconf <-
function(x,conf,confconf=NULL){
# Quantile value of BPA x at level conf, optional with Wilk's bounds confconf
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

#Function for calculating quantile
myconf<-function(x,conf){
down=cbind(x[,1],x[,3]);
down<-down[order(down[,1]),];
down[,2]=cumsum(down[,2]);
up=cbind(x[,2],x[,3]);
up<-up[order(up[,1]),];
up[,2]=cumsum(up[,2]);
y=numeric();
y[1]=down[min(which(down[,2]>=conf-1E-12)),1];
y[2]=up[min(which(up[,2]>=conf-1E-12)),1];
dsconf=y;
}

#Function for calculating order values for Wilk's bounds
calcwilks<-function(level,n,confidence){
  na <- floor(n*level)
  tmpup <- cumsum(dbinom(0:(n-na),n,1-level))
  tmpdown <- cumsum(dbinom(0:(na-1),n,level))
  indxup <- length(tmpup [tmpup <=(1-confidence)])
  indxdown <- length(tmpdown [tmpdown <=(1-confidence)])
  if(indxup == 0 || indxdown == 0){
print("Warning, too little samples for calculating upper conf");
    return(c(1/n,1))
  }
  rup <- (1+n-indxup )/n
rdown<- indxdown /n 
r=c(rdown, rup)
}

if (is.matrix(x)==FALSE){
x=matrix(x,ncol=3);
}
#Quantile value
y=myconf(x,conf)
if(!is.null(confconf)){
n=dim(x)[1];
#Quantile value for Wilk's bounds
confwilks=calcwilks(conf,n,confconf);
confdown=myconf(x,confwilks[1]);
confup=myconf(x,confwilks[2]);
y=rbind(y,confup,confdown)
}
y

}

