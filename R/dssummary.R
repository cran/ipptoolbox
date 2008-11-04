`dssummary` <-
function(x){
# Summary stats of a BPA
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
#Function for calculating quantile
erg=list()
erg$Q01=dsconf(x,0.01)
erg$Q05=dsconf(x,0.05)
erg$Q25=dsconf(x,0.25)
erg$Med=dsconf(x,0.5)
erg$Q75=dsconf(x,0.75)
erg$Q95=dsconf(x,0.95)
erg$Q99=dsconf(x,0.99)
erg$E=dsexpect(x)
erg$Sd=sqrt(dsvariance(x))
erg=t(as.data.frame(erg))
erg
}

