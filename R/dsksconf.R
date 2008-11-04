`dsksconf` <-
function(x,conf=0.95,lims=c(-Inf,Inf)){
# Generates BPA from Kolmogorov-Smirnov bounds on data set x (points, intervals)
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================

# Function to determine KS alpha from confidence level
getalpha2=function(n,alpha){
if(n<200){
f=function(ks,n,conf){.C("pkolmogorov2x", p = as.double(ks), as.integer(n), PACKAGE = "stats")$p-conf}
} else {
f=function(ks,n,conf){pkstwo(sqrt(n)*ks)-conf}
}
alpha=uniroot(f,c(0,1),n,alpha)$root

print(alpha)
alpha
}

getalpha=function(n,alpha){
# Matrix for alpha values depending on sample size and confidence
amat=matrix(c(
1	,	.900	,	.925	,	.950	,	.975	,	.995	,
2	,	.684	,	.726	,	.776	,	.842	,	.929	,
3	,	.565	,	.597	,	.642	,	.708	,	.828	,
4	,	.494	,	.525	,	.564	,	.624	,	.733	,
5	,	.446	,	.474	,	.510	,	.565	,	.669	,
6	,	.410	,	.436	,	.470	,	.521	,	.618	,
7	,	.381	,	.405	,	.438	,	.486	,	.577	,
8	,	.358	,	.381	,	.411	,	.457	,	.543	,
9	,	.339	,	.360	,	.388	,	.432	,	.514	,
10	,	.322	,	.342	,	.368	,	.410	,	.490	,
11	,	.307	,	.326	,	.352	,	.391	,	.468	,
12	,	.295	,	.313	,	.338	,	.375	,	.450	,
13	,	.284	,	.302	,	.325	,	.361	,	.433	,
14	,	.274	,	.292	,	.314	,	.349	,	.418	,
15	,	.266	,	.283	,	.304	,	.338	,	.404	,
16	,	.258	,	.274	,	.295	,	.328	,	.392	,
17	,	.250	,	.266	,	.286	,	.318	,	.381	,
18	,	.244	,	.259	,	.278	,	.309	,	.371	,
19	,	.237	,	.252	,	.272	,	.301	,	.363	,
20	,	.231	,	.246	,	.264	,	.294	,	.356	,
25	,	.210	,	.220	,	.240	,	.270	,	.320	,
30	,	.190	,	.200	,	.220	,	.240	,	.290	,
35	,	.180	,	.190	,	.210	,	.230	,	.270	
),ncol=6,byrow=TRUE)
alim=c(1.07,1.14,1.22,1.36,1.63)
# Vector for confidence levels
aconf=c(.80, .85, .90, .95, .99 )
if(length(which(aconf>=alpha))==0){
print("Confidence interval too large");
}
confind=min(which(aconf>=alpha))
# If sample number exceeds entries in table, use analytical approximation for large n
if(length(which(amat[,1]>=n))==0){
alpha=alim[confind]/sqrt(n);
}else{
# Else table lookup
tabind=min(which(amat[,1]>=n))
alpha=amat[tabind,confind+1]
}
alpha
}

pkstwo <- function(x, tol = 1e-06) {
        if (is.numeric(x)) 
            x <- as.vector(x)
        else stop("argument 'x' must be numeric")
        p <- rep(0, length(x))
        p[is.na(x)] <- NA
        IND <- which(!is.na(x) & (x > 0))
        if (length(IND) > 0) {
            p[IND] <- .C("pkstwo", as.integer(length(x[IND])), 
                p = as.double(x[IND]), as.double(tol), PACKAGE = "stats")$p
        }
        return(p)
    }


# Function that converts the KS bounds for the data to a BPA
convert2ds=function(x){
x=unique(x);
x=rbind(c(-Inf,0,0),x);
xs=sort(c(x[,2],x[,3]));
xs=unique(xs);
erg=matrix(NA,length(xs)-1,3);
for (i in 2:length(xs)){
erg[i-1,1]=x[1+max(which(x[,2]<xs[i])),1]
erg[i-1,2]=x[1+max(which(x[,3]<xs[i])),1]
erg[i-1,3]=xs[i]-xs[i-1]
}
erg
}


if (is.matrix(x)==FALSE){
x=matrix(x);
}
# If x is not in interval form, create zero-width intervals
if (dim(x)[2]==1){
x=cbind(x,x);
}
x[,1]=sort(x[,1]);
x[,2]=sort(x[,2]);
# First value starts with -Inf, last ends with Inf
xcor=sort(c(-Inf,x));
erg=cbind(xcor,xcor,xcor);
n=dim(x)[1];

# Get alpha value from confidence level
a=getalpha2(n,conf);
# Generate upper and lover bounds for data values
for (i in 1:length(xcor)){
erg[i,2]=min(1,(sum(x[,1]<=xcor[i])/n)+a);
erg[i,3]=max(0,(sum(x[,2]<=xcor[i])/n)-a);
}
erg[length(xcor),]=c(Inf,1,1)
# Generate BPA from these bounds
erg=convert2ds(erg)

# Replace infinite values with bounds
erg[erg==-Inf]=lims[1]
erg[erg==Inf]=lims[2]
erg
}

