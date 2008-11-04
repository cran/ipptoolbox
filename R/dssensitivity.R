dssensitivity <-
function(x, parnums, fhandle, uncfn, mcIT, pinch_samples, pinch_type=NULL,optimizer=dsmonotonous,fnoptions=NULL){
# Sensitivity measure routine based on Monte-Carlo sampling.
#=========================================================================   
# [euncquot]=dssensitivity(x, parnums, fhandle, uncfn, mcIT, pinch_samples, pinch_type) 
# 
# Sensitivity measure routine based on Monte-Carlo sampling. Performs a
# sensitivity analysis on the function f(x) using an uncertainty measure "uncfun". Sensitivity is
# measured by reducing the input variable uncertainty and detecting the
# change in the overall uncertainty.
# Input uncertainty is reduced by iterative, random "pinching", replacing an input x_i by a
# random point/interval/distribution enclosed by x_i.
#
# Input:
# x: Input vector of Demster-Shafer structures
# parnums: Indices of elements in x that are the objective of the
# sensitivity analysis
# fhandle: Function handle of f(x)
# uncfun: Uncertainty measure (i. e. dsaggunc, dsavgwidth, ...)
# mcIT: Number of MC samples per pinch_sample input
# pinch_samples: Number of pinch_samples
# pinch_type (optional): Type of pinching, pinch_type='point' for sampling points in
# x_i (default), pinch_type='interval' for sampling intervals annd
# pinch_type='distribution' for distributions.
#
# Output:
# euncquot: Quotient of uncertainty in f(x) between pinched and unpinched
# input x_i
#
# Usage - Determining the sensitivity of the average width and of the
# dissonance measure:
# myfn=inline('sqrt(x(:,1).^2.+abs(x(:,2)))','x')
# mu=dsstruct([2,3,1])
# dss1=dsodf('norminv',200,mu)
# dss2=dsodf('norminv',500,mu)
# [a]=dssensitivity([dss1,dss2],[1,2],myfn,'dsavgwidth',200,100);
# bar(100*(a));
# xlabel('Input');
# ylabel('S in # - dsavgwidth');
# [a2]=dssensitivity([dss1,dss2],[1,2],myfn,'dsdissonance',200,100,'interval');
# figure
# bar(100*(a2));
# xlabel('Input');
# ylabel('S in # - dsdissonance');
#=========================================================================
# Reference : Ferson S, Tucker WT: Sensitivity in Risk Analyses with
# Uncertain Numbers. Sandia National Laboratories, Albuquerque (2006).
# Link      : http://citeseer.ist.psu.edu/660030.html
# Reference : Limbourg P, Savic R, Petersen J, Kochs H-D: Fault Tree
# Analysis in an Early Design Stage using the Dempster-Shafer Theory of
# Evidence. In: Terje Aven JEV (ed) European Conference on Safety and
# Reliability – ESREL 2007, pp. 713-722. Taylor and Francis, Stavanger, Norway (2007).
# Link      :
# http://iit.uni-duisburg.de/forschung/veroeffentlichungsdateien/2007/LSKP07.pdf
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
uncy=0;
uncgds=0;
uncraw=0;
euncpar=matrix(NA,length(x),1);
euncraw=matrix(NA,length(x),1);
euncnorm=matrix(NA,length(x),1);
euncquot=matrix(NA,length(x),1);
if (is.null(pinch_type)){
    pinch_type='point';
}
for (parnum in parnums){
dss=x;
dsst=x[[parnum]];
mcsel=matrix(0,dim(dsst)[1],1);
mcsel[1]=dsst[1,3];
unc=matrix(0,pinch_samples,1);
if (dim(dsst)[1]>1){
for (i in 2:dim(dsst)[1]){
	mcsel[i]=mcsel[i-1]+dsst[i,3];
}
}
for (i in 1:pinch_samples){
    if (pinch_type == 'point' | pinch_type == 'interval'){
    numsamp=1;
ind1=matrix(NA,numsamp,1);
    for (j in 1:numsamp) {
    ind1[j]=sum(mcsel>=runif(1,0,1));
    }
foc=matrix(NA,numsamp,3);
    tempfoc=dsst[ind1,];
    if (pinch_type=='point'){
    foc[,1]=runif(numsamp,tempfoc[1],tempfoc[2]);
    foc[,2]=foc[,1];
    }
    if (pinch_type=='interval'){
    foc[,1]=tempfoc[1];
    foc[,2]=tempfoc[2];
    }
  foc[,3]=tempfoc[3]/numsamp;
    }
    if (pinch_type =='distribution'){    
    foc=dsst;        
    foc[,1]=runif(dim(dsst)[1],dsst[,1],dsst[,2]);
    foc[,2]=foc[,1];        
    }
    dss[[parnum]]=foc;
    ytemp=dsevalmc(fhandle,list(x,dss),mcIT,optimizer,fnoptions=fnoptions)[[1]];
    gds=ytemp[[1]];
    y=ytemp[[2]];
    uncy[i]=uncfn(y);
    uncgds[i]=uncfn(gds);
    uncraw[i]=uncgds[i]-uncy[i];
#print(c(uncy[i],uncgds[i],uncy[i]/uncgds[i]));
}
euncpar[parnum]=uncfn(dsst);
uncnorm=uncraw/euncpar[parnum];
euncraw[parnum]=mean(uncraw);
euncnorm[parnum]=mean(uncnorm);
euncquot[parnum]=1-mean(uncy[i])/mean(uncgds[i]);
#uncquot=100*(1-median(uncquot));
if(is.na(euncquot[parnum])){
print("Sensitivity index not computable, set to zero for parameter");
print(parnum);
euncquot[parnum]=0;
}
}
euncquot;
}

