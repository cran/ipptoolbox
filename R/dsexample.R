`dsexample` <-
function(){
# Example file calculating the reliability of a three-component series-parallel system
#=========================================================================   
# y=dsexample
# This is an exemplary use of the toolbox.
# Calculate the reliability of a simple three-component
# series-parallel system with dependent inputs.
# Please use "type dsexample" or open with an editor
# 
#=========================================================================
# Copyright (c) Philipp Limbourg, University of Duisburg-Essen
# www.uni-duisburg-essen.de/informationslogistik/
#=========================================================================
print('This is an example.');
print('Calculate the reliability of a simple two-component series system with dependent inputs.');
print('1. Define system function. Components 2 and 3 are in parallel:');

myfn=(function(x){ z <- pmin(x[,1],pmax(x[,2],x[,3]));z});

print('2. Define component reliability.');
print('Component 1 fails with an exponential cdf:');
print('Define lambdas');
lambda1=dsstruct(c(1/12000,1/10000,1));

print('Sample cdf using outer discretization, 1000 samples');

component1=dsodf('qexp',1000,lambda1);

print('Component 2 fails with a Weibull cdf:');

a2<-dsstruct(c(5000,10000,1));

print('The b value has two focal elements:');

b2<-dsstruct(matrix(c(1.2,1.5,0.1,1.1,1.8,0.9),2,byrow=TRUE));
component2=dsodf('qweibull',1000,b2,a2);

print('Component 3 fails with a Weibull cdf. Unfortunately we have two estimates');

a3_1=dsstruct(c(2000,5000,1));b3_1=dsstruct(matrix(c(1.2,1.5,0.8,1.1,1.8,0.2),2,byrow=TRUE));
a3_2=dsstruct(c(1200,2500,1));b3_2=dsstruct(c(1.7,1.9,1));
component3_1=dsodf('qweibull',1000,b3_1,a3_1);
component3_2=dsodf('qweibull',1000,b3_2,a3_2);

print('Now we aggregate the two estimates using weighted mixing. Expert 1 is twice as competent as Expert 2');

w=c(2,1);
component3=dswmix(component3_1,component3_2,w=w);

print('3. Now we propagate the uncertain variables through our system function');

print('3. The propagation using independence assumption');
y=dsevalmc(myfn,list(component1,component2,component3),10000,optimizer=dsbound)[[1]];

#'Setting up the correlation matrix. Component 2 and 3 have a positive correlation';

corr=c(0,0,0.5);

#'Now start the propagation. The optimizer should be dsbound as the function is monotonously increasing. 1000 iterations, time to get a coffee (-;';

ycorr=dsevalmc(myfn,list(component1,component2,component3),10000,optimizer=dsbound,corr)[[1]];


print('Now plot the system function');

dscdf(y);

print('Deriving some statistics');
statsy=list()
statsy$exp=dsexpect(y);
statsy$med=dsconf(y,0.5);
statsy$conf95=dsconf(y,0.95);
statsy$aggWidth=dsaggwidth(y);
print(statsy);
print('Sensitivity');
a=dssensitivity(list(component1,component2,component3),c(1,2,3),myfn,dsaggwidth,20,20);
print("s1,s2,s3:")
print(a)
}

