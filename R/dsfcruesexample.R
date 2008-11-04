dsfcruesexample <-
function(){
#Load necessary packages (please install if missing)
library(evd)
library(triangle)

#Data on Q:
Qdonnes=c(3854, 1256, 1649, 1605, 341, 1149, 868, 1148, 1227, 1991, 1255,
1366, 1100, 1837, 351, 1084, 1924, 843, 2647, 1248, 2417, 1125, 
903, 1462, 378, 1230, 1149, 1400, 2078, 1433, 917, 1530, 2442, 
2151, 1909, 630, 2435, 1920, 1512, 1377, 3330, 1858, 1359, 714, 
1528, 1035, 1026, 1127, 1839, 771, 1730, 1889, 3320, 352, 885, 
759, 731, 1711, 1906, 1543, 1307, 1275, 2706, 582, 1260, 1331, 
1283, 1348, 1048, 1348, 383, 1526, 789, 811, 1073, 965, 619, 
3361, 523, 493, 424, 2017, 1958, 3192, 1556, 1169, 1511, 1515, 
2491, 881, 846, 856, 1036, 1830, 1391, 1334, 1512, 1792, 136, 
891, 635, 733, 758, 1368, 935, 1173, 547, 669, 331, 227, 2037, 
3224, 1525, 766, 1575, 1695, 1235, 1454, 2595, 706, 1837, 1629, 
1421, 2204, 956, 971, 1383, 541, 703, 2090, 800, 651, 1153, 704, 
1771, 1433, 238, 122, 1306, 733, 793, 856, 1903, 1594, 740, 3044, 
1128, 522, 642)

Zmdonnes=c(55.09, 55, 54.87, 54.28, 54.74, 55.48, 55.36, 55.39, 54.8, 
55.18, 54.94, 54.42, 55.34, 55.3, 54.31, 55.57, 54.28, 55.49, 
54.49, 55.11, 55.15, 54.4, 55.87, 55.63, 54.93, 55.61, 54.95, 
55.38, 54.57)

Zvdonnes=c(50.39, 50.28, 50.23, 49.92, 50.51, 50.42, 50.16, 50.16, 49.76, 
50.17, 50.71, 50.08, 49.95, 50.63, 49.51, 50.77, 49.98, 50.3, 
50.1, 50.12, 50.54, 49.21, 50.55, 50.67, 50, 50.7, 50.27, 50.06, 
49.49)


Qlimits=c(10,10000)
Kslimits=c(5,60)
Zvlimits=c(48,52.8)
Zmlimits=c(53.2,58)

#Define Pbox with imprecise Gumbel parameter a.
a=dsstruct(c(1000,1040,1));
b=558
#Sample 10000 points from the modified Gumbel distribution
Q=dsadf('qgumbel',10000,a,b)
#Cutting can be performed using Intersection rule
QlimitsBPA=dsstruct(c(Qlimits,1))
Q=dsdempstersrule(Q,QlimitsBPA)
#plot Q
dscdf(Q,xlab="Q",ylab="")


#plot Q versus data in a qq plot
dev.new()
dsqqplot(Q,Qdonnes)
#What is the Kolmogorov-Smirnov probability for the best distribution in Q?
p=dskstest(Qdonnes,Q)
print("K-S probability is:")
print(p)


#Quantify Ks, no distribution, no data. Estimations
#Mean: 30, Min: 5, Max: 60
#5% quantile at 20, 95% quantile at 40
Ks1=dsminmeanmax(1000,5,30,60)
Ks2=dsstruct(rbind(c(5,20,0.05),c(20,40,0.9),c(40,60,0.05)));
dev.new()
dscdf(Ks1,xlab="Ks")
dev.new()
dscdf(Ks2,xlab="Ks")
dev.new()
Ks=dsintersect(Ks1,Ks2)
dscdf(Ks,xlab="Ks",ylab="")

#Quantify Zv
#Laplace fitting from data
Zv=dslapconf(Zvdonnes,Zvlimits)
dev.new()
dscdf(Zv,xlab="Zv")

#Quantify Zm
#Kolmogorov-Smirnov 50% fitting from data
Zm=dsksconf(Zmdonnes,conf=0.5,Zmlimits)
dev.new()
dscdf(Zm,xlab="Zm")

#Define the function fcrue
fcrue=function(x,...){
Q=x[,1];
Ks=x[,2];
Zm=x[,3];
Zv=x[,4];
const=1;
length=5000;
B=300;
result=Zv+(pmax(Q,0)*const/(pmax(Ks,0)*B*sqrt((Zm - Zv)/length)))^0.6;
}

#Evaluate the function, propagate the uncertainty
#If function is nonmonotonous, choose 'dsopt' as internal optimizer (slow)
#If function is monotonously increasing, choose 'dsbound' (faster)
#If function is monotonous, but you don't know if increasing or decreasing, choose 'dsmonotonous' (medium)
#10000 samples
temp=dsevalmc(fcrue,list(Q,Ks,Zm,Zv),10000,dsmonotonous);
Zc=temp[[1]]
dev.new()
dscdf(Zc,xlab="Zc",ylab="",xrange=c(48,65))
print("Median")
print(dsconf(Zc,0.5))
print("Q99 with 95% Wilks bounds")
print(dsconf(Zc,0.99,0.95))
print("Bel/Pl(Zc<=55.5)")
print(dsbelpl(Zc,c(-Inf,55.5)))
print("Bel/Pl(Zc>=55.5)")
print(dsbelpl(Zc,c(55.5,Inf)))
print("Exp. value")
print(dsexpect(Zc))
print("Variance, standard deviation")
print(dsvariance(Zc))
print(sqrt(dsvariance(Zc)))
print("Aggregated width")
print(dsaggwidth(Zc))

#Do sensitivity analysis
sens=(dssensitivity(list(Q,Ks,Zm,Zv),c(1,2,3,4),fcrue,dsaggwidth,mcIT=20,pinch_samples=20,pinch_type='distribution'));
dev.new()
barplot(sens,beside=TRUE,names=list("Q","Ks","Zm","Zv"))
title("Sensitivity on aggregated width")
#Study with correlation between Zm and Zv
temp=dsevalmc(fcrue,list(Q,Ks,Zm,Zv),10000,dsmonotonous,corr=c(0,0,0,0,0,0.66));
Zccorr=temp[[1]]
dev.new()
dscdf(Zccorr,xlab="Zccorr",ylab="",xrange=c(48,65))
}

