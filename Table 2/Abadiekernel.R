install.packages("readstata13")
library(readstata13)

install.packages("np")
library(np)
#Raw Data
data <- read.dta13 ("/Users/jack1021/Documents/2017C/Semi DID/empirical example/Sequeira 2016/Replication_Final/Bribes_Regression.dta")

Y=data$lba
D=data$tariff_change_2008
T=data$post_2008

clear_agent2=(data$clear_agent==2)
clear_agent3=(data$clear_agent==3)
clear_agent4=(data$clear_agent==4)
clear_agent5=(data$clear_agent==5)
clear_agent6=(data$clear_agent==6)
clear_agent7=(data$clear_agent==7)
clear_agent8=(data$clear_agent==8)

hc_group2=(data$hc_group==2)
hc_group3=(data$hc_group==3)
hc_group4=(data$hc_group==4)
hc_group5=(data$hc_group==5)
hc_group6=(data$hc_group==6)
hc_group7=(data$hc_group==7)
hc_group8=(data$hc_group==8)
hc_group9=(data$hc_group==9)
hc_group10=(data$hc_group==10)
hc_group11=(data$hc_group==11)
hc_group12=(data$hc_group==12)
hc_group13=(data$hc_group==13)
hc_group14=(data$hc_group==14)
hc_group15=(data$hc_group==15)


X=cbind(data$lvalue_tonnage, data$tariff2007,data$tariff_change_2008, data$differentiated, data$agri, data$perishable,data$dfs, data$day_w_arrival, data$psi, data$monitor, data$rsa, data$term, data$post_2008,clear_agent2,clear_agent3,clear_agent4,clear_agent5,clear_agent6,clear_agent7,clear_agent8,hc_group2,hc_group3,hc_group4,hc_group5,hc_group6,hc_group7,hc_group8,hc_group9,hc_group10,hc_group11,hc_group12,hc_group13,hc_group14,hc_group15,data$hc_4digit)

#Working data
Working_data=cbind(Y,D,T,X)

Working_data=Working_data[complete.cases(Working_data), ]

Y=Working_data[,1]
D=Working_data[,2]
T=Working_data[,3]
X=Working_data[,4]


###Estimation Function
AbadieDIDkernel=function(Y,D,X,T){
  ##Trimming
  N=length(Y)
  h=npregbw(D~X)$bw
  g=function(x){
    w=dnorm((X-x)/h)
    w=w/sum(w)
    return(sum(D*w))
  }
  
  ghat=c(0)
  for (j in 1:N){
    ghat[j]=g(X[j])
  }
  
  index=which(ghat<0.97 & ghat>0.03)
  
  Y=Y[index]
  D=D[index]
  X=X[index]
  T=T[index]
  
  ##Estimation
  N=length(Y)
  lambda=mean(T)
  h=npregbw(D~X)$bw
  thetahat=c(0)
  
  g=function(x){
    w=dnorm((X-x)/h)
    w=w/sum(w)
    return(sum(D*w))
  }
  
  ghat=c(0)
  for (j in 1:N){
    ghat[j]=g(X[j])
  }
  
  thetahat=mean((T-lambda)*Y*(D-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D))
  
  #Standard Deviation
  G_2_lambda=-(1-2*lambda)*thetahat/(lambda*(1-lambda))-mean(Y*(D-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D))
  
  ell2=function(x){
    w=dnorm((X*(1-D)-x)/h)
    w=w/sum(w)
    return(sum((T*(1-D)-lambda)*Y*(1-D)*w))
  }
  
  ellhat2=c(0)
  for (j in 1:N){
    ellhat2[j]=ell2(X[j])
  }
  
  var=mean((((T-lambda)*Y-ellhat2)*(D-ghat)/(1-ghat)/mean(D)/(lambda*(1-lambda))-D*thetahat/mean(D)+G_2_lambda*(T-lambda))^2)
  sd=sqrt(var)/sqrt(N)
  return(c(thetahat, sd))
}


AbadieDIDkernel(Y,D,X,T)
