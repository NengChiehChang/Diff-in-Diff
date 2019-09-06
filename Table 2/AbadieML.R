install.packages("readstata13")
library(readstata13)

install.packages("glmnet")
library(glmnet)

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


X=cbind(data$post_2008, data$tariff2007,data$lvalue_tonnage, data$differentiated, data$agri, data$perishable,data$dfs, data$day_w_arrival, data$monitor, data$psi, data$rsa, data$post_2008,clear_agent2,clear_agent3,clear_agent4,clear_agent5,clear_agent6,clear_agent7,clear_agent8,hc_group2,hc_group3,hc_group4,hc_group5,hc_group6,hc_group7,hc_group8,hc_group9,hc_group10,hc_group11,hc_group12,hc_group13,hc_group14,hc_group15,data$hc_4digit)

#Working data
Working_data=cbind(Y,D,T,X)

Working_data=Working_data[complete.cases(Working_data), ]

Y=Working_data[,1]
D=Working_data[,2]
T=Working_data[,3]
X=Working_data[,4:12]

index=t(combn(length(X[1,]),2))
p=X

for (i in 1:length(index[,1])){
  p=cbind(p,X[,index[i]]*X[,index[i+length(index[,1])]])
}

###Estimation Function
AbadieML=function(Y,D,X,T){
  ##Trimming
  #CV=cv.glmnet(p,D,family="binomial",alpha=1)
  fit=glmnet(X,D,family="binomial",alpha=1,lambda=0.01638)
  beta1hat=fit$beta
  beta1hat <- as.numeric(as.character(beta1hat))
  
  ghat=1/(1+exp(-X%*%beta1hat))
  
  index=which(ghat<0.98 & ghat>0.02)
  
  Y=Y[index]
  D=D[index]
  X=X[index,]
  T=T[index]
  
  
  ##Estimation
  N=length(Y)
  #CV=cv.glmnet(X,D,family="binomial",alpha=1)
  fit=glmnet(X,D,family="binomial",alpha=1,lambda=0.00645)
  
  beta1hat=fit$beta
  beta1hat <- as.numeric(as.character(beta1hat))
  
  ghat=1/(1+exp(-X%*%beta1hat))
  lambda=mean(T)
  s=(T-lambda)*Y*(D-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D)
  
  thetahat=mean(s)
  
  ##Variance Estimation
  G_2_lambda=-(1-2*lambda)*thetahat/(lambda*(1-lambda))-mean(Y*(D-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D))
  
  index=which(D==0)
  YY=Y[index]
  
  TT=T[index]
  
  XX=X[index,]
  
  #CV=cv.glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1)
  fit=glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1,lambda=0.0433)
  beta2hat=fit$beta
  beta2hat <- as.numeric(as.character(beta2hat))
  
  ellhat2=X%*%beta2hat
  
  var=mean((((T-lambda)*Y-ellhat2)*(D-ghat)/(1-ghat)/mean(D)/(lambda*(1-lambda))-D*thetahat/mean(D)+G_2_lambda*(T-lambda))^2)
  sd=sqrt(var)/sqrt(N)
  
  return(c(thetahat, sd))
}

AbadieML(Y,D,p,T)