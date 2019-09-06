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

####

DMLML=function(Y,D,p,T){
  N=length(Y)
  B=100
  set.seed(123)
  random=sample(1:1000,B)
  
  thetabar=c(0)
  for (l in 1:B){
    k=2
    samplesplit=function(k,N){
      c1=1:N
      smp_size <- floor((1/k) * length(c1))
      
      ## set the seed to make your partition reproducible
      set.seed(random[l])
      train_ind <- sample(seq_len(length(c1)), size = smp_size)
      
      k1 <- c1[train_ind]
      k2 <- c1[-train_ind]
      return(rbind(k1,k2))
    }
    K=samplesplit(k,N)
    
    thetaDML=c(0)
    
    for (q in 1:k){
      ##Trimming
      set.seed(333)
      CV=cv.glmnet(p[-K[q,],],D[-K[q,]],family="binomial",alpha=1)
      fit=glmnet(p[-K[q,],],D[-K[q,]],family="binomial",alpha=1,lambda=CV$lambda.1se)
      beta1hat=fit$beta
      beta1hat <- as.numeric(as.character(beta1hat))
      
      ghat=1/(1+exp(-p[K[q,],]%*%beta1hat))
      
      index1=K[q,][which(ghat<0.97 & ghat>0.03)]
      
      ##Estimation
      ghat=1/(1+exp(-p[index1,]%*%beta1hat))
      
      lambda=mean(T[-K[q,]])
      
      index=which(D[-K[q,]]==0)
      YY=Y[-K[q,]]
      YY=YY[index]
      TT=T[-K[q,]]
      TT=TT[index]
      XX=p[-K[q,],]
      XX=XX[index,]
      
      set.seed(333)
      CV=cv.glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1)
      fit=glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1,lambda=CV$lambda.1se)
      beta2hat=fit$beta
      beta2hat <- as.numeric(as.character(beta2hat))
      
      ellhat2=p[index1,]%*%beta2hat
      
      s=((T[index1]-lambda)*Y[index1]-ellhat2)*(D[index1]-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D[index1])
      s=s[which(s<abs(min(s)))]
      
      thetaDML[q]=mean(s)
      
    }
    
    thetabar[l]=mean(thetaDML)
    
    
  }
  finaltheta=mean(thetabar)
  finaltheta
  
  ##Variance
  var=c(0)
  for (m in 1:B){
    k=2
    samplesplit=function(k,N){
      c1=1:N
      smp_size <- floor((1/k) * length(c1))
      
      ## set the seed to make your partition reproducible
      set.seed(random[m])
      train_ind <- sample(seq_len(length(c1)), size = smp_size)
      
      k1 <- c1[train_ind]
      k2 <- c1[-train_ind]
      return(rbind(k1,k2))
    }
    K=samplesplit(k,N)
    
    varDML=c(0)
    for (q in 1:k){
      ##Trimming
      set.seed(333)
      CV=cv.glmnet(p[-K[q,],],D[-K[q,]],family="binomial",alpha=1)
      fit=glmnet(p[-K[q,],],D[-K[q,]],family="binomial",alpha=1,lambda=CV$lambda.1se)
      beta1hat=fit$beta
      beta1hat <- as.numeric(as.character(beta1hat))
      
      ghat=1/(1+exp(-p[K[q,],]%*%beta1hat))
      
      index1=K[q,][which(ghat<0.97 & ghat>0.03)]
      
      ##Estimation
      ghat=1/(1+exp(-p[index1,]%*%beta1hat))
      lambda=mean(T[-K[q,]])
      
      index=which(D[-K[q,]]==0)
      YY=Y[-K[q,]]
      YY=YY[index]
      TT=T[-K[q,]]
      TT=TT[index]
      XX=p[-K[q,],]
      XX=XX[index,]
      
      set.seed(333)
      CV=cv.glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1)
      fit=glmnet(XX,(TT-lambda)*YY,family="gaussian",alpha=1,lambda=CV$lambda.1se)
      beta2hat=fit$beta
      beta2hat <- as.numeric(as.character(beta2hat))
      
      ellhat2=p[index1,]%*%beta2hat
      
      
      
      G=-(1-2*lambda)*finaltheta/(lambda*(1-lambda))-mean(Y[index1]*(D[index1]-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D[index1]))
      
      s=((T[index1]-lambda)*Y[index1]-ellhat2)*(D[index1]-ghat)/(1-ghat)/mean(D[index1])/(lambda*(1-lambda))-D[index1]*finaltheta/mean(D[index1])+G*(T[index1]-lambda)
      s=s[which(s<abs(min(s)))]
      
      varDML[q]=mean(s^2)
    }
    
    var[m]=mean(varDML)
  }
  
  sd=sqrt(mean(var))/sqrt(N)
  sd
  return(c(finaltheta,sd))
}

DMLML(Y,D,p,T)
