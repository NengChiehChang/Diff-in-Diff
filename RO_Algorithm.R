#The function DMLDiD returns the DMLDiD estimator and its estimated variance. 
#Data is randomly splitted into K=2 parts.
#To obtain a robust result, I repeat B=100 times and return the average of the 100 DMLDiD estimators. 
#Inputs:
#Y1: the outcome variable at t=1
#Y0: the outcome variable at t=0
#D: the treatment indicator
#p: the basis functions of control variables
#Trimming: I trim out the observations with propensity score value less than 0.05 and greater than 0.95

#Packages
install.packages("glmnet")
library(glmnet)

install.packages("randomForest")
library(randomForest)


#Algorithm

DMLDiD=function(Y1,Y0,D,p){
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
      
      index1=K[q,][which(ghat<0.95 & ghat>0.05)]
      
      ##Estimation
      ghat=1/(1+exp(-p[index1,]%*%beta1hat))
      
      index=which(D[-K[q,]]==0)
      y=Y1[-K[q,]]-Y0[-K[q,]]
      y=y[index]
      XX=X[-K[q,],]
      XX=XX[index,]
      
      model=randomForest(XX,y)
      ellhat=predict(model,X[index1,])
      
      thetaDML[q]=mean((Y1[index1]-Y0[index1])/mean(D[index1])*(D[index1]-ghat)/(1-ghat)-(D[index1]-ghat)/mean(D[index1])/(1-ghat)*ellhat)
      
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
      
      index=which(D[-K[q,]]==0)
      y=Y1[-K[q,]]-Y0[-K[q,]]
      y=y[index]
      XX=X[-K[q,],]
      XX=XX[index,]
      
      model=randomForest(XX,y)
      ellhat=predict(model,X[index1,])
      
      G=-finaltheta/mean(D[index1])
      
      s=(Y1[index1]-Y0[index1])/mean(D[index1])*(D[index1]-ghat)/(1-ghat)-(D[index1]-ghat)/mean(D[index1])/(1-ghat)*ellhat-finaltheta+G*(D-mean(D[index1]))
      
      varDML[q]=mean(s^2)
    }
    
    var[m]=mean(varDML)
  }
  
  sd=sqrt(mean(var))/sqrt(N)
  sd
  return(c(finaltheta,sd))
}

DMLDiD(Y1,Y0,D,p)