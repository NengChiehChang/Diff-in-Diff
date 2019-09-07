install.packages('e1071') 
library(e1071) 

install.packages("randomForest")
library(randomForest)

#Data generating
set.seed(666)
B=500
N=200
p=100
s=5
X=array(rnorm(N*B*p,0.3,0.1),dim=c(B,N,p))
gamma=c(s:1,rep(0,(p-s)))/s

D=matrix(0,B,N)
z=matrix(0,B,N)
pr=matrix(0,B,N)

for (i in 1:B){
  z[i,] = X[i,,]%*%gamma        # linear combination 
  pr[i,] = 1/(1+exp(-z[i,]))         # P(D=1 given X) Probit
  D[i,] = rbinom(N,1,pr[i,]) 
}


beta1=gamma+0.5
beta2=gamma+1
theta=3

e1=matrix(rnorm(B*N,0,0.1),B,N)
e2=matrix(rnorm(B*N,0,0.1),B,N)
e3=matrix(rnorm(B*N,0,0.1),B,N)

Y00=matrix(0,B,N)
Y01=matrix(0,B,N)
Y11=matrix(0,B,N)

for (i in 1:B){
  
  Y00[i,]=1+e1[i,]
  Y01[i,]=Y00[i,]+1+e2[i,]
  Y11[i,]=theta+Y01[i,]+e3[i,]
}

Y0=Y00
Y1=Y01*(1-D)+Y11*D

T=matrix(rbinom(B*N,1,0.5),B,N)

Y=Y0+T*(Y1-Y0)




#####################################################################################
#Abadie's Estimator
thetahat=c(0)
for (i in 1:B){
  treatment=ifelse(D[i,]==1,'yes','no')
  data1=data.frame(treatment, X[i,,])
  
  svm_model <- svm(treatment ~ .
                   , data = data1
                   , probability = TRUE
  )
  class <- predict(svm_model, data1, probability = TRUE)
  p1 <- attr(class, "probabilities")[,1]
  
  ghat=p1
  
  lambda=mean(T[i,])
  thetahat[i]=mean((T[i,]-lambda)*Y[i,]*(D[i,]-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D[i,]))
  
}

hist(thetahat,breaks=100,main="Abadie",xlab="",ylab="")




#HD
###Sample splitting parameters
k=2
k1=c(1:(N/k))
k2=c((N/k+1):(2*N/k))
K=rbind(k1,k2)

#DML
thetabar=c(0)
for (i in 1:B){
  thetaDML=c(0)
  for (q in 1:k){
    treatment=ifelse(D[i,-K[q,]]==1,'yes','no')
    data1=data.frame(treatment, X[i,-K[q,],])
    
    svm_model <- svm(treatment ~ .
                     , data = data1
                     , probability = TRUE
    )
    treatment=ifelse(D[i,K[q,]]==1,'yes','no')
    data2=data.frame(treatment, X[i,K[q,],])
    
    class <- predict(svm_model, data2, probability = TRUE)
    p1 <- attr(class, "probabilities")[,1]
    
    ghat=p1
    
    lambda=mean(T[i,])
    
    index=which(D[i,-K[q,]]==0)
    YY=Y[i,-K[q,]]
    YY=YY[index]
    TT=T[i,-K[q,]]
    TT=TT[index]
    XX=X[i,-K[q,],1:s]
    XX=XX[index,]
    
    
    model2=randomForest(XX,YY*TT)
    model3=randomForest(XX,YY)
    
    ellhat2=predict(model2,X[i,K[q,],1:s])
    ellhat3=predict(model3,X[i,K[q,],1:s])
    
    
    
    thetaDML[q]=mean((T[i,K[q,]]-lambda)*Y[i,K[q,]]*(D[i,K[q,]]-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D[i,K[q,]])-(D[i,K[q,]]-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D[i,K[q,]])*(ellhat2-lambda*ellhat3))
  }
  thetabar[i]=mean(thetaDML)
}




par(mfrow=c(1,2))
hist(thetahat[which(thetahat>0)],breaks=100,main="Abadie",xlab="",ylab="")
hist(thetabar,breaks=100,main="DMLDiD",xlab="",ylab="")
