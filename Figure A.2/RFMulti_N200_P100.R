install.packages("randomForest")
library(randomForest)

#Data generating #Multilevel
set.seed(666)
B=500
N=200
p=100
s=5
X=array(rnorm(N*B*p,0.3,1),dim=c(B,N,p))
gamma0=c(s:1,rep(0,(p-s)))/s
gamma1=gamma0+0.5
gamma2=gamma0+0.6


pr0=matrix(0.3,B,N)
pr1=matrix(0.3,B,N)
pr2=matrix(0.4,B,N)

p <- c(0.3,0.3,0.4)
x <- runif(B*N)
W <- numeric(B*N)
for (i in 0:length(p)) { 
  W <- W + (x >= sum(p[0:i])) 
}

W=matrix(W,B,N)

W=W-1

beta1=gamma0+0.5
theta1=3
theta2=6

e1=matrix(rnorm(B*N,0,0.1),B,N)
e2=matrix(rnorm(B*N,0,0.1),B,N)
e3=matrix(rnorm(B*N,0,0.1),B,N)
e4=matrix(rnorm(B*N,0,0.1),B,N)

Y00=matrix(0,B,N)
Y01=matrix(0,B,N)
Y11=matrix(0,B,N)
Y21=matrix(0,B,N)

for (i in 1:B){
  
  Y00[i,]=X[i,,]%*%beta1+e1[i,]
  Y01[i,]=Y00[i,]+1+e2[i,]
  Y11[i,]=theta1+Y01[i,]+e3[i,]
  Y21[i,]=theta2+Y01[i,]+e4[i,]
}

Y0=Y00

Y1=Y01*(W==0)+Y11*(W==1)+Y21*(W==2)



#####################################################################################
#Abadie's Estimator
thetahat=c(0)
for (i in 1:B){
  treatment=as.factor(W[i,])
  data1=data.frame(treatment, X[i,,])
  
  rf = randomForest(treatment~., data = data1, norm.votes = TRUE, proximity = TRUE)
  p1 = predict(rf, data1, type = "prob")
  
  ghat0=p1[,1]
  ghat2=p1[,3]  
  thetahat[i]=mean(((W[i,]==2)*ghat0-(W[i,]==0)*ghat2)*(Y1[i,]-Y0[i,])/ghat0/mean(W[i,]==2))
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
    treatment=as.factor(W[i,-K[q,]])
    data1=data.frame(treatment, X[i,-K[q,],])
    
    rf = randomForest(treatment~., data = data1, norm.votes = TRUE, proximity = TRUE)
    treatment=as.factor(W[i,K[q,]])
    data2=data.frame(treatment, X[i,K[q,],])
    
    p1 = predict(rf, data2, type = "prob")
    
    ghat0=p1[,1]
    ghat2=p1[,3]  
    
    index=which(W[i,-K[q,]]==0)
    y=Y1[i,-K[q,]]-Y0[i,-K[q,]]
    y=y[index]
    XX=X[i,-K[q,],1:s]
    XX=XX[index,]
    
    
    model=randomForest(XX,y)
    
    
    ellhat_w=predict(model,X[i,K[q,],1:s])
    
    thetaDML[q]=mean(((W[i,K[q,]]==2)*ghat0-(W[i,K[q,]]==0)*ghat2)*(Y1[i,K[q,]]-Y0[i,K[q,]])/ghat0/mean(W[i,K[q,]]==2)-ellhat_w*(((W[i,K[q,]]==2)-ghat2)/mean(W[i,K[q,]]==2)-((W[i,K[q,]]==0)-ghat0)*ghat2/mean(W[i,K[q,]]==2)/ghat0))
  }
  thetabar[i]=mean(thetaDML)
}




par(mfrow=c(1,2))
hist(thetahat,breaks=100,main="Abadie",xlab="",ylab="")
hist(thetabar,breaks=50,main="HD",xlab="",ylab="")



