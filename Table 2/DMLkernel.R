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



###Estimation


DMLkernel=function(Y,D,X,T){
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
    h=npregbw(D~X)$bw
    for (q in 1:k){
      #Trimming
      g=function(x){
        w=dnorm((X[-K[q,]]-x)/h)
        w=w/sum(w)
        return(sum(D[-K[q,]]*w))
      }
      
      ghat=c(0)
      for (j in 1:length(K[q,])){
        ghat[j]=g(X[K[q,][j]])
      }
      
      index1=K[q,][which(ghat<0.97 & ghat>0.03)]
      
      ##Estimation
      ghat=c(0)
      for (j in 1:length(index1)){
        ghat[j]=g(X[index1[j]])
      }
      
      lambda=mean(T)
      
      index=which(D[-K[q,]]==0)
      YY=Y[-K[q,]]
      YY=YY[index]
      TT=T[-K[q,]]
      TT=TT[index]
      XX=X[-K[q,]]
      XX=XX[index]
      
      
      ell2=function(x){
        w=dnorm((XX-x)/h)
        w=w/sum(w)
        return(sum(YY*(TT-lambda)*w))
      }
      
      ellhat2=c(0)
      for (j in 1:length(index1)){
        ellhat2[j]=ell2(X[index1[j]])
      }
      s=(D[index1]-ghat)*((T[index1]-lambda)*Y[index1]-ellhat2)/(lambda*(1-lambda))/mean(D[index1])/(1-ghat)
      thetaDML[q]=mean(s)
      
    }
    thetabar[l]=mean(thetaDML)
  }
  finalthetabar=mean(thetabar)
  
  ###Standard Deviation
  #Standard Deviation
  var=c()
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
      #Trimming
      g=function(x){
        w=dnorm((X[-K[q,]]-x)/h)
        w=w/sum(w)
        return(sum(D[-K[q,]]*w))
      }
      
      ghat=c(0)
      for (j in 1:length(K[q,])){
        ghat[j]=g(X[K[q,][j]])
      }
      
      index1=K[q,][which(ghat<0.97 & ghat>0.03)]
      
      ##Estimation
      g=function(x){
        w=dnorm((X[-K[q,]]-x)/h)
        w=w/sum(w)
        return(sum(D[-K[q,]]*w))
      }
      
      ghat=c(0)
      for (j in 1:length(index1)){
        ghat[j]=g(X[index1[j]])
      }
      
      lambda=mean(T[-K[q,]])
      
      index=which(D[-K[q,]]==0)
      YY=Y[-K[q,]]
      YY=YY[index]
      TT=T[-K[q,]]
      TT=TT[index]
      XX=X[-K[q,]]
      XX=XX[index]
      
      
      ell2=function(x){
        w=dnorm((XX-x)/h)
        w=w/sum(w)
        return(sum((TT-lambda)*YY*w))
      }
      
      ellhat2=c(0)
      for (j in 1:length(index1)){
        ellhat2[j]=ell2(X[index1[j]])
      }
      
      G=-(1-2*lambda)*finalthetabar/(lambda*(1-lambda))-mean(Y[index1]*(D[index1]-ghat)/(1-ghat)/(lambda*(1-lambda))/mean(D[index1]))
      
      s=((T[index1]-lambda)*Y[index1]-ellhat2)*(D[index1]-ghat)/(1-ghat)/mean(D[index1])/(lambda*(1-lambda))-D[index1]*finalthetabar/mean(D[index1])+G*(T[index1]-lambda)
      varDML[q]=mean(s^2)
    }
    
    var[m]=mean(varDML)
  }
  finalvar=mean(var)
  
  SD=sqrt(finalvar)/sqrt(N)
  return(c(finalthetabar,SD))
}

DMLkernel(Y,D,X,T)