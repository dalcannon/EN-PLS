#EN,Ridge,Lasso
enetreg <- function(X,y,avec,lamvec,cal,val,test,scale=F,yconv){
  #This function is a wrapper function to optimize the elastic net regression
  #Uses function glmnet (by Friedman et. al) to calculate EN Regression
  #Inputs:
  #X = matrix of data
  #y = vector of property
  #avec = vector of alpha values. alpha values must be between 0 and 1.
  #lamvec = vector of lambda indices in sequence to use for EN. By default, a sequence
  #of 100 lambda values are calculated (see glmnet package documentation). The values of lamvec must be between 1 and 100
  #cal, val, test: indices for calibration, optimization, and test sets
  #scale = Is X autoscaled when EN is calculated? If FALSE, X is mean-centered
  #yconv = Converts units of y for calculation of RMSEP
  #Outputs: 
  #RMSEopt= RMSEopt for EN
  #RMSEPtest=RMSEP for test set using optimal EN model
  #nvars=Number of variables selected by EN
  #vars=vecter of variable indices (channel numbers) for the variables selected by EN
  #vars.pos=vector of variable indices for selected variables with positive EN regression coefficients
  #vars.neg=vector of variable indices for selected variables with negative EN regression coefficients
  #optnlv=optimal number of latent variables chosen 
  #alflam= Values of alpha and lambda for best EN model
  #PRESSopt=Predicted Residual Error Sum of Squares for optimization set
  #PRESStest=Predicted Residual Error Sum of Squares for test set
  #encoef=Elastic Net Regression coefficients for best EN model
  #R code by Cannon Giglio, 2016-2018, Brown Research Group, University of Delaware
  require(glmnet)
  #varexpl
  if(min(avec)<0 | max(avec)>1){ stop("alpha values must be between 0 and 1")}
  if(min(lamvec)<=0 | max(lamvec)>100){stop("lambda indices (lamvec) must be greater than 0 and less than/equal to 100")}
  if(missing(yconv)){
    yconv <- rep(1,length(y))
  }
  glpval <- array(100,dim=c(length(val),length(avec),length(lamvec)),dimnames=list(val,avec,lamvec))
  RMSEPval <- matrix(100,nrow=length(avec),ncol=length(lamvec),dimnames=list(avec,lamvec))
  lammat <- matrix(100,nrow=length(avec),ncol=length(lamvec),dimnames=list(avec,lamvec))
  #valbias <- matrix(100,nrow=length(avec),ncol=length(lamvec),dimnames=list(avec,lamvec))
  Nvarmat <- matrix(100,nrow=length(avec),ncol=length(lamvec),dimnames=list(avec,lamvec))
  #varexpl <- array(100,dim=c(length(avec),length(lamvec),2),dimnames=list(avec,lamvec,c("X","y")))
  X0 <- scale(X,center=colMeans(X[cal,]),scale=F)
  y <- scale(y,center=mean(y[cal]),scale=F)
  if(scale==T){    X <- scale(X,center=colMeans(X[cal,]),scale=apply(X[cal,],2,sd))  }
  else{X <- X0}
  for(i in 1:length(avec)){
    enp <- glmnet(X[cal,],y[cal],family="gaussian",standardize=F,type.gaussian="naive",alpha=avec[i])
    lammat[i,] <- enp$lambda[lamvec]
    for(j in 1:length(lamvec)){
      if(lamvec[j]>length(enp$lambda)){ next }
      enpcoef <- as.matrix(coef(enp, s=enp$lambda[lamvec[j]])[2:(ncol(X)+1),])
      enpcs <- as.vector(apply(enpcoef, 2, function(z)sum(z!=0))) #Number of active variables
      Nvarmat[i,j] <- enpcs 
      glpval[,i,j] <- predict(enp,newx=X[val,],s=enp$lambda[lamvec[j]])
      RMSEPval[i,j] <- sqrt(mean(yconv[val]^2*(glpval[,i,j]-y[val])^2))
      #RMSEPval[i,j] <- sqrt(sum(abs(glpval[,i,j]-y[val])^2)/length(y[val]))
      #valbias[i,j] <- abs(sum(y[val] - glpval[,i,j])/length(val))
    }
  }  
  bestENrmsep <- min(RMSEPval)
  ENmin <-  which(RMSEPval==min(RMSEPval),arr.ind=T)[1,1:2]
  glm.en <- glmnet(X[cal,],y[cal],family="gaussian",standardize=F,type.gaussian="naive",alpha=avec[ENmin[1]])
  encoef <- as.matrix(coef(glm.en, s=glm.en$lambda[lamvec[ENmin[2]]])[2:(ncol(X)+1),])
  vars <- which(encoef!=0)
  nvars <- sum(encoef!=0)
  enalflam <- c(avec[ENmin[1]],glm.en$lambda[lamvec[ENmin[2]]])
  #ENbias <- c(valbias[ENmin[1],ENmin[2]])
  #names(ENbias) <- c("Val Bias")
  RMSEopt <- bestENrmsep[1]
  valpreds <- predict(glm.en,newx=X[val,],s=glm.en$lambda[lamvec[ENmin[2]]])
  PRESSopt= sum((valpreds-y[val])^2) #sum(yconv[val]^2*(valpreds-y[val])^2)
  testpreds <- predict(glm.en,newx=X[test,],s=glm.en$lambda[lamvec[ENmin[2]]])
  RMSEPtest <- sqrt(mean(yconv[test]^2*(testpreds-y[test])^2))
  PRESStest <- sum(yconv[test]^2*(testpreds-y[test])^2) #sum(yconv[test]^2*(testpreds-df$y[test])^2)
  return(list(RMSEopt=RMSEopt,RMSEPtest=RMSEPtest,nvars=nvars,vars=vars,
              PRESSopt=PRESSopt,PRESStest=PRESStest,enalflam=enalflam,encoef=encoef))
  #return(list(RMSEPval=RMSEPval,nvars=Nvarmat,ENRMSEP=bestENrmsep,ENnvars=ENnvars,ENalflam=enalflam, ENcoef=encoef,
  #            ENvars=envars,valbias=valbias,ENbias=ENbias,lambdas=lammat,RMSEPtest=RMSEPtest))
}