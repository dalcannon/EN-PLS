enpls <- function(X,y,avec,lamvec,cal,val,test,Ncomp,scale=F,yconv,PLSsc=F){
  #This function uses the elastic net regression to select variables, and PLS regression to model the data 
  #using the selected variables. 
  #Inputs:
  #X = matrix of data
  #y = vector of property
  #avec = vector of alpha values. alpha must be between 0 and 1. PLS-based outputs assume that avec[1]=0
  #lamvec = vector of lambda's in sequence to use for EN-PLS. By default, 100 are calculated. 
  #cal, val, test: indices for calibration, optimization, and test sets
  #Ncomp = Number of PLS components
  #scale = Is X autoscaled when EN is calculated?
  #PLSsc = If TRUE, Autoscales X when PLS regression is calculated
  #yconv = Converts units of y for calculation of RMSEP
  #Outputs: 
  #RMSEopt=Optimal RMSEopt for EN-PLS 
  #RMSEPtest=RMSEP for test set using optimal EN-PLS model
  #nvars=Number of variables selected by EN-PLS
  #vars=vecter of variable indices (channel numbers) for the variables selected by EN-PLS
  #vars.pos=vector of variable indices for selected variables with positive EN regression coefficients
  #vars.neg=vector of variable indices for selected variables with negative EN regression coefficients
  #optnlv=optimal number of latent variables chosen 
  #alflam= Values of alpha and lambda for best ENPLS model
  #PRESSopt=Predicted Residual Error Sum of Squares for optimization set
  #PRESStest=Predicted Residual Error Sum of Squares for test set
  #coefs=Elastic Net Regression coefficients for best EN-PLS model
  #For details on the EN-PLS Regression see "Using elastic net regression to perform spectrally relevant
  #variable selection" by C. Giglio and S.D. Brown, J. Chemometrics, 2018. 
  #R code by Cannon Giglio, 2016-2018, Brown Research Group, University of Delaware
  if(min(avec)<0 | max(avec)>1){ stop("alpha values must be between 0 and 1")}
  if(min(lamvec)<=0 | max(lamvec)>100){stop("lambda indices must be greater than 0 and less than/equal to 100")}
  y <- as.vector(y)
  if(length(y)!=nrow(X)){stop("Length of y does not Equal Number of Samples in X")}
  require(pls)
  require(glmnet)
  if(missing(yconv)){
    yconv <- rep(1,length(y))
  }
  lammat <- matrix(0, nrow=length(avec),ncol=length(lamvec),dimnames=list(avec,lamvec))
  Nvarmat <- matrix(0, nrow=length(avec),ncol=length(lamvec),dimnames=list(avec,lamvec))
  NLVmat <- matrix(0, nrow=length(avec),ncol=length(lamvec),dimnames=list(avec,lamvec))
  RMSEPval <- array(100, dim=c(length(avec),length(lamvec),Ncomp),dimnames=list(avec,lamvec,1:Ncomp))
  #biasarr <- array(100, dim=c(length(avec),length(lamvec),Ncomp),dimnames=list(avec,lamvec,1:Ncomp))
  #varexpl <- array(100, dim=c(length(avec),length(lamvec),2),dimnames=list(avec,lamvec,c("X","y")))
  X <- scale(X,center=colMeans(X[cal,]),scale=F)
  y <- scale(y,center=mean(y[cal]),scale=F)
  if(scale==T){    X0 <- scale(X,center=colMeans(X[cal,]),scale=apply(X[cal,],2,sd))  }
  else{X0 <- X}
  for(i in 1:length(avec)){
    enp <- glmnet(X0[cal,],y[cal],alpha=avec[i],family="gaussian",standardize=F,type.gaussian="naive",nlambda=100)
    for(j in 1:length(lamvec)){
      if(lamvec[j]>length(enp$lambda)){ next }
      lammat[i,j] <- enp$lambda[lamvec[j]]
      enpcoef <- as.matrix(coef(enp, s=enp$lambda[lamvec[j]])[2:(ncol(X)+1),])
      enpcs <- as.vector(apply(enpcoef, 2, function(z)sum(z!=0))) #number of active variables
      if(enpcs == 0){ next  }
      else{
        colist <- rep(0,(enpcs))
        colist <- which(enpcoef != 0) #vector/list of nonzero coefficients
      }
      Nvarmat[i,j] <- enpcs
      if(j==1){
        whatcoef0 <- 0
      }
      whatcoef = which(enpcoef!=0)
      if(j>1 & length(whatcoef)==length(whatcoef0)){
        if(isTRUE(any((whatcoef-whatcoef0)!=0))==F){
          RMSEPval[i,j,] <- RMSEPval[i,(j-1),]
          NLVmat[i,j] <- NLVmat[i,(j-1)]
          #biasarr[i,j,] <- biasarr[i,(j-1),]
          #varexpl[i,j,] <- varexpl[i,(j-1),]
          next
        }
      }
      whatcoef0 <- whatcoef
      if(PLSsc==T){
        enpdf.cal <- data.frame(y = y[cal], X = I(X0[cal,colist]))
        df.val <- data.frame(y=y[val], X=I(X0[val,colist]))
      }
      else{
        enpdf.cal <- data.frame(y = y[cal], X = I(X[cal,colist]))
        df.val <- data.frame(y=y[val], X=I(X[val,colist]))
      }
      enpls <- plsr(y ~ X, data=enpdf.cal,subset=cal, ncomp=min(enpcs,Ncomp))
      enRMSEP <- rep(0,enpls$ncomp)
      valpreds <- predict(enpls, newdata = df.val, ncomp=1:min(enpcs,Ncomp))
      for(L in 1:enpls$ncomp){
        enRMSEP[L] <- sqrt(mean(yconv[val]^2*(valpreds[,,L]-df.val$y)^2))
        #biasarr[i,j,L] <- abs(sum(df.val[,1] - valpreds[,,L])/nrow(df.val))
      }
      if(length(enRMSEP) < Ncomp){
        enlen <- length(enRMSEP)
        fill <- rep(100,(Ncomp-enlen))
        RMSEPval[i,j,] <- c(enRMSEP,fill)
      }
      else { RMSEPval[i,j,] <- enRMSEP }
      optmin <- which.min(enRMSEP)
      NLVmat[i,j] <- optmin
    }
  }
  BVval <- min(RMSEPval) #Best Validation RMSEP for Elastic Net
  BVlamc <- which(RMSEPval==BVval,arr.ind=T)[1,1:3] #Indices 
  nvars <- Nvarmat[BVlamc[1],BVlamc[2]]
  RMSEopt <- RMSEPval[BVlamc[1],BVlamc[2],BVlamc[3]]
  alflam <- c(avec[BVlamc[1]],lammat[BVlamc[1],BVlamc[2]])
  BVenp <- glmnet(X0[cal,],y[cal],alpha=avec[BVlamc[1]],family="gaussian",standardize=F,type.gaussian="naive",nlambda=100)
  ENcoef <- as.matrix(coef(BVenp, s=BVenp$lambda[lamvec[BVlamc[2]]])[2:(ncol(X)+1),])
  vars <- which(ENcoef != 0) 
  vars.pos <- which(ENcoef > 0) 
  vars.neg <- which(ENcoef < 0) 
  optnlv <- NLVmat[BVlamc[1],BVlamc[2]]
  #BVbias <- biasarr[BVlamc[1],BVlamc[2],BVlamc[3]]#,optnlv]
  if(PLSsc==T){df <- data.frame(y = y, X = I(X0))}
  else{    df <- data.frame(y = y, X = I(X))  }
  enpls.t <- plsr(y ~ X[,vars], data=df,subset=cal, ncomp=optnlv)
  valpreds <- predict(enpls.t, newdata = df[val,], ncomp=optnlv)
  PRESSopt= sum((valpreds-df$y[val])^2) #sum(yconv[val]^2*(valpreds-y[val])^2)
  testpreds <- predict(enpls.t, newdata = df[test,], ncomp=optnlv)
  PRESStest= sum(yconv[test]^2*(testpreds-df$y[test])^2) #sum((testpreds-df$y[test])^2)
  RMSEPtest <- sqrt(mean(yconv[test]^2*(testpreds[,,1]-df$y[test])^2))
  #testbias <- abs(sum(df$y[test] - testpreds[,,1])/length(test))
  #testSEP <- sqrt(RMSEPtest^2 - testbias^2)
  #ENSEP <- c(sqrt(BVval^2-BVbias^2),sqrt(RMSEPtest^2 - testbias^2))
  #return(list(RMSEPopt=RMSEPopt,RMSEPval=RMSEPval,nvars=Nvarmat,lamvals=lammat,NLV=NLVmat,
  #            vars=vars,vars.pos=vars.pos,vars.neg=vars.neg,alflam=alflam,nvars=nvars,
  #            optnlv=optnlv,ENcoef=ENcoef,RMSEPtest=RMSEPtest,testpreds=testpreds))
  return(list(RMSEopt=RMSEopt,RMSEPtest=RMSEPtest,nvars=nvars,vars=vars,vars.pos=vars.pos,vars.neg=vars.neg,optnlv=optnlv,
              alflam=alflam,PRESSopt=PRESSopt,PRESStest=PRESStest,coefs=ENcoef))
}