#VIP
VIPpls <- function(X,y,cal,val,test,Ncomp,yconv,sc=F){
  #This function is an implementation of the Variable Importance in Projection (VIP)
  #method  (see citation below) for variable slection in Partial Least Squares (PLS) regression. 
  #Wold S, Johansson A. In: Cocchi M, ed. PLS-Partial Least Squares Projections to Latent Structures. 
  #Leiden, Netherlands: ESCOM Science Publishers; 1993:523-550.
  #Inputs:
  #X = matrix of data
  #y = vector of property
  #cal, val: vectors for calibration, validation sets
  #Ncomp = Number of PLS components
  #yconv = conversion factor for calculating RMSEPs for y
  #Outputs:
  #RMSEopt=Optimal RMSEopt for VIP
  #RMSEPtest=RMSEP for test set using optimal VIP model
  #nvars=Number of variables selected by VIP
  #vars=vecter of variable indices (channel numbers) for the variables selected by VIP
  #optnlv=optimal number of latent variables chosen 
  #PRESSopt=Predicted Residual Error Sum of Squares for optimization set
  #PRESStest=Predicted Residual Error Sum of Squares for test set
  #VIP = VIP scores vector
  #cutoff = VIP threshold used for optimal VIP-PLS model
  # R code by Cannon Giglio, 2016-2018, Brown Research Group, University of Delaware
  require(pls)
  if(missing(yconv)){
    yconv <- rep(1,length(y))
  }
  X <- scale(X,center=colMeans(X[cal,]),scale=F)
  y <- scale(y,center=mean(y[cal]),scale=F)
  df.cal <- data.frame(y = y[cal], X = I(X[cal,]))
  df.val <- data.frame(y = y[val], X = I(X[val,]))
  mvr <- plsr(y ~ X, ncomp = Ncomp, data = df.cal,y=TRUE,scale=sc)
  rmsep_res <- RMSEP(mvr, newdata = df.val, intercept = FALSE)
  rmsepval <- c(1:Ncomp)
  for(i in 1:Ncomp){rmsepval[i] <- unlist(rmsep_res)[[i]] }
  optLVs <- which.min(rmsepval)
  Wts <- loading.weights(mvr)[,1:Ncomp]
  Scores <- scores(mvr)[,1:Ncomp]
  J <- nrow(Wts)
  Vipmat <- matrix(0,nrow=Ncomp,ncol=J)
  Viprankmat <- matrix(0,nrow=Ncomp,ncol=J)
  for(i in 1:Ncomp){
    for(j in 1:J){
      res <- residuals(mvr)[,,i]
      Sco <- Scores[,1:i]
      d <- t(t(Sco)%*%(mvr$y-res))%*%solve(t(Sco)%*%Sco) 
      Vipmat[i,j] <- sqrt((t(Wts[j,1:i]^2)%*%t((d^2)%*%(t(Sco)%*%Sco))*J)/sum((d^2)%*%(t(Sco)%*%Sco))) 
    }
  }
  for(i in 1:Ncomp){
    Viprankmat[i,] <- order(Vipmat[i,],decreasing=T)
  }
  VIP <- Vipmat[optLVs,]
  Viprank <- Viprankmat[optLVs,]
  RMSEPval <- matrix(100, nrow=J, ncol=Ncomp)
  valbias <- matrix(100, nrow=J, ncol=Ncomp)
  VipoptLVs <- rep(100,J)
  for(j in 1:J){ 
    whichVip <- Viprank[1:j]
    Vipdf.cal <- data.frame(y = y[cal], X = I(X[cal,whichVip]))
    Vipdf.val <- data.frame(y = y[val], X = I(X[val,whichVip]))
    Vippls <- plsr(y ~ X, data=Vipdf.cal,subset=cal,ncomp=min(j,Ncomp))
    valpreds <- predict(Vippls, newdata = Vipdf.val, ncomp=1:min(j,Ncomp))
    Viprmsepval <- c(1:Vippls$ncomp)
    for(L in 1:min(j,Ncomp)){ 
      Viprmsepval[L] <- sqrt(mean(yconv[val]^2*(valpreds[,,L]-Vipdf.val$y)^2))
      valbias[j,L] <- abs(sum(Vipdf.val[,1] - valpreds[,,L])/nrow(Vipdf.val))
    }
    #vip.rmsepval <- RMSEP(Vippls, newdata = Vipdf.val, intercept = FALSE,ncomp=1:min(j,Ncomp))
    #Viprmsepval <- c(1:min(j,Ncomp))
    #for(k in 1:min(j,Ncomp)){ Viprmsepval[k] <- unlist(vip.rmsepval)[[k]]}
    if(anyNA(Viprmsepval==T)){
      watNaNs <- which(is.nan(Viprmsepval))
      Viprmsepval[watNaNs] <- 222
    }
    if(length(Viprmsepval) < Ncomp){
      Viplen <- length(Viprmsepval)
      fill <- rep(100,(Ncomp-Viplen))
      RMSEPval[j,] <- c(Viprmsepval,fill)
    }
    else { RMSEPval[j,] <- Viprmsepval }
    VipoptLVs[j] <- which.min(Viprmsepval)
  }
  minRMSEP <- min(RMSEPval)
  nvars <- which(RMSEPval==min(RMSEPval),arr.ind=T)[1,1]
  if(nvars==J){
    minRMSEP <- min(RMSEPval[-J,])
    nvars <- which(RMSEPval==minRMSEP,arr.ind=T)[1,1]
  }
  optnlv <- VipoptLVs[nvars]
  vars <- sort(Viprank[1:nvars])
  cutoff <- VIP[Viprank[nvars]]
  RMSEopt=minRMSEP[1]
  #bestbias <- c(valbias[nvars,optnlv])
  #SEP <- sqrt(minRMSEP^2-bestbias[1]^2)
  #names(bestbias) <- c("Validation")
  #row.names(RMSEPval) <- c(1:J) 
  df <- data.frame(y = y, X = I(X))
  vippls.t <- plsr(y ~ X[,vars], data=df,subset=cal, ncomp=optnlv)
  valpreds <- predict(vippls.t, newdata = df[val,], ncomp=optnlv)
  PRESSopt= sum((valpreds-df$y[val])^2) #sum(yconv[val]^2*(valpreds-y[val])^2)
  testpreds <- predict(vippls.t, newdata = df[test,], ncomp=optnlv)
  RMSEPtest <- sqrt(mean(yconv[test]^2*(testpreds[,,1]-df$y[test])^2))
  PRESStest= sum(yconv[test]^2*(testpreds-df$y[test])^2) #sum(yconv[test]^2*(testpreds-df$y[test])^2)
  return(list(RMSEopt=RMSEopt,RMSEPtest=RMSEPtest,nvars=nvars,vars=vars,optnlv=optnlv,
              PRESSopt=PRESSopt,PRESStest=PRESStest,VIP=VIP,cutoff=cutoff))
  #return(list(RMSEPval=RMSEPval,VIP=VIP,Viprank=Viprank,nvars=nvars,optnlv=optnlv,
  #            Vipvars=bestVip,cutoff=cutoff,
  #            rmsepval=rmsepval,optLVs=VipoptLVs,minRMSEP=minRMSEP,bias=bestbias,SEP=SEP,RMSEPtest=RMSEPtest))
}