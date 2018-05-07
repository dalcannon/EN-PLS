SRpls <- function(X,y,cal,val,test,Ncomp,yconv,sc=F){#,centerX=T){
  #This function is an implementation of the Selectivity Ratio (SR)
  #method (see citation below) for variable slection in Partial Least Squares (PLS) regression. 
  #Chun H, Keles S. Sparse partial least squares regression for simultaneous dimension 
  #reduction and variable selection. J R Stat Soc Series B. 2010;72(1):3-25.
  #Inputs:
  #X = matrix of data
  #y = vector of property
  #cal, val,test: vectors for calibration, validation, and test sets
  #Ncomp = Number of PLS components
  #PLSsc = If TRUE, Autoscales X when PLS regression is calculated
  #yconv = Converts units of y for calculation of RMSEP
  #Outputs:
  #RMSEopt=Optimal RMSEopt for SR
  #RMSEPtest=RMSEP for test set using optimal SR model
  #nvars=Number of variables selected by SR
  #vars=vecter of variable indices (channel numbers) for the variables selected by SR
  #optnlv=optimal number of latent variables chosen 
  #PRESSopt=Predicted Residual Error Sum of Squares for optimization set
  #PRESStest=Predicted Residual Error Sum of Squares for test set
  #SR = Selectivity Ratio vector
  #cutoff = SR threshold used for optimal SR-PLS model
  # R code by Cannon Giglio, 2016-2018, Brown Research Group, University of Delaware
  if(missing(yconv)){
    yconv <- rep(1,length(y))
  }
  require(pls)
  df.cal <-  data.frame(y = scale(y[cal], scale = FALSE), X = I(scale(X[cal,], scale = FALSE)))
  df.val <- data.frame(y=scale(y[val], center=mean(y[cal]),scale=F), X=I(scale(X[val,], center=colMeans(X[cal,]), scale=F)))
  mvr <- plsr(y ~ X, ncomp = Ncomp, data = df.cal,y=T,scale=sc)
  rmsep_val <- RMSEP(mvr, newdata = df.val, intercept = FALSE)
  rmsepval <- c(1:Ncomp)
  for(i in 1:Ncomp){rmsepval[i] <- unlist(rmsep_val)[[i]] }
  optLVs <- which.min(rmsepval)
  J <- ncol(X)
  ## if(optLVs==1){bpls <- loading.weights(mvr)[,1]} #Previous version, NOT equal to loading weights
  ###else{bpls <- apply(loading.weights(mvr)[,1:optLVs],1,sum)} #Previous version, NOT equal to loading weights
  bpls <- as.vector(coef(mvr,ncomp=optLVs))
  wtp <- bpls/base:::norm(bpls,type="2") 
  ttp <- df.cal$X%*%wtp 
  ptp <- (t(df.cal$X)%*%ttp)/as.numeric(t(ttp)%*%ttp)
  Xtp <- ttp%*%t(ptp)
  Etp <- df.cal$X - Xtp
  varexpl <- apply(Xtp,2,var)
  varres <- apply(Etp,2,var)
  SR <- varexpl/varres 
  SRrank <- order(SR,decreasing=T)
  RMSEPval <- matrix(100, nrow=J, ncol=Ncomp)
  SRoptLVs <- rep(100,J)
  for(j in 1:J){ 
    whichSR <- SRrank[1:j]
    SRdf.cal <- data.frame(y = scale(y[cal], scale = F), X = I(scale(X[cal,whichSR], scale = F)))
    if(j==1){
      SRdf.val <- data.frame(y=scale(y[val], center=mean(y[cal]), scale=F), X=I(scale(X[val,whichSR], center=mean(X[cal,whichSR]), scale=F)))
    }
    else{
      SRdf.val <- data.frame(y=scale(y[val], center=mean(y[cal]),scale=F), X=I(scale(X[val,whichSR], center=colMeans(X[cal,whichSR]),scale=F)))
    }
    SRpls <- plsr(y ~ X, data=SRdf.cal,subset=cal,ncomp=min(j,Ncomp))
    valpreds <- predict(SRpls, newdata = SRdf.val, ncomp=1:min(j,Ncomp))
    SRrmsepval <- c(1:min(j,Ncomp))
    for(L in 1:min(j,Ncomp)){
      SRrmsepval[L] <- sqrt(mean(yconv[val]^2*(valpreds[,,L]-SRdf.val$y)^2))
    }
    if(anyNA(SRrmsepval==T)){
      watNaNs <- which(is.nan(SRrmsepval))
      SRrmsepval[watNaNs] <- 222
    }
    if(length(SRrmsepval) < Ncomp){
      SRlen <- length(SRrmsepval)
      fill <- rep(100,(Ncomp-SRlen))
      RMSEPval[j,] <- c(SRrmsepval,fill)
    }
    else { RMSEPval[j,] <- SRrmsepval }
    SRoptLVs[j] <- which.min(SRrmsepval)
  }
  minRMSEP <- min(RMSEPval)
  nvars <- which(RMSEPval==min(RMSEPval),arr.ind=T)[1,1]
  if(nvars==J){
    minRMSEP <- min(RMSEPval[-J,])
    nvars <- which(RMSEPval==minRMSEP,arr.ind=T)[1,1]
  }
  RMSEopt <- minRMSEP[1]
  optnlv <- SRoptLVs[nvars]
  vars <- sort(SRrank[1:nvars])
  cutoff <- SR[SRrank[nvars]]
  df <- data.frame(y = y, X = I(X))
  SRpls.t <- plsr(y ~ X[,vars], data=df,subset=cal, ncomp=optnlv)
  valpreds <- predict(SRpls.t, newdata = df[val,], ncomp=optnlv)
  PRESSopt= sum((valpreds-df$y[val])^2) #sum(yconv[val]^2*(valpreds-y[val])^2)
  testpreds <- predict(SRpls.t, newdata = df[test,], ncomp=optnlv)
  RMSEPtest <- sqrt(mean(yconv[test]^2*(testpreds[,,1]-df$y[test])^2))
  PRESStest= sum(yconv[test]^2*(testpreds-df$y[test])^2) #sum(yconv[test]^2*(testpreds-df$y[test])^2)
  return(list(RMSEopt=RMSEopt,RMSEPtest=RMSEPtest,nvars=nvars,vars=vars,optnlv=optnlv,
              PRESSopt=PRESSopt,PRESStest=PRESStest,SR=SR,cutoff=cutoff))
}