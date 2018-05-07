splsreg <- function(X,y,cal,val,test,Ncomp,etavec,type,yconv,sc=F){#,coefs=F,nvar=F){
  #This function is a wrapper function to optimize the Sparse Partial Least Squares Regression (SPLS)
  #method (see citation below)
  #Chun H, Keles S. Sparse partial least squares regression for simultaneous dimension 
  #reduction and variable selection. J R Stat Soc Series B. 2010;72(1):3-25.
  #Inputs:
  #X = matrix of data
  #y = vector of property
  #cal, val, test: vectors for calibration, validation, test sets
  #etavec= vector of values for the SPLS thresholding parameter eta. Values must be between 0 and 1.
  #Ncomp = Number of PLS components
  #type = Which version of SPLS to use? Choices are "simpls"(default) and "nipals" (See spls package documentation)
  #yconv = conversion factor for calculating RMSEPs for y
  #sc= Is X autoscales?
  #Outputs:
  #RMSEopt=Optimal RMSEopt for SPLS
  #RMSEPtest=RMSEP for test set using optimal SPLS model
  #nvars=Number of variables selected by SPLS
  #vars=vecter of variable indices (channel numbers) for the variables selected by SPLS
  #optnlv=optimal number of latent variables chosen 
  #PRESSopt=Predicted Residual Error Sum of Squares for optimization set
  #PRESStest=Predicted Residual Error Sum of Squares for test set
  #eta = Value of eta used for optimal SPLS model
  #coefs = SPLS regression coefficients
  #R code by Cannon Giglio, 2016-2018, Brown Research Group, University of Delaware
  #type=c("simpls","nipals")
  if(missing(type) | type=="simpls"){splstype <- "simpls"}
  if(type=="nipals"){splstype <- "pls2"}
  if(missing(yconv)){
    yconv <- rep(1,length(y))
  }
  else{ splstype <- type}
  require(pls)
  require(spls)
  Kvec <- as.vector(1:Ncomp)
  Yval <- scale(y[val], center=mean(y[cal]),scale=F)
  RMSEPval <- matrix(100, nrow=length(etavec),ncol=Ncomp,dimnames=list(etavec,Kvec))
  Nvarmat <- matrix(100, nrow=length(etavec),ncol=Ncomp,dimnames=list(etavec,Kvec))
  for(j in 1:length(etavec)){
    for(k in 1:Ncomp){
      Spls <- spls:::spls(X[cal,],y[cal],K=Kvec[k],eta=etavec[j],kappa=0.5,select=splstype,scale.x=sc,scale.y=F)
      Splscoef <- coef(Spls)
      Nvarmat[j,k] <- as.vector(apply(Splscoef, 2, function(z)sum(z!=0)))
      spredval <- predict(Spls,newx=X[val,],type="fit")
      RMSEPval[j,k] <- sqrt(mean(yconv[val]^2*(y[val]-spredval)^2))
    }
  }
  RMSEopt <- min(RMSEPval)
  bestspls <- which(RMSEPval==min(RMSEPval),arr.ind=T)[1,1:2] #bestspls[1:2] #1=eta, 2=Ncomps
  nvars <- Nvarmat[bestspls[1],bestspls[2]]
  eta <- etavec[bestspls[1]]
  optnlv <- Kvec[bestspls[2]]
  SPLSbest <- spls(X[cal,],y[cal],K=Kvec[bestspls[2]],eta=etavec[bestspls[1]],kappa=0.5,select="simpls",scale.x=sc,scale.y=F)
  coefs <- coef(SPLSbest)
  vars <- which(coefs!=0)
  df <- data.frame(y = y, X = I(X))
  spls.t <- plsr(y ~ X[,vars], data=df,subset=cal, ncomp=optnlv)
  valpreds <- predict(spls.t, newdata = df[val,], ncomp=optnlv)
  PRESSopt= sum((valpreds-df$y[val])^2) #sum(yconv[val]^2*(valpreds-y[val])^2)
  testpreds <- predict(spls.t, newdata = df[test,], ncomp=optnlv)
  RMSEPtest <- sqrt(mean(yconv[test]^2*(testpreds[,,1]-df$y[test])^2))
  PRESStest= sum(yconv[test]^2*(testpreds-df$y[test])^2) #sum(yconv[test]^2*(testpreds-df$y[test])^2)
  return(list(RMSEopt=RMSEopt,RMSEPtest=RMSEPtest,nvars=nvars,vars=vars,optnlv=optnlv,
              PRESSopt=PRESSopt,PRESStest=PRESStest,coefs=coefs,eta=eta))
  #return(list(Nvars=Nvars,RMSEopt=RMSEopt,nvars=nvars,eta=eta,optnlv=optnlv,
  #            bestcoef=bestcoef,vars=vars,RMSEPtest=RMSEPtest))
}