#Reproducible code from the paper "Using elastic net regression to perform spectrally relevant
#variable selection" by C. Giglio and S.D. Brown, J. Chemometrics, 2018
#R code by Cannon Giglio, 2016-2018, Brown Research Group, University of Delaware

#install the following packages:
install.packages("pls")
install.packages("glmnet")
install.packages("scales")
install.packages("spls")
install.packages("R.matlab")

#load packages
library(pls)
library(glmnet)
library(spls)
library(scales)
library(R.matlab)

#Source the following functions:
source('enpls.R')
source('enetreg.R')
source('VIPpls.R')
source('SRpls.R')
source('splsreg.R')
source('rectbands.R')

#Move the Data and functions into your current workspace directory
#Beer dataset
#Download iPLS toolbox to get dataset: http://www.models.life.ku.dk/iToolbox
#load the file 'nirbeer.mat' into current workspace directory
beerdat <-readMat('nirbeer.mat') #if nirbeer.mat is not in your current workspace directory,
#use readMat('/yourdirectory/nirbeer.mat') to read the file
cbergX0 <- rbind(beerdat$Xcal,beerdat$Xtest)
cbergY0 <- rbind(beerdat$ycal,beerdat$ytest)
beerxax <- as.vector(beerdat$xaxis)

cbcal <- c(seq(1,23,2),seq(27,40,2),40) #calibration set: odd samples + sample 40
cbopt <- c(seq(2,24,2),25,seq(26,38,2)) #optimization set: even samples + sample 25 
#test set: samples 41:60
cbs <- c(cbcal,cbopt,41:60) #store indices of cal,opt, test sets
cbergX <- cbergX0[cbs,]
cbergY <- cbergY0[cbs] 

#Shootout Dataset
#The Shootout data was formerly accessible via the package "ChemometricsWithRData"
#However, the package has been removed from CRAN (as of May 1 2018). 
#The file 'shootout.RData' was originally obtained from the ChemometricsWithRData package before its removal
load("shootout.RData")
ShootoutX1 <- rbind(shootout$calibrate.1,shootout$validate.1,shootout$test.1)
#ShootoutX2 <- rbind(shootout$calibrate.2,shootout$validate.2,shootout$test.2) #not used for present work
ShootoutY <- rbind(shootout$calibrate.Y,shootout$validate.Y,shootout$test.Y)
tramadolX <- msc(ShootoutX1,reference=colMeans(ShootoutX1[1:155,]))
tramadolY <- (ShootoutY[,3]/ShootoutY[,1])*100 #Y= mass % of API
tyconv <- ShootoutY[,1]/100 #used to convert mass% to mass
tramxax <- seq(600,1898,2)
trcal <-  c(1:155)
tropt <-  c(196:445)
trtest <- c(156:195,446:655) 


#Test the variable selection methods
############################################
#Beer Data
enavec <- c(0,0.04,seq(0.08,1,0.08),1) #alpha values to use for EN and EN-PLS
enlamvec <- seq(5,100,5) #lambda indices to use for EN and EN-PLS
Cbpls <- enpls(cbergX,cbergY,avec=0,lamvec=seq(50),1:20,21:40,41:60,10,scale=F) #PLS (EN-PLS w/alpha=0 is the same result as usual PLS)
Cbenpls.mc <- enpls(cbergX,cbergY,avec=enavec,lamvec=enlamvec,1:20,21:40,41:60,10,scale=F) #EN-PLS (MC)
Cbenpls.sc <- enpls(cbergX,cbergY,avec=enavec,lamvec=enlamvec,1:20,21:40,41:60,10,scale=T) #EN-PLS s/c
Cbenet.sc <- enetreg(cbergX,cbergY,avec=enavec,lamvec=enlamvec,1:20,21:40,41:60,scale=T) #EN
CbVIP <- VIPpls(cbergX,cbergY,1:20,21:40,41:60,15,sc=F) 
CbSR <- SRpls(cbergX,cbergY,1:20,21:40,41:60,10,sc=F)
etaseq <- c(seq(0,.8,.05),seq(0.81,0.99,0.01)) 
Cbspls <- splsreg(cbergX,cbergY,1:20,21:40,41:60,10,etaseq,type="simpls",sc=F)

#RMSEopt
CbRMSEopt<- round(c(Cbpls$RMSEopt,Cbenpls.mc$RMSEopt,Cbenpls.sc$RMSEopt,Cbenet.sc$RMSEopt,CbVIP$RMSEopt,Cbspls$RMSEopt,CbSR$RMSEopt),3)
#RMSEPtest
CbRMSEPtest <- round(c(Cbpls$RMSEPtest,Cbenpls.mc$RMSEPtest,Cbenpls.sc$RMSEPtest,Cbenet.sc$RMSEPtest,CbVIP$RMSEPtest,Cbspls$RMSEPtest,CbSR$RMSEPtest),3)
#Number of Variables Selected
Cbnvars<- c(Cbpls$nvars,Cbenpls.mc$nvars,Cbenpls.sc$nvars,Cbenet.sc$nvars,CbVIP$nvars,Cbspls$nvars,CbSR$nvars)
#Number of Latent Variables
Cbvars<- c(Cbpls$vars,Cbenpls.mc$vars,Cbenpls.sc$vars,Cbenet.sc$vars,CbVIP$vars,Cbspls$vars,CbSR$vars)
#Tuning Parameters
Cbenpls.mc$alflam
Cbenpls.sc$alflam
Cbenet.sc$enalflam
CbVIP$cutoff
CbSR$cutoff
Cbspls$eta
#F-tests
Cbpresstest <- c(Cbpls$PRESStest,Cbenpls.mc$PRESStest,Cbenpls.sc$PRESStest,Cbenet.sc$PRESStest,CbVIP$PRESStest,
                Cbspls$PRESStest,CbSR$PRESStest)
CbFmat <- matrix(0,nrow=7,ncol=7)
rownames(CbFmat) <- c("PLS","EN-PLSmc","EN-PLSsc","EN","VIP","SPLS","SR")
colnames(CbFmat) <- c("PLS","EN-PLSmc","EN-PLSsc","EN","VIP","SPLS","SR")
for(i in 1:7){
  for(j in 1:7){
    if(i!=j){
      if(Cbpresstest[i]>Cbpresstest[j]){CbFmat[i,j] <- Cbpresstest[i]/Cbpresstest[j]}
      else{ CbFmat[i,j] <- Cbpresstest[j]/Cbpresstest[i]}
    }
    else{ CbFmat[i,j] <- 0} 
  }
}
round(CbFmat,3)

#Generate Figures
#par(mfrow=c(3,2))
matplot(beerxax,t(cbergX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Cberg (scaled) - ENPLS")
rectbands(Cbenpls.sc$vars,waveseq=beerxax,ymin=-2,ymax=10,width=0.5,plotcol="firebrick",alfa=0.5)
matplot(beerxax,t(cbergX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Cberg (mean-centered) - ENPLS")
rectbands(Cbenpls.mc$vars,waveseq=beerxax,ymin=-2,ymax=10,width=0.5,plotcol="darkorange2",alfa=0.5)
matplot(beerxax,t(cbergX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Cberg - EN sc")
rectbands(Cbenet.sc$vars,waveseq=beerxax,ymin=-2,ymax=10,width=0.5,plotcol="darkgoldenrod2",alfa=0.5) #purple
matplot(beerxax,t(cbergX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Cberg - VIP")
rectbands(sort(CbVIP$vars),waveseq=beerxax,ymin=-2,ymax=10,width=0.5,plotcol="forestgreen",alfa=0.5)
matplot(beerxax,t(cbergX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Cberg - SPLS CK")
rectbands(sort(Cbspls$vars),waveseq=beerxax,ymin=-2,ymax=10,width=0.5,plotcol="royalblue",alfa=0.5)
matplot(beerxax,t(cbergX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Cberg - SR")
rectbands(sort(CbSR$vars),waveseq=beerxax,ymin=-2,ymax=10,width=0.5,plotcol="purple3",alfa=0.5)

######################################
#Tablet Data
Trpls <-enpls(tramadolX,tramadolY,avec=0,lamvec=seq(50),trcal,tropt,trtest,10,scale=F,yconv=tyconv)
enavec <- c(0,0.04,seq(0.08,1,0.08),1)  #alpha values to use for EN and EN-PLS
enlamvec <- seq(5,100,5) #lambda indices to use for EN and EN-PLS
Trenpls.mc <-enpls(tramadolX,tramadolY,avec=enavec,lamvec=enlamvec,trcal,tropt,trtest,10,scale=F,yconv=tyconv)
Trenpls.sc <-enpls(tramadolX,tramadolY,avec=enavec,lamvec=enlamvec,trcal,tropt,trtest,10,scale=T,yconv=tyconv)
Trenet.sc <- enetreg(tramadolX,tramadolY,avec=enavec,lamvec=enlamvec,trcal,tropt,trtest,scale=T,yconv=tyconv)
TrVIP <- VIPpls(tramadolX,tramadolY,trcal,tropt,trtest,10,yconv=tyconv,sc=F)
etaseq <- c(seq(0,.8,.05),seq(0.81,0.99,0.01))
Trspls <- splsreg(tramadolX,tramadolY,trcal,tropt,trtest,10,etaseq,type="simpls",yconv=tyconv,sc=F)
TrSR <- SRpls(tramadolX,tramadolY,trcal,tropt,trtest,10,yconv=tyconv,sc=F)

#RMSEopt:
TrRMSEopt<- round(c(Trpls$RMSEopt,Trenpls.mc$RMSEopt,Trenpls.sc$RMSEopt,Trenet.sc$RMSEopt,TrVIP$RMSEopt,Trspls$RMSEopt,TrSR$RMSEopt),3)
#RMSEPtest:
TrRMSEPtest <- round(c(Trpls$RMSEPtest,Trenpls.mc$RMSEPtest,Trenpls.sc$RMSEPtest,Trenet.sc$RMSEPtest,TrVIP$RMSEPtest,Trspls$RMSEPtest,TrSR$RMSEPtest),3)
#Number of Variables Selected:
Trnvars<- c(Trpls$nvars,Trenpls.mc$nvars,Trenpls.sc$nvars,Trenet.sc$nvars,TrVIP$nvars,Trspls$nvars,TrSR$nvars)
#Number of Latent Variables:
Trnlv<- c(Trpls$optnlv,Trenpls.mc$optnlv,Trenpls.sc$optnlv,NA,TrVIP$optnlv,Trspls$optnlv,TrSR$optnlv)
#Tuning Parameters:
Trenpls.mc$alflam
Trenpls.sc$alflam
Trenet.sc$enalflam
TrVIP$cutoff
Trspls$eta
TrSR$cutoff
#F-tests to compare PRESS values:
Trpresstest <- c(Trpls$PRESStest,Trenpls.mc$PRESStest,Trenpls.sc$PRESStest,Trenet.sc$PRESStest,TrVIP$PRESStest,
                 Trspls$PRESStest,TrSR$PRESStest)
TrFmat <- matrix(0,nrow=7,ncol=7)
rownames(TrFmat) <- c("PLS","EN-PLSmc","EN-PLSsc","EN","VIP","SPLS","SR")
colnames(TrFmat) <- c("PLS","EN-PLSmc","EN-PLSsc","EN","VIP","SPLS","SR")
for(i in 1:7){
  for(j in 1:7){
    if(i!=j){
      if(Trpresstest[i]>Trpresstest[j]){TrFmat[i,j] <- Trpresstest[i]/Trpresstest[j]}
      else{ TrFmat[i,j] <- Trpresstest[j]/Trpresstest[i]}
    }
    else{ TrFmat[i,j] <- 0} 
  }
}
round(TrFmat,3)

#Generate Figures
#par(mfrow=c(3,2))
matplot(tramxax,t(tramadolX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Tramadol (scaled) - ENPLS")
rectbands(Trenpls.sc$vars,waveseq=tramxax,ymin=-2,ymax=10,width=0.5,plotcol="firebrick",alfa=0.5)
matplot(tramxax,t(tramadolX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Tramadol (mean-centered) - ENPLS")
rectbands(Trenpls.mc$vars,waveseq=tramxax,ymin=-2,ymax=10,width=0.5,plotcol="darkorange2",alfa=0.5)
matplot(tramxax,t(tramadolX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Tramadol - EN sc")
rectbands(Trenet.sc$vars,waveseq=tramxax,ymin=-2,ymax=10,width=0.5,plotcol="darkgoldenrod2",alfa=0.5) #purple
matplot(tramxax,t(tramadolX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Tramadol - VIP")
rectbands(sort(TrVIP$vars),waveseq=tramxax,ymin=-2,ymax=10,width=0.5,plotcol="forestgreen",alfa=0.5)
matplot(tramxax,t(tramadolX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Tramadol - SPLS")
rectbands(sort(Trspls$vars),waveseq=tramxax,ymin=-2,ymax=10,width=0.5,plotcol="royalblue",alfa=0.5)
matplot(tramxax,t(tramadolX),type="l",lty=1,col="gray50",xlab="Wavelength (nm)",ylab="Absorbance")
title(main="Tramadol - SR")
rectbands(sort(TrSR$vars),waveseq=tramxax,ymin=-2,ymax=10,width=0.5,plotcol="purple3",alfa=0.5)