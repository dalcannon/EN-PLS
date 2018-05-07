
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
source('enpls Function.R')
source('enetreg Function.R')
source('VIPpls Function.R')
source('SRpls Function.R')
source('splsreg function.R')
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