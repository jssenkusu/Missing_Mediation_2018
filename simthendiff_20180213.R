
################################################################################
# By John Ssenkusu (ssenk001@umn.edu, jssenkusu@gmail.com)
# Purpose: Simulation study code for missing data in mediation analysis in
#          longitudinal studies
# Created: February 10, 2017
# Modified: February 13, 2018
################################################################################

library(foreach)
library(doParallel)
library(utils)

source('C:/Users/name/R code/datagen_mix_mar_20180213.R')

# SOBEL TEST

sobel.test <- function(a, b, var.a, var.b)
{
  sobel.var <- a*a*var.b + b*b*var.a
  sobel.var
}

True.ACME <- 2.4*-0.9
True.ADE <- -1.5
  
nsim <- 2000
prop.miss <- c(0.1, 0.2, 0.3, 0.4)
n.miss <- length(prop.miss)

cores=detectCores()
cl <- makeCluster(cores-1) # not to overload your computer
# cl <- makeCluster(cores)
registerDoParallel(cl)


strt <- Sys.time() # start time


df <- foreach(i=1:nsim, .combine = "rbind", .export=ls(envir=globalenv()))  %:%  
    foreach(j=1:n.miss, .combine = "cbind", .multicombine=TRUE)  %dopar%
    {      
    library(MASS)
    library(pan)
    library(mice)
    library(mitools)
    library(lme4)
    library(mlmmm)
   
    library(reshape)
    library(simglm)
      

    data.complete <- as.data.frame(simdata(prop.miss[j]))

    data.complete_wide <- reshape(data.complete, idvar = "id", timevar = "time", 
                                   v.names=c("M","Y", "M.mcar","Y.mcar","M.mar","Y.mar"), 
                                  direction = "wide")
    
    data.complete_wide$M6_M1 <- data.complete_wide$M.6 - data.complete_wide$M.1
    
    ############################################################################
    
    ### ESTIMATES WHEN DATA IS COMPLETE
    
    nomiss.M <- lm(M6_M1 ~ group + X1 + X2, data = data.complete_wide)
    nomiss.Y <- lm(Y.6 ~ group + M6_M1 + X1 + X2, data = data.complete_wide)
    
    ACME.nomiss <- coefficients(nomiss.M)[2] * coefficients(nomiss.Y)[3]
    ADE.nomiss <- coefficients(nomiss.Y)[2]    
    

    ADE_SE_nomiss <- sqrt(diag(vcov(nomiss.Y))[2])
    
    ACME_nomiss_SE <- sqrt(sobel.test(a=coefficients(nomiss.M)[2], b=coefficients(nomiss.Y)[3], 
                                 var.a=diag(vcov(nomiss.M))[2], var.b=diag(vcov(nomiss.Y))[3]))
    ACME_cov_nomiss <- ACME.nomiss-ACME_nomiss_SE*1.96 <=True.ACME & 
                       ACME.nomiss+ACME_nomiss_SE*1.96 >= True.ACME # Does CI include true estimate?
    
    ADE_cov_nomiss <- ADE.nomiss-ADE_SE_nomiss*1.96 <=True.ADE & 
                      ADE.nomiss+ADE_SE_nomiss*1.96 >= True.ADE 
    
    ############################################################################
    
    ### ESTIMATES WHEN NO IMPUTATION IS CONDUCTED
    
    ### missing cases are excluded from models 
    
    data.complete_wide$mcar_M6_M1 <- data.complete_wide$M.mcar.6 - data.complete_wide$M.mcar.1 
    data.complete_wide$mar_M6_M1 <- data.complete_wide$M.mar.6 - data.complete_wide$M.mar.1
    
    ### MCAR
    
    noimp.mcar_M <- lm(mcar_M6_M1 ~ group + X1 + X2, data = data.complete_wide)
    noimp.mcar_Y <- lm(Y.mcar.6 ~ group + mcar_M6_M1 + X1 + X2, data = data.complete_wide)

    ACME.noimp.mcar <- coefficients(noimp.mcar_M)[2] * coefficients(noimp.mcar_Y)[3]
    ADE.noimp.mcar <- coefficients(noimp.mcar_Y)[2]    
 
    ADE_SE_noimp.mcar <- sqrt(diag(vcov(noimp.mcar_Y))[2])
    
    ACME_noimp.mcar_SE <- sqrt(sobel.test(a=coefficients(noimp.mcar_M)[2], b=coefficients(noimp.mcar_Y)[3], 
                                     var.a=diag(vcov(noimp.mcar_M))[2], var.b=diag(vcov(noimp.mcar_Y))[3]))
    ACME_cov_noimp.mcar <- ACME.noimp.mcar-ACME_noimp.mcar_SE*1.96 <=True.ACME & 
                            ACME.noimp.mcar+ACME_noimp.mcar_SE*1.96 >= True.ACME # Does CI include true estimate?
    
    ADE_cov_noimp.mcar <- ADE.noimp.mcar-ADE_SE_noimp.mcar*1.96 <=True.ADE & 
                          ADE.noimp.mcar+ADE_SE_noimp.mcar*1.96 >= True.ADE 
    
    ### MAR
    
    noimp.mar_M <- lm(mar_M6_M1 ~ group + X1 + X2, data = data.complete_wide)
    noimp.mar_Y <- lm(Y.mar.6 ~ group + mar_M6_M1 + X1 + X2, data = data.complete_wide)

    ACME.noimp.mar <- coefficients(noimp.mar_M)[2] * coefficients(noimp.mar_Y)[3]
    ADE.noimp.mar <- coefficients(noimp.mar_Y)[2]    
    
    
    ADE_SE_noimp.mar <- sqrt(diag(vcov(noimp.mar_Y))[2])
    
    ACME_noimp.mar_SE <- sqrt(sobel.test(a=coefficients(noimp.mar_M)[2], b=coefficients(noimp.mar_Y)[3], 
                                    var.a=diag(vcov(noimp.mar_M))[2], var.b=diag(vcov(noimp.mar_Y))[3]))
    ACME_cov_noimp.mar <- ACME.noimp.mar-ACME_noimp.mar_SE*1.96 <=True.ACME & 
                          ACME.noimp.mar+ACME_noimp.mar_SE*1.96 >= True.ACME # Does CI include true estimate?
    
    ADE_cov_noimp.mar <- ADE.noimp.mar-ADE_SE_noimp.mar*1.96 <=True.ADE & 
                         ADE.noimp.mar+ADE_SE_noimp.mar*1.96 >= True.ADE  
    
    ############################################################################
    
    ### MULTIPLE IMPUTATION CONSIDERING CLUSTERING USING 'PAN'

  
############## Missing completely at random (MCAR)  #############
  
mcar.data <- data.complete[,c("group","X1","X2")] 
mcar.X <- cbind(mcar.data[,0, drop=FALSE], 1, mcar.data[,1:3] ) 
names(mcar.X)[1] <- "int"

time <- model.matrix( ~ as.factor(time), data=data.complete)[,2:6]

pred <- cbind(mcar.X$int, mcar.X$group, time, mcar.X$X1, mcar.X$X2, mcar.X$group*time)

xcol <- 1:14
zcol <- 1
Y.mcar <- as.vector(data.complete$Y.mcar)
M.mcar <- as.vector(data.complete$M.mcar)
YM.mcar <- cbind(Y.mcar, M.mcar)
subj <- as.vector(data.complete$id)

# Conducting the E-M algorithmn to get priors for the error variance and REs variance
# Choice of priors is guided by Joseph Schafer (2001) - Multiple imputation with PAN

model.em_mcar <- mlmmm.em(y=YM.mcar, subj, pred, xcol, zcol, maxits = 400, eps = 0.0001)
a = r = c = 2 
Binv_mcar <- a * model.em_mcar$sigma
Dinv_mcar <- c * model.em_mcar$psi

prior_mcar <- list(a=a, Binv=Binv_mcar, c=c, Dinv=Dinv_mcar) # setting priors

# # Conducting diagnostics to determine the burn-in
# 
# mcar.YM_pan <- pan(y=YM.mcar, subj, pred, xcol, zcol, prior_mcar, seed=194,iter=1000)
# 
# plot(1:500,log(mcar.YM_pan$sigma[1,1,]),type="l", main="Variance for Y")
# acf(log((mcar.YM_pan$sigma[1,1,])))
# 
# plot(1:500,log(mcar.YM_pan$sigma[2,2,]),type="l", main="Variance for M")
# acf(log((mcar.YM_pan$sigma[2,2,])))
# 
# par(mar = rep(2, 4))
# par(mfrow=c(5,3))
# for(i in 1:14) plot(1:1000,mcar.YM_pan$beta[i,1,],type="l")
# par(mfrow=c(5,3))
# for(i in 1:14) acf(mcar.YM_pan$beta[i,1,])

# dev.off()

# Imputing 10 data sets and storing values for the outcome and the mediator

mcar.YM_pan <- pan(y=YM.mcar,subj, pred, xcol, zcol, prior_mcar,seed=i+194,iter=5000)
mcar.YM1 <- mcar.YM_pan$y
mcar.YM_pan <- pan(y=YM.mcar,subj, pred, xcol, zcol, prior_mcar,seed=i+257,iter=100,start=mcar.YM_pan$last)
mcar.YM2 <- mcar.YM_pan$y
mcar.YM_pan <- pan(y=YM.mcar,subj,pred, xcol, zcol, prior_mcar,seed=i+483,iter=100,start=mcar.YM_pan$last)
mcar.YM3 <- mcar.YM_pan$y
mcar.YM_pan <- pan(y=YM.mcar,subj,pred, xcol, zcol, prior_mcar,seed=i+106,iter=100,start=mcar.YM_pan$last)
mcar.YM4 <- mcar.YM_pan$y
mcar.YM_pan <- pan(y=YM.mcar,subj,pred, xcol, zcol, prior_mcar,seed=i+835,iter=100,start=mcar.YM_pan$last)
mcar.YM5 <- mcar.YM_pan$y
mcar.YM_pan <- pan(y=YM.mcar,subj,pred, xcol, zcol, prior_mcar,seed=i+876,iter=100,start=mcar.YM_pan$last)
mcar.YM6 <- mcar.YM_pan$y
mcar.YM_pan <- pan(y=YM.mcar,subj,pred, xcol, zcol, prior_mcar,seed=i+661,iter=100,start=mcar.YM_pan$last)
mcar.YM7 <- mcar.YM_pan$y
mcar.YM_pan <- pan(y=YM.mcar,subj,pred, xcol, zcol, prior_mcar,seed=i+163,iter=100,start=mcar.YM_pan$last)
mcar.YM8 <- mcar.YM_pan$y
mcar.YM_pan <- pan(y=YM.mcar,subj,pred, xcol, zcol, prior_mcar,seed=i+989,iter=100,start=mcar.YM_pan$last)
mcar.YM9 <- mcar.YM_pan$y
mcar.YM_pan <- pan(y=YM.mcar,subj,pred, xcol, zcol, prior_mcar,seed=i+979,iter=100,start=mcar.YM_pan$last)
mcar.YM10 <- mcar.YM_pan$y
          

d.mcar.Y1 <- data.frame(y=mcar.YM1[,1],subj,pred,m=mcar.YM1[,2])
d.mcar.Y2 <- data.frame(y=mcar.YM2[,1],subj,pred,m=mcar.YM2[,2])
d.mcar.Y3 <- data.frame(y=mcar.YM3[,1],subj,pred,m=mcar.YM3[,2])
d.mcar.Y4 <- data.frame(y=mcar.YM4[,1],subj,pred,m=mcar.YM4[,2])
d.mcar.Y5 <- data.frame(y=mcar.YM5[,1],subj,pred,m=mcar.YM5[,2])
d.mcar.Y6 <- data.frame(y=mcar.YM6[,1],subj,pred,m=mcar.YM6[,2])
d.mcar.Y7 <- data.frame(y=mcar.YM7[,1],subj,pred,m=mcar.YM7[,2])
d.mcar.Y8 <- data.frame(y=mcar.YM8[,1],subj,pred,m=mcar.YM8[,2])
d.mcar.Y9 <- data.frame(y=mcar.YM9[,1],subj,pred,m=mcar.YM9[,2])
d.mcar.Y10 <- data.frame(y=mcar.YM10[,1],subj,pred,m=mcar.YM10[,2])

# combining all 10 imputed data sets

d.mcar_Y <- imputationList(list(d.mcar.Y1, d.mcar.Y2, d.mcar.Y3,d.mcar.Y4, d.mcar.Y5, d.mcar.Y6, d.mcar.Y7, 
                                d.mcar.Y8, d.mcar.Y9,d.mcar.Y10)) 
# removing data sets from memory

rm(mcar.YM1,mcar.YM2,mcar.YM3,mcar.YM4,mcar.YM5,mcar.YM6,mcar.YM7,mcar.YM8,mcar.YM9,mcar.YM10,
   d.mcar.Y1,d.mcar.Y2,d.mcar.Y3,d.mcar.Y4,d.mcar.Y5,d.mcar.Y6,d.mcar.Y7,d.mcar.Y8,d.mcar.Y9,d.mcar.Y10)

# Function that renames the variable names in imputed data sets

rename_Y_dataset <- function(x) {
  names(x) <- c("Y", "id", "int","group","time.2","time.3","time.4","time.5","time.6","X1","X2",
                "grp_time.2","grp_time.3","grp_time.4","grp_time.5","grp_time.6","M")
  return(x)
}

# Fitting LM for 10 data sets and combining the results from the analysis 

ACME.mcar_pan <- rep(NA, 10)
ADE.mcar_pan <- rep(NA, 10)
sobel_mcar_pan <- rep(NA, 10)
ADE.var_mcar_pan <- rep(NA, 10)

for(k in 1:10)
{
  data.mcar <- as.data.frame(d.mcar_Y$imputations[k])
  data.mcar <- rename_Y_dataset(data.mcar)
  data.mcar <- data.mcar[,c("Y", "id", "group","X1","X2","M")]
  time <- data.complete$time
  data.mcar <- cbind(data.mcar, time)
  
  data.mcar_wide <- reshape(data.mcar, idvar = "id", timevar = "time", v.names=c("M","Y"), 
                            direction = "wide")
  
  data.mcar_wide$M6_M1 <- data.mcar_wide$M.6 - data.mcar_wide$M.1
  
  lm_mcar.M <- lm(M6_M1 ~ group + X1 + X2, data = data.mcar_wide)
  lm_mcar.Y <- lm(Y.6 ~ group + M6_M1 + X1 + X2, data = data.mcar_wide)
  
  ACME.mcar_pan[k] <- coefficients(lm_mcar.M)[2] * coefficients(lm_mcar.Y)[3]
  ADE.mcar_pan[k] <- coefficients(lm_mcar.Y)[2]  
  
  sobel_mcar_pan[k] <- sobel.test(a=coefficients(lm_mcar.M)[2], b=coefficients(lm_mcar.Y)[3], 
                                  var.a=diag(vcov(lm_mcar.M))[2], var.b=diag(vcov(lm_mcar.Y))[3])
  ADE.var_mcar_pan[k] <- diag(vcov(lm_mcar.Y))[2]
}

  ACME_mcar_pan <- mean(ACME.mcar_pan)
  ADE_mcar_pan <- mean(ADE.mcar_pan)
  ACME_mcar_pan_SE <- sqrt(mean(sobel_mcar_pan) + (1+(1/10))*var(ACME.mcar_pan))
  
  ACME_cov_mcar_pan <- ACME_mcar_pan-ACME_mcar_pan_SE*1.96 <= True.ACME & 
                       ACME_mcar_pan+ACME_mcar_pan_SE*1.96 >= True.ACME # Does CI include true estimate?
  
  ADE_SE_mcar_pan <- sqrt(mean(ADE.var_mcar_pan) + (1+(1/10))*var(ADE.mcar_pan))
  
  ADE_cov_mcar_pan <- ADE_mcar_pan-ADE_SE_mcar_pan*1.96 <= True.ADE & 
                      ADE_mcar_pan+ADE_SE_mcar_pan*1.96 >= True.ADE 

############## Missing at random (MAR)  ##############

Y.mar <- as.vector(data.complete$Y.mar)
M.mar <- as.vector(data.complete$M.mar)
YM.mar <- cbind(Y.mar, M.mar)

# Using E-M algorithmn to get priors for error variance and RE variance

model.em_mar <- mlmmm.em(y=YM.mar, subj, pred, xcol, zcol, maxits = 400, eps = 0.0001)
a = r = c = 2 
Binv_mar <- a * model.em_mar$sigma
Dinv_mar <- c * model.em_mar$psi

prior_mar <- list(a=a, Binv=Binv_mar, c=c, Dinv=Dinv_mar) # setting priors

# Imputing 10 data sets and storing values for the outcome and the mediator

mar.YM_pan <- pan(y=YM.mar,subj,pred, xcol, zcol, prior_mar,seed=i+387,iter=5000)
mar.YM1 <- mar.YM_pan$y
mar.YM_pan <- pan(y=YM.mar,subj,pred, xcol, zcol, prior_mar,seed=i+579,iter=100,start=mar.YM_pan$last)
mar.YM2 <- mar.YM_pan$y
mar.YM_pan <- pan(y=YM.mar,subj,pred, xcol, zcol, prior_mar,seed=i+621,iter=100,start=mar.YM_pan$last)
mar.YM3 <- mar.YM_pan$y
mar.YM_pan <- pan(y=YM.mar,subj,pred, xcol, zcol, prior_mar,seed=i+467,iter=100,start=mar.YM_pan$last)
mar.YM4 <- mar.YM_pan$y
mar.YM_pan <- pan(y=YM.mar,subj,pred, xcol, zcol, prior_mar,seed=i+612,iter=100,start=mar.YM_pan$last)
mar.YM5 <- mar.YM_pan$y
mar.YM_pan <- pan(y=YM.mar,subj,pred, xcol, zcol, prior_mar,seed=i+787,iter=100,start=mar.YM_pan$last)
mar.YM6 <- mar.YM_pan$y
mar.YM_pan <- pan(y=YM.mar,subj,pred, xcol, zcol, prior_mar,seed=i+489,iter=100,start=mar.YM_pan$last)
mar.YM7 <- mar.YM_pan$y
mar.YM_pan <- pan(y=YM.mar,subj,pred, xcol, zcol, prior_mar,seed=i+825,iter=100,start=mar.YM_pan$last)
mar.YM8 <- mar.YM_pan$y
mar.YM_pan <- pan(y=YM.mar,subj,pred, xcol, zcol, prior_mar,seed=i+751,iter=100,start=mar.YM_pan$last)
mar.YM9 <- mar.YM_pan$y
mar.YM_pan <- pan(y=YM.mar,subj,pred, xcol, zcol, prior_mar,seed=i+822,iter=100,start=mar.YM_pan$last)
mar.YM10 <- mar.YM_pan$y

d.mar.Y1 <- data.frame(y=mar.YM1[,1],subj,pred,m=mar.YM1[,2])
d.mar.Y2 <- data.frame(y=mar.YM2[,1],subj,pred,m=mar.YM2[,2])
d.mar.Y3 <- data.frame(y=mar.YM3[,1],subj,pred,m=mar.YM3[,2])
d.mar.Y4 <- data.frame(y=mar.YM4[,1],subj,pred,m=mar.YM4[,2])
d.mar.Y5 <- data.frame(y=mar.YM5[,1],subj,pred,m=mar.YM5[,2])
d.mar.Y6 <- data.frame(y=mar.YM6[,1],subj,pred,m=mar.YM6[,2])
d.mar.Y7 <- data.frame(y=mar.YM7[,1],subj,pred,m=mar.YM7[,2])
d.mar.Y8 <- data.frame(y=mar.YM8[,1],subj,pred,m=mar.YM8[,2])
d.mar.Y9 <- data.frame(y=mar.YM9[,1],subj,pred,m=mar.YM9[,2])
d.mar.Y10 <- data.frame(y=mar.YM10[,1],subj,pred,m=mar.YM10[,2])

# combining all 10 imputed data sets

d.mar_Y <- imputationList(list(d.mar.Y1, d.mar.Y2, d.mar.Y3,d.mar.Y4, d.mar.Y5, d.mar.Y6, d.mar.Y7, 
                                d.mar.Y8, d.mar.Y9,d.mar.Y10)) 

# removing data sets from memory

rm(mar.YM1,mar.YM2,mar.YM3,mar.YM4,mar.YM5,mar.YM6,mar.YM7,mar.YM8,mar.YM9,mar.YM10,
   d.mar.Y1,d.mar.Y2,d.mar.Y3,d.mar.Y4,d.mar.Y5,d.mar.Y6,d.mar.Y7,d.mar.Y8,d.mar.Y9,d.mar.Y10)

# Fitting LM for 10 data sets and combining the results from the analysis 

ACME.mar_pan <- rep(NA, 10)
ADE.mar_pan <- rep(NA, 10)
sobel_mar_pan <- rep(NA, 10)
ADE.var_mar_pan <- rep(NA, 10)

for(m in 1:10)
{
  data.mar <- as.data.frame(d.mar_Y$imputations[m])
  data.mar <- rename_Y_dataset(data.mar)
  data.mar <- data.mar[,c("Y", "id", "group","X1","X2","M")]
  time <- data.complete$time
  data.mar <- cbind(data.mar, time)
  
  data.mar_wide <- reshape(data.mar, idvar = "id", timevar = "time", v.names=c("M","Y"), 
                            direction = "wide")
  
  data.mar_wide$M6_M1 <- data.mar_wide$M.6 - data.mar_wide$M.1
  
  lm_mar.M <- lm(M6_M1 ~ group + X1 + X2, data = data.mar_wide)
  lm_mar.Y <- lm(Y.6 ~ group + M6_M1 + X1 + X2, data = data.mar_wide)
  
  ACME.mar_pan[m] <- coefficients(lm_mar.M)[2] * coefficients(lm_mar.Y)[3]
  ADE.mar_pan[m] <- coefficients(lm_mar.Y)[2]  
  
  sobel_mar_pan[m] <- sobel.test(a=coefficients(lm_mar.M)[2], b=coefficients(lm_mar.Y)[3], 
                                  var.a=diag(vcov(lm_mar.M))[2], var.b=diag(vcov(lm_mar.Y))[3])
  ADE.var_mar_pan[m] <- diag(vcov(lm_mar.Y))[2]
}

ACME_mar_pan <- mean(ACME.mar_pan)
ADE_mar_pan <- mean(ADE.mar_pan)
ACME_mar_pan_SE <- sqrt(mean(sobel_mar_pan) + (1+(1/10))*var(ACME.mar_pan))

ACME_cov_mar_pan <- ACME_mar_pan-ACME_mar_pan_SE*1.96 <= True.ACME & 
                    ACME_mar_pan+ACME_mar_pan_SE*1.96 >= True.ACME # Does CI include true estimate?

ADE_SE_mar_pan <- sqrt(mean(ADE.var_mar_pan) + (1+(1/10))*var(ADE.mar_pan))

ADE_cov_mar_pan <- ADE_mar_pan-ADE_SE_mar_pan*1.96 <=True.ADE & 
                   ADE_mar_pan+ADE_SE_mar_pan*1.96 >= True.ADE 

###############################################################

#### IGNORING CLUSTERING

numb_imput <- 10

################################################################

#### MULTIPLE IMPUTATION USING THE LINEAR MODEL METHOD
## 'MICE' is used in this imputation

data.complete_wide <- as.data.frame(data.complete_wide)

mcar_wide_mice_lm <- data.complete_wide[,c("group","X1","X2","M.mcar.1","M.mcar.6","Y.mcar.6")]
mar_wide_mice_lm <- data.complete_wide[,c("group","X1","X2","M.mar.1","M.mar.6","Y.mar.6")]

mcar_wide_mice_lm$mcar_M6_M1 <- mcar_wide_mice_lm$M.mcar.6 - mcar_wide_mice_lm$M.mcar.1
mar_wide_mice_lm$mar_M6_M1 <- mar_wide_mice_lm$M.mar.6 - mar_wide_mice_lm$M.mar.1

mcar_wide_mice_lm <- mcar_wide_mice_lm[,c("group","X1","X2","mcar_M6_M1","Y.mcar.6")]
mar_wide_mice_lm <- mar_wide_mice_lm[,c("group","X1","X2","mar_M6_M1","Y.mar.6")]

mcar_mice_lm <- mice(mcar_wide_mice_lm, m=numb_imput, maxit = 15,printFlag=FALSE, method = 'pmm', seed = 4+i)
mar_mice_lm <- mice(mar_wide_mice_lm, m=numb_imput, maxit = 15,printFlag=FALSE, method = 'pmm', seed = 14+i)

long.mcar_mice_lm_imp <- mice::complete(mcar_mice_lm, action='long', include=FALSE)
long.mar_mice_lm_imp <- mice::complete(mar_mice_lm, action='long', include=FALSE)

ACME.mcar_mice_lm <- rep(NA, 10)
ADE.mcar_mice_lm <- rep(NA, 10)
sobel_mcar_mice_lm <- rep(NA, 10)
ADE.var_mcar_mice_lm <- rep(NA, 10)

ACME.mar_mice_lm <- rep(NA, 10)
ADE.mar_mice_lm <- rep(NA, 10)
sobel_mar_mice_lm <- rep(NA, 10)
ADE.var_mar_mice_lm <- rep(NA, 10)

for (q in 1:10)
{
  data.mcar <- long.mcar_mice_lm_imp[which(long.mcar_mice_lm_imp$.imp==q), ]
  data.mar <- long.mar_mice_lm_imp[which(long.mar_mice_lm_imp$.imp==q), ]
  
  lm_mcar.M <- lm(mcar_M6_M1 ~ group + X1 + X2, data = data.mcar)
  lm_mcar.Y <- lm(Y.mcar.6 ~ group + mcar_M6_M1 + X1 + X2, data = data.mcar)
  
  lm_mar.M <- lm(mar_M6_M1 ~ group + X1 + X2, data = data.mar)
  lm_mar.Y <- lm(Y.mar.6 ~ group + mar_M6_M1 + X1 + X2, data = data.mar)
  
  ACME.mcar_mice_lm[q] <- coefficients(lm_mcar.M)[2] * coefficients(lm_mcar.Y)[3]
  ADE.mcar_mice_lm[q] <- coefficients(lm_mcar.Y)[2]  
  
  sobel_mcar_mice_lm[q] <- sobel.test(a=coefficients(lm_mcar.M)[2], b=coefficients(lm_mcar.Y)[3], 
                                   var.a=diag(vcov(lm_mcar.M))[2], var.b=diag(vcov(lm_mcar.Y))[3])
  ADE.var_mcar_mice_lm[q] <- diag(vcov(lm_mcar.Y))[2]
  
  ACME.mar_mice_lm[q] <- coefficients(lm_mar.M)[2] * coefficients(lm_mar.Y)[3]
  ADE.mar_mice_lm[q] <- coefficients(lm_mar.Y)[2]  
  
  sobel_mar_mice_lm[q] <- sobel.test(a=coefficients(lm_mar.M)[2], b=coefficients(lm_mar.Y)[3], 
                                      var.a=diag(vcov(lm_mar.M))[2], var.b=diag(vcov(lm_mar.Y))[3])
  ADE.var_mar_mice_lm[q] <- diag(vcov(lm_mar.Y))[2]  
}

ACME_mcar_mice_lm <- mean(ACME.mcar_mice_lm)
ADE_mcar_mice_lm <- mean(ADE.mcar_mice_lm)
ACME_mcar_mice_SE_lm <- sqrt(mean(sobel_mcar_mice_lm) + (1+(1/10))*var(ACME.mcar_mice_lm))

ACME_cov_mcar_mice_lm <- ACME_mcar_mice_lm-ACME_mcar_mice_SE_lm*1.96 <= True.ACME & 
                         ACME_mcar_mice_lm+ACME_mcar_mice_SE_lm*1.96 >= True.ACME 

ADE_SE_mcar_mice_lm <- sqrt(mean(ADE.var_mcar_mice_lm) + (1+(1/10))*var(ADE.mcar_mice_lm))

ADE_cov_mcar_mice_lm <- ADE_mcar_mice_lm-ADE_SE_mcar_mice_lm*1.96 <= True.ADE & 
                        ADE_mcar_mice_lm+ADE_SE_mcar_mice_lm*1.96 >= True.ADE 


ACME_mar_mice_lm <- mean(ACME.mar_mice_lm)
ADE_mar_mice_lm <- mean(ADE.mar_mice_lm)
ACME_mar_mice_SE_lm <- sqrt(mean(sobel_mar_mice_lm) + (1+(1/10))*var(ACME.mar_mice_lm))

ACME_cov_mar_mice_lm <- ACME_mar_mice_lm-ACME_mar_mice_SE_lm*1.96 <= True.ACME & 
                        ACME_mar_mice_lm+ACME_mar_mice_SE_lm*1.96 >= True.ACME 

ADE_SE_mar_mice_lm <- sqrt(mean(ADE.var_mar_mice_lm) + (1+(1/10))*var(ADE.mar_mice_lm))

ADE_cov_mar_mice_lm <- ADE_mar_mice_lm-ADE_SE_mar_mice_lm*1.96 <= True.ADE & 
                       ADE_mar_mice_lm+ADE_SE_mar_mice_lm*1.96 >= True.ADE 


################################################################

### MULTIPLE IMPUTATION IGNORING CLUSTERING
## 'MICE' is used in this imputation
## Data is arranged in long format

mcar.data_mice0 <- data.complete[,c("group","time","X1","X2","M.mcar","Y.mcar")]
mar.data_mice0 <- data.complete[,c("group","time","X1","X2","M.mar","Y.mar")]

mcar.data_mice0$group_time <- mcar.data_mice0$group*mcar.data_mice0$time
mar.data_mice0$group_time <- mar.data_mice0$group*mar.data_mice0$time

factor_time <- model.matrix( ~ as.factor(time), data=mcar.data_mice0)[,2:6]
factor_group.time <- model.matrix( ~ as.factor(group_time), data=mcar.data_mice0)[,2:6]

time_factors <- as.data.frame(cbind(factor_time,factor_group.time))
names(time_factors) <- c("time.2","time.3","time.4","time.5","time.6",
              "grp_time.2","grp_time.3","grp_time.4","grp_time.5","grp_time.6")

mcar.data_mice1 <- as.data.frame(cbind(mcar.data_mice0, time_factors))
mar.data_mice1 <- as.data.frame(cbind(mar.data_mice0, time_factors))

drops <- c("time","group_time")
mcar.data_mice <- mcar.data_mice1[ , !(names(mcar.data_mice1) %in% drops)] # Drops time and its interaction with group
mar.data_mice <- mar.data_mice1[ , !(names(mar.data_mice1) %in% drops)] # Drops time and its interaction with group

### Imputing missing for both MCAR and MAR

mcar.mice_imputed <- mice(mcar.data_mice, m=numb_imput, maxit = 15,printFlag=FALSE, method = 'pmm', seed = 339+i)
mar.mice_imputed <- mice(mar.data_mice, m=numb_imput, maxit = 15, printFlag=FALSE, method = 'pmm', seed = 753+i)

### MCAR

long.mcar_mice <- mice::complete(mcar.mice_imputed, action='long', include=FALSE)

long.mcar_mice <- cbind(long.mcar_mice, id2=data.complete$id, time=data.complete$time)

# Create an identifier id incoporating the data set imputation id

for (b in 1:nrow(long.mcar_mice))
{
  long.mcar_mice$id[b] <- paste(c(as.character(long.mcar_mice$id2[b]),
                                  as.character(long.mcar_mice$.imp[b])),
                             collapse=".")
}

long.mcar_mice2 <- long.mcar_mice[,c(".imp","group","time","X1","X2","M.mcar","Y.mcar","id")] 
wide.mcar_mice <- reshape(long.mcar_mice2, idvar = "id", timevar = "time", 
                              v.names=c("M.mcar","Y.mcar"), direction = "wide")
ACME.mcar_mice <- rep(NA, 10)
ADE.mcar_mice <- rep(NA, 10)
sobel_mcar_mice <- rep(NA, 10)
ADE.var_mcar_mice <- rep(NA, 10)

for (c in 1:10)
{
  data.mcar <- wide.mcar_mice[which(wide.mcar_mice$.imp==c), ]
  data.mcar$M6_M1 <- data.mcar$M.mcar.6 - data.mcar$M.mcar.1
  
  lm_mcar.M <- lm(M6_M1 ~ group + X1 + X2, data = data.mcar)
  lm_mcar.Y <- lm(Y.mcar.6 ~ group + M6_M1 + X1 + X2, data = data.mcar)

  ACME.mcar_mice[c] <- coefficients(lm_mcar.M)[2] * coefficients(lm_mcar.Y)[3]
  ADE.mcar_mice[c] <- coefficients(lm_mcar.Y)[2]  
  
  sobel_mcar_mice[c] <- sobel.test(a=coefficients(lm_mcar.M)[2], b=coefficients(lm_mcar.Y)[3], 
                                  var.a=diag(vcov(lm_mcar.M))[2], var.b=diag(vcov(lm_mcar.Y))[3])
  ADE.var_mcar_mice[c] <- diag(vcov(lm_mcar.Y))[2]
}

ACME_mcar_mice <- mean(ACME.mcar_mice)
ADE_mcar_mice <- mean(ADE.mcar_mice)
ACME_mcar_mice_SE <- sqrt(mean(sobel_mcar_mice) + (1+(1/10))*var(ACME.mcar_mice))

ACME_cov_mcar_mice <- ACME_mcar_mice-ACME_mcar_mice_SE*1.96 <= True.ACME & 
                      ACME_mcar_mice+ACME_mcar_mice_SE*1.96 >= True.ACME 

ADE_SE_mcar_mice <- sqrt(mean(ADE.var_mcar_mice) + (1+(1/10))*var(ADE.mcar_mice))

ADE_cov_mcar_mice <- ADE_mcar_mice-ADE_SE_mcar_mice*1.96 <= True.ADE & 
                     ADE_mcar_mice+ADE_SE_mcar_mice*1.96 >= True.ADE 


### MAR

long.mar_mice <- mice::complete(mar.mice_imputed, action='long', include=FALSE)
long.mar_mice <- cbind(long.mar_mice, id2=data.complete$id, time=data.complete$time)

# Create an identifier id incoporating the data set imputation id

for (h in 1:nrow(long.mar_mice))
{
  long.mar_mice$id[h] <- paste(c(as.character(long.mar_mice$id2[h]), 
                                  as.character(long.mar_mice$.imp[h])), 
                                collapse="")
}

long.mar_mice2 <- long.mar_mice[,c(".imp","group","time","X1","X2","M.mar","Y.mar","id")] 
wide.mar_mice <- reshape(long.mar_mice2, idvar = "id", timevar = "time", 
                          v.names=c("M.mar","Y.mar"), direction = "wide")
ACME.mar_mice <- rep(NA, 10)
ADE.mar_mice <- rep(NA, 10)
sobel_mar_mice <- rep(NA, 10)
ADE.var_mar_mice <- rep(NA, 10)

for (d in 1:10)
{
  data.mar <- wide.mar_mice[which(wide.mar_mice$.imp==d), ]
  data.mar$M6_M1 <- data.mar$M.mar.6 - data.mar$M.mar.1
  
  lm_mar.M <- lm(M6_M1 ~ group + X1 + X2, data = data.mar)
  lm_mar.Y <- lm(Y.mar.6 ~ group + M6_M1 + X1 + X2, data = data.mar)
  
  ACME.mar_mice[d] <- coefficients(lm_mar.M)[2] * coefficients(lm_mcar.Y)[3]
  ADE.mar_mice[d] <- coefficients(lm_mar.Y)[2]  
  
  sobel_mar_mice[d] <- sobel.test(a=coefficients(lm_mar.M)[2], b=coefficients(lm_mar.Y)[3], 
                                   var.a=diag(vcov(lm_mar.M))[2], var.b=diag(vcov(lm_mar.Y))[3])
  ADE.var_mar_mice[d] <- diag(vcov(lm_mar.Y))[2]
}

ACME_mar_mice <- mean(ACME.mar_mice)
ADE_mar_mice <- mean(ADE.mar_mice)
ACME_mar_mice_SE <- sqrt(mean(sobel_mar_mice) + (1+(1/10))*var(ACME.mar_mice))

ACME_cov_mar_mice <- ACME_mar_mice-ACME_mar_mice_SE*1.96 <= True.ACME & 
                     ACME_mar_mice+ACME_mar_mice_SE*1.96 >= True.ACME 

ADE_SE_mar_mice <- sqrt(mean(ADE.var_mar_mice) + (1+(1/10))*var(ADE.mar_mice))

ADE_cov_mar_mice <- ADE_mar_mice-ADE_SE_mar_mice*1.96 <= True.ADE & 
                     ADE_mar_mice+ADE_SE_mar_mice*1.96 >= True.ADE 

##########

propmiss <- prop.miss[j]


 unlist(list(propmiss, ACME.nomiss, ADE.nomiss, ACME_mcar_pan, ADE_mcar_pan,
             ACME_mar_pan, ADE_mar_pan, ACME_mcar_mice, ADE_mcar_mice, ACME_mar_mice, ADE_mar_mice,  
             ACME_mcar_pan_SE, ACME_mar_pan_SE, ACME_cov_mcar_pan, ACME_cov_mar_pan, 
             ACME_mcar_mice_SE, ACME_mar_mice_SE, ACME_cov_mcar_mice, 
             ACME_cov_mar_mice, ACME_nomiss_SE, ADE_SE_nomiss, ACME_cov_nomiss, ADE_cov_nomiss,
             ADE_SE_mcar_pan, ADE_cov_mcar_pan, ADE_SE_mar_pan, ADE_cov_mar_pan, ADE_SE_mcar_mice, 
             ADE_cov_mcar_mice, ADE_SE_mar_mice, ADE_cov_mar_mice,
             ACME_mcar_mice_lm, ACME_mcar_mice_SE_lm, ADE_mcar_mice_lm, ADE_SE_mcar_mice_lm,
             ACME_cov_mcar_mice_lm, ADE_cov_mcar_mice_lm,
             ACME_mar_mice_lm, ACME_mar_mice_SE_lm, ADE_mar_mice_lm, ADE_SE_mar_mice_lm,
             ACME_cov_mar_mice_lm, ADE_cov_mar_mice_lm,
             
             ACME.noimp.mcar, ADE.noimp.mcar, ADE_SE_noimp.mcar,  ACME_noimp.mcar_SE, 
             ACME_cov_noimp.mcar, ADE_cov_noimp.mcar, 
             ACME.noimp.mar, ADE.noimp.mar,  ADE_SE_noimp.mar, ACME_noimp.mar_SE, 
             ACME_cov_noimp.mar, ADE_cov_noimp.mar), 
        recursive = TRUE, use.names = FALSE)

}

df <- as.data.frame(df)

print(Sys.time()-strt)
stopCluster(cl)

######################################################################################################

est <- c(df$result.1, df$result.2, df$result.3, df$result.4)
Estimates <- matrix(est, ncol=55, byrow=TRUE)

colnames(Estimates) <- c("propmiss", "ACME.nomiss", "ADE.nomiss", "ACME_mcar_pan", "ADE_mcar_pan",
                         "ACME_mar_pan", "ADE_mar_pan", "ACME_mcar_mice", "ADE_mcar_mice", "ACME_mar_mice", "ADE_mar_mice",  
                         "ACME_mcar_pan_SE", "ACME_mar_pan_SE", "ACME_cov_mcar_pan", "ACME_cov_mar_pan", 
                         "ACME_mcar_mice_SE", "ACME_mar_mice_SE", "ACME_cov_mcar_mice", 
                         "ACME_cov_mar_mice", "ACME_nomiss_SE", "ADE_SE_nomiss", "ACME_cov_nomiss", "ADE_cov_nomiss",
                         "ADE_SE_mcar_pan", "ADE_cov_mcar_pan", "ADE_SE_mar_pan", "ADE_cov_mar_pan", "ADE_SE_mcar_mice", 
                         "ADE_cov_mcar_mice", "ADE_SE_mar_mice", "ADE_cov_mar_mice",
                         "ACME_mcar_mice_lm", "ACME_mcar_mice_SE_lm", "ADE_mcar_mice_lm", "ADE_SE_mcar_mice_lm",
                         "ACME_cov_mcar_mice_lm", "ADE_cov_mcar_mice_lm",
                         "ACME_mar_mice_lm", "ACME_mar_mice_SE_lm", "ADE_mar_mice_lm", "ADE_SE_mar_mice_lm",
                         "ACME_cov_mar_mice_lm", "ADE_cov_mar_mice_lm",
                         
                         "ACME.noimp.mcar", "ADE.noimp.mcar", "ADE_SE_noimp.mcar",  "ACME_noimp.mcar_SE", 
                         "ACME_cov_noimp.mcar", "ADE_cov_noimp.mcar", 
                         "ACME.noimp.mar", "ADE.noimp.mar",  "ADE_SE_noimp.mar", "ACME_noimp.mar_SE", 
                         "ACME_cov_noimp.mar", "ADE_cov_noimp.mar")

est_miss_20180213 <- as.data.frame(Estimates)

save(est_miss_20180213, file = "C:/Users/John/Google Drive/PhD/Thesis/Paper 2/BMC Medical Research Methodology/R code/simthendiff_20180213.RData")

######################################################################################################





