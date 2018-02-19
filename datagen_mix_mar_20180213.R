
################################################################################
# By John Ssenkusu (ssenk001@umn.edu, jssenkusu@gmail.com)
# Purpose: Simulation study on missing data and mediation analysis in 
#          longitudinal data. Data generation code.
# Created: December 27, 2017
# Modified: February 13, 2018
################################################################################

# Generates data whose MAR mechanisms depends on the outcome and the other
# half it depends on the mediator

library(MASS)
library(reshape)
library(simglm)

### DATA GENERATION

simdata <- function(prop.miss){
  
  samplesize = 300
  numb_obs <- 6   # Number of observations per subject
  sigma_err <- matrix(c(1,0,0,1),2,2) # Error variance for the outcome and mediator
  mu_err <- c(0, 0) # Error mean the outcome and mediator
  
  sigma_bi <- matrix(c(4,0,0,4),2,2) # Random intercept variance for the outcome and mediator
  mu_bi <- c(0, 0)  # Random effect mean 
  
  errors <- mvrnorm(n = samplesize*numb_obs, mu=mu_err, Sigma=sigma_err, tol = 1e-6, empirical = FALSE)
  id_err <- sort(rep(1:samplesize, numb_obs))
  time <- rep(1:numb_obs, samplesize)
  group <- sample(c(1,0), size = samplesize, replace=TRUE)
  
  bi <- mvrnorm(n = samplesize, mu=mu_bi, Sigma=sigma_bi, tol = 1e-6, empirical = FALSE)
  
  X1 <- rnorm(samplesize, mean=10, sd = 1.2)
  X2 <- rbinom(samplesize, size=1 , prob = 0.5)
  
  dt1 <- cbind(rep(1:samplesize), group, X1, X2, bi)
  colnames(dt1) <- c("id","group","X1","X2" ,"bi_y","bi_m")
  dt1 <- as.data.frame(dt1)
  
  dt2 <- cbind(id_err, time, errors)
  colnames(dt2) <- c("id","time","e_y","e_m")
  dt2 <- as.data.frame(dt2)
  
  data <- merge(dt2, dt1, by="id")
  
  X_m <- model.matrix(~ group + as.factor(time) + X1 + X2 + group*as.factor(time), data = data)
  
  # Generate mediator 
  
  M_interac_terms <- c(-0.2,-0.6,-0.5,-0.7,-0.9)
  
  m_beta <- matrix(c(8, -1.7, 3.3, 3.8, 4.3, 4.8, 5.3, -0.30, 1.9, M_interac_terms), ncol=1)
  
  data$M <- X_m %*% m_beta + data$bi_m + data$e_m
  
  data_wide <- reshape(data, idvar = "id", timevar = "time", 
                       v.names=c("e_y","e_m","M"), direction = "wide")
  data_wide$M6_M1 <- data_wide$M.6 - data_wide$M.1
  data_wide$M5_M1 <- data_wide$M.5 - data_wide$M.1
  data_wide$M4_M1 <- data_wide$M.4 - data_wide$M.1
  data_wide$M3_M1 <- data_wide$M.3 - data_wide$M.1
  data_wide$M2_M1 <- data_wide$M.2 - data_wide$M.1
  data_wide$M1_M1 <- data_wide$M.1 - data_wide$M.1
  
  # Generate outcome
  
  int <- c(0,1.2,1.25,1.45,1.6,1.75)-1  # Intercept
  d <- c(2.4,2.4,2.4,2.4,2.4,2.4)       # Treatment/exposure effect
  
  data_wide$Y.1 <- int[1]-1.5*data_wide$group+d[1]*data_wide$M1_M1+0.22*data_wide$X1-1.55*data_wide$X2+data_wide$bi_y+data_wide$e_y.1
  data_wide$Y.2 <- int[2]-1.5*data_wide$group+d[2]*data_wide$M2_M1+0.22*data_wide$X1-1.55*data_wide$X2+data_wide$bi_y+data_wide$e_y.2
  data_wide$Y.3 <- int[3]-1.5*data_wide$group+d[3]*data_wide$M3_M1+0.22*data_wide$X1-1.55*data_wide$X2+data_wide$bi_y+data_wide$e_y.3
  data_wide$Y.4 <- int[4]-1.5*data_wide$group+d[4]*data_wide$M4_M1+0.22*data_wide$X1-1.55*data_wide$X2+data_wide$bi_y+data_wide$e_y.4
  data_wide$Y.5 <- int[5]-1.5*data_wide$group+d[5]*data_wide$M5_M1+0.22*data_wide$X1-1.55*data_wide$X2+data_wide$bi_y+data_wide$e_y.5
  data_wide$Y.6 <- int[6]-1.5*data_wide$group+d[6]*data_wide$M6_M1+0.22*data_wide$X1-1.55*data_wide$X2+data_wide$bi_y+data_wide$e_y.6
  
  data_long <- reshape(data_wide, idvar = "id", varying = c("e_y.1","e_m.1","M.1" ,"e_y.2","e_m.2","M.2","e_y.3","e_m.3","M.3","e_y.4",
                                                            "e_m.4","M.4","e_y.5","e_m.5","M.5","e_y.6","e_m.6","M.6","Y.1","Y.2","Y.3",
                                                            "Y.4","Y.5","Y.6"),
                       timevar="time", times=c(1,2,3,4,5,6), direction = "long")
  
  data_long <- data_long[order(data_long$id,data_long$time),]
  
  data.complete <- data_long[,c("id","group","time","M","Y","X1","X2")]
  
  # ### Ckecking on the parameter estimates
  # 
  # library(lme4)
  # summary(lmer(M~group+as.factor(time)+X1+X2+group*as.factor(time)+(1|id), data=data.complete))
  # summary(lmer(Y~group+M+as.factor(time)+X1+X2+(1|id), data=data.complete))
  # 
  # library(nlme)
  # model.M <- lme(M~group+as.factor(time)+X1+X2+group*as.factor(time), random=~1|id, data=data.complete)
  # model.Y <- lme(Y~group+M+as.factor(time)+X1+X2, random=~1|id, data=data.complete)
  # 
  # summary(model.M)
  # summary(model.Y)
  ######################################################################################
  
  ### Creating data that is MCAR and MAR from a complete data set
  
  ### MCAR FOR THE OUTCOME AND MEDIATOR
  ### The baseline time observations is always observed (never missing).
  ### When the mediator is missing, the outcome is also missing 
  ### The percentage of missingness is NOT based on subjects, rather on measurements for 
  ### the outcome or mediator. Since the process is random, the proportion of the difference M6-M1 MCAR is approximately as expected.
  ### This is different from how MAR is simulated, it is based on number of subjects missing atleast one observation.
  ### When the outcome is missing, the mediator is set to missing also and viceversa.
  
  sim_data <- data.complete
  sim_data$M <- as.numeric(sim_data$M)
  sim_data$Y <- as.numeric(sim_data$Y)
  
  names(sim_data)[names(sim_data)=="M"] <- "sim_data"
  simdata1 <- sim_data[ which(sim_data$time > 1),] # Data for time 2 to 6
  # Generating MCAR from time 2 to time 6
  sim_mcar <- random_missing(sim_data=simdata1, resp_var = "sim_data", 
                             miss_prop=prop.miss, clust_var = NULL,
                             within_id = NULL)
  
  sim_mcar$Y.mcar <- sim_mcar$Y
  sim_mcar$Y.mcar[is.na(sim_mcar$sim_data2) ] <- NA
  names(sim_mcar)[names(sim_mcar)=="sim_data2"] <- "M.mcar"
  names(sim_mcar)[names(sim_mcar)=="sim_data"] <- "M"
  sim_mcar <- sim_mcar[,c("id","group","time","M","Y","X1","X2","Y.mcar","M.mcar")]
  
  simdata2 <- sim_data[ which(sim_data$time == 1),] # Data for time 1
  simdata2$M.mcar <- simdata2$sim_data  
  simdata2$Y.mcar <- simdata2$Y 
  names(simdata2)[names(simdata2)=="sim_data"] <- "M"
  
  data.complete <- rbind(simdata2, sim_mcar)
  data.complete <- data.complete[order(data.complete$id, data.complete$time),]

  ### MAR FOR THE OUTCOME AND MEDIATOR 
  ## The MAR mechanism is mediator depended half of the times and the other half,
  ## it is dependent on the outcome.
  
  ## Spliting complete data set into two equal data sets
  
  data.complete_y <- data.complete[ which(data.complete$id <= (samplesize/2)), ]
  data.complete_m <- data.complete[ which(data.complete$id > (samplesize/2)), ]
  
### MAR based on the outcome
  
  y.mar_y = matrix(data.complete_y$Y, ncol=numb_obs, nrow=(samplesize/2), byrow=TRUE)
  m.mar_y = matrix(data.complete_y$M, ncol=numb_obs, nrow=(samplesize/2), byrow=TRUE)
  
  y.mar.1 <- y.mar_y[,1]
  
  dif_y <- matrix(0, nrow=(samplesize/2), ncol=numb_obs)
  for(i in 1:(samplesize/2)){
    for(j in 2:numb_obs){
      dif_y[i,j] = y.mar_y[i,j]-y.mar_y[i,j-1]
    }
  }
  
  min.row_y <- apply(dif_y[1:(samplesize/2), 2:(numb_obs-1)], MARGIN=1, FUN = min) # select minimum for each subject row excluding first visit
  
  ymar1.min <- cbind(y.mar.1, min.row_y)
  ymar1.min_order <- ymar1.min[order(ymar1.min[,1],ymar1.min[,2],decreasing=FALSE),]
  cutoffs.miss_y <- ymar1.min_order[ceiling(prop.miss*(samplesize/2)), ]
  
  max.diff_y <- max(ymar1.min_order[1:ceiling(prop.miss*(samplesize/2)),2])
  
  for(p in 1:(samplesize/2)){
    for(q in 2:numb_obs){
      dif2 = y.mar_y[p,q]-y.mar_y[p,q-1]
      dif2.round <- round(dif2, digits = 5)
      if(y.mar.1[p] <= cutoffs.miss_y[1])
      {
        # if(dif2.round <= cutoffs.miss[2])
        if(dif2.round <= max.diff_y)
        {  
          k = q + 1 
          if(q != numb_obs) { y.mar_y[p,k:numb_obs] = NA;  m.mar_y[p,k:numb_obs] = NA; break }  # if score goes down by more than cutoff.miss, drop subsequent obs. The information leading to the drop is observed and remains in data set.
        } }
    }
  }
  
  data.complete_y$Y.mar = as.vector(t(y.mar_y))
  data.complete_y$M.mar = as.vector(t(m.mar_y))
  
### MAR based on the mediator
  
  y.mar_m = matrix(data.complete_m$Y, ncol=numb_obs, nrow=(samplesize/2), byrow=TRUE)
  m.mar_m = matrix(data.complete_m$M, ncol=numb_obs, nrow=(samplesize/2), byrow=TRUE)
  
  m.mar.1 <- m.mar_m[,1]
  
  dif_m <- matrix(0, nrow=(samplesize/2), ncol=numb_obs)
  for(r in 1:(samplesize/2)){
    for(s in 2:numb_obs){
      dif_m[r,s] = m.mar_m[r,s]-m.mar_m[r,s-1]
    }
  }
  
  min.row_m <- apply(dif_m[1:(samplesize/2), 2:(numb_obs-1)], MARGIN=1, FUN = min) # select minimum for each subject row excluding first visit
  
  mmar1.min <- cbind(m.mar.1, min.row_m)
  mmar1.min_order <- mmar1.min[order(mmar1.min[,1],mmar1.min[,2],decreasing=FALSE),]
  cutoffs.miss_m <- mmar1.min_order[ceiling(prop.miss*(samplesize/2)), ]
  
  max.diff_m <- max(mmar1.min_order[1:ceiling(prop.miss*(samplesize/2)),2])
  
  for(c in 1:(samplesize/2)){
    for(d in 2:numb_obs){
      dif3 = m.mar_m[c,d]-m.mar_m[c,d-1]
      dif3.round <- round(dif3, digits = 5)
      if(m.mar.1[c] <= cutoffs.miss_m[1])
      {
        if(dif3.round <= max.diff_m)
        {  
          f = d + 1 
          if(d != numb_obs) { m.mar_m[c,f:numb_obs] = NA;  y.mar_m[c,f:numb_obs] = NA; break }  # if score goes down by more than cutoff.miss, drop subsequent obs. The information leading to the drop is observed and remains in data set.
        } }
    }
  }
  
  data.complete_m$Y.mar = as.vector(t(y.mar_m))
  data.complete_m$M.mar = as.vector(t(m.mar_m))
  
  ## Recombining the data
  
  data.complete <- rbind(data.complete_y, data.complete_m)
 
  data.complete <- data.complete[order(data.complete$id, data.complete$time),]
  
  return(list(data.complete))
  
}

#################################################################################################

 
  
  
  
  
  