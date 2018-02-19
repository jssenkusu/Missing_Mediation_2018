
################################################################################
# By John Ssenkusu (ssenk001@umn.edu, jssenkusu@gmail.com)
# Purpose: Generate latex tables for the missing data Paper 
# Created: February 25, 2017
# Modified: February 13, 2018
################################################################################

library(xtable)

est <- get(load("C:/Users/name/R code/simthendiff_20180213.rdata"))
# est <- get(load("C:/Users/name/R code/diffthensim_20180213.rdata"))

### Computing bias

true_ACME <- 2.4*-0.9
true_ADE <- -1.5

est$bias_ACME_nomiss <- est$ACME.nomiss - true_ACME
est$bias_ACME_mcar_pan <- est$ACME_mcar_pan - true_ACME
est$bias_ACME_mar_pan <- est$ACME_mar_pan - true_ACME
est$bias_ACME_mcar_mice <- est$ACME_mcar_mice - true_ACME
est$bias_ACME_mar_mice <- est$ACME_mar_mice - true_ACME
est$bias_ACME_mcar_mice_lm <- est$ACME_mcar_mice_lm - true_ACME
est$bias_ACME_mar_mice_lm <- est$ACME_mar_mice_lm - true_ACME

est$bias_ADE_nomiss <- est$ADE.nomiss - true_ADE
est$bias_ADE_mcar_pan <- est$ADE_mcar_pan - true_ADE
est$bias_ADE_mar_pan <- est$ADE_mar_pan - true_ADE
est$bias_ADE_mcar_mice <- est$ADE_mcar_mice - true_ADE
est$bias_ADE_mar_mice <- est$ADE_mar_mice - true_ADE
est$bias_ADE_mcar_mice_lm <- est$ADE_mcar_mice_lm - true_ADE
est$bias_ADE_mar_mice_lm <- est$ADE_mar_mice_lm - true_ADE

### computing square error

est$SE_ACME_nomiss <- est$bias_ACME_nomiss^2
est$SE_ACME_mcar_pan <- est$bias_ACME_mcar_pan^2
est$SE_ACME_mar_pan <- est$bias_ACME_mar_pan^2
est$SE_ACME_mcar_mice <- est$bias_ACME_mcar_mice^2
est$SE_ACME_mar_mice <- est$bias_ACME_mar_mice^2
est$SE_ACME_mcar_mice_lm <- est$bias_ACME_mcar_mice_lm^2
est$SE_ACME_mar_mice_lm <- est$bias_ACME_mar_mice_lm^2

est$SE_ADE_nomiss <- est$bias_ADE_nomiss^2
est$SE_ADE_mcar_pan <- est$bias_ADE_mcar_pan^2
est$SE_ADE_mar_pan <- est$bias_ADE_mar_pan^2
est$SE_ADE_mcar_mice <- est$bias_ADE_mcar_mice^2
est$SE_ADE_mar_mice <- est$bias_ADE_mar_mice^2
est$SE_ADE_mcar_mice_lm <- est$bias_ADE_mcar_mice_lm^2
est$SE_ADE_mar_mice_lm <- est$bias_ADE_mar_mice_lm^2

### Bias and square error for 'no imputation'

est$bias_ACME_noimp_mcar <- est$ACME.noimp.mcar - true_ACME
est$bias_ACME_noimp_mar <- est$ACME.noimp.mar - true_ACME

est$bias_ADE_noimp_mcar <- est$ADE.noimp.mcar - true_ADE
est$bias_ADE_noimp_mar <- est$ADE.noimp.mar - true_ADE

est$SE_ACME_noimp_mcar <- est$bias_ACME_noimp_mcar^2
est$SE_ACME_noimp_mar <- est$bias_ACME_noimp_mar^2

est$SE_ADE_noimp_mcar <- est$bias_ADE_noimp_mcar^2
est$SE_ADE_noimp_mar <- est$bias_ADE_noimp_mar^2

#################################################################

### SUMMARIZING SIMULATION RESULTS INTO TABLES

prop_miss <- unique(est$propmiss)

ACME_estimates <- matrix(0, nrow = 1, ncol = 12)
ADE_estimates <- matrix(0, nrow = 1, ncol = 12)


for(k in 1:length(prop_miss))
{
  dt <- est[est$propmiss==prop_miss[k],]
  
  estM <- colMeans(dt)
  SD <- apply(dt[,c(2:11,32,34,38,40,44:45,50:51)], 2, sd)
  
  ACME_est <- matrix(c(estM[44], estM[84], estM[88], SD[15], estM[47], estM[48],
                       estM[50], estM[85], estM[89], SD[17], estM[53], estM[54],
                       estM[8], estM[59], estM[73], SD[7], estM[16], estM[18],
                       estM[10], estM[60], estM[74], SD[9], estM[17], estM[19],
                       estM[32], estM[61], estM[75], SD[11], estM[33], estM[36],
                       estM[38], estM[62], estM[76], SD[13], estM[39], estM[42],
                       estM[4], estM[57], estM[71], SD[3], estM[12], estM[14],
                       estM[6], estM[58], estM[72], SD[5], estM[13], estM[15]),
                     ncol = 12, nrow = 4, byrow = TRUE)

  ADE_est <- matrix(c(estM[45], estM[86], estM[90], SD[16], estM[46], estM[49],
                      estM[51], estM[87], estM[91], SD[18], estM[52], estM[55],
                      estM[9], estM[66], estM[80], SD[8], estM[28], estM[29],
                      estM[11], estM[67], estM[81], SD[10], estM[30], estM[31],
                      estM[34], estM[68], estM[82], SD[12], estM[35], estM[37],
                      estM[40], estM[69], estM[83], SD[14], estM[41], estM[43],
                      estM[5], estM[64], estM[78], SD[4], estM[24], estM[25],
                      estM[7], estM[65], estM[79], SD[6], estM[26], estM[27]),
                    ncol = 12, nrow = 4, byrow = TRUE)
  

  ACME_estimates <- rbind(ACME_estimates, ACME_est)
  ADE_estimates <- rbind(ADE_estimates, ADE_est)

}


ACME_results <- ACME_estimates[2:17,1:12]
ADE_results <- ADE_estimates[2:17,1:12]

misprop <- sort(rep(prop_miss, times=4))
misprop <- c(0, misprop)

# Adding estimates when data is complete

nomiss_ACME_est <- cbind(estM[2], estM[56], estM[70], SD[1], estM[20], estM[22],
                         estM[2], estM[56], estM[70], SD[1], estM[20], estM[22])
nomiss_ADE_est <- cbind(estM[3], estM[63], estM[77], SD[2], estM[21], estM[23],
                        estM[3], estM[63], estM[77], SD[2], estM[21], estM[23])

ACME_Table <- rbind(nomiss_ACME_est, ACME_results)
ADE_Table <- rbind(nomiss_ADE_est, ADE_results)

Table_ACME <- cbind(misprop, ACME_Table )
Table_ADE <- cbind(misprop, ADE_Table)

###############################################################################

colnames(Table_ACME) <- c("Proportion","Estimate","Bias","MSE","MC-SD","MS-SE","CR",
                          "Estimate","Bias","MSE","MC-SD","MS-SE","CR")

colnames(Table_ADE) <- c("Proportion","Estimate","Bias","MSE","MC-SD","Avg-SE","CR",
                              "Estimate","Bias","MSE","MC-SD","Avg-SE","CR")


rownames(Table_ACME) <- c("Complete","No imputation","IC (long)","LM","CC",
                          "No imputation2","IC (long)2","LM2","CC2",
                          "No imputation3","IC (long)3","LM3","CC3",
                          "No imputation4","IC (long)4","LM4","CC4")

rownames(Table_ADE) <- c("Complete","No imputation","IC (long)","LM","CC",
                         "No imputation2","IC (long)2","LM2","CC2",
                         "No imputation3","IC (long)3","LM3","CC3",
                         "No imputation4","IC (long)4","LM4","CC4")


print(xtable(Table_ACME, digits = c(0,2,2,2,2,2,2,2,2,2,2,2,2,2), align = "lccccccccccccc",
             caption = "ACME estimates when no imputation is conducted and after multiple imputation using different approaches in longitudinal data at differing proportions of 
                        of missing data. True ACME = -2.16, Sample size = 300 and number of Monte-Carlo 
             data sets = 2000"), caption.placement = "top")

print(xtable(Table_ADE, digits = c(0,2,2,2,2,2,2,2,2,2,2,2,2,2), align = "lccccccccccccc",
             caption = "ADE estimates when no imputation is conducted and after multiple imputation using different approaches in longitudinal data at differing proportions of 
                        of missing data. True ADE = -1.50, Sample size = 300 and number of Monte-Carlo 
             data sets = 2000"), caption.placement = "top")


###############################################################################


