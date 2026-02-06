
#####################################################################################################################
####                                        RUNNING LOGISTIC REGRESSION MODEL                                    ####
#####################################################################################################################

library(readxl)
data <- read_excel("mtb_wh_extended.v1.with_DR_trt_data.metadata.gyrAB_mutations.xlsx", sheet = 1)
dim(data)
# [1] 550  17

# checking right variable format of variables to be used in model
data$gyrAB_mutated = as.factor(data$gyrAB_mutated)
data$Fq = as.factor(data$Fq)
data$time_gap = as.numeric(data$time_gap)
data$main_lineage = as.factor(data$main_lineage)
data$tb_profiler_baseline_drug_resistance = as.factor(data$tb_profiler_baseline_drug_resistance)
data$tb_profiler_baseline_drug_resistance2 = as.factor(data$tb_profiler_baseline_drug_resistance2)

model_vars = c("gyrAB_mutated","Fq","time_gap","main_lineage","tb_profiler_baseline_drug_resistance")
data_clean <- data[complete.cases(data[, model_vars]), ]
dim(data_clean)
# [1]  405 17
data_clean = data_clean[-which(is.infinite(data_clean$time_gap)),]
dim(data_clean)
# [1]  400 17

model1 <- glm(gyrAB_mutated ~ Fq,
             data = data_clean,
             family = binomial)
summary(model1)
# Call:
#   glm(formula = gyrAB_mutated ~ Fq, family = binomial, data = data_clean)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -4.1972     0.7121  -5.894 3.77e-09 ***
#   Fq1           2.4098     0.7334   3.286  0.00102 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 260.07  on 399  degrees of freedom
# Residual deviance: 238.69  on 398  degrees of freedom
# AIC: 242.69
# 
# Number of Fisher Scoring iterations: 6

# Calculating CI for ORs
coef_est <- coef(model1)
se <- sqrt(diag(vcov(model1)))
lower <- coef_est - 1.96 * se
upper <- coef_est + 1.96 * se
OR <- exp(coef_est)
lower_OR <- exp(lower)
upper_OR <- exp(upper)

model2 <- glm(gyrAB_mutated ~ tb_profiler_baseline_drug_resistance2,
              data = data_clean,
              family = binomial)

summary(model2)
# Call:
#   glm(formula = gyrAB_mutated ~ tb_profiler_baseline_drug_resistance2, 
#       family = binomial, data = data_clean)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                     -1.9384     0.1735 -11.172  < 2e-16 ***
#   tb_profiler_baseline_drug_resistance2Sensitive  -1.9328     0.7352  -2.629  0.00856 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 260.07  on 399  degrees of freedom
# Residual deviance: 248.07  on 398  degrees of freedom
# AIC: 252.07
# 
# Number of Fisher Scoring iterations: 6


coef_est <- coef(model2)
se <- sqrt(diag(vcov(model2)))
lower <- coef_est - 1.96 * se
upper <- coef_est + 1.96 * se
OR <- exp(coef_est)
lower_OR <- exp(lower)
upper_OR <- exp(upper)


model3 <- glm(gyrAB_mutated ~ time_gap,
              data = data_clean,
              family = binomial)
summary(model3)
# Call:
#   glm(formula = gyrAB_mutated ~ time_gap, family = binomial, data = data_clean)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -2.193e+00  1.817e-01 -12.070   <2e-16 ***
#   time_gap    -5.265e-05  8.470e-04  -0.062     0.95    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 260.07  on 399  degrees of freedom
# Residual deviance: 260.06  on 398  degrees of freedom
# AIC: 264.06
# 
# Number of Fisher Scoring iterations: 4

coef_est <- coef(model3)
se <- sqrt(diag(vcov(model3)))
lower <- coef_est - 1.96 * se
upper <- coef_est + 1.96 * se
OR <- exp(coef_est)
lower_OR <- exp(lower)
upper_OR <- exp(upper)


model4 <- glm(gyrAB_mutated ~ main_lineage,
              data = data_clean,
              family = binomial)
summary(model4)
# Call:
#   glm(formula = gyrAB_mutated ~ main_lineage, family = binomial, 
#       data = data_clean)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          -1.657e+01  1.697e+03  -0.010    0.992
# main_lineagelineage2  1.436e+01  1.697e+03   0.008    0.993
# main_lineagelineage3 -4.731e-09  2.008e+03   0.000    1.000
# main_lineagelineage4  1.443e+01  1.697e+03   0.009    0.993
# main_lineagelineage5 -4.739e-09  2.400e+03   0.000    1.000
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 260.07  on 399  degrees of freedom
# Residual deviance: 258.10  on 395  degrees of freedom
# AIC: 268.1
# 
# Number of Fisher Scoring iterations: 15

coef_est <- coef(model4)
se <- sqrt(diag(vcov(model4)))
lower <- coef_est - 1.96 * se
upper <- coef_est + 1.96 * se
OR <- exp(coef_est)
lower_OR <- exp(lower)
upper_OR <- exp(upper)


model5 <- glm(gyrAB_mutated ~ Fq + time_gap + tb_profiler_baseline_drug_resistance2,
              data = data_clean,
              family = binomial)
summary(model5)
# Call:
#   glm(formula = gyrAB_mutated ~ Fq + time_gap + tb_profiler_baseline_drug_resistance2, 
#       family = binomial, data = data_clean)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                    -3.821e+00  8.094e-01  -4.721 2.35e-06 ***
#   Fq1                                             2.060e+00  8.133e-01   2.532   0.0113 *  
#   time_gap                                        4.157e-05  7.442e-04   0.056   0.9555    
# tb_profiler_baseline_drug_resistance2Sensitive -7.052e-01  8.257e-01  -0.854   0.3930    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 260.07  on 399  degrees of freedom
# Residual deviance: 237.86  on 396  degrees of freedom
# AIC: 245.86
# 
# Number of Fisher Scoring iterations: 7

coef_est <- coef(model5)
se <- sqrt(diag(vcov(model5)))
lower <- coef_est - 1.96 * se
upper <- coef_est + 1.96 * se
OR <- exp(coef_est)
lower_OR <- exp(lower)
upper_OR <- exp(upper)


