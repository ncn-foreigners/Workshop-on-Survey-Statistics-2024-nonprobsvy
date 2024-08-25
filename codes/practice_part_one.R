###########################
### Packages Installation
###########################

# remove.packages("nonprobsvy")
# NONPROBSVY
# Install from CRAN
install.packages("nonprobsvy")
# Install a development version
devtools::install_github("ncn-foreigners/nonprobsvy@dev", force = TRUE)

# Other required packages for workshops
install.packages("sampling")
install.packages("survey")
install.packages("ggplot")

library(sampling)
library(survey)
library(nonprobsvy)
library(ggplot2)

###########################
### Simulation
###########################

seed_for_sim <- 2024-8-27
set.seed(seed_for_sim)
N <- 10000
n_A <- 500
p <- 50
alpha_vec1 <- c(-2, 1, 1, 1,1, rep(0, p-5))
alpha_vec2 <- c(0,0,0,3,3,3,3, rep(0, p-7))
beta_vec <- c(1,0,0,1,1,1,1, rep(0, p-7))
X <- cbind("(Intercept)"=1, matrix(rnorm(N*(p-1)), nrow=N, byrow=T, dimnames = list(NULL, paste0("X",1:(p-1)))))
Y <- 1 + as.numeric(X %*% beta_vec) +   rnorm(N) ## linear model
#Y <- 1 + exp(3*sin(as.numeric(X %*% beta_vec))) + X[, "X5"] + X[, "X6"] + rnorm(N) ## nonlinear model
pi_B <- plogis(as.numeric(X %*% alpha_vec1)) ## linear probability
#pi_B <- plogis(3.5 + as.numeric(log(X^2) %*% alpha_vec2) - sin(X[, "X3"] + X[, "X4"]) - X[,"X5"] + X[, "X6"]) ## nonlinear probability
flag_B <- rbinom(N, 1, prob = pi_B)
pi_A <- inclusionprobabilities(0.25 + abs(X[, "X1"]) + 0.03*abs(Y), n_A)
flag_A <- UPpoisson(pik = pi_A)
pop_data <- data.frame(pi_A, pi_B, flag_A, flag_B, Y, X[, 2:p])

X_totals <- colSums(X)
X_means <- colMeans(X[,-1])
sample_A_svy <- svydesign(ids = ~ 1, 
                          probs = ~ pi_A, 
                          pps = "brewer", 
                          data = pop_data[pop_data$flag_A == 1, ])

sample_B <- pop_data[pop_data$flag_B == 1, ]

y_true <- mean(Y)
y_true

y_naive <- mean(sample_B$Y)


###########################
### Inverse Probability Weighting
###########################


### Maximum Likelihood Approach
ipw_logit_sample <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                            target = ~ Y,
                            data = sample_B,
                            svydesign = sample_A_svy,
                            method_selection = "logit")

cbind(ipw_logit_sample$output,ipw_logit_sample$confidence_interval)

### Calibation constrains with GEE aproach
ipw_logit_sample_gee <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                                target = ~ Y,
                                data = sample_B,
                                svydesign = sample_A_svy,
                                method_selection = "logit",
                                control_selection = controlSel(est_method_sel = "gee", h = 2)) 

cbind(ipw_logit_sample_gee$output,ipw_logit_sample_gee$confidence_interval)


###########################
### Mass Imputation
###########################

### Nearest Neighbour
mi_sample_nn <- nonprob(outcome = Y ~ X3 + X4 + X5 + X6,
                        data = sample_B,
                        svydesign = sample_A_svy,
                        method_outcome = "nn",
                        family_outcome = "gaussian",
                        control_outcome = controlOut(k = 3))

cbind(mi_sample_nn$output,mi_sample_nn$confidence_interval)

### Predictive Mean Matching y_hat - y_hat
mi_sample_pmm <- nonprob(outcome = Y ~ X3 + X4 + X5 + X6,
                        data = sample_B,
                        svydesign = sample_A_svy,
                        method_outcome = "pmm",
                        family_outcome = "gaussian")

cbind(mi_sample_pmm$output,mi_sample_pmm$confidence_interval)

### GLM
mi_sample <- nonprob(outcome = Y ~ X3 + X4 + X5 + X6,
                     data = sample_B,
                     svydesign = sample_A_svy,
                     method_outcome = "glm",
                     family_outcome = "gaussian")

cbind(mi_sample$output,mi_sample$confidence_interval)

###########################
### Doubly Robust Method
###########################

### Maximum Likelihood and GLM
dr_logit_sample_mle <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                               outcome = Y ~ X1 + X2 + X3 + X4,
                               data = sample_B,
                               svydesign = sample_A_svy,
                               method_selection = "logit",
                               verbose = TRUE)

cbind(dr_logit_sample_mle$output,dr_logit_sample_mle$confidence_interval)

### GEE and GLM
dr_logit_sample_gee <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                               outcome = Y ~ X1 + X2 + X3 + X4,
                               data = sample_B,
                               svydesign = sample_A_svy,
                               method_selection = "logit",
                               control_selection = controlSel(est_method_sel = "gee", h = 1),
                               verbose = TRUE)

cbind(dr_logit_sample_gee$output,dr_logit_sample_gee$confidence_interval)


### Bias Minimization - must be run on dev version
dr_logit_sample_mm <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                              outcome = Y ~ X1 + X2 + X3 + X4,
                              data = sample_B,
                              svydesign = sample_A_svy,
                              method_selection = "logit",
                              verbose = TRUE,
                              control_inference = controlInf(bias_correction = TRUE))

cbind(dr_logit_sample_mm$output,dr_logit_sample_mm$confidence_interval)

###########################
### Variable Selection
###########################

### SCAD penalty for IPW 
ipw_logit_sample_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                                 target = ~ Y,
                                 data = sample_B,
                                 svydesign = sample_A_svy,
                                 method_selection = "logit",
                                 control_selection = controlSel(penalty = "SCAD"),
                                 control_inference = controlInf(vars_selection = TRUE,
                                                                num_boot = 50),
                                 verbose = TRUE)

cbind(ipw_logit_sample_scad$output,ipw_logit_sample_scad$confidence_interval)

### SCAD penalty for Mass Imputation
mi_sample_scad <- nonprob(outcome = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                          data = sample_B,
                          svydesign = sample_A_svy,
                          method_outcome = "glm",
                          family_outcome = "gaussian",
                          control_outcome = controlOut(penalty = "SCAD"),
                          control_inference = controlInf(vars_selection = TRUE),
                          verbose = TRUE)

cbind(mi_sample_scad$output,mi_sample_scad$confidence_interval)

### SCAD penalty for DR method
dr_logit_sample_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                                outcome = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                                data = sample_B,
                                svydesign = sample_A_svy,
                                method_selection = "logit",
                                control_selection = controlSel(penalty = "SCAD"),
                                control_inference = controlInf(vars_selection = TRUE),
                                verbose = TRUE)

cbind(dr_logit_sample_scad$output, dr_logit_sample_scad$confidence_interval)

### SCAD penalty for DR with bias minimization approach - must be run on dev version

dr_logit_sample_mm_scad <- nonprob(selection = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                                   outcome = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                                   data = sample_B,
                                   svydesign = sample_A_svy,
                                   method_selection = "logit",
                                   control_selection = controlSel(penalty = "SCAD"),
                                   control_inference = controlInf(vars_selection = TRUE, 
                                                                  bias_correction = TRUE),
                                   verbose = TRUE)
cbind(dr_logit_sample_mm_scad$output, dr_logit_sample_mm_scad$confidence_interval)

###########################
### Variable Selection
###########################

### IPW
ipw_logit_sample_boot <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                                     target = ~ Y,
                                     data = sample_B,
                                     svydesign = sample_A_svy,
                                     method_selection = "logit",
                                     control_inference = controlInf(var_method = "bootstrap", num_boot = 100),
                                     verbose = TRUE)

cbind(ipw_logit_sample_boot$output,ipw_logit_sample_boot$confidence_interval)

# MI
mi_sample_nn_boot <- nonprob(outcome = Y ~ X3 + X4 + X5 + X6,
                        data = sample_B,
                        svydesign = sample_A_svy,
                        method_outcome = "glm",
                        family_outcome = "gaussian",
                        control_inference = controlInf(var_method = "bootstrap", num_boot = 100),
                        verbose = TRUE)

cbind(mi_sample_nn$output,mi_sample_nn$confidence_interval)

# DR
dr_logit_sample_mle_boot <- nonprob(selection = ~ X1 + X2 + X3 + X4,
                                    outcome = Y ~ X1 + X2 + X3 + X4,
                                    data = sample_B,
                                    svydesign = sample_A_svy,
                                    method_selection = "logit",
                                    control_inference = controlInf(var_method = "bootstrap", num_boot = 50),
                                    verbose = TRUE)
cbind(dr_logit_sample_mle_boot$output,dr_logit_sample_mle_boot$confidence_interval)

tab_res <- rbind(
  cbind(est="ipw with sample", ipw_logit_sample$output,ipw_logit_sample$confidence_interval),
  cbind(est="ipw with sample (gee h=1)" , ipw_logit_sample_gee$output,ipw_logit_sample_gee$confidence_interval),
  cbind(est="ipw with sample (mle bootstrap)" , ipw_logit_sample_boot$output,ipw_logit_sample_boot$confidence_interval),
  cbind(est="ipw with sample SCAD" , ipw_logit_sample_scad$output,ipw_logit_sample_scad$confidence_interval),
  cbind(est="mi with sample (glm)" , mi_sample$output,mi_sample$confidence_interval),
  cbind(est="mi with sample (nn)" , mi_sample_nn$output,mi_sample_nn$confidence_interval),
  cbind(est="mi with sample (pmm)" , mi_sample_pmm$output,mi_sample_pmm$confidence_interval),
  cbind(est="mi with sample (glm) SCAD" , mi_sample_scad$output, mi_sample_scad$confidence_interval),
  cbind(est="mi with sample (nn bootstrap)" , mi_sample_nn_boot$output, mi_sample_nn_boot$confidence_interval),
  cbind(est="dr with sample (glm)" , dr_logit_sample_mle$output,dr_logit_sample_mle$confidence_interval),
  cbind(est="dr with sample (gee h=1)" , dr_logit_sample_gee$output,dr_logit_sample_gee$confidence_interval),
  cbind(est="dr with sample (bias min)" , dr_logit_sample_mm$output,dr_logit_sample_mm$confidence_interval),
  cbind(est="dr with sample (mle) SCAD" , dr_logit_sample_scad$output,dr_logit_sample_scad$confidence_interval),
  cbind(est="dr with sample (bias min) SCAD", dr_logit_sample_mm_scad$output, dr_logit_sample_mm_scad$confidence_interval),
  cbind(est="dr with sample (mle bootstrap) ", dr_logit_sample_mle_boot$output, dr_logit_sample_mle_boot$confidence_interval)
)

tab_res

ggplot(data = tab_res, aes(x = est, y = mean, ymin = lower_bound, ymax = upper_bound)) +
  geom_point() +
  geom_pointrange() +
  geom_hline(yintercept = y_true, linetype = "dashed", color = "darkgreen") + 
  geom_hline(yintercept = y_naive, linetype = "dashed", color = "red") + 
  coord_flip()
