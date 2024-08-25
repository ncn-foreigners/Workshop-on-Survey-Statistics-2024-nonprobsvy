## ----install, eval=FALSE---------------------------------------------------------------------------------------------
## install.packages("ggplot2") ## for visualisation
## install.packages("survey") ## for probabilistic survey
## install.packages("nonprobsvy") ## lastest version from CRAN
## remove.packages("nonprobsvy")

## ----packages, message = FALSE---------------------------------------------------------------------------------------
library(nonprobsvy)
library(survey)
library(ggplot2)


## ----reading-jvs-----------------------------------------------------------------------------------------------------
jvs <- read.csv("data-raw/jvs.csv",
                  colClasses = c("character", "numeric", rep("character", 3), "numeric"))
head(jvs)


## ----svydesign-------------------------------------------------------------------------------------------------------
jvs_svy <- svydesign(ids = ~ 1,  ## stratified sampling (no clustering)
                     weights = ~ weight,  ## final weight
                     strata = ~ size + nace + region, ## stratification variables
                     data = jvs
                     #,calibrate.formula = ~ size + nace + region ## if weights are calibrated to known population sizes in these variables see: https://cran.r-project.org/web/packages/survey/vignettes/precalibrated.pdf
                     )

svytotal(~size, jvs_svy)


## ----data-reading------------------------------------------------------------------------------------------------------
admin <- read.csv("data-raw/admin.csv",
                 colClasses = c("character", "numeric", rep("character", 3), "logical")
                 )
head(admin)


## ----ipw-stadnard----------------------------------------------------------------------------------------------------
est1_logit <- nonprob(
  selection = ~ region + private + nace + size,
  target = ~ single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit"
)

est1_probit <- nonprob(
  selection = ~ region + private + nace + size,
  target = ~ single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "probit"
)

est1_cloglog <- nonprob(
  selection = ~ region + private + nace + size,
  target = ~ single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "cloglog"
)

est_ipw_standard <- rbind(cbind(est1_logit$output, est1_logit$confidence_interval),
                          cbind(est1_probit$output, est1_probit$confidence_interval),
                          cbind(est1_cloglog$output, est1_cloglog$confidence_interval))
est_ipw_standard$est <- "ipw"
rownames(est_ipw_standard) <- NULL
est_ipw_standard


## ----ipw-structure---------------------------------------------------------------------------------------------------
str(est1_logit,1)


## ----ipw-structure-sel-----------------------------------------------------------------------------------------------
str(est1_logit$selection,1)


## ----ipw-summary-----------------------------------------------------------------------------------------------------
summary(est1_logit)


## ----ipw-cal---------------------------------------------------------------------------------------------------------
est2_logit <- nonprob(
  selection = ~ region + private + nace + size,
  target = ~ single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_selection = controlSel(h = 1, est_method_sel = "gee")
)

est2_probit <- nonprob(
  selection = ~ region + private + nace + size,
  target = ~ single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "probit",
  control_selection = controlSel(h = 1, est_method_sel = "gee")
)

est2_cloglog <- nonprob(
  selection = ~ region + private + nace + size,
  target = ~ single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "cloglog",
  control_selection = controlSel(h = 1, est_method_sel = "gee")
)

est_ipw_calib <- rbind(cbind(est2_logit$output, est2_logit$confidence_interval),
                       cbind(est2_probit$output, est2_probit$confidence_interval),
                       cbind(est2_cloglog$output, est2_cloglog$confidence_interval))
est_ipw_calib$est <- "ipw calib"
rownames(est_ipw_calib) <- NULL
est_ipw_calib


## ----ipw-cal-weights----------------------------------------------------------------------------------------------------
admin$ipw1_weight <- est1_logit$weights
admin$ipw2_weight <- est2_logit$weights


## ----ipw-cal-totals--------------------------------------------------------------------------------------------------
c(jvs=sum(weights(jvs_svy)), ipw1_mle=sum(admin$ipw1_weight), ipw2_gee=sum(admin$ipw2_weight))


## ----ipw-cal-klasa1--------------------------------------------------------------------------------------------------
svytotal(~size, jvs_svy)


## ----ipw-cal-klasa2--------------------------------------------------------------------------------------------------
xtabs(ipw1_weight ~ size, admin)
xtabs(ipw2_weight ~ size, admin)


## ----ipw-bootstrap---------------------------------------------------------------------------------------------------
set.seed(2024-08-27)
est3_logit <- nonprob(
  selection = ~ region + private + nace + size,
  target = ~ single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_inference = controlInf(var_method = "bootstrap", num_boot = 50),
  verbose = T,
)

## ----ipw-bootstrap-summary-------------------------------------------------------------------------------------------
summary(est3_logit)


## ----ipw-scad--------------------------------------------------------------------------------------------------------
set.seed(2024-08-27)
est4_logit <- nonprob(
  selection = ~ region + private + nace + size,
  target = ~ single_shift,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  control_selection = controlSel(nfolds = 5, nlambda = 10),
  control_inference = controlInf(vars_selection = TRUE),
  verbose = TRUE
)


## ----ipw-scad-summary------------------------------------------------------------------------------------------------
summary(est4_logit)


## ----ipw-comparison---------------------------------------------------------------------------------------------------
ipw_summary <- rbind(cbind(est1_logit$output, est1_logit$confidence_interval),
                     cbind(est2_logit$output, est2_logit$confidence_interval),
                     cbind(est3_logit$output, est3_logit$confidence_interval),
                     cbind(est4_logit$output, est4_logit$confidence_interval))
rownames(ipw_summary) <- NULL
ipw_summary$est <- c("ipw (st)", "ipw (cal)", "ipw (boot)", "ipw (scad)")
ipw_summary


## ----mi-glm-lp-------------------------------------------------------------------------------------------------------
est5_glm <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "glm",
  family_outcome = "gaussian"
)

cbind(est5_glm$output,est5_glm$confidence_interval)


## ----mi-glm-binom----------------------------------------------------------------------------------------------------
est5_glm_biom <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "glm",
  family_outcome = "binomial"
)
cbind(est5_glm_biom$output,est5_glm_biom$confidence_interval)

## ----mi-glm-binom-summary--------------------------------------------------------------------------------------------
summary(est5_glm_biom)


## ----mi-glm-binom-structure------------------------------------------------------------------------------------------
str(est5_glm_biom$outcome,1)


## ----mi-glm-nn-------------------------------------------------------------------------------------------------------
est6_glm_nn <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "nn",
  control_outcome = controlOut(k=5)
)
cbind(est6_glm_nn$output,est6_glm_nn$confidence_interval)


## ----mi-glm-pmm-1----------------------------------------------------------------------------------------------------
set.seed(2023-12-8)
est6_glm_pmm1 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "pmm",
  control_outcome = controlOut(k=5, predictive_match=1),
  verbose = TRUE
)
cbind(est6_glm_pmm1$output, est6_glm_pmm1$confidence_interval)


## ----mi-glm-pmm-2----------------------------------------------------------------------------------------------------
set.seed(2023-12-8)
est6_glm_pmm2 <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "pmm",
  control_outcome = controlOut(k=5, predictive_match=2),
  verbose = TRUE
)
cbind(est6_glm_pmm2$output, est6_glm_pmm2$confidence_interval)


## ----mi-glm-scad-----------------------------------------------------------------------------------------------------
set.seed(2024-4-24)
est7_glm_sel <- nonprob(
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_outcome = "glm",
  family_outcome = "binomial",
  control_outcome = controlOut(nfolds = 5, nlambda = 10),
  control_inference = controlInf(vars_selection = TRUE),
  verbose = TRUE
)

## ----mi-glm-scad-result-----------------------------------------------------------------------------------------------
cbind(est7_glm_sel$output,est7_glm_sel$confidence_interval)

## ----mi-glm-scad-summary---------------------------------------------------------------------------------------------
summary(est7_glm_sel)


## ----mi-summary--------------------------------------------------------------------------------------------------
mi_summary <- rbind(cbind(est5_glm$output, est5_glm$confidence_interval),
                     cbind(est5_glm_biom$output, est5_glm_biom$confidence_interval),
                     cbind(est6_glm_nn$output, est6_glm_nn$confidence_interval),
                     cbind(est6_glm_pmm1$output, est6_glm_pmm1$confidence_interval),
                    cbind(est6_glm_pmm2$output, est6_glm_pmm2$confidence_interval),
                     cbind(est7_glm_sel$output, est7_glm_sel$confidence_interval))
rownames(mi_summary) <- NULL
mi_summary$est <- c("mi (lm)", "mi (glm)", "mi (nn)", "mi (pmm1)", "mi (pmm2)", "mi (glm, scad)")
mi_summary


## ----dr-glm-binom----------------------------------------------------------------------------------------------------
est8_dr1 <- nonprob(
  selection = ~ region + private + nace + size,
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  method_outcome = "glm",
  family_outcome = "binomial",
  pop_size = sum(weights(jvs_svy))
)

cbind(est8_dr1$output,est8_dr1$confidence_interval)

## ----dr-glm-binom-summary--------------------------------------------------------------------------------------------
summary(est8_dr1)


## ----dr-glm-binom-structure------------------------------------------------------------------------------------------
str(est8_dr1,1)


## ----dr-glm-calib----------------------------------------------------------------------------------------------------
est8_dr2 <- nonprob(
  selection = ~ region + private + nace + size,
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  method_outcome = "glm",
  family_outcome = "binomial",
  control_selection = controlSel(h = 1, est_method_sel = "gee")
)

cbind(est8_dr2$output,est8_dr2$confidence_interval)


## ----dr-glm-bootstrap------------------------------------------------------------------------------------------------
set.seed(2024-4-24)
est8_dr3 <- nonprob(
  selection = ~ region + private + nace + size,
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  method_outcome = "glm",
  family_outcome = "binomial",
  control_inference = controlInf(var_method = "bootstrap", num_boot = 50),
  verbose = TRUE
)

cbind(est8_dr3$output,est8_dr3$confidence_interval)


## ----dr-glm-scad-----------------------------------------------------------------------------------------------------
set.seed(2024-4-24)
est9_dr1 <- nonprob(
  selection = ~ region + private + nace + size,
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  method_outcome = "glm",
  family_outcome = "binomial",
  control_selection = controlSel(nfolds = 5, nlambda = 10),
  control_outcome = controlOut(nfolds = 5, nlambda = 10),
  control_inference = controlInf(vars_selection = TRUE),
  verbose = TRUE
)

cbind(est9_dr1$output,est9_dr1$confidence_interval)


## ----dr-glm-scad-bias-min--------------------------------------------------------------------------------------------
set.seed(2024-4-24)
est9_dr2 <- nonprob(
  selection = ~ region + private + nace + size,
  outcome = single_shift ~ region + private + nace + size,
  svydesign = jvs_svy,
  data = admin,
  method_selection = "logit",
  method_outcome = "glm",
  family_outcome = "binomial",
  control_selection = controlSel(nfolds = 5, nlambda = 10),
  control_outcome = controlOut(nfolds = 5, nlambda = 10),
  control_inference = controlInf(vars_selection = TRUE, bias_correction = TRUE),
  verbose = TRUE
)

cbind(est9_dr2$output,est9_dr2$confidence_interval)


## ----dr-summary-------------------------------------------------------------------------------------------------
dr_summary <- rbind(cbind(est8_dr1$output, est8_dr1$confidence_interval),
                    cbind(est8_dr2$output, est8_dr2$confidence_interval),
                    cbind(est8_dr3$output, est8_dr3$confidence_interval),
                    cbind(est9_dr1$output, est9_dr1$confidence_interval),
                    cbind(est9_dr2$output, est9_dr2$confidence_interval))
rownames(dr_summary) <- NULL
dr_summary$est <- c("dr (ipw)", "dr (ipw cal)", "dr (ipw, boot)", "dr (scad)", "dr (scad, min)")
dr_summary


## ----summary-plot------------------------------------------------------------------------------------------------
wyniki <- rbind(ipw_summary, mi_summary, dr_summary)

ggplot(data = wyniki, aes(y = est, x = mean, xmin = lower_bound, xmax = upper_bound)) +
  geom_point() +
  geom_vline(xintercept = mean(admin$single_shift), linetype = "dotted", color = "red") +
  geom_errorbar() +
  labs(x = "Point estimator and confidence interval", y = "estimators")

