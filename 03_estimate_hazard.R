


rm(list=ls())

# lod pacakges
library(tidyverse)
library(data.table)
library(survival)
library(skimr)
library(scales)
library(origami)
#library(future)
library(survtmle)
options(sl3.verbose = TRUE)
library(sl3)

# set dplyr functions b/c of server issues 
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

################################################################################
# load data
################################################################################

# data 
db <- load(file = "data/clean_processed_data.RData")
db
dim(dat)
dat %>% select(1:6, ipcw_ftype) %>% head

# sl3 output
db <- load("~/240b_final_proj/data/sl3_fit.RData")
db <- load("~/240b_final_proj/data/sl3_cv_fit.RData") # cv output
db
CV_sl_fit

# remove gene expression from dat
# and make trt variables
set.seed(12345)
dat <- dat %>% select(-starts_with("ILMN")) %>% 
  mutate(trt = as.numeric(runif(nrow(.)) > 0.6))
table(dat$trt)

# create long format dataset
t0 <- 120 # we want the survival up to t0 = t
long_dat <- dat %>% 
  rename(trt = trt) %>%
  mutate(
    id = row_number()
  ) %>%
  makeDataList(
    # can ignore these options or see ?makeDataList for more info, but briefly
    J = 1,                 # the failure type of interest is ftype==1
    ntrt = 2,              # two treatment options: vaccine or placebo arms
    uniqtrt = c(0, 1),     # unique treatment types: vaccine or placebo arms
    t0 = t0                # we're interested in timepoints 1-4
  )
dat_long_data <- as.data.table(long_dat[[1]])
dim(dat_long_data)
head(dat_long_data)
stop()

dat_long_data[id == 1]  
# NOTE: there are only as many rows as time points for each observation
dat_long_data[id == 5]
# and that N1 is an indicator of an infection (failure of interest) at time t
dat_long_data[id == 137]


dat_long_data %>% head

# fit a single Super Learner for all hazard estimation regressions on {A, W, t}

# create tasks
dat_long_data_task <- make_sl3_Task(
  data = dat_long_data,
  covariates = c("t", "trt", "cancer_stage", "largest_diameter"),
  outcome = "N1",
  # make sure data splitting for cross-validation is done by person, not by row
  id = "id",
  weights = "ipcw_ftype"
)

# specify the Super Learner ensemblel
lrnr_mean <- Lrnr_mean$new()
lrnr_glm <- Lrnr_glm_fast$new()
sl_mod <- Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_glm),
  #slscreener_glment, slscreener_rf),
  metalearner = Lrnr_nnls$new()
)
sl_hazard_all_times <- sl_mod$train(dat_long_data_task)

# make datatable 
dat_dt <- as.data.table(dat)
dat_dt[, id := 1:.N]
dat_dt[, y := as.numeric(ftype == 1)]
dat_dt <- dat_dt[rep(1:.N, t0)]
dat_dt[, t := 1:t0, by = "id"]
setorder(dat_dt, id)
dim(dat_dt)
dat_dt %>% head
table(dat_dt$id)
ftype#setnames(dat_dt, "vax", "trt")

length(dat_dt$id)
length(sl_hazard_all_times$predict(dat_long_data_task))
surv_est_long <- data.table(
  t = dat_long_data$t,
  id = dat_long_data$id,
  one_minus_hazard = 1 - sl_hazard_all_times$predict(dat_long_data_task)
)
surv_est_long %>% head
surv_est <- surv_est_long[, list(surv = prod(one_minus_hazard), freq = .N),
                          by = "id"]
surv_est[, freq := NULL]
surv_est[, cum_inc := (1 - surv)]
surv_est %>% head

surv_est %>% head
surv_est_long %>% head
surv_est_long %>% tail
surv_est_long %>% ggplot() + 
  geom_step(aes(x = t, y = one_minus_hazard, color = id))


# survtmle  --------

dat <- dat %>% mutate(A = trt, 
                      W1 = cancer_stage, W2 = largest_diameter)
t0 <- 120
fit <- survtmle(ftime = dat$ftime, ftype = dat$ftype,
                trt = dat$trt, adjustVars = select(dat, W1, W2),
                glm.trt = "W1 + W2",
                glm.ftime = "trt + W1 + W2 + t",
                glm.ctime = "W1 + W2 + t",
                method = "hazard",
                t0 = t0)







# extract cumulative incidence at each timepoint
tpfit <- timepoints(fit, times = seq(10, 120, 10))
tpfit
fit$ic
plot(tpfit)

# construct simultaneous confidence intervals --------
times <- seq(10, 120, 10)
fit_iter <- lapply(times, function(t0) {
  tmp <- survtmle(ftime = dat$ftime, ftype = dat$ftype,
                  trt = dat$trt, adjustVars = select(dat, W1, W2),
                  glm.trt = "W1 + W2",
                  glm.ftime = "trt + W1 + W2 + t",
                  glm.ctime = "W1 + W2 + t",
                  method = "hazard",
                  t0 = t0)
  return(tmp)
})


# collect eICs
ic_all_t <- lapply(1:length(fit_iter), function(x) {
  tmp <- data.frame(fit_iter[[x]]$ic)
  colnames(tmp) <- paste0(colnames(tmp), "_", times[x])
  return(tmp)
}) %>% reduce(bind_cols) %>% head
ic_z0 <- ic_all_t %>% select(contains("z0"))
ic_z1 <- ic_all_t %>% select(contains("z1"))

corr_ic_est <- cor(ic_z0)
ci_level <- 0.95
library(mvtnorm)
mvn_ic_mult <- qmvnorm(ci_level, tail = "both", corr = corr_ic_est)
ci_mult <- c(-1, 1) * mvn_ic_mult$quantile

length(tpfit)

haz_ci <- lapply(1:length(tpfit), function(x) {
  trt0 <- tpfit[[x]]$est[1] + ci_mult * var(ic_z0[, x]) # trt=0
  trt1 <- tpfit[[x]]$est[2] + ci_mult * var(ic_z1[, x]) # trt=1
  return(c(tpfit[[x]]$est[1], trt0, tpfit[[x]]$est[2], trt1))
}) %>% reduce(rbind) %>% data.frame %>% 
  rename(haz_trt0 = X1, ll_trt0 = X2, ul_trt0 = X3, 
         haz_trt1 = X4, ll_trt1 = X5, ul_trt1 = X6) %>% 
  mutate(t = times) %>% select(7, 1:6)
haz_ci %>% select(1:4) %>% bind

t0 <- select(haz_ci, 1:4) %>% mutate(trt="0")
t1 <- select(haz_ci, 1,5:7) %>% mutate(trt="1")

names(t1) <- names(t0) <- c("t", "haz", "ll", "ul", "trt")
rbind(t1, t0) %>% ggplot(aes(x=t, y = haz, color = trt)) + 
  #geom_step() + 
  #geom_step(aes(x=t, y = ll, color = trt), linetype = 2) + 
  #geom_step(aes(x=t, y = ul, color = trt), linetype = 2)
  geom_smooth(aes(x=t, y = haz, color = trt), method = "loess", linetype = 1) + 
  geom_smooth(aes(x=t, y = ul, color = trt), method = "loess", linetype = 2) + 
  geom_smooth(aes(x=t, y = ll, color = trt), method = "loess", linetype = 2) + 
  scale_x_discrete(name = "Time", limits = as.character(times)) + theme_bw()

