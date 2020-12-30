
rm(list=ls())


# lod pacakges
library(tidyverse)
library(data.table)
library(survival)
library(skimr)
library(scales)
library(origami)
library(future)
library(survtmle)
#options(sl3.verbose = TRUE)
library(sl3)

# set dplyr functions b/c of server issues 
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# load data
db <- load(file = paste0(#"~/Google Drive/01_semester_archive/11_Fall_2020/ph240b/", 
  "final_proj.RData"))
db

dim(cov)
cov %>% head
dim(gene_exp)


#-------------------------------------------------------------------------------
# Set up data
#-------------------------------------------------------------------------------

#cov <- cov %>%
#  mutate(disease_free_survival_months = round(scales::rescale(disease_free_survival_months, c(1, 150)), 0),
#    survival_120_months = ifelse(disease_free_survival_months > 120, 120, 
#                                 disease_free_survival_months),
#    status_120_months = ifelse(disease_free_survival_months < 120 & status_of_dfs == 2, 
#                          1, 0))
#summary(cov$disease_free_survival_months)
#with(cov, table(status_of_dfs, status_120_months))                                 


dat <- cov %>% mutate(id = paste0("ID", 1:nrow(.)), 
                      status = as.integer(status_of_dfs-1),
                      #ftype = status_of_dfs-1, # 1 = recurrence, 0 no reccurence
                      ftype = status, # 1 = recurrence, 0 no reccurence
                      time = disease_free_survival_months,
                      #time = survival_120_months,
                      ftime = time, 
                      cancer_stage = as.factor(cancer_stage), 
                      largest_diameter = as.integer(largest_diameter)) %>% 
  select(ftype, ftime, status, cancer_stage, largest_diameter) %>%  
  bind_cols(gene_exp)

dat %>% select(1:10) %>% head
table(dat$ftype, dat$status) %>% addmargins



#-------------------------------------------------------------------------------
# Estimate P(C>t|W) via KM
#-------------------------------------------------------------------------------

# view survival time by outcome
with(dat, tapply(ftime, ftype, summary))

# define CV function for KM
cv_km <- function(fold) {
  
  train_data <- training(dat)
  valid_data <- validation(dat)
  
  #train_data <- dat %>% sample_n(nrow(.)*.9)
  #valid_data <- dat %>% sample_frac(0.1)
  #data <- dat
  
  # fit KM
  mod <- survfit(Surv(ftime, ftype==0) ~ 0, data = train_data)
  
  # get survival at timepoints for validation
  km_summary <- summary(mod, time = unique(valid_data$ftime))
  
  ## correct estimates to be P(C >= t | A, W)
  #km_fix <- c(1, km_summary$surv)
  km_fix_df <- data.frame(km_fix_probs = c(1, km_summary$surv), 
                          km_fix_times = c(0, km_summary$time))
  
  # extract G_t for validation data
  G_t <- rep(NA, nrow(valid_data))
  for (this_time in km_fix_df$km_fix_times[-1]) { 
    n0 <- sum(valid_data$ftime == this_time)
    G_t[valid_data$ftime == this_time] <- 
      rep(km_fix_df$km_fix_probs[km_fix_df$km_fix_times == this_time], n0)
  }
  #out <- G_t
  return(list(G_t = G_t))
}

# 10-fold CV for KM predictions
dat_folds <- make_folds(dat)
results <- cross_validate(cv_fun = cv_km, 
                          folds = dat_folds,
                          use_future = FALSE, .combine = FALSE, 
                          .combine_control = list())

# add G_t to dataframe indexing by the folds
dat$G_t <- NA
for (i in 1:length(results$G_t)) {
  dat$G_t[dat_folds[[i]]$validation_set] <- results$G_t[[i]]
}

# replace NAs
dat <- dat %>% mutate(G_t = replace(G_t, is.na(G_t), min(G_t, na.rm = T)))
dat %>% arrange(G_t) %>% select(ftime, ftype, G_t) %>% head


# add in weight and create failure indicator
dim(dat)
dat <- dat %>%
  mutate(
    ipcw_ftype = ifelse(ftype == 1, 
                        as.numeric(ftype == 1) / G_t,
                        0),
    # create failure time indicator
    y = as.numeric(ftime > 120)
  )
table(dat$y)
with(dat, tapply(ipcw_ftype, status, summary))
with(dat, table(y, status))

#-------------------------------------------------------------------------------
# sl3
#-------------------------------------------------------------------------------

# First task is to perform variable selection for ------------------------------
# gene expression values

#options(sl3.verbose = TRUE)
# set sl3 task
genes <- grep("ILMN", colnames(dat), value = T)#[1:300]
length(genes)
head(genes)
dat_task <- make_sl3_Task(
  data = dat,
  #covariates = c("cancer_stage", "largest_diameter"),
  covariates = c("cancer_stage", "largest_diameter", genes),
  outcome = "y",
  weights = "ipcw_ftype",
  folds = dat_folds # same folds as used for P(C>t)
)


# set learners --------

# define xgboost learner
grid_params <- list(max_depth = c(4, 6, 8),
                    eta = c(0.001, 0.01, 0.1, 0.2, 0.3),
                    nrounds = c(20, 50))
grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)
params_default <- list(nthread = getOption("sl.cores.learners", 1))
xgb_learners <- apply(grid, MARGIN = 1, function(params_tune) {
  do.call(Lrnr_xgboost$new, c(params_default, as.list(params_tune)))})
  
# caret learner
#lrnr_caret_xgboost <- make_learner(Lrnr_caret, algorithm = "xgboost", 
#                                   tuneLength = 5)

# other learners
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
lrnr_ridge <- make_learner(Lrnr_glmnet, alpha = 0)
lrnr_elasticnet <- make_learner(Lrnr_glmnet, alpha = .5)
lrnr_gam <- make_learner(Lrnr_gam)
#lrnr_polspline <- make_learner(Lrnr_polspline)

# set stack
stack <- make_learner(
  Stack, 
  lrnr_glm, lrnr_mean, lrnr_ridge, lrnr_lasso, lrnr_elasticnet,
  #lrnr_gam, lrnr_polspline,
  xgb_learners[[1]], xgb_learners[[2]], xgb_learners[[3]], xgb_learners[[4]], xgb_learners[[5]]
  ,
  xgb_learners[[6]], xgb_learners[[7]], xgb_learners[[8]], xgb_learners[[9]], xgb_learners[[10]],
  xgb_learners[[11]], xgb_learners[[12]], xgb_learners[[13]], xgb_learners[[14]], xgb_learners[[15]],
  xgb_learners[[16]], xgb_learners[[17]], xgb_learners[[18]], xgb_learners[[19]], xgb_learners[[20]],
  xgb_learners[[21]], xgb_learners[[22]], xgb_learners[[23]], xgb_learners[[24]], xgb_learners[[25]], 
  xgb_learners[[26]], xgb_learners[[27]], xgb_learners[[28]], xgb_learners[[29]], xgb_learners[[30]]
  )

stack <- make_learner(
  Stack, 
  lrnr_glm, lrnr_mean, lrnr_ridge, lrnr_lasso, lrnr_elasticnet
)


# screeners --------
screen_rf <- make_learner(Lrnr_screener_randomForest, #nVar = 20, ntree = 20)
                          nVar = ncol(dat)-1, 
                          #nVar = 20,
                          ntree = 200)
screen_lasso <- make_learner(Lrnr_glmnet)

# add screener to pipline
screener_pipeline <- make_learner(Pipeline, 
                                  screen_lasso,
                                  #screen_rf, 
                                  stack)

# complete stack
fancy_stack <- make_learner(Stack, screener_pipeline, stack)

# set sl
sl <- Lrnr_sl$new(
  learners = fancy_stack,
  #learners = dt_stack,
  loss = loss_loglik_binomial
)

availableCores() 
options(future.fork.enable = FALSE) ## automatically set in RStudio
#supportsMulticore()
#plan(multiprocess)
#f <- future(1)
#class(f)
#f <- future(Sys.getpid())
#value(f)
Sys.getpid()
plan(multiprocess, workers = 2)
options(sl3.verbose = TRUE)
sl_fit <- sl$train(dat_task)
print("#### initial SL fit complete ####")

sl_fit$predict() %>% summary()
sl_fit$learner_fits$`Pipeline(Lrnr_screener_randomForest_20_200->Stack)`


 #save(sl_fit, file = "/my/home/cadams/240b_final_proj/data/sl3_parallel_testing.RData")

save.image(file = "/my/home/cadams/240b_final_proj/data/sl3_parallel_testing.RData")
print("#### save 1 ####")

# cross validate
CV_sl_fit <- CV_lrnr_sl(sl_fit, dat_task, loss_loglik_binomial)
print("#### CV sl fit complete ####")

save.image(file = "/my/home/cadams/240b_final_proj/data/sl3_parallel_testing.RData")

print("#### save 2 ####")


stop()


dt_sl <- delayed_learner_train(sl, dat_task)
dt_sl <- delayed_learner_train(sl_fit, dat_task)

plot(dt_sl, color = FALSE, height = "400px", width = "90%")
dt_sl$compute()

# Variable importance for sl_fit
varimp(sl_fit)



# Seccond task is to use ----------------------------

# nohup R CMD BATCH --vanilla sl3_parallel_testing.R & 
