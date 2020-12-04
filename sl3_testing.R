


library(tidyverse)
library(data.table)
library(SuperLearner)
library(origami)
library(sl3)


################################################################################
# sl3: https://tlverse.org/sl3/
################################################################################

# load example data set
data(cpp)
cpp <- cpp %>%
  dplyr::filter(!is.na(haz)) %>%
  mutate_all(~ replace(., is.na(.), 0))

# use covariates of intest and the outcome to build a task object
covars <- c("apgar1", "apgar5", "parity", "gagebrth", "mage", "meducyrs",
            "sexn")
task <- sl3_Task$new(cpp, covariates = covars, outcome = "haz")

# set up screeners and learners via built-in functions and pipelines
slscreener <- Lrnr_pkg_SuperLearner_screener$new("screen.glmnet")
glm_learner <- Lrnr_glm$new()
screen_and_glm <- Pipeline$new(slscreener, glm_learner)
SL.glmnet_learner <- Lrnr_pkg_SuperLearner$new(SL_wrapper = "SL.glmnet")

# stack learners into a model (including screeners and pipelines)
learner_stack <- Stack$new(SL.glmnet_learner, glm_learner, screen_and_glm)
stack_fit <- learner_stack$train(task)
preds <- stack_fit$predict()
head(preds)

################################################################################
# 
################################################################################

db <- load(file = paste0("~/Google Drive/01_semester_archive/11_Fall_2020/ph240b/", 
                         "final_proj/final_proj.RData"))
db

dim(cov)
dim(gene_exp)
colnames(gene_exp) %>% head

dat <- cov %>% mutate(status = as.integer(status_of_dfs-1), 
                      time = disease_free_survival_months) %>% select(time, status) %>% 
  bind_cols(gene_exp[, 1:1000])

# use covariates of intest and the outcome to build a task object
covars <- grep("ILMN", colnames(dat), value = T)
task <- sl3_Task$new(dat, covariates = covars, outcome = "status")

# set up screeners and learners via built-in functions and pipelines
slscreener <- Lrnr_pkg_SuperLearner_screener$new("screen.glmnet")
glm_learner <- Lrnr_glm$new()
screen_and_glm <- Pipeline$new(slscreener, glm_learner)
SL.glmnet_learner <- Lrnr_pkg_SuperLearner$new(SL_wrapper = "SL.glmnet")


# stack learners into a model (including screeners and pipelines)
learner_stack <- Stack$new(SL.glmnet_learner, glm_learner, screen_and_glm)
stack_fit <- learner_stack$train(task)
preds <- stack_fit$predict()
dim(preds)
head(preds)



#######

dat <- dat %>% select(-ID_REF)
outcome <- "status"
covars <- colnames(dat)[-which(names(dat) == outcome)]

# create the sl3 task
dat_task <- make_sl3_Task(
  data = dat,
  covariates = covars,
  outcome = outcome
)
dat_task


sl3_list_properties()


# choose base learners
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_lasso <- make_learner(Lrnr_glmnet) # alpha default is 1
lrnr_ridge <- make_learner(Lrnr_glmnet, alpha = 0)
lrnr_elasticnet <- make_learner(Lrnr_glmnet, alpha = .5)

# make stack
stack <- make_learner(
  Stack,
  lrnr_glm, lrnr_mean, lrnr_ridge, lrnr_lasso, lrnr_elasticnet
)

# screeners
screen_glmnet <- make_learner(Lrnr_glmnet, alpha = 1)

# update stack with screener
fancy_stack <- make_learner(Stack, screen_glmnet, stack)

# we can visualize the stack
dt_stack <- delayed_learner_train(fancy_stack, dat_task)
plot(dt_stack, color = FALSE, height = "400px", width = "90%")


# make the superlearner
sl <- make_learner(Lrnr_sl, learners = fancy_stack)

# train
sl_fit <- sl$train(dat_task)

#  Get
sl_preds <- sl_fit$predict()
head(sl_preds)

# which covariates are selected on the full data?
s#creen_glmnet$train(dat_task)

