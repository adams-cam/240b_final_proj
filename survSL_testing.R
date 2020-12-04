

rm(list=ls())

db <- load(file = paste0("~/Google Drive/01_semester_archive/11_Fall_2020/ph240b/", 
                   "final_proj/final_proj.RData"))
db

dim(cov)
dim(gene_exp)

#install.packages("remotes")
#remotes::install_github("tlverse/sl3")
#install.packages("skimr")
#install.packages("~/Downloads/SurvSL", repos = NULL, type = "source")

library(sl3)
library(tidyverse)
library(skimr)
library(survival)
library(survSL)
library(broom)


# get data together
cov %>% head
gene_exp_only <- gene_exp %>% select(-1)

dat <- cov %>% mutate(status = as.integer(status_of_dfs-1), 
                      time = disease_free_survival_months) %>% select(time, status) 
dat %>% head


# identify variables missing > 5% of their data, and remove them
# median imputation for the rest
gene_exp_only_miss_mat <- gene_exp_only %>% is.na() %>% colMeans
gene_exp_only %>% na.omit %>% dim
summary(gene_exp_only_miss_mat)
table(gene_exp_only_miss_mat)
gene_exp_only_miss_mat_na <- gene_exp_only_miss_mat %>% 
  select(!gene_exp_only_miss_mat)
gene_exp_only_miss_mat_na <- 
  apply(gene_exp_only_miss_mat_na, 2, function(x) {
  x[is.na(x)] <- median(x, na.rm = T)
})

dim(gene_exp_only)
X <- data.frame(t(na.omit(t(gene_exp_only[, 1:100]))))
X <- gene_exp_only[, 10:75]
X <- gene_exp_only[, c(10:75, 77:100)]
gene_exp_only_log <- log(gene_exp_only)
X <- gene_exp_only_log[, 1:10000]
#X <- gene_exp_only_log
dim(X)
design <- model.matrix(~., X)

# create surv object 
class(dat$time)
class(dat$status)
table(dat$status)

Y <- with(dat, Surv(time, status))
class(Y)
as.matrix(Y)

# define sl_library
survSL::listWrappers()
sl_library <- c("survSL.glmnet", "survSL.coxph")
#sl_library <- c("survSL.coxph")
#sl_library <- c("survSL.glmnet")

# sl
sl_fit1 <- survSuperLearner(Y = Y, X = design, 
                        SL.library = sl_library, 
                        method = "method.convex1", 
                        verbose = T, 
                        cvControl = list(V = 5L))

sl_fit1
sl_fit1$coef
#sl3::varimp(sl_fit1)
sl_fit1$fitLibrary
sl_fit1$library.predict
sl_fit1$SL.predict

# variable importance

#Y <- dat
library(glmnet)
library(glmnetUtils)
dim(X)
class(X)
class(Y)
options(expressions = 5e5)
Cstack_info()
design <- model.matrix(~., X)
glmnet_fit1 <- cv.glmnet(x = design, y = as.matrix(Y), family = "cox", alpha = 1)
#glmnet_fit1 <- cv.glmnet(x = design, y = as.matrix(Y), family = "cox", alpha = 0.1)

top_hits <- c("ILMN_2339377", "ILMN_3243185", "ILMN_1774974", "ILMN_2120575", 
              "ILMN_1710495", "ILMN_2396198", "ILMN_2178775", "ILMN_1813544",
              "ILMN_2123665", "ILMN_1687840", "ILMN_1743319", "ILMN_2388746")
grep(paste0(top_hits, collapse = "|"), colnames(gene_exp_only))
X <- gene_exp_only %>% select(all_of(top_hits)) %>% log
design <- model.matrix(~., X)
glmnet_fit1 <- cv.glmnet(x = design, y = as.matrix(Y), family = "cox", alpha = 1)
coef(glmnet_fit1)
plot(glmnet_fit1)

coef(glmnet_fit1, s = glmnet_fit1$lambda.1se)
coef(glmnet_fit1, s = glmnet_fit1$lambda)
coef.min = coef(glmnet_fit1, s = "lambda.min")
active.min = which(coef.min != 0)
index.min = coef.min[active.min]
index.min
plot(glmnet_fit1)

stop()
# testing
set.seed(1)
N <- 200
X <- matrix(rnorm(N*10), N, 10)
event_time <- rexp(N, rate = exp(0.2*X[,1] + -0.9*X[,2] - 0.4*X[,3] + 0.1*X[,1]*X[,4] + rnorm(N, 0, .5)))
censor_time <- rexp(N, rate = 1/2)
time <- pmin(event_time, censor_time)  # observed follow-up time
status <- as.numeric(event_time <= censor_time) # event indicator
X <- as.data.frame(X)
Y <- Surv(time, status)
Y

sl_library <- c("survSL.glmnet", "survSL.coxph")
# sl_library <- c("survSL.glmnet", "survSL.coxph", "survSL.CoxBoost", "survSL.gbm", "survSL.glmboost")
fit <- survSuperLearner(Y = Surv(time, status), X = X, SL.library = sl_library, method = 'method.convex1')
fit
