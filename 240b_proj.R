

rm(list=ls())

library(data.table)
library(tidyverse)





# read in raw data
dat <- fread("data/GSE44001_series_matrix.txt", fill = T) %>% 
  t() %>% data.frame()

dim(dat)
class(dat)
dat$X1

dat# fix column names
dat[1:5, 1:10] 
dat[1:5, 11:20] 
dat[1:5, 21:30] 
dat[1:5, 31:40] 
dat[1:5, 41:50]
dat[1:5, 51:60]
dat[1:5, 61:70]
dat[, (ncol(dat)-10):ncol(dat)]
table(dat$X38)



# get subject data --------

cov <- dat %>% select(31, 33:38, 46) %>% slice(2:n())
colnames(cov) <- c("Sample_type", "Sample_source_name", "Sample_organism_ch1", 
                   "cancer_stage", "largest_diameter", "disease_free_survival_months", 
                   "status_of_dfs", "sample_description")

# fix data
cov$cancer_stage <- with(cov, str_split(cancer_stage, ": ")) %>% 
  do.call(rbind, .) %>% data.frame %>% select(2) %>% unlist() %>% as.character()
cov$disease_free_survival_months <- 
  with(cov, str_split(disease_free_survival_months, ": ")) %>% 
  do.call(rbind, .) %>% data.frame %>% select(2) %>% unlist() %>% as.numeric() %>% 
  rescale(., c(0, 120)) 
cov$largest_diameter <- 
  with(cov, str_split(largest_diameter, ": ")) %>% 
  do.call(rbind, .) %>% data.frame %>% select(2) %>% unlist() %>% as.numeric()
cov$status_of_dfs <- 
  with(cov, str_split(status_of_dfs, ": ")) %>% 
  do.call(rbind, .) %>% data.frame %>% select(2) %>% unlist() %>% as.numeric()

with(cov, table(cancer_stage))
with(cov, table(status_of_dfs))
with(cov, summary(disease_free_survival_months))
with(cov, summary(largest_diameter))

with(cov, tapply(disease_free_survival_months, cancer_stage, summary)) %>% do.call(cbind,.)
with(cov, tapply(largest_diameter, cancer_stage, summary)) %>% do.call(cbind,.)
with(cov, tapply(disease_free_survival_months, status_of_dfs, summary))
table(cov$status_of_dfs)

library(scales)
rescale(cov$disease_free_survival_months, c(0, 120)) %>% 
  tapply(., cov$status_of_dfs, summary)
tapply(cov$disease_free_survival_months/2.5, cov$status_of_dfs, summary)
tapply(cov$disease_free_survival_months/30, cov$status_of_dfs, summary)
mean(cov$disease_free_survival_months/30)
# make summary table of cov data

library(tableone)
library(xtable)
tbl <- cov %>% #rename(Recurrance = status_of_dfs,)
  CreateTableOne(data = ., 
                       vars = c("disease_free_survival_months", 
                                "cancer_stage", 
                                "largest_diameter"), 
                       strata = "status_of_dfs", test = F) 

tabAsStringMatrix <- print(tbl, printToggle = FALSE, noSpaces = TRUE)
xtable(tabAsStringMatrix)

# get RNA data --------

# select ID, and RNA expresssion measurements
gene_exp <- dat %>% select(60:ncol(.)) %>% slice(2:n()) 

# get ID and gene names
var_names <- dat %>% select(60:ncol(.)) %>% slice(1) %>% 
  unlist() %>% unname() %>% as.character()

# apply to data
colnames(gene_exp) <- var_names
gene_exp %>% select((ncol(.)-5):ncol(.)) %>% head
gene_exp <- gene_exp %>% select(-ncol(.)) # drop last column

# make gene epxression levels numeric
# ID_REF is character
gene_exp <- gene_exp %>% mutate_if(is.factor, as.numeric) %>% 
  mutate(ID_REF = as.character(ID_REF))

# combine data and save --------

dim(cov)
dim(gene_exp)


missing_pp <- colMeans(is.na(gene_exp))
hist(missing_pp)
table(missing_pp)

#save(list = c("cov", "gene_exp"), file = "~/Google Drive/01_semester_archive/11_Fall_2020/ph240b/final_project/final_proj.RData")
load(file = paste0("~/Google Drive/01_semester_archive/11_Fall_2020/ph240b/", 
                   "final_project/final_proj.RData"))

# time plot --------
cov$id <- factor(1:nrow(cov), 
                 levels = order(cov$disease_free_survival_months,
                                cov$status_of_dfs,
                                decreasing = T))
cov %>% #sample_frac(0.25) %>% 
  ggplot(aes(x = id)) + 
  geom_linerange(aes(ymin = 0, ymax = disease_free_survival_months)) + 
  geom_point(aes(y = disease_free_survival_months, 
                 shape = factor(status_of_dfs)), 
             size = 2) + 
  coord_flip() + facet_wrap(~status_of_dfs)

p1 <- cov %>% #sample_frac(0.25) %>% 
  mutate(Outcome = ifelse(status_of_dfs == 2, "Recurrence", "No Recurrence")) %>% 
  ggplot(aes(x = disease_free_survival_months, 
             fill = Outcome)) + 
  geom_density(alpha = 0.25) + 
  labs(x = "Disease free survival") + theme_bw()


dat <- structure(list(ID = 1:5, eventA = c(0L, 1L, 1L, 0L, 1L), 
               eventB = c(1L, 0L, 0L, 1L, 0L), t1 = c(7, 5, 10, 4.5, 2),
               t2 = c(7, 5, 10, 4.5, 8), censored = c(0, 0, 0, 0, 1)), 
               .Names = c("ID", "eventA", "eventB", "t1", "t2", "censored"), 
               class = "data.frame", row.names = c(NA, -5L))
dat$event <- with(dat, ifelse(eventA, "A", "B"))
dat$id.ordered <- factor(x = dat$ID, levels = order(dat$t2, decreasing = T))
dat %>% head


# clustering --------

# hclust
h_complete <- gene_exp %>% select(-1) %>% dist() %>% hclust(., method = "complete")
h_avg <- gene_exp %>% select(-1) %>% dist() %>% hclust(., method = "average")
h_med <- gene_exp %>% select(-1) %>% dist() %>% hclust(., method = "median")
h_cent <- gene_exp %>% select(-1) %>% dist() %>% hclust(., method = "centroid")

plot(h_complete)
plot(h_avg)
#plot(h_med)
#plot(h_cent)

# complete works best
plot(h_complete)

# variable importantance


stop()

################################################################################
# get gene annotation for illumina probes
################################################################################

probeID=c("ILMN_1690170", "ILMN_2410826", "ILMN_1675640", "ILMN_1801246",
          "ILMN_1658247", "ILMN_1740938", "ILMN_1657871", "ILMN_1769520",
          "ILMN_1778401")
library(illuminaHumanv4.db) #Get this library if you don't have

data.frame(Gene=unlist(mget(x = probeID,envir = illuminaHumanv4SYMBOL)))

