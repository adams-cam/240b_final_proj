# 01_get_data.R


# publication:
# Genetic profiling to predict recurrence of early cervical cancer
# GEO asscension: GSE44001
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44001
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44001/matrix/


library(tidyverse)
library(data.table)
library(utils)
library(R.utils) 
library(scales)
library(tableone)
library(xtable)

################################
# Extract data from internet
################################

# download normalized expression data and available covariates from GEO
url_link <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44001/", 
                   "matrix/GSE44001_series_matrix.txt.gz")
download.file(url = url_link , 
                     destfile = "data/GSE44001_series_matrix.txt.gz", 
                     method = "auto")

url_link <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44001/", 
                   "suppl/GSE44001_non-normalzied.txt.gz")
download.file(url = url_link , 
                     destfile = "data/GSE44001_non-normalzied.txt.gz", 
                     method = "auto")


# extract .gz
gunzip("data/GSE44001_series_matrix.txt.gz")
gunzip("data/GSE44001_non-normalzied.txt.gz")


################################
# read in data and create dataset
################################

# read in data -------

# normalized gene exprssion values and covariates
dat <- fread("data/GSE44001_series_matrix.txt", fill = T) %>% t() %>% data.frame()

# non-normalized data with QA/QC information
raw_qc <- fread("data/GSE44001_non-normalzied.txt.gz", fill = T) 


# get subject data --------
dat %>% select(31, 33:38, 46) %>% head
cov <- dat %>% select(31, 33:38, 46) %>% slice(2:n())
colnames(cov) <- c("Sample_type", "Sample_source_name", "Sample_organism_ch1", 
                   "cancer_stage", "largest_diameter", "disease_free_survival_months", 
                   "status_of_dfs", "sample_description")

# fix covarites
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


# check covariates
with(cov, table(cancer_stage))
with(cov, table(status_of_dfs))
with(cov, summary(disease_free_survival_months))
with(cov, summary(largest_diameter))

with(cov, tapply(disease_free_survival_months, cancer_stage, summary)) %>% do.call(cbind,.)
with(cov, tapply(largest_diameter, cancer_stage, summary)) %>% do.call(cbind,.)
with(cov, tapply(disease_free_survival_months, status_of_dfs, summary))
table(cov$status_of_dfs)

# make summary table of cov data 
tbl <- cov %>% #rename(Recurrance = status_of_dfs,)
  CreateTableOne(data = ., 
                 vars = c("disease_free_survival_months", 
                          "cancer_stage", 
                          "largest_diameter"), 
                 strata = "status_of_dfs", test = F) 

tabAsStringMatrix <- print(tbl, printToggle = FALSE, noSpaces = TRUE)
#xtable(tabAsStringMatrix) # latex for presentations


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
tbl <- table(gene_exp < 0.001)

# identify variables that failed detection p-value threshold (p<0.05)
dim(raw_qc)
raw_qc %>% select(590:601)

# fix var names
new_var_names <- raw_qc %>% slice(5) %>% unlist
new_var_names <- sub(" ", "_", new_var_names) # remove spaces
new_var_names[grep("Detection_Pval", new_var_names)] <- 
  paste0("det_pval_SAMPLE_", 1:300) # make unique det P varnames

# subset data only
raw_qc <- raw_qc %>% slice(6:n())
colnames(raw_qc) <- new_var_names
raw_qc %>% select(1:10) %>% headcoul

# identify bad reads
det_pval_mat <- raw_qc %>% select(contains("det")) # remove measurements
table(det_pval_mat < 0.05) %>% prop.table()# identify which are bad
pval_cutoff_pp <- rowMeans(det_pval_mat < 0.05) # get % of p<0.05 for each sample
hist(pval_cutoff_pp)

det_pval_mat <- apply(., 2, function(gene) { gene < 0.05 })
raw_qc %>% select(1:10) %>% head
# combine data and save --------

dim(cov)
dim(gene_exp)


missing_pp <- colMeans(is.na(gene_exp))
hist(missing_pp)
table(missing_pp)

#save(list = c("cov", "gene_exp"), file = "~/Google Drive/01_semester_archive/11_Fall_2020/ph240b/final_project/final_proj.RData")
stop()

# use for 
tmp <- fread("~/Downloads/GSE44001_non-normalzied.txt", fill = T) #%>% 
dim(tmp)
tmp[, 590:601]
tmp %>% head
  t() %>% data.frame()
tmp %>% select(31, 33:38, 46) %>% head


