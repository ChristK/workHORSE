## project population to 2071 ( 2021 + 50)
## adapting file ./ONS_data/pop_size/transform_pops.R
## Vincy Huang
## Jun2021

## set-up ####
if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("Rcpp", "dqrng", "fst", "qs", "gamlss", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multicore)

# load data 
population_before_2041 <- read.fst("./ONS_data/pop_size/pop_proj.fst", as.data.table = TRUE)

# train model ####
population_train <- copy(population_before_2041)
population_train[, sex := factor(sex, levels = c("men", "women"), labels = c("men", "women"))]

model_population<- glm(pops ~ age + sex + year + LAD17CD , family = gaussian, data = population_train)
qsave(model_population, "./preparatory_work/marginal_distr_smok_quit_yrs.qs")

# pedict till 2071 ####
population_predict <- CJ(age = unique(population_train$age), sex = unique(population_train$sex), 
              year = 2042: 2071)
