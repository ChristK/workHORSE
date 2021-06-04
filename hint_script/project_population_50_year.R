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

# fit model to pedict till 2071
