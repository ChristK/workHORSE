## generate .fst for age 5 - 89
## Vincy Huang
## 02-Jun-2021

## Set-up 
# load packages
if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("Rcpp", "dqrng", "fst", "qs", "gamlss", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multicore)
smok_cess_model <- qread("./lifecourse_models/smok_cess_model.qs")
dt <- read_fst("./lifecourse_models/smok_cess_table.fst", as.data.table = TRUE)

trms <- all.vars(formula(smok_cess_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 15:89, sex = unique(dt$sex), qimd = unique(dt$qimd),
              ethnicity = unique(dt$ethnicity), sha = unique(dt$sha))
newdata[, age := scale(age_int, 43.8, 15.8)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(smok_cess_model, .SD, data = smok_cess_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]
write_fst(newdata, "./lifecourse_models/smok_cess_table.fst", 100L)

print("Table saved")