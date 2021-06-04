## generate .fst for age 5 - 89
## Vincy Huang
## 02-Jun-2021

## Set-up ##
# load packages
if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("Rcpp", "dqrng", "fst", "qs", "gamlss", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multicore)

## Process file ##
# find all .qs files
all_file_fit_model <- list.files("./preparatory_work/", pattern = ".R", full.names = TRUE)
all_file_qs <- list.files("./lifecourse_models", pattern = "qs", full.names = TRUE) # find all qs file
file_name_fst <- gsub("model.qs", "table.fst", all_file_qs) # change name to fst file name


# read file
# for (i in 2:length(all_file_qs)) {
# file_qs <- all_file_qs[i]
# file_fst <- file_name_fst[i]

# cess_model ####
file_qs <- all_file_qs[22]
file_fst <- file_name_fst[22]
model <- qread(file_qs)
origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
              sex = unique(origin_table$sex), 
              qimd = unique(origin_table$qimd),
              ethnicity = unique(origin_table$ethnicity), 
              sha = unique(origin_table$sha)
              ) # create completed dataset
newdata[, age := scale(age_int, 43.8, 15.8)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <-
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(model, .SD, data = model$data), .SDcols = trms]) # predict mu
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
#}

## inci_model ####
file_qs <- all_file_qs[27]
file_fst <- file_name_fst[27]
model <- qread(file_qs)
origin_table <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
              sex = unique(origin_table$sex), 
              qimd = unique(origin_table$qimd),
              ethnicity = unique(origin_table$ethnicity), 
              sha = unique(origin_table$sha)
) # create completed dataset
newdata[, age := scale(age_int, 51.7, 17.4)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <-
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(model, .SD, data = model$data), .SDcols = trms]) # predict mu
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))

## 

# Problems 
# scale(age_int, XX) is diff
# CJ has smoke_status: alcohol[3],

# quit_years ####
# read file
file_qs <- all_file_qs[28]
file_fst <- file_name_fst[28]
model <- qread(file_qs)
origin_table <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
              sex = unique(origin_table$sex), 
              qimd = unique(origin_table$qimd),
              ethnicity = unique(origin_table$ethnicity), 
              sha = unique(origin_table$sha)
) # create completed dataset
newdata[, age := scale(age_int, 51.7, 17.4)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <-
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(model, .SD, data = model$data), .SDcols = trms]) # predict mu
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))