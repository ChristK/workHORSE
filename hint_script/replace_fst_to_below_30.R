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

## activey days 1 ####
file_qs <- all_file_qs[1]
file_fst <- file_name_fst[1]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset
newdata[, age := scale(age_int, 49.6, 17)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("pa", 0:7)) := data.table(rowCumsums(predict(model, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, c("age_int", "pa7") := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()
#}

## af dgn 2 ####
file_qs <- all_file_qs[2]
file_fst <- file_name_fst[2]
model <- qread(file_qs)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset
newdata[, age := scale(age_int, 50.9, 17.6)]
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
remove(newdata)
gc()
#}

## alcohol 3 ####
file_qs <- all_file_qs[3]
file_fst <- file_name_fst[3]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha),
              smok_status = unique(model$data$smok_status)
) # create completed dataset
newdata[, age := scale(age_int, 51.3, 17.5)]
newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu") := predictAll(model, .SD, data = model$data), .SDcols = trms])
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()
#}

## bmi 4 ####
file_qs <- all_file_qs[4]
file_fst <- file_name_fst[4]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha),
              smok_status = unique(model$data$smok_status)
) # create completed dataset

newdata[, age := scale(age_int, 50.8, 16.8)]
newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu", "tau") := predictAll(model, .SD, data = model$data), .SDcols = trms])
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()
#}

## bmi med 5 ####
file_qs <- all_file_qs[5]
file_fst <- file_name_fst[5]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha),
              sbp = seq(110, 200, 10)
) # create completed dataset
newdata[, age := scale(age_int, 53.2, 17.2)]
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
remove(newdata)
gc()
#}

## chd duration 6

## ckd 7 ####
file_qs <- all_file_qs[7]
file_fst <- file_name_fst[7]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(age_int = 5:89) #
newdata[, age := scale(age_int, 52.2, 16.4)]
newdata[, c(paste0("ckd", 0:4)) := data.table(matrixStats::rowCumsums(predict(model, type = "p", newdata = .SD))), .SDcols = trms]

newdata[, age := age_int]
newdata[, c("age_int", "ckd4") := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()
#}


## dm_dgn 9####
file_qs <- all_file_qs[9]
file_fst <- file_name_fst[9]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha),
              bmi = unique(round(clamp(model$data$bmi, 18, 50), -1))
) # create completed dataset
newdata[, age := scale(age_int, 62.5, 13.3)]
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
remove(newdata)
gc()

## dm_dur 10 ####
file_qs <- all_file_qs[10]
file_fst <- file_name_fst[10]
model <- qread(file_qs)
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(age_int = 5:89, sex = unique(model$data$sex))
newdata[, age := scale(age_int, 62.5, 13.3)]
newdata <- split(newdata, by = "sex")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma") := predictAll(model, .SD, data = model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
## dm 11 with warning ####
file_qs <- all_file_qs[11]
file_fst <- file_name_fst[11]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              bmi = unique(round(clamp(model$data$bmi, 18, 50)))
) # create completed dataset
newdata[, age := scale(age_int, 51.6, 16.6)]
newdata <- split(newdata, by = "ethnicity")
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
remove(newdata)
gc()

## education imputation 12 - issue with file####
file_qs <- all_file_qs[12]
file_fst <- file_name_fst[12]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset
newdata[, age := scale(age_int, 54.52, 15.28)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("ed", 1:7)) := data.table(rowCumsums(predict(model, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))
newdata <- rbindlist(newdata)
newdata[, age := age_int]
newdata[, c("age_int", "ed7") := NULL]
# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()
#}

## education 13 ####
file_qs <- all_file_qs[13]
file_fst <- file_name_fst[13]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset
newdata[, age := scale(age_int, 54.52, 15.28)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("ed", 1:7)) := data.table(rowCumsums(predict(model, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))
newdata <- rbindlist(newdata)
newdata[, age := age_int]
newdata[, c("age_int", "ed7") := NULL]
# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()
#}

## ethnicity 14 ####
file_qs <- all_file_qs[14]
file_fst <- file_name_fst[14]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age = 5:89, sex = unique(model$data$sex), qimd = unique(model$data$qimd),
              sha = unique(model$data$sha))
newdata[, (levels(model$data$ethnicity)) := data.table(matrixStats::rowCumsums(predict(model, type = "p", newdata = .SD))), .SDcols = trms]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()
#}

## ets 15 ####
file_qs <- all_file_qs[15]
file_fst <- file_name_fst[15]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha),
              smok_status = unique(model$data$smok_status)
) # create completed dataset
newdata[, age := scale(age_int, 50.8, 17.4)]
newdata <- split(newdata, by = "year")
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(model, .SD, data = model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
newdata[, age := age_int]
newdata[, c("age_int") := NULL]

write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()


## famcvd 16 ####
file_qs <- all_file_qs[16]
file_fst <- file_name_fst[16]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset
newdata[, age := scale(age_int, 47.3, 15.6)]
newdata <- split(newdata, by = "sha") # process by year to make small job
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
remove(newdata)
gc()
#}
## frt 17 ####
file_qs <- all_file_qs[17]
file_fst <- file_name_fst[17]
model <- qread(file_qs)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset
newdata[, age := scale(age_int, 50.6, 17.4)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu", "tau") := predictAll(model, .SD, data = model$data), .SDcols = trms]) # predict mu, sigma, nu, tau
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]
# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
remove(newdata)
gc()
## hdl 18 ####
file_qs <- all_file_qs[18]
file_fst <- file_name_fst[18]
model <- qread(file_qs)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset
newdata[, age := scale(age_int, 52, 16.6)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu", "tau") := predictAll(model, .SD, data = model$data), .SDcols = trms]) # predict mu, sigma, nu, tau
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]
# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()

## income 19 ####
file_qs <- all_file_qs[19]
file_fst <- file_name_fst[19]
model <- qread(file_qs)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha),
              education = unique(model$data$education)
) # create completed dataset
newdata[, age := scale(age_int, 49.96, 16.98)]
newdata <- split(newdata, by = "education")
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("inc", 1:5)) := data.table(rowCumsums(predict(model, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))
newdata <- rbindlist(newdata)
newdata[, age := age_int]
newdata[, c("age_int", "inc5") := NULL]

write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))

remove(newdata)
gc()


## predm 20
## sbp 21 ####
file_qs <- all_file_qs[21]
file_fst <- file_name_fst[21]
model <- qread(file_qs)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha),
              smok_status = unique(model$data$smok_status)
) # create completed dataset
newdata[, age := scale(age_int, 52.1, 17.1)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu", "tau") := predictAll(model, .SD, data = model$data), .SDcols = trms]) # predict mu, sigma, nu, tau
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]
# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()

## cess_model 22 ####
file_qs <- all_file_qs[22]
file_fst <- file_name_fst[22]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
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
remove(newdata)
gc()


## cig_curr 23 ####
file_qs <- all_file_qs[23]
file_fst <- file_name_fst[23]
model <- qread(file_qs)
#model$data <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset

newdata[, age := scale(age_int, 45, 15.3)]
newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu") := predictAll(model, .SD, data = model$data), .SDcols = trms])
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()

## cig_ex 24 ####
file_qs <- all_file_qs[24]
file_fst <- file_name_fst[24]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset
newdata[, age := scale(age_int, 57.4, 16.6)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu", "tau") := predictAll(model, .SD, data = model$data), .SDcols = trms]) # predict mu, sigma, nu, tau
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()

## dur_curr 25####
file_qs <- all_file_qs[25]
file_fst <- file_name_fst[25]
model <- qread(file_qs)
#model$data <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset
newdata[, age := scale(age_int, 45, 15.3)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <-
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma") := predictAll(model, .SD, data = model$data), .SDcols = trms]) # predict mu, sigma
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()

## dur_ex 26  ####
file_qs <- all_file_qs[26]
file_fst <- file_name_fst[26]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha),
              smok_status = 2:3
) # create completed dataset
newdata[, age := scale(age_int, 56.2, 17.0)]
newdata[, smok_status := factor(smok_status, 1:4)]
newdata <- split(newdata, by = "year")
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma") := predictAll(model, .SD, data = model$data), .SDcols = trms])
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()

## inci_model 27 ####
file_qs <- all_file_qs[27]
file_fst <- file_name_fst[27]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
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
remove(newdata)
gc()

## 
# Problems 
# scale(age_int, XX) is diff

## quit_years 28 ####
# read file
file_qs <- all_file_qs[28]
file_fst <- file_name_fst[28]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset
newdata[, age := scale(age_int, 56.2, 17.0)]
newdata <- split(newdata, by = "year") # process by year to make small job
newdata <-
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma") := predictAll(model, .SD, data = model$data), .SDcols = trms]) # predict mu, sigma
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()

## smok_status 29 ####
file_qs <- all_file_qs[29]
file_fst <- file_name_fst[29]
model <- qread(file_qs)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset

newdata[, age := scale(age_int, 50.7, 17.4)]
newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu") := predictAll(model, .SD, data = model$data), .SDcols = trms])
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]
# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()

## statin 30 ####
file_qs <- all_file_qs[30]
file_fst <- file_name_fst[30]
model <- qread(file_qs)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha),
              tchol = seq(2, 12, 1)
) # create completed dataset
newdata[, age := scale(age_int, 52.4, 16.6)]
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

remove(newdata)
gc()

## tchol 32 ####
file_qs <- all_file_qs[32]
file_fst <- file_name_fst[32]
model <- qread(file_qs)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset

newdata[, age := scale(age_int, 52, 16.6)]
newdata <- split(newdata, by = "year")
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu", "tau") := predictAll(model, .SD, data = model$data), .SDcols = trms])
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]
# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()

## veg 33 ####
file_qs <- all_file_qs[33]
file_fst <- file_name_fst[33]
model <- qread(file_qs)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 5:89, 
              sex = unique(model$data$sex), 
              qimd = unique(model$data$qimd),
              ethnicity = unique(model$data$ethnicity), 
              sha = unique(model$data$sha)
) # create completed dataset

newdata[, age := scale(age_int, 50.6, 17.4)]
newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu") := predictAll(model, .SD, data = model$data), .SDcols = trms])
newdata <- rbindlist(newdata) # bind all years
#setattr(newdata, "distribution", distr_nam) # ISSUE here: object 'distr_nam' not found
newdata[, age := age_int]
newdata[, age_int := NULL]
# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
remove(newdata)
gc()

# epi disease ####
# #disease_after_20 <- read.csv("./disease_epidemiology/disease_epi.csv")
# disease_after_20 <- fread("./disease_epidemiology/disease_epi.csv", stringsAsFactors = TRUE)
# subset_disease <- disease_after_20[, -c("sex", "age", "qimd", "disease")]
# # subset_disease <- disease_after_20[!names(disease_after_20) %in% c("sex", "age", "qimd", "disease")]
# col_nam <- colnames(subset_disease)
# disease_before_20 <- CJ(age = 5:19,
#                     sex = unique(disease_after_20$sex),
#                     qimd = unique(disease_after_20$qimd),
#                     disease = unique(disease_after_20$disease))
# disease_before_20[, col_nam] <- 0L 
# setcolorder(disease_after_20, c("age", "sex", "qimd", "disease"))
# disease_epi_l <- rbind(disease_before_20, disease_after_20)
# disease_epi_l[is.na(disease_epi_l$in_case_fatality_rates), c("in_case_fatality_rates")] <- 0L
# 
# write_fst(disease_epi_l, "./disease_epidemiology/disease_epi_l.fst", 100L)
# # create a table with row numbers for each mc -- QUESTION: do I need?
# # disease_epi_l[, rn := .I]
# # tt <- disease_epi_l[, .(from = min(rn), to = max(rn)), keyby = mc]
# # write_fst(tt, "./disease_epidemiology/disease_epi_indx.fst", 100L)
# print(paste0("Table saved ", "disease_epi"))
# remove(newdata)
# gc()


disease_epi_after20 <- read.fst("./disease_epidemiology/disease_epi_l.fst", as.data.table = TRUE)
disease_epi_20 <- disease_epi_after20[age == 20L]
disease_epi_before_20 <- data.frame()
for (i in 5:19) {
  subset <- disease_epi_20
  subset[, age := i]
  disease_epi_before_20 <- rbind(disease_epi_before_20, subset)
}
disease_epi <- rbind(disease_epi_before_20, disease_epi_after20)
summary(disease_epi)
write_fst(disease_epi, "./lifecourse_models/disease_epi_l.fst", 100L)
print(paste0("Table saved ", "disease_epi"))
remove(disease_epi)
gc()


# cst_prvl ####
cst_prvl_after_20 <- fread("./lifecourse_models/corticosteroids_prvl.csv", stringsAsFactors = TRUE)

cst_prvl_after_20 <- as.data.table(cst_prvl_after_20)
cst_prvl_20 <- cst_prvl_after_20[age == 20L]
cst_prvl_before_20 <- data.frame()
for (i in 5:19) {
  subset <- cst_prvl_20
  subset[, age := i]
  cst_prvl_before_20 <- rbind(cst_prvl_before_20, subset)
}
cst_prvl <- rbind(cst_prvl_before_20, cst_prvl_after_20)
write.csv(cst_prvl, "./lifecourse_models/corticosteroids_prvl.csv", row.names = FALSE)
print(paste0("Table saved ", "cst_prvl"))
remove(cst_prvl)
gc()

# chd duration ####
chd_duration_after20 <- read.fst("./lifecourse_models/chd_duration_table.fst", as.data.table = TRUE)
chd_duration_20 <- chd_duration_after20[age == 20L]
chd_duration_before_20 <- data.frame()
for (i in 5:19) {
  subset <- chd_duration_20
  subset[, age := i]
  chd_duration_before_20 <- rbind(chd_duration_before_20, subset)
}
chd_duration <- rbind(chd_duration_before_20, chd_duration_after20)
summary(chd_duration)
write_fst(chd_duration, "./lifecourse_models/chd_duration_table.fst", 100L)
print(paste0("Table saved ", "chd"))
remove(chd_duration)
gc()

# stroke duration ####
#stroke_duration_model <- qread("./lifecourse_models/stroke_duration_model.qs")
stroke_duration_after20 <- read.fst("./lifecourse_models/stroke_duration_table.fst", as.data.table = TRUE)
stroke_duration_20 <- stroke_duration_after20[age == 20L]
stroke_duration_before_20 <- data.frame()
for (i in 5:19) {
  subset <- stroke_duration_20
  subset[, age := i]
  stroke_duration_before_20 <- rbind(stroke_duration_before_20, subset)
}
stroke_duration <- rbind(stroke_duration_before_20, stroke_duration_after20)
write_fst(stroke_duration, "./lifecourse_models/stroke_duration_table.fst", 100L)
print(paste0("Table saved ", "cst_prvl"))
remove(newdata)
gc()
# copd duration ####
copd_duration_after20 <- read.fst("./lifecourse_models/copd_duration_table.fst", as.data.table = TRUE)
copd_duration_20 <- copd_duration_after20[age == 20L]
copd_duration_before_20 <- data.frame()
for (i in 5:19) {
  subset <- copd_duration_20
  subset[, age := i]
  copd_duration_before_20 <- rbind(copd_duration_before_20, subset)
}
copd_duration <- rbind(copd_duration_before_20, copd_duration_after20)
summary(copd_duration)
write_fst(copd_duration, "./lifecourse_models/copd_duration_table.fst", 100L)
print(paste0("Table saved ", "copd"))
remove(copd_duration)
gc()
