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
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
#}

## af dgn 2 ####
file_qs <- all_file_qs[2]
file_fst <- file_name_fst[2]
model <- qread(file_qs)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
#}

## alcohol 3 ####
file_qs <- all_file_qs[3]
file_fst <- file_name_fst[3]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
#}

## bmi 4 ####
file_qs <- all_file_qs[4]
file_fst <- file_name_fst[4]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
#}

## bmi med 5 ####
file_qs <- all_file_qs[5]
file_fst <- file_name_fst[5]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
#}


## ckd 7 ####
file_qs <- all_file_qs[7]
file_fst <- file_name_fst[7]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(age_int = 15:89) #
newdata[, age := scale(age_int, 52.2, 16.4)]
newdata[, c(paste0("ckd", 0:4)) := data.table(rowCumsums(predict(model, type = "p", newdata = .SD))), .SDcols = trms]

newdata[, age := age_int]
newdata[, c("age_int", "ckd4") := NULL]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
#}







## education imputation 12 ####
file_qs <- all_file_qs[12]
file_fst <- file_name_fst[12]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
#}

## education 13 ####
file_qs <- all_file_qs[13]
file_fst <- file_name_fst[13]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
#}

## ethnicity 14 (can't run) ####
file_qs <- all_file_qs[14]
file_fst <- file_name_fst[14]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, sex = unique(model$data$sex), qimd = unique(model$data$qimd),
              sha = unique(model$data$sha))
newdata[, (levels(model$data$ethnicity)) := data.table(matrixStats::rowCumsums(predict(model, type = "p", newdata = .SD))), .SDcols = trms]

# write .fst file
write_fst(newdata, file_fst , 100L) # export .fst
print(paste0("Table saved ", file_fst))
#}

## ets 15 ####
file_qs <- all_file_qs[15]
file_fst <- file_name_fst[15]
model <- qread(file_qs)
origin_table <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
              sex = unique(origin_table$sex), 
              qimd = unique(origin_table$qimd),
              ethnicity = unique(origin_table$ethnicity), 
              sha = unique(origin_table$sha),
              smok_status = unique(origin_table$smok_status)
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


## famcvd 16 ####
file_qs <- all_file_qs[16]
file_fst <- file_name_fst[16]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
#}

## hdl 18 ####
file_qs <- all_file_qs[18]
file_fst <- file_name_fst[18]
model <- qread(file_qs)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
## income 19 ####
file_qs <- all_file_qs[19]
file_fst <- file_name_fst[19]
model <- qread(file_qs)
origin_table <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(age_int = 15:89, 
              sex = unique(origin_table$sex), 
              qimd = unique(origin_table$qimd),
              ethnicity = unique(origin_table$ethnicity), 
              sha = unique(origin_table$sha),
              education = unique(origin_table$education)
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



## sbp 21 ####
file_qs <- all_file_qs[21]
file_fst <- file_name_fst[21]
model <- qread(file_qs)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
## cess_model 22 ####
file_qs <- all_file_qs[22]
file_fst <- file_name_fst[22]
model <- qread(file_qs)
#origin_table <- read_fst(path = file_fst, as.data.table = TRUE)
# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
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
#}

## cig_curr 23 ####
file_qs <- all_file_qs[23]
file_fst <- file_name_fst[23]
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

## cig_ex 24 ####
file_qs <- all_file_qs[24]
file_fst <- file_name_fst[24]
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

## dur_curr 25####
file_qs <- all_file_qs[25]
file_fst <- file_name_fst[25]
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

## dur_ex 26  ####
file_qs <- all_file_qs[26]
file_fst <- file_name_fst[26]
model <- qread(file_qs)
origin_table <- read_fst(path = file_fst, as.data.table = TRUE)

# generate fst
trms <- all.vars(formula(model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:73, age_int = 15:89, 
              sex = unique(origin_table$sex), 
              qimd = unique(origin_table$qimd),
              ethnicity = unique(origin_table$ethnicity), 
              sha = unique(origin_table$sha),
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
## inci_model 27 ####
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

## quit_years 28 ####
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

## smok_status 29 ####
file_qs <- all_file_qs[29]
file_fst <- file_name_fst[29]
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


