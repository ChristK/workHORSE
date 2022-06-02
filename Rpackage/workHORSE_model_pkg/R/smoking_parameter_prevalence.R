# # TODO: 
# check if need this script
# this script is used for a stand alone app now

unlink(dir(path = "Rpackage/IMPACThint_model_package/src/", pattern = ".o$", 
           full.names = TRUE)) # remove trouble files in src folder

library(future) # parallel the process 
library(IMPACThintCalib)
source("./global.R")

# structure for calculating prevalence
# remove HSE calculation & validation 
# use calib parameter for smoking parameter numbers

mc_iter <- 100L # modify to create a different synth pop
interation_index <- mc_iter
design$get_lags(mc_iter)
dt <- SynthPop$new(mc_iter, design) #run until ncc , with problems on 5(frt, veg, smok all, hdl, alcohol)

# parallel process ####
options(future.fork.enable = TRUE)
options(future.rng.onMisuse = "ignore")
plan(multicore, workers = 10) # cores: 10

# population table ####
population_table = CJ(
  age = 16:89, 
  year = 15:69, 
  sex = unique(dt$pop$sex), 
  qimd = unique(dt$pop$qimd)
)

#  function ####
## mortality ####
mortality <- function(dt, mc_, design) {
  ## Tobacco impacted mortality 
  tobacco_model(scenario_nam = "",
                mc = mc_, 
                dt = dt, 
                design_ = design,
                timing = FALSE)
  ## identify actually population life history
  if (all(key(dt) != c("pid", "year"))) {setkey(dt, pid, year)} # set key for C++
  dt[, dead_ref := mortality_all_cause(year,
                                       pid,
                                       prob = revise_mortality_prb,
                                       random_number,
                                       init_year = design$sim_prm$init_year)] # true is dead
  
  dt[, dead_remove := identify_longdead(dead_ref, pid_mrk)] # mark individual history the year after dead
  return(dt)
}

## status prevalence ####


### never
get_ns_prvl <- function(dt) { # dt a synthpop object
  never_smkr_numer <- dt$pop[!(dead_remove) & year >= 3 & smok_status == "1", .(never_smkr_numer = sum(wt)), keyby = .(year, age, sex, qimd)] # Note that not all age/sex/year/qimd combinations are present
  # Never smoker prevalence denominator
  total  <- dt$pop[!(dead_remove) & year >= 3, .(smkr_denom = sum(wt)), keyby = .(year, age, sex, qimd)] # Note that not all age/sex/year/qimd combinations are present
  absorb_dt(total, never_smkr_numer)
  total[, `:=` (
    never_smkr_prvl = never_smkr_numer/smkr_denom,
    never_smkr_numer = NULL,
    smkr_denom = NULL
  )] 
  absorb_dt(population_table, total)
  setnafill(population_table, "const", 0, cols = "never_smkr_prvl")
  return(population_table)
}


get_cs_prvl <- function(dt) { # dt a synthpop object
  current_smoker <- dt$pop[!(dead_remove) & year >= 3 & smok_status == "4", .(current_smkr_numer = sum(wt)), keyby = .(year, age, sex, qimd)] # Note that not all age/sex/year/qimd combinations are present
  # # Never smoker prevalence denominator
  total  <- dt$pop[!(dead_remove) & year >= 3, .(smkr_denom = sum(wt)), keyby = .(year, age, sex, qimd)] # Note that not all age/sex/year/qimd combinations are present

  absorb_dt(total, current_smoker)
  total[, `:=` (
    current_smkr_prvl = current_smkr_numer/smkr_denom,
    current_smkr_numer = NULL,
    smkr_denom = NULL
  )] 
  absorb_dt(population_table, total)
  setnafill(population_table, "const", 0, cols = "current_smkr_prvl")
  return(population_table)
}

get_fs_prvl <- function(dt) { # dt a synthpop object
  former_smoker <- dt$pop[!(dead_remove) & year >= 3 & (smok_status == "2" | smok_status == "3") , .(former_smkr_numer = sum(wt)), keyby = .(year, age, sex, qimd)] # Note that not all age/sex/year/qimd combinations are present
  # # Never smoker prevalence denominator
  total  <- dt$pop[!(dead_remove) & year >= 3, .(smkr_denom = sum(wt)), keyby = .(year, age, sex, qimd)] # Note that not all age/sex/year/qimd combinations are present
  
  absorb_dt(total, former_smoker)
  total[, `:=` (
    former_smkr_prvl = former_smkr_numer/smkr_denom,
    former_smkr_numer = NULL, 
    smkr_denom = NULL
  )] 
  absorb_dt(population_table, total)
  setnafill(population_table, "const", 0, cols = "former_smkr_prvl")
  return(population_table)
}

# look-up tables ####
# Assign smok_incid probabilities
tbl <-
  read_fst("./lifecourse_model/smoke_initiation_table.fst", as.data.table = TRUE)
# calibration
tbl[age < 8L, mu := 0]
lookup_dt(dt$pop, tbl, check_lookup_tbl_validity = TRUE)
setnames(dt$pop, "mu", "prb_smok_incid_orig")

## Assign smok_cessation probabilities
tbl <-
  read_fst("./lifecourse_model/smoke_cessation_table.fst",
           as.data.table = TRUE)
tbl[age < 15L, mu := 0]
lookup_dt(dt$pop, tbl, check_lookup_tbl_validity = TRUE)
setnames(dt$pop, "mu", "prb_smok_cess_orig")

## Handle smok_relapse probabilities
tbl <-
  read_fst("./lifecourse_model/smok_relapse_table.fst",
           as.data.table = TRUE)
tbl <-
  dcast(tbl, sex + qimd ~ smok_quit_yrs, value.var = "pr")
nam <- tbl[, paste0(sex, " ", qimd)]
tbl <-
  as.matrix(tbl[, mget(paste0(1:15))], rownames = nam)

# TODO: add calibrated parameters

# Initiation
message("Multiply parameters")
tbl_relapse_origin = copy(tbl)

## plot baseline prevalence ####
dt$pop[, `:=`(
  prb_smok_incid = prb_smok_incid_orig, # revise: UI input
  smok_status = 1L)] # define a smoking status
dt$pop[, prb_smok_cess := prb_smok_cess_orig] # revise: UI input
dt$pop[(pid_mrk), smok_status := smok_status_ref] # TODO: keep or no?

simsmok(dt$pop, tbl, design$sim_prm$smoking_relapse_limit) # apply initiation, cessation value to synth pop
dt$pop <- mortality(dt$pop, mc_iter, design) # mortality comes here

if ("never_smkr_prvl" %in% names(population_table)){
  population_table[, never_smkr_prvl := NULL]}
if ("former_smkr_prvl" %in% names(population_table)){
  population_table[, former_smkr_prvl := NULL]}
if ("current_smkr_prvl" %in% names(population_table)){
  population_table[, current_smkr_prvl := NULL]}

smoke_prevalence_original = get_ns_prvl(dt)
smoke_prevalence_original = get_cs_prvl(dt) # 
smoke_prevalence_original = get_fs_prvl(dt)# 

# plot 
# aggregate by year
smoke_prevalence_original <- smoke_prevalence_original[, 
                                              .(current_smoker = mean(current_smkr_prvl),
                                                 never_smoker = mean(never_smkr_prvl),
                                                 former_smoker = mean(never_smkr_prvl)),
                                                  keyby = .(year)]

# plot 
smoke_prevalence_original[,{plot(x = year, y = current_smoker, type = "l",  ylim = c(0, 1))}]

## modify the parameter ####
multiply_initiation = scenario_parms$smoking_initiation_sc1
multiply_cessation = scenario_parms$smoking_cessation_sc1
multiply_relapse = scenario_parms$smoking_relapse_sc1

dt$pop[, `:=`(
        prb_smok_incid = prb_smok_incid_orig * multiply_initiation, # revise: UI input
        smok_status = 1L)] # define a smoking status
dt$pop[, prb_smok_cess := prb_smok_cess_orig * multiply_cessation] # revise: UI input
tbl = tbl_relapse_origin * multiply_relapse # revise: UI input

dt$pop[(pid_mrk), smok_status := smok_status_ref] # TODO: keep or no?
simsmok(dt$pop, tbl, design$sim_prm$smoking_relapse_limit) # apply initiation, cessation value to synth pop
dt$pop <- mortality(dt$pop, mc_iter, design) # mortality comes here
      
if ("never_smkr_prvl" %in% names(population_table)){
        population_table[, never_smkr_prvl := NULL]}
if ("former_smkr_prvl" %in% names(population_table)){
        population_table[, former_smkr_prvl := NULL]}
if ("current_smkr_prvl" %in% names(population_table)){
        population_table[, current_smkr_prvl := NULL]}

smoke_prevalence = get_ns_prvl(dt)
smoke_prevalence = get_cs_prvl(dt) # 
smoke_prevalence = get_fs_prvl(dt)# 

# aggregate by year
smoke_prevalence <- smoke_prevalence[, .(current_smoker = mean(current_smkr_prvl),
                    never_smoker = mean(never_smkr_prvl),
                    former_smoker = mean(never_smkr_prvl)),
                    keyby = .(year)]

# plot 
smoke_prevalence[,{
  plot(x = year, y = current_smoker, col = "red", type = "l",  ylim = c(0, 1))}]

