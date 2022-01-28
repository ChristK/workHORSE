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
# smoke_prevalence_loop = data.table()
# smoke_prevalence_relapse = data.table()
# smoke_prevalence_cessation = data.table()
# smoke_prevalence = data.table()


# plot of the before any multiplication and return back to UI


# for (multiply_initiation in seq(0, 1, 0.5)) {
#   dt$pop[, `:=`(
#     prb_smok_incid = prb_smok_incid_orig * multiply_initiation,
#     smok_status = 1L # define a smoking status
#   )]
#   for (multiply_cessation in seq(0, 1, 0.5)) { 
#     dt$pop[, prb_smok_cess := prb_smok_cess_orig * multiply_cessation]
#   # relapse
#     for (multiply_relapse in seq(0,1,0.5)) {
#       print(paste0("initiation is ",multiply_initiation, "and cessation is", multiply_cessation, "and relapse is",multiply_relapse , collapse = " "))
#      
#       

# input UI parameter 


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
smoke_prevalence_original = get_cs_prvl(dt) # add row with cs and result
smoke_prevalence_original = get_fs_prvl(dt)# add row with cs and result

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
      # calculate prevalence
      # smoke_prevalence_loop[,  `:=`(initiation = multiply_initiation, 
      #                               cessation = multiply_cessation,
      #                               relapse = multiply_relapse)]
smoke_prevalence = get_ns_prvl(dt)
smoke_prevalence = get_cs_prvl(dt) # add row with cs and result
smoke_prevalence = get_fs_prvl(dt)# add row with cs and result

# aggregate by year
smoke_prevalence <- smoke_prevalence[, .(current_smoker = mean(current_smkr_prvl),
                    never_smoker = mean(never_smkr_prvl),
                    former_smoker = mean(never_smkr_prvl)),
                    keyby = .(year)]

# plot 
smoke_prevalence[,{
  plot(x = year, y = current_smoker, col = "red", type = "l",  ylim = c(0, 1))}]

      #smoke_prevalence_relapse = rbind(smoke_prevalence_loop, smoke_prevalence_relapse)# add the matrix to the next one
   
#        }
#   smoke_prevalence_cessation = rbind(smoke_prevalence_relapse, smoke_prevalence_cessation)# add the matrix to the next one
#   }
# smoke_prevalence = rbind(smoke_prevalence, smoke_prevalence_cessation)# add the matrix to the next one
# }



# reformat data to wide format (year tag along each age, ini, cess, relap combination)
never_smoker_prevalence = smoke_prevalence[, -c('current_smkr_prvl', 'former_smkr_prvl')]
current_smoker_prevalence = smoke_prevalence[, -c('never_smkr_prvl', 'former_smkr_prvl')]
former_smoker_prevalence = smoke_prevalence[, -c('current_smkr_prvl', 'never_smkr_prvl')]

write.fst(never_smoker_prevalence, path = "./never_smoker_prevalence.fst" , compress = 90)
write.fst(current_smoker_prevalence, path = "./current_smoker_prevalence.fst", compress = 90)
write.fst(former_smoker_prevalence, path = "./former_smoker_prevalence.fst", compress = 90)

write.csv(never_smoker_prevalence[year < 35], file = "./never_smoker_prevalence.csv" )
write.csv(current_smoker_prevalence[year < 35], file = "./current_smoker_prevalence.csv")
write.csv(former_smoker_prevalence[year < 35], file = "./former_smoker_prevalence.csv")

never_smoker_prevalence_wide = dcast(never_smoker_prevalence, age + sex + qimd + initiation + cessation + relapse ~ year, 
                                     fun = sum, 
                                     value.var = c("never_smkr_prvl")) # TODO: result is not correct
current_smoker_prevalence_wide = dcast(current_smoker_prevalence, age + sex + qimd + initiation + cessation + relapse ~ year, 
                                     fun = sum, 
                                     value.var = c("current_smkr_prvl")) # TODO: result is not correct

former_smoker_prevalence_wide = dcast(former_smoker_prevalence, age + sex + qimd + initiation + cessation + relapse ~ year, 
                                     fun = sum, 
                                     value.var = c("former_smkr_prvl")) # TODO: result is not correct

remove(smoke_prevalence)

# plot of the changed prevalence and return back to UI


# # plot_path <- paste0("./HSE_data/plot/prevalence_after_calib/", interation_index, "_never_smoker.jpg", collapse = "")
# # jpeg(plot_path)
# par(mfrow=c(1,2))
# par(mar=c(2,2,1,1))
# 
# na.omit(truth_tbl_gamlss)[, .(obs = mean(never_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(never_smkr_prvl_obs),
#                               gamlss = mean(never_smkr_prvl_gamlss, na.rm = TRUE),
#                               mdl = mean(never_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(year)][, {
#                                 plot(x = year, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = year, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = year, y = gamlss, lty = 3, col = "black")
#                                 lines(x = year, y = mdl, lty = 2, col = "blue")}]
# na.omit(truth_tbl_gamlss)[, .(obs = mean(never_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(never_smkr_prvl_obs),
#                               gamlss = mean(never_smkr_prvl_gamlss, na.rm = TRUE),
#                               mdl = mean(never_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(age)][, {
#                                 plot(x = age, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = age, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = age, y = gamlss, lty = 3, col = "black")
#                                 lines(x = age, y = mdl, lty = 2, col = "blue")}]
# 
# par(mfrow=c(1,1))
# par(mar=c(2,2,1,1))
# na.omit(truth_tbl_gamlss)[, .(difference = (mean(never_smkr_prvl_mdl, na.rm = TRUE) - mean(never_smkr_prvl_obs, na.rm = TRUE))),
#                           keyby = .(age)][, {
#                             plot(x = age, y = difference,
#                                  type = "n",  col = "blue",
#                                  #ylim = c(-1, 1),
#                                  main = "model - observed",
#                                  ylab = "difference %")
#                             axis(side = 1, seq(min(truth_tbl_gamlss$age), max(truth_tbl_gamlss$age), 10))
#                             axis(side = 2, seq(-1, 1 , 0.1))
#                             lines(x = age, y = difference, lty = 2, col = "blue")
#                           }]
# 
# 
# #dev.off()
# 
# # Relapse
# message("Calibrate relapse - below 65")
# tbl_temp <- copy(tbl)
# tbl_above_65_temp <- copy(tbl)
# tbl_above_65_temp <- tbl_above_65_temp * 1.7
# 
# sub1_tbl_calibrated %<-% calibrate_relapse_below_65(tbl_below_65_temp = tbl_temp, tbl_above_65 = tbl_above_65_temp, tbl_relapse_origin = tbl,
#                                                     row_index = 1, age_start = 16, age_end = 65)
# sub2_tbl_calibrated %<-% calibrate_relapse_below_65(tbl_below_65_temp = tbl_temp, tbl_above_65 = tbl_above_65_temp, tbl_relapse_origin = tbl,
#                                                     row_index = 2, age_start = 16, age_end = 65)
# sub3_tbl_calibrated %<-% calibrate_relapse_below_65(tbl_below_65_temp = tbl_temp, tbl_above_65 = tbl_above_65_temp, tbl_relapse_origin = tbl,
#                                                     row_index = 3, age_start = 16, age_end = 65)
# sub4_tbl_calibrated %<-% calibrate_relapse_below_65(tbl_below_65_temp = tbl_temp, tbl_above_65 = tbl_above_65_temp, tbl_relapse_origin = tbl,
#                                                     row_index = 4, age_start = 16, age_end = 65)
# sub5_tbl_calibrated %<-% calibrate_relapse_below_65(tbl_below_65_temp = tbl_temp, tbl_above_65 = tbl_above_65_temp, tbl_relapse_origin = tbl,
#                                                     row_index = 5, age_start = 16, age_end = 65)
# sub6_tbl_calibrated %<-% calibrate_relapse_below_65(tbl_below_65_temp = tbl_temp, tbl_above_65 = tbl_above_65_temp, tbl_relapse_origin = tbl,
#                                                     row_index = 6, age_start = 16, age_end = 65)
# sub7_tbl_calibrated %<-% calibrate_relapse_below_65(tbl_below_65_temp = tbl_temp, tbl_above_65 = tbl_above_65_temp, tbl_relapse_origin = tbl,
#                                                     row_index = 7, age_start = 16, age_end = 65)
# sub8_tbl_calibrated %<-% calibrate_relapse_below_65(tbl_below_65_temp = tbl_temp, tbl_above_65 = tbl_above_65_temp, tbl_relapse_origin = tbl,
#                                                     row_index = 8, age_start = 16, age_end = 65)
# sub9_tbl_calibrated %<-% calibrate_relapse_below_65(tbl_below_65_temp = tbl_temp, tbl_above_65 = tbl_above_65_temp, tbl_relapse_origin = tbl,
#                                                     row_index = 9, age_start = 16, age_end = 65)
# sub10_tbl_calibrated %<-% calibrate_relapse_below_65(tbl_below_65_temp = tbl_temp, tbl_above_65 = tbl_above_65_temp, tbl_relapse_origin = tbl,
#                                                      row_index = 10, age_start = 16, age_end = 65)
# tbl_calibrated_below_65 <- rbind(sub1_tbl_calibrated, sub2_tbl_calibrated,
#                                  sub3_tbl_calibrated, sub4_tbl_calibrated,
#                                  sub5_tbl_calibrated, sub6_tbl_calibrated,
#                                  sub7_tbl_calibrated, sub8_tbl_calibrated,
#                                  sub9_tbl_calibrated, sub10_tbl_calibrated) 
# remove(sub1_tbl_calibrated, sub2_tbl_calibrated,
#        sub3_tbl_calibrated, sub4_tbl_calibrated,
#        sub5_tbl_calibrated, sub6_tbl_calibrated,
#        sub7_tbl_calibrated, sub8_tbl_calibrated,
#        sub9_tbl_calibrated, sub10_tbl_calibrated)
# 
# tbl_calibrated_below_65 <- as.matrix(tbl_calibrated_below_65)
# rownames(tbl_calibrated_below_65) <- nam
# colnames(tbl_calibrated_below_65) <- seq(1:15)
# relapse_modifier_below_65 <- as.data.table(tbl_calibrated_below_65/tbl)   # relapse modifer tbl
# path <- paste0("lifecourse_model/relapse_modifier/below_65", mc_iter, ".fst", collapse = "")
# message("Write relapse below_65")
# # write.fst(relapse_modifier_below_65, path, compress = 100L)
# 
# ## Remove for loop: Apply modified relapse to population
# simsmok_relapse_calibrate(dt$pop, tbl_calibrated_below_65, tbl_above_65_temp, design$sim_prm$smoking_relapse_limit) # apply initiation, cessation value to synth pop
# dt$pop <- mortality(dt$pop, mc_iter, design)
# model_tbl <- get_fs_prvl(dt) # get never_smoker prevaelence
# if ("former_smkr_prvl_mdl" %in% names(truth_tbl_gamlss)){
#   truth_tbl_gamlss[, former_smkr_prvl_mdl := NULL]}
# absorb_dt(truth_tbl_gamlss, model_tbl)
# truth_tbl_gamlss[is.na(truth_tbl_gamlss$former_smkr_prvl_mdl), former_smkr_prvl_mdl := former_smoker_gamlss]
# 
# #Plot current smoker prevalence by year and by age
# par(mfrow=c(1,2))
# par(mar=c(2,2,1,1))
# 
# na.omit(truth_tbl_gamlss)[, .(obs = mean(former_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(former_smkr_prvl_obs),
#                               gamlss = mean(former_smoker_gamlss, na.rm = TRUE),
#                               mdl = mean(former_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(year)][, {
#                                 plot(x = year, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = year, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = year, y = gamlss, lty = 3, col = "black")
#                                 lines(x = year, y = mdl, lty = 2, col = "blue")}]
# na.omit(truth_tbl_gamlss)[, .(obs = mean(former_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(former_smkr_prvl_obs),
#                               gamlss = mean(former_smoker_gamlss, na.rm = TRUE),
#                               mdl = mean(former_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(age)][, {
#                                 plot(x = age, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = age, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = age, y = gamlss, lty = 3, col = "black")
#                                 lines(x = age, y = mdl, lty = 2, col = "blue")}]
# 
# par(mfrow=c(1,1))
# na.omit(truth_tbl_gamlss)[, .(difference_former = (mean(former_smkr_prvl_mdl, na.rm = TRUE) - mean(former_smoker_gamlss, na.rm = TRUE))),
#                           keyby = .(age)][, {
#                             plot(x = age, y = difference_former,
#                                  type = "n",  col = "blue",
#                                  # ylim = c(-20, 50),
#                                  main = "model - observed",
#                                  ylab = "difference")
#                             axis(side = 1, seq(min(truth_tbl_gamlss$age), max(truth_tbl_gamlss$age), 10))
#                             axis(side = 2, seq(-1, 1 ))
#                             lines(x = age, y = difference_former, lty = 2, col = "blue")
#                           }]
# 
# 
# # Relapse above 65 ####
# message("Calibrate relapse - above 65")
# tbl_above_65_origin = copy(tbl)
# tbl_above_65_origin = tbl_above_65_origin * 1.7
# sub1_tbl_calibrated %<-% calibrate_relapse_above_65(tbl_calibrated_below_65, tbl_above_65_temp, tbl_above_65_origin = tbl_above_65_origin,
#                                                     row_index = 1, age_start = 66, age_end = max(dt$pop$age))
# sub2_tbl_calibrated %<-% calibrate_relapse_above_65(tbl_calibrated_below_65, tbl_above_65_temp, tbl_above_65_origin = tbl_above_65_origin,
#                                                     row_index = 2, age_start = 66, age_end = max(dt$pop$age))
# sub3_tbl_calibrated %<-% calibrate_relapse_above_65(tbl_calibrated_below_65, tbl_above_65_temp, tbl_above_65_origin = tbl_above_65_origin,
#                                                     row_index = 3, age_start = 66, age_end = max(dt$pop$age))
# sub4_tbl_calibrated %<-% calibrate_relapse_above_65(tbl_calibrated_below_65, tbl_above_65_temp, tbl_above_65_origin = tbl_above_65_origin,
#                                                     row_index = 4, age_start = 66, age_end = max(dt$pop$age))
# sub5_tbl_calibrated %<-% calibrate_relapse_above_65(tbl_calibrated_below_65, tbl_above_65_temp, tbl_above_65_origin = tbl_above_65_origin,
#                                                     row_index = 5, age_start = 66, age_end = max(dt$pop$age))
# sub6_tbl_calibrated %<-% calibrate_relapse_above_65(tbl_calibrated_below_65, tbl_above_65_temp, tbl_above_65_origin = tbl_above_65_origin,
#                                                     row_index = 6, age_start = 66, age_end = max(dt$pop$age))
# sub7_tbl_calibrated %<-% calibrate_relapse_above_65(tbl_calibrated_below_65, tbl_above_65_temp, tbl_above_65_origin = tbl_above_65_origin,
#                                                     row_index = 7, age_start = 66, age_end = max(dt$pop$age))
# sub8_tbl_calibrated %<-% calibrate_relapse_above_65(tbl_calibrated_below_65, tbl_above_65_temp, tbl_above_65_origin = tbl_above_65_origin,
#                                                     row_index = 8, age_start = 66, age_end = max(dt$pop$age))
# sub9_tbl_calibrated %<-% calibrate_relapse_above_65(tbl_calibrated_below_65, tbl_above_65_temp, tbl_above_65_origin = tbl_above_65_origin,
#                                                     row_index = 9, age_start = 66, age_end = max(dt$pop$age))
# sub10_tbl_calibrated %<-% calibrate_relapse_above_65(tbl_calibrated_below_65, tbl_above_65_temp, tbl_above_65_origin = tbl_above_65_origin,
#                                                      row_index = 10, age_start = 66, age_end = max(dt$pop$age))
# tbl_calibrated_above_65 <- rbind(sub1_tbl_calibrated, sub2_tbl_calibrated,
#                                  sub3_tbl_calibrated, sub4_tbl_calibrated,
#                                  sub5_tbl_calibrated, sub6_tbl_calibrated,
#                                  sub7_tbl_calibrated, sub8_tbl_calibrated,
#                                  sub9_tbl_calibrated, sub10_tbl_calibrated)
# remove(sub1_tbl_calibrated, sub2_tbl_calibrated,
#        sub3_tbl_calibrated, sub4_tbl_calibrated,
#        sub5_tbl_calibrated, sub6_tbl_calibrated,
#        sub7_tbl_calibrated, sub8_tbl_calibrated,
#        sub9_tbl_calibrated, sub10_tbl_calibrated)
# 
# tbl_calibrated_above_65 <- as.matrix(tbl_calibrated_above_65)
# rownames(tbl_calibrated_above_65) <- nam
# colnames(tbl_calibrated_above_65) <- seq(1:15)
# relapse_modifier_above_65 <- as.data.table(tbl_calibrated_above_65/tbl_above_65_origin)   # relapse modifer tbl
# path <- paste0("lifecourse_model/relapse_modifier/", mc_iter, ".fst", collapse = "")
# message("Write relapse")
# # write.fst(relapse_modifier_above_65, path, compress = 100L)
# 
# ## Remove for loop: Apply modified relapse to population
# simsmok_relapse_calibrate(dt$pop, tbl_calibrated_below_65, tbl_calibrated_above_65, design$sim_prm$smoking_relapse_limit) # apply initiation, cessation value to synth pop
# dt$pop <- mortality(dt$pop, mc_iter, design)
# model_tbl <- get_fs_prvl(dt) # get never_smoker prevaelence
# if ("former_smkr_prvl_mdl" %in% names(truth_tbl_gamlss)){
#   truth_tbl_gamlss[, former_smkr_prvl_mdl := NULL]}
# absorb_dt(truth_tbl_gamlss, model_tbl)
# truth_tbl_gamlss[is.na(truth_tbl_gamlss$former_smkr_prvl_mdl), former_smkr_prvl_mdl := former_smoker_gamlss]
# 
# #Plot current smoker prevalence by year and by age
# par(mfrow=c(1,2))
# par(mar=c(2,2,1,1))
# 
# na.omit(truth_tbl_gamlss)[, .(obs = mean(former_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(former_smkr_prvl_obs),
#                               gamlss = mean(former_smoker_gamlss, na.rm = TRUE),
#                               mdl = mean(former_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(year)][, {
#                                 plot(x = year, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = year, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = year, y = gamlss, lty = 3, col = "black")
#                                 lines(x = year, y = mdl, lty = 2, col = "blue")}]
# na.omit(truth_tbl_gamlss)[, .(obs = mean(former_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(former_smkr_prvl_obs),
#                               gamlss = mean(former_smoker_gamlss, na.rm = TRUE),
#                               mdl = mean(former_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(age)][, {
#                                 plot(x = age, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = age, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = age, y = gamlss, lty = 3, col = "black")
#                                 lines(x = age, y = mdl, lty = 2, col = "blue")}]
# 
# # Cessation
# message("Calibrate Cessation")
# age_seq_cessation <- c(16:91)
# 
# cessation_1m_dt %<-% calibrate_cessation(dt, tbl_calibrated_below_65, tbl_calibrated_above_65, qimd_input = "1 most deprived",  sex_input = "men", mc_iter, design)
# cessation_2m_dt %<-% calibrate_cessation(dt, tbl_calibrated_below_65, tbl_calibrated_above_65, qimd_input = "2",  sex_input = "men", mc_iter, design)
# cessation_3m_dt %<-% calibrate_cessation(dt, tbl_calibrated_below_65, tbl_calibrated_above_65, qimd_input = "3",  sex_input = "men", mc_iter, design)
# cessation_4m_dt %<-% calibrate_cessation(dt, tbl_calibrated_below_65, tbl_calibrated_above_65, qimd_input = "4",  sex_input = "men", mc_iter, design)
# cessation_5m_dt %<-% calibrate_cessation(dt, tbl_calibrated_below_65, tbl_calibrated_above_65, qimd_input = "5 least deprived",  sex_input = "men", mc_iter, design)
# cessation_1f_dt %<-% calibrate_cessation(dt, tbl_calibrated_below_65, tbl_calibrated_above_65, qimd_input = "1 most deprived",  sex_input = "women", mc_iter, design)
# cessation_2f_dt %<-% calibrate_cessation(dt, tbl_calibrated_below_65, tbl_calibrated_above_65, qimd_input = "2",  sex_input = "women", mc_iter, design)
# cessation_3f_dt %<-% calibrate_cessation(dt, tbl_calibrated_below_65, tbl_calibrated_above_65, qimd_input = "3",  sex_input = "women", mc_iter, design)
# cessation_4f_dt %<-% calibrate_cessation(dt, tbl_calibrated_below_65, tbl_calibrated_above_65, qimd_input = "4",  sex_input = "women", mc_iter, design)
# cessation_5f_dt %<-% calibrate_cessation(dt, tbl_calibrated_below_65, tbl_calibrated_above_65, qimd_input = "5 least deprived",  sex_input = "women", mc_iter, design)
# dt$pop <- rbind(cessation_1m_dt, cessation_2m_dt, cessation_3m_dt, cessation_4m_dt, cessation_5m_dt,
#                 cessation_1f_dt, cessation_2f_dt, cessation_3f_dt, cessation_4f_dt, cessation_5f_dt )
# remove(cessation_1m_dt, cessation_2m_dt, cessation_3m_dt, cessation_4m_dt, cessation_5m_dt,
#        cessation_1f_dt, cessation_2f_dt, cessation_3f_dt, cessation_4f_dt, cessation_5f_dt)
# 
# cessation_modifier <- dt$pop[age %in% age_seq_cessation, .(age, sex, qimd, smok_status, prb_smok_cess_orig, prb_smok_cess)]
# cessation_modifier <- cessation_modifier[, modifier := prb_smok_cess - prb_smok_cess_orig]
# cessation_modifier[, c("prb_smok_cess_orig", "prb_smok_cess") := NULL]
# cessation_modifier <- unique(cessation_modifier)
# path <- paste0("lifecourse_model/cessation_modifier/", mc_iter, ".fst", collapse = "")
# message("Write cessation")
# # write.fst(cessation_modifier, path, compress = 100L)
# 
# 
# ## Remove for loop: Apply modified relapse to population
# simsmok_relapse_calibrate(dt$pop, tbl_calibrated_below_65, tbl_calibrated_above_65, design$sim_prm$smoking_relapse_limit) # apply initiation, cessation value to synth pop
# dt$pop <- mortality(dt$pop, mc_iter, design)
# model_tbl <- get_fs_prvl(dt) # get never_smoker prevaelence
# if ("former_smkr_prvl_mdl" %in% names(truth_tbl_gamlss)){
#   truth_tbl_gamlss[, former_smkr_prvl_mdl := NULL]}
# absorb_dt(truth_tbl_gamlss, model_tbl)
# truth_tbl_gamlss[is.na(truth_tbl_gamlss$former_smkr_prvl_mdl), former_smkr_prvl_mdl := former_smoker_gamlss]
# 
# model_tbl <- get_cs_prvl(dt) # get never_smoker prevaelence
# if ("current_smkr_prvl_mdl" %in% names(truth_tbl_gamlss)){
#   truth_tbl_gamlss[, current_smkr_prvl_mdl := NULL]}
# absorb_dt(truth_tbl_gamlss, model_tbl)
# truth_tbl_gamlss[is.na(truth_tbl_gamlss$current_smkr_prvl_mdl), current_smkr_prvl_mdl := current_smkr_prvl_gamlss]
# 
# ### Plot former smoker prevalence by year and by age
# # former_smoker_two_plots_calib_by_age
# par(mfrow=c(1,2))
# par(mar=c(2,2,1,1))
# 
# na.omit(truth_tbl_gamlss)[, .(obs = mean(former_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(former_smkr_prvl_obs),
#                               gamlss = mean(former_smoker_gamlss, na.rm = TRUE),
#                               mdl = mean(former_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(year)][, {
#                                 plot(x = year, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = year, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = year, y = gamlss, lty = 3, col = "black")
#                                 lines(x = year, y = mdl, lty = 2, col = "blue")}]
# na.omit(truth_tbl_gamlss)[, .(obs = mean(former_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(former_smkr_prvl_obs),
#                               gamlss = mean(former_smoker_gamlss, na.rm = TRUE),
#                               mdl = mean(former_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(age)][, {
#                                 plot(x = age, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = age, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = age, y = gamlss, lty = 3, col = "black")
#                                 lines(x = age, y = mdl, lty = 2, col = "blue")}]
# 
# # plot current smoker prevalence by year and by age
# # current_smoker_two_plots_calib_by_age
# 
# 
# na.omit(truth_tbl_gamlss)[, .(obs = mean(current_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(current_smkr_prvl_obs),
#                               gamlss = mean(current_smkr_prvl_gamlss, na.rm = TRUE),
#                               mdl = mean(current_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(year)][, {
#                                 plot(x = year, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = year, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = year, y = gamlss, lty = 3, col = "black")
#                                 lines(x = year, y = mdl, lty = 2, col = "blue")}]
# na.omit(truth_tbl_gamlss)[, .(obs = mean(current_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(current_smkr_prvl_obs),
#                               gamlss = mean(current_smkr_prvl_gamlss, na.rm = TRUE),
#                               mdl = mean(current_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(age)][, {
#                                 plot(x = age, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = age, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = age, y = gamlss, lty = 3, col = "black")
#                                 lines(x = age, y = mdl, lty = 2, col = "blue")}]
# 
# 
# # Recalibe Relapse AGAIN ####
# message("Re-calibrate relapse - below 65")
# tbl_relapse_below_65_temp <- copy(tbl_calibrated_below_65)
# tbl_relapse_above_65_temp <- copy(tbl_calibrated_above_65)
# 
# sub1_tbl_calibrated %<-% calibrate_recalibrate_relapse_below_65(tbl_relapse_below_65_temp, tbl_relapse_above_65_temp, tbl_relapse_origin = tbl_calibrated_below_65,
#                                                                 row_index = 1, age_start = 16, age_end = 65)
# sub2_tbl_calibrated %<-% calibrate_recalibrate_relapse_below_65(tbl_relapse_below_65_temp, tbl_relapse_above_65_temp, tbl_relapse_origin = tbl_calibrated_below_65,
#                                                                 row_index = 2, age_start = 16, age_end = 65)
# sub3_tbl_calibrated %<-% calibrate_recalibrate_relapse_below_65(tbl_relapse_below_65_temp, tbl_relapse_above_65_temp, tbl_relapse_origin = tbl_calibrated_below_65,
#                                                                 row_index = 3, age_start = 16, age_end = 65)
# sub4_tbl_calibrated %<-% calibrate_recalibrate_relapse_below_65(tbl_relapse_below_65_temp, tbl_relapse_above_65_temp, tbl_relapse_origin = tbl_calibrated_below_65,
#                                                                 row_index = 4, age_start = 16, age_end = 65)
# sub5_tbl_calibrated %<-% calibrate_recalibrate_relapse_below_65(tbl_relapse_below_65_temp, tbl_relapse_above_65_temp, tbl_relapse_origin = tbl_calibrated_below_65,
#                                                                 row_index = 5, age_start = 16, age_end = 65)
# sub6_tbl_calibrated %<-% calibrate_recalibrate_relapse_below_65(tbl_relapse_below_65_temp, tbl_relapse_above_65_temp, tbl_relapse_origin = tbl_calibrated_below_65,
#                                                                 row_index = 6, age_start = 16, age_end = 65)
# sub7_tbl_calibrated %<-% calibrate_recalibrate_relapse_below_65(tbl_relapse_below_65_temp, tbl_relapse_above_65_temp, tbl_relapse_origin = tbl_calibrated_below_65,
#                                                                 row_index = 7, age_start = 16, age_end = 65)
# sub8_tbl_calibrated %<-% calibrate_recalibrate_relapse_below_65(tbl_relapse_below_65_temp, tbl_relapse_above_65_temp, tbl_relapse_origin = tbl_calibrated_below_65,
#                                                                 row_index = 8, age_start = 16, age_end = 65)
# sub9_tbl_calibrated %<-% calibrate_recalibrate_relapse_below_65(tbl_relapse_below_65_temp, tbl_relapse_above_65_temp, tbl_relapse_origin = tbl_calibrated_below_65,
#                                                                 row_index = 9, age_start = 16, age_end = 65)
# sub10_tbl_calibrated %<-% calibrate_recalibrate_relapse_below_65(tbl_relapse_below_65_temp, tbl_relapse_above_65_temp, tbl_relapse_origin = tbl_calibrated_below_65,
#                                                                  row_index = 10, age_start = 16, age_end = 65)
# tbl_recalibrate_relapse_below_65 <- rbind(sub1_tbl_calibrated, sub2_tbl_calibrated,
#                                           sub3_tbl_calibrated, sub4_tbl_calibrated,
#                                           sub5_tbl_calibrated, sub6_tbl_calibrated,
#                                           sub7_tbl_calibrated, sub8_tbl_calibrated,
#                                           sub9_tbl_calibrated, sub10_tbl_calibrated) 
# remove(sub1_tbl_calibrated, sub2_tbl_calibrated,
#        sub3_tbl_calibrated, sub4_tbl_calibrated,
#        sub5_tbl_calibrated, sub6_tbl_calibrated,
#        sub7_tbl_calibrated, sub8_tbl_calibrated,
#        sub9_tbl_calibrated, sub10_tbl_calibrated)
# 
# tbl_recalibrate_relapse_below_65 <- as.matrix(tbl_recalibrate_relapse_below_65)
# rownames(tbl_recalibrate_relapse_below_65) <- nam
# colnames(tbl_recalibrate_relapse_below_65) <- seq(1:15)
# recalibrate_relapse_modifier_below_65 <- as.data.table(tbl_recalibrate_relapse_below_65/tbl_calibrated_below_65)   # relapse modifer tbl
# setnafill(data.table(recalibrate_relapse_modifier_below_65), type = 'const', fill = 0, nan = NaN)
# 
# path <- paste0("lifecourse_model/recalibrate_relapse_modifier/below_65", mc_iter, ".fst", collapse = "")
# message("Write relapse below_65")
# # write.fst(recalibrate_relapse_modifier_below_65, path, compress = 100L)
# 
# ## Remove for loop: Apply modified relapse to population
# simsmok_relapse_calibrate(dt$pop, tbl_recalibrate_relapse_below_65, tbl_relapse_above_65_temp, design$sim_prm$smoking_relapse_limit) # apply initiation, cessation value to synth pop
# dt$pop <- mortality(dt$pop, mc_iter, design)
# model_tbl <- get_fs_prvl(dt) # get never_smoker prevaelence
# if ("former_smkr_prvl_mdl" %in% names(truth_tbl_gamlss)){
#   truth_tbl_gamlss[, former_smkr_prvl_mdl := NULL]}
# absorb_dt(truth_tbl_gamlss, model_tbl)
# truth_tbl_gamlss[is.na(truth_tbl_gamlss$former_smkr_prvl_mdl), former_smkr_prvl_mdl := former_smoker_gamlss]
# 
# 
# #Plot current smoker prevalence by year and by age
# par(mfrow=c(1,2))
# par(mar=c(2,2,1,1))
# 
# na.omit(truth_tbl_gamlss)[, .(obs = mean(former_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(former_smkr_prvl_obs),
#                               gamlss = mean(former_smoker_gamlss, na.rm = TRUE),
#                               mdl = mean(former_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(year)][, {
#                                 plot(x = year, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = year, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = year, y = gamlss, lty = 3, col = "black")
#                                 lines(x = year, y = mdl, lty = 2, col = "blue")}]
# na.omit(truth_tbl_gamlss)[, .(obs = mean(former_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(former_smkr_prvl_obs),
#                               gamlss = mean(former_smoker_gamlss, na.rm = TRUE),
#                               mdl = mean(former_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(age)][, {
#                                 plot(x = age, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = age, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = age, y = gamlss, lty = 3, col = "black")
#                                 lines(x = age, y = mdl, lty = 2, col = "blue")}]
# 
# par(mfrow=c(1,1))
# na.omit(truth_tbl_gamlss)[, .(difference_former = (mean(former_smkr_prvl_mdl, na.rm = TRUE) - mean(former_smoker_gamlss, na.rm = TRUE))),
#                           keyby = .(age)][, {
#                             plot(x = age, y = difference_former,
#                                  type = "n",  col = "blue",
#                                  # ylim = c(-20, 50),
#                                  main = "model - observed",
#                                  ylab = "difference")
#                             axis(side = 1, seq(min(truth_tbl_gamlss$age), max(truth_tbl_gamlss$age), 10))
#                             axis(side = 2, seq(-1, 1 ))
#                             lines(x = age, y = difference_former, lty = 2, col = "blue")
#                           }]
# 
# 
# # Recalibrate Relapse above 65 ####
# message("Re-calibrate relapse - above 65")
# tbl_recalibrate_relapse_above_65 = copy(tbl_calibrated_above_65)
# 
# sub1_tbl_calibrated %<-% calibrate_recalibrate_relapse_above_65(tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, tbl_calibrated_above_65,
#                                                                 row_index = 1, age_start = 66, age_end = max(dt$pop$age))
# sub2_tbl_calibrated %<-% calibrate_recalibrate_relapse_above_65(tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, tbl_calibrated_above_65,
#                                                                 row_index = 2, age_start = 66, age_end = max(dt$pop$age))
# sub3_tbl_calibrated %<-% calibrate_recalibrate_relapse_above_65(tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, tbl_calibrated_above_65,
#                                                                 row_index = 3, age_start = 66, age_end = max(dt$pop$age))
# sub4_tbl_calibrated %<-% calibrate_recalibrate_relapse_above_65(tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, tbl_calibrated_above_65,
#                                                                 row_index = 4, age_start = 66, age_end = max(dt$pop$age))
# sub5_tbl_calibrated %<-% calibrate_recalibrate_relapse_above_65(tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, tbl_calibrated_above_65,
#                                                                 row_index = 5, age_start = 66, age_end = max(dt$pop$age))
# sub6_tbl_calibrated %<-% calibrate_recalibrate_relapse_above_65(tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, tbl_calibrated_above_65,
#                                                                 row_index = 6, age_start = 66, age_end = max(dt$pop$age))
# sub7_tbl_calibrated %<-% calibrate_recalibrate_relapse_above_65(tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, tbl_calibrated_above_65,
#                                                                 row_index = 7, age_start = 66, age_end = max(dt$pop$age))
# sub8_tbl_calibrated %<-% calibrate_recalibrate_relapse_above_65(tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, tbl_calibrated_above_65,
#                                                                 row_index = 8, age_start = 66, age_end = max(dt$pop$age))
# sub9_tbl_calibrated %<-% calibrate_recalibrate_relapse_above_65(tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, tbl_calibrated_above_65,
#                                                                 row_index = 9, age_start = 66, age_end = max(dt$pop$age))
# sub10_tbl_calibrated %<-% calibrate_recalibrate_relapse_above_65(tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, tbl_calibrated_above_65,
#                                                                  row_index = 10, age_start = 66, age_end = max(dt$pop$age))
# tbl_recalibrate_relapse_above_65 <- rbind(sub1_tbl_calibrated, sub2_tbl_calibrated,
#                                           sub3_tbl_calibrated, sub4_tbl_calibrated,
#                                           sub5_tbl_calibrated, sub6_tbl_calibrated,
#                                           sub7_tbl_calibrated, sub8_tbl_calibrated,
#                                           sub9_tbl_calibrated, sub10_tbl_calibrated)
# remove(sub1_tbl_calibrated, sub2_tbl_calibrated,
#        sub3_tbl_calibrated, sub4_tbl_calibrated,
#        sub5_tbl_calibrated, sub6_tbl_calibrated,
#        sub7_tbl_calibrated, sub8_tbl_calibrated,
#        sub9_tbl_calibrated, sub10_tbl_calibrated)
# 
# tbl_recalibrate_relapse_above_65 <- as.matrix(tbl_recalibrate_relapse_above_65)
# rownames(tbl_recalibrate_relapse_above_65) <- nam
# colnames(tbl_recalibrate_relapse_above_65) <- seq(1:15)
# recalibrate_relapse_modifier_above_65 <- as.data.table(tbl_recalibrate_relapse_above_65/tbl_calibrated_above_65)   # relapse modifer tbl
# setnafill(data.table(recalibrate_relapse_modifier_above_65), type = 'const', fill = 0, nan = NaN)
# 
# path <- paste0("lifecourse_model/relapse_modifier/", mc_iter, ".fst", collapse = "")
# message("Write relapse")
# # write.fst(recalibrate_relapse_modifier_above_65, path, compress = 100L)
# 
# ## Remove for loop: Apply modified relapse to population
# simsmok_relapse_calibrate(dt$pop, tbl_recalibrate_relapse_below_65, tbl_recalibrate_relapse_above_65, design$sim_prm$smoking_relapse_limit) # apply initiation, cessation value to synth pop
# dt$pop <- mortality(dt$pop, mc_iter, design)
# 
# model_tbl <- get_fs_prvl(dt) # get never_smoker prevaelence
# if ("former_smkr_prvl_mdl" %in% names(truth_tbl_gamlss)){
#   truth_tbl_gamlss[, former_smkr_prvl_mdl := NULL]}
# absorb_dt(truth_tbl_gamlss, model_tbl)
# truth_tbl_gamlss[is.na(truth_tbl_gamlss$former_smkr_prvl_mdl), former_smkr_prvl_mdl := former_smoker_gamlss]
# 
# model_tbl <- get_cs_prvl(dt) # get never_smoker prevaelence
# if ("current_smkr_prvl_mdl" %in% names(truth_tbl_gamlss)){
#   truth_tbl_gamlss[, current_smkr_prvl_mdl := NULL]}
# absorb_dt(truth_tbl_gamlss, model_tbl)
# truth_tbl_gamlss[is.na(truth_tbl_gamlss$current_smkr_prvl_mdl), current_smkr_prvl_mdl := current_smkr_prvl_gamlss]
# 
# model_tbl <- get_ns_prvl(dt) # get never_smoker prevaelence
# if ("never_smkr_prvl_mdl" %in% names(truth_tbl_gamlss)){
#   truth_tbl_gamlss[, never_smkr_prvl_mdl := NULL]}
# absorb_dt(truth_tbl_gamlss, model_tbl)
# truth_tbl_gamlss[is.na(truth_tbl_gamlss$never_smkr_prvl_mdl), never_smkr_prvl_mdl := never_smkr_prvl_gamlss]
# 
# ### Plot former smoker prevalence by year and by age ####
# # former_smoker_two_plots_calib_by_age
# # plot_path <- paste0("./HSE_data/plot/prevalence_after_calib/", interation_index, "_former_smoker.jpg", collapse = "")
# # jpeg(plot_path)
# par(mfrow=c(1,2))
# par(mar=c(2,2,1,1))
# 
# na.omit(truth_tbl_gamlss)[, .(obs = mean(former_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(former_smkr_prvl_obs),
#                               gamlss = mean(former_smoker_gamlss, na.rm = TRUE),
#                               mdl = mean(former_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(year)][, {
#                                 plot(x = year, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = year, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = year, y = gamlss, lty = 3, col = "black")
#                                 lines(x = year, y = mdl, lty = 2, col = "blue")}]
# na.omit(truth_tbl_gamlss)[, .(obs = mean(former_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(former_smkr_prvl_obs),
#                               gamlss = mean(former_smoker_gamlss, na.rm = TRUE),
#                               mdl = mean(former_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(age)][, {
#                                 plot(x = age, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = age, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = age, y = gamlss, lty = 3, col = "black")
#                                 lines(x = age, y = mdl, lty = 2, col = "blue")}]
# # #dev.off()
# par(mfrow=c(1,1))
# na.omit(truth_tbl_gamlss)[, .(difference = mean(former_smkr_prvl_mdl, na.rm = TRUE) - mean(former_smoker_gamlss, na.rm = TRUE)),
#                           keyby = .(age)][, {
#                             plot(x = age, y = difference,
#                                  type = "n",  col = "blue",
#                                  # ylim = c(-20, 50),
#                                  main = "model - observed",
#                                  ylab = "difference")
#                             axis(side = 1, seq(min(truth_tbl_gamlss$age), max(truth_tbl_gamlss$age), 10))
#                             # axis(side = 2, seq(-15, 30 , 1))
#                             lines(x = age, y = difference, lty = 2, col = "blue")
#                           }]
# 
# ### plot current smoker prevalence by year and by age ####
# # current_smoker_two_plots_calib_by_age
# # plot_path <- paste0("./HSE_data/plot/prevalence_after_calib/", interation_index, "_current_smoker.jpg", collapse = "")
# # jpeg(plot_path)
# 
# 
# na.omit(truth_tbl_gamlss)[, .(obs = mean(current_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(current_smkr_prvl_obs),
#                               gamlss = mean(current_smkr_prvl_gamlss, na.rm = TRUE),
#                               mdl = mean(current_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(year)][, {
#                                 plot(x = year, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = year, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = year, y = gamlss, lty = 3, col = "black")
#                                 lines(x = year, y = mdl, lty = 2, col = "blue")}]
# na.omit(truth_tbl_gamlss)[, .(obs = mean(current_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(current_smkr_prvl_obs),
#                               gamlss = mean(current_smkr_prvl_gamlss, na.rm = TRUE),
#                               mdl = mean(current_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(age)][, {
#                                 plot(x = age, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = age, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = age, y = gamlss, lty = 3, col = "black")
#                                 lines(x = age, y = mdl, lty = 2, col = "blue")}]
# 
# #dev.off()
# 
# par(mfrow=c(1,1))
# na.omit(truth_tbl_gamlss)[, .(difference = mean(current_smkr_prvl_mdl, na.rm = TRUE) - mean(current_smkr_prvl_obs, na.rm = TRUE)),
#                           keyby = .(age)][, {
#                             plot(x = age, y = difference,
#                                  type = "n",  col = "blue",
#                                  # ylim = c(-20, 50),
#                                  main = "model - observed",
#                                  ylab = "difference")
#                             axis(side = 1, seq(min(truth_tbl_gamlss$age), max(truth_tbl_gamlss$age), 10))
#                             # axis(side = 2, seq(-15, 30 , 1))
#                             lines(x = age, y = difference, lty = 2, col = "blue")
#                           }]
# 
# # plot never smoker prevalence by year and by age
# # current_smoker_two_plots_calib_by_age
# # plot_path <- paste0("./HSE_data/plot/prevalence_after_calib/", interation_index, "_never_smoker.jpg", collapse = "")
# # jpeg(plot_path)
# par(mfrow=c(1,2))
# par(mar=c(2,2,1,1))
# 
# na.omit(truth_tbl_gamlss)[, .(obs = mean(never_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(never_smkr_prvl_obs),
#                               gamlss = mean(never_smkr_prvl_gamlss, na.rm = TRUE),
#                               mdl = mean(never_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(year)][, {
#                                 plot(x = year, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = year, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = year, y = gamlss, lty = 3, col = "black")
#                                 lines(x = year, y = mdl, lty = 2, col = "blue")}]
# na.omit(truth_tbl_gamlss)[, .(obs = mean(never_smkr_prvl_obs_org, na.rm = TRUE),
#                               obs_revised = mean(never_smkr_prvl_obs),
#                               gamlss = mean(never_smkr_prvl_gamlss, na.rm = TRUE),
#                               mdl = mean(never_smkr_prvl_mdl, na.rm = TRUE)), keyby = .(age)][, {
#                                 plot(x = age, y = obs, type = "l",  ylim = c(0, 1))
#                                 lines(x = age, y = obs_revised, lty = 3, col = "red")
#                                 lines(x = age, y = gamlss, lty = 3, col = "black")
#                                 lines(x = age, y = mdl, lty = 2, col = "blue")}]
# 
# par(mfrow=c(1,1))
# na.omit(truth_tbl_gamlss)[, .(difference = mean(never_smkr_prvl_mdl, na.rm = TRUE) - mean(never_smkr_prvl_obs, na.rm = TRUE)),
#                           keyby = .(age)][, {
#                             plot(x = age, y = difference,
#                                  type = "n",  col = "blue",
#                                  # ylim = c(-20, 50),
#                                  main = "model - observed",
#                                  ylab = "difference")
#                             axis(side = 1, seq(min(truth_tbl_gamlss$age), max(truth_tbl_gamlss$age), 10))
#                             lines(x = age, y = difference, lty = 2, col = "blue")
#                           }]
# 
# 
# 
# # dev.off()
# # } # for loop
