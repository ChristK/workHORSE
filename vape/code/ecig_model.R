## Script for method 2 ##
## Continute form ecig_model.r to generate ecig realise_prob

source("global.R")

smoke_dt = fread("./vape/data/smk_output.csv", data.table = TRUE) # add the .N for this file
setnames(smoke_dt,  old = c("agegrp10"),"age_group")
smoke_dt[, year := year - 2000L]
smoke_dt = smoke_dt[! (age_group == "All" | sex == "All"| qimd == "All") ]
smoke_dt[, aux := as.numeric(substr(age_group, 1,2)) + 5 ]
smoke_dt[, `:=`(cohort = year - aux,
                aux = NULL)]

# AP model ####
fit_ap_current = fread("./vape/data/model_current_ap.csv", data.table = TRUE)
fit_ap_current = clone_dt(fit_ap_current, 30) # clone_dt for 30 times
setnames(fit_ap_current, ".id", "mc")
fit_ap_current[, prob_realised := rnorm(.N, mean = ecig, sd = se)] 
fit_ap_current[, prob_realised := pmin(pmax(prob_realised, 0), 1)] #limit to between (0,1)

fit_ap_never = fread("./vape/data/model_never_ap.csv", data.table = TRUE)
fit_ap_never = clone_dt(fit_ap_never, 30)
setnames(fit_ap_never, ".id", "mc")
fit_ap_never[, prob_realised := rnorm(.N, mean = ecig, sd = se)] 
fit_ap_never[, prob_realised := pmin(pmax(prob_realised, 0), 1)] #limit to between (0,1)

fit_ap_former = fread("./vape/data/model_former_ap.csv", data.table = TRUE)
fit_ap_former = clone_dt(fit_ap_former, 30)
setnames(fit_ap_former, ".id", "mc")
fit_ap_former[, prob_realised := rnorm(.N, mean = ecig, sd = se)] 
fit_ap_former[, prob_realised := pmin(pmax(prob_realised, 0), 1)] #limit to between (0,1)

smoke_dt[, smoke_former_prob := 1 - smok_never_smok_xps - smok_active_smok_xps]
current_smoker_ap = fit_ap_current[smoke_dt, on = .(year, age_group, sex, qimd, mc)] # TODO: add mc
current_smoker_ap[, `:=`(dual_total = round(prob_realised * smok_active_smok_xps * N),
                       smoke_exclusive_total = round((1-prob_realised) * smok_active_smok_xps * N)) ]
current_smoker_ap[scenario == "sc1" & sex == "men" & qimd == 2 & age_group == "30-39" & mc == 5][, plot(x = year, y = dual_total)]   

never_smoker_ap = fit_ap_never[smoke_dt, on = .(year, age_group, sex, qimd, mc)]
never_smoker_ap[, `:=`(ecig_exclusive_total = prob_realised * smok_never_smok_xps * N,
                       neither_total = (1-prob_realised) * smok_never_smok_xps * N)]
never_smoker_ap[scenario == "sc1" & sex == "men" & qimd == 2 & age_group == "30-39" & mc == 5][, plot(x = year, y = ecig_exclusive_total)]   

former_smoker_ap = fit_ap_former[smoke_dt, on = .(year, age_group, sex, qimd, mc)]
former_smoker_ap[, `:=`(former_smoker_vaper_total = prob_realised * smoke_former_prob * N,
                       former_smoker_nonvaper_total = (1-prob_realised) * smoke_former_prob * N)]
former_smoker_ap[scenario == "sc1" & sex == "men" & qimd == 2 & age_group == "30-39" & mc == 5][, plot(x = year, y = former_smoker_vaper_total)]   

# PC model ####
fit_pc_current = fread("./vape/data/model_current_pc.csv", data.table = TRUE)
fit_pc_current = clone_dt(fit_pc_current, 30)
setnames(fit_pc_current, ".id", "mc")
fit_pc_current[, prob_realised := rnorm(.N, mean = ecig, sd = se)] 
fit_pc_current[, prob_realised := pmin(pmax(prob_realised, 0), 1)] #limit to between (0,1)


fit_pc_former = fread("./vape/data/model_former_pc.csv", data.table = TRUE)
fit_pc_former = clone_dt(fit_pc_former, 30)
setnames(fit_pc_former, ".id", "mc")
fit_pc_former[, prob_realised := rnorm(.N, mean = ecig, sd = se)] 
fit_pc_former[, prob_realised := pmin(pmax(prob_realised, 0), 1)] #limit to between (0,1)

fit_pc_never = fread("./vape/data/model_never_pc.csv", data.table = TRUE)
fit_pc_never = clone_dt(fit_pc_never, 30)
setnames(fit_pc_never, ".id", "mc")
fit_pc_never[, prob_realised := rnorm(.N, mean = ecig, sd = se)] 
fit_pc_never[, prob_realised := pmin(pmax(prob_realised, 0), 1)] #limit to between (0,1)

current_smoker_pc = fit_pc_current[smoke_dt, on = .(year, cohort, sex, qimd, mc)]
current_smoker_pc[, `:=`(dual_total = prob_realised * smok_active_smok_xps * N,
                         smoke_exclusive_total = (1-prob_realised) * smok_active_smok_xps * N)]
current_smoker_pc[scenario == "sc1" & sex == "men" & qimd == 2 & age_group == "30-39" & mc == 5][, plot(x = year, y = dual_total)]   

never_smoker_pc = fit_pc_never[smoke_dt, on = .(year, cohort, sex, qimd, mc)]
never_smoker_pc[, `:=`(ecig_exclusive_total = prob_realised * smok_never_smok_xps * N,
                       neither_total = (1-prob_realised) * smok_never_smok_xps * N)]
never_smoker_pc[scenario == "sc1" & sex == "men" & qimd == 2 & age_group == "30-39" & mc == 5][, plot(x = year, y = ecig_exclusive_total)]   

former_smoker_pc = fit_pc_former[smoke_dt, on = .(year, cohort, sex, qimd, mc)]
former_smoker_pc[, `:=`(former_smoker_vaper_total = prob_realised * smoke_former_prob * N,
                        former_smoker_nonvaper_total = (1-prob_realised) * smoke_former_prob * N)]
former_smoker_pc[scenario == "sc1" & sex == "men" & qimd == 2 & age_group == "30-39" & mc == 5][, plot(x = year, y = former_smoker_vaper_total)]   

# summarise results ####
## AP ####
tt_never = never_smoker_ap[, .(mc, scenario, year, age_group, cohort, sex, qimd, N, ecig_exclusive_total, neither_total)]
tt_former = former_smoker_ap[, .(mc, scenario, year, age_group, cohort, sex, qimd, N, former_smoker_vaper_total, former_smoker_nonvaper_total)]
tt_current = current_smoker_ap[, .(mc, scenario, year, age_group, cohort, sex, qimd, N, dual_total, smoke_exclusive_total)]

ttt = copy(tt_never)
ttt = tt_former[ttt, on = .(mc, scenario, year, age_group, cohort, sex, qimd, N), allow.cartesian=TRUE]
ttt = tt_current[ttt, on = .(mc, scenario, year, age_group, cohort, sex, qimd, N), allow.cartesian=TRUE]

ttt = ttt[, `:=`(vape = dual_total + former_smoker_vaper_total + ecig_exclusive_total,
                 non_vape = former_smoker_nonvaper_total + smoke_exclusive_total + neither_total,
                 smoke = dual_total + smoke_exclusive_total,
                 non_smoke = ecig_exclusive_total + neither_total)]
ttt = ttt[, `:=`(vape_prevalence = vape/N,
                 smoke_prevalence = smoke/N,
                 non_smoke_prevalence = non_smoke/N,
                 non_vape_prevalence =  non_vape/N,
                 dual_prevalence = dual_total/N, 
                 ecig_exclusive_prevalence = ecig_exclusive_total/N,
                 smoke_exclusive_prevalence = smoke_exclusive_total/N,
                 neither_prevalence = neither_total/N
                 )]
ap_dt = ttt[,  .(vape_prevalence = weighted.mean(vape_prevalence, N), 
                  smoke_prevalence = weighted.mean(smoke_prevalence, N), 
                  non_smoke_prevalence = weighted.mean(non_smoke_prevalence, N), 
                  non_vape_prevalence = weighted.mean(non_vape_prevalence, N),
                 dual_prevalence = weighted.mean(dual_prevalence, N), 
                 ecig_exclusive_prevalence = weighted.mean(ecig_exclusive_prevalence, N), 
                 smoke_exclusive_prevalence = weighted.mean(smoke_exclusive_prevalence, N), 
                 neither_prevalence = weighted.mean(neither_prevalence, N)
                 ), 
              by = .(scenario, year)]
fwrite(ttt, "./vape/data/result_ap.csv")
remove(tt_never, tt_former, tt_current, ttt)

# 
# p = ggplot(ap_dt[scenario == "sc1"], 
#        aes(x = year, y = smoke_prevalence)) +
#   geom_line()
# 
# p = ggplot(ap_dt[scenario == "sc1"], 
#            aes(x = year, y = dual_prevalence)) +
#   geom_line()
# 
# p = ggplot(ap_dt[scenario == "sc1"], 
#            aes(x = year, y = ecig_exclusive_prevalence)) +
#   geom_line()
# 
# p = ggplot(ap_dt[scenario == "sc1"], 
#            aes(x = year, y = smoke_exclusive_prevalence)) +
#   geom_line()
# 
# p = ggplot(ap_dt[scenario == "sc1"], 
#            aes(x = year, y = neither_prevalence)) +
#   geom_line()
# 
# p = ggplot(ap_dt[scenario == "sc1"], 
#            aes(x = year, y = vape_prevalence, color = age_group)) +
#   geom_line(aes(group = age_group))
# 
# ggplot(ap_dt[scenario == "sc1"], 
#        aes(x = year, y = non_vape_prevalence, color = age_group)) +
#   geom_line(aes(group = age_group))
# 
# ggsave(filename = "./output/ecig_exclusive_total.jpg", p)

# PC

tt_never = never_smoker_pc[, .(mc, scenario, year, age_group, cohort, sex, qimd, N, ecig_exclusive_total, neither_total)]
tt_former = former_smoker_pc[, .(mc, scenario, year, age_group, cohort, sex, qimd, N, former_smoker_vaper_total, former_smoker_nonvaper_total)]
tt_current = current_smoker_pc[, .(mc, scenario, year, age_group, cohort, sex, qimd, N, dual_total, smoke_exclusive_total)]

ttt = copy(tt_never)
ttt = tt_former[ttt, on = .(mc, scenario, year, age_group, cohort, sex, qimd, N), allow.cartesian=TRUE]
ttt = tt_current[ttt, on = .(mc, scenario, year, age_group, cohort, sex, qimd, N), allow.cartesian=TRUE]
ttt = ttt[, `:=`(vape = dual_total + former_smoker_vaper_total + ecig_exclusive_total,
                     non_vape = former_smoker_nonvaper_total + smoke_exclusive_total + neither_total,
                     smoke = dual_total + smoke_exclusive_total,
                     non_smoke = ecig_exclusive_total + neither_total)]
ttt = ttt[, `:=`(vape_prevalence = vape/N,
                     smoke_prevalence = smoke/N,
                     non_smoke_prevalence = non_smoke/N,
                     non_vape_prevalence =  non_vape/N,
                    dual_prevalence = dual_total/N, 
                    ecig_exclusive_prevalence = ecig_exclusive_total/N,
                    smoke_exclusive_prevalence = smoke_exclusive_total/N,
                    neither_prevalence = neither_total/N
                 )]

fwrite(ttt, "./vape/data/results_pc.csv")

pc_dt = ttt[, .(vape_prevalence = weighted.mean(vape_prevalence, N),
             smoke_prevalence = weighted.mean(smoke_prevalence, N), 
             non_smoke_prevalence = weighted.mean(non_smoke_prevalence, N), 
             non_vape_prevalence = weighted.mean(non_vape_prevalence, N),
             dual_prevalence = weighted.mean(dual_prevalence, N), 
             ecig_exclusive_prevalence = weighted.mean(ecig_exclusive_prevalence, N), 
             smoke_exclusive_prevalence = weighted.mean(smoke_exclusive_prevalence, N), 
             neither_prevalence = weighted.mean(neither_prevalence, N)),
      by = .(scenario, year)]
remove(tt_never, tt_former, tt_current, ttt)

# pc_dt[, cohort := cohort + 2000]
# 
# ggplot(pc_dt[scenario == "sc1"], 
#        aes(x = year, y = vape_prevalence)) +
#   geom_line()
# 
# ggplot(pc_dt[scenario == "sc1"], 
#        aes(x = year, y = dual_prevalence)) +
#   geom_line()
# 
# ggplot(pc_dt[scenario == "sc1"], 
#        aes(x = year, y = ecig_exclusive_prevalence)) +
#   geom_line()
# 
# ggplot(pc_dt[scenario == "sc1"], 
#        aes(x = year, y = smoke_exclusive_prevalence)) +
#   geom_line()
# 
# ggplot(pc_dt[scenario == "sc1"], 
#        aes(x = year, y = vape_prevalence)) +
#   geom_line()
# 
# pc_dt[scenario == "sc1" & year  == 22 & sex == "men" & qimd == "3" & age_group == "40-49", .(scenario, year, age_group, vape_prevalence, smoke_prevalence, non_smoke_prevalence, non_vape_prevalence)]

# Continue with analysis_ecig.r to plot graphs
