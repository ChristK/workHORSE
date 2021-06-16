library(data.table)
library(fst)
sp <- fread("./output/smk_output.csv", stringsAsFactors = TRUE)

# reduce to 1 entry per mc
sp <- sp[, lapply(.SD, mean), keyby = .(year, sex, agegrp10, qimd, scenario, mc)]

# overall
tt <- sp[sex == "All" & agegrp10 == "All" & qimd == "All", ]
tt <- tt[, lapply(.SD, quantile, probs = 0.5), 
  keyby = .(year, sex, agegrp10, qimd, scenario), 
  .SDcols = !"mc"]
plot(tt$year, tt$smok_active_smok_xps, col = tt$scenario)
plot(tt$year, tt$smok_never_smok_xps, col = tt$scenario)

tt <- sp[sex == "All" & agegrp10 == "All" & qimd == "All" & year == 2070, ]
tt <- split(tt, tt$scenario)
signif(quantile(tt$sc1$smok_active_smok_xps - tt$sc2$smok_active_smok_xps, probs = c(0.5, 0.025, 0.975)), 2)

# by qimd
tt <- sp[sex == "All" & agegrp10 == "All" & qimd == "1 most deprived", ]
tt <- tt[, lapply(.SD, quantile, probs = 0.5), 
  keyby = .(year, sex, agegrp10, qimd, scenario), 
  .SDcols = !"mc"]
plot(tt$year, tt$smok_active_smok_xps, col = tt$scenario)
plot(tt$year, tt$smok_never_smok_xps, col = tt$scenario)

tt <- sp[sex == "All" & agegrp10 == "All" & qimd == "1 most deprived" & year == 2070, ]
tt <- split(tt, tt$scenario)
signif(quantile(tt$sc1$smok_active_smok_xps  - tt$sc2$smok_active_smok_xps, probs = c(0.5, 0.025, 0.975)), 2)

tt <- sp[sex == "All" & agegrp10 == "All" & qimd == "5 least deprived" & year == 2070, ]
tt <- split(tt, tt$scenario)
signif(quantile(tt$sc1$smok_active_smok_xps  - tt$sc2$smok_active_smok_xps, probs = c(0.5, 0.025, 0.975)), 2)

tt <- sp[sex == "All" & agegrp10 == "All" & qimd == "4" & year == 2070, ]
tt <- split(tt, tt$scenario)
signif(quantile(tt$sc1$smok_active_smok_xps  - tt$sc2$smok_active_smok_xps, probs = c(0.5, 0.025, 0.975)), 2)

# 15-24
tt <- sp[sex == "All" & agegrp10 == "15-24" & qimd == "All", ]
tt <- tt[, lapply(.SD, quantile, probs = 0.5), 
  keyby = .(year, sex, agegrp10, qimd, scenario), 
  .SDcols = !"mc"]
plot(tt$year, tt$smok_active_smok_xps, col = tt$scenario)
plot(tt$year, tt$smok_never_smok_xps, col = tt$scenario)



tt <- fread("~/pCloudDrive/My Papers/Lancet PH conf/Vincy/workHORSE_raw_2021-06-13.csv")
tt[, unique(scenario)] # "sc2" "sc3" "sc4"
tt[, unique(friendly_name)] # "levelup"          "levelup_fat"      "levelup_fat_mort"

tt[year == 2070, .(cpp_cml = sum(cpp_cvd_cml + cpp_t2dm_cml + cpp_copd_cml + cpp_poststroke_dementia_cml +
    cpp_lung_ca_cml + cpp_breast_ca_cml + cpp_colon_ca_cml),
  cypp_cml = sum(cypp_cvd_cml + cypp_t2dm_cml + cypp_copd_cml + cypp_poststroke_dementia_cml +
      cypp_lung_ca_cml + cypp_breast_ca_cml + cypp_colon_ca_cml),
  dpp_cml = sum(dpp_chd_cml + dpp_stroke_cml + dpp_nonmodelled_cml + dpp_copd_cml +
      dpp_lung_ca_cml + dpp_breast_ca_cml + dpp_colon_ca_cml)
),
  keyby = .(scenario, mc)
][, lapply(.SD,
  function(x) {
    signif(quantile(
      x, prob = c(0.5, 0.025, 0.975), na.rm = TRUE
    ), 2)
  }), keyby = scenario, .SDcols = !"mc"]



tt[year == 2070, .(qalys = sum(net_utility_cml)),
  keyby = .(scenario, mc)
][, as.list(signif(quantile(
  qalys,
  prob = c(0.5, 0.025, 0.975),
  na.rm = TRUE
), 2)), keyby = scenario]

# healthcare
tt[year == 2070, .(costs = sum(total_hcp_cost_cml)),
  keyby = .(scenario, mc)
][, as.list(signif(quantile(
  costs/1e6,
  prob = c(0.5, 0.025, 0.975),
  na.rm = TRUE
), 2)), keyby = scenario]

# healthcare + social care
tt[year == 2070, .(costs = sum(total_hscp_cost_cml)),
  keyby = .(scenario, mc)
][, as.list(signif(quantile(
  costs/1e6,
  prob = c(0.5, 0.025, 0.975),
  na.rm = TRUE
), 2)), keyby = scenario]

# societal
tt[year == 2070, .(costs = sum(societal_cost_cml)),
  keyby = .(scenario, mc)
][, as.list(signif(quantile(
  costs/1e9,
  prob = c(0.5, 0.025, 0.975),
  na.rm = TRUE
), 2)), keyby = scenario]


# equity
tt[year == 2070, .(cpp_cml = sum(cpp_cvd_cml + cpp_t2dm_cml + cpp_copd_cml + cpp_poststroke_dementia_cml +
    cpp_lung_ca_cml + cpp_breast_ca_cml + cpp_colon_ca_cml),
  cypp_cml = sum(cypp_cvd_cml + cypp_t2dm_cml + cypp_copd_cml + cypp_poststroke_dementia_cml +
      cypp_lung_ca_cml + cypp_breast_ca_cml + cypp_colon_ca_cml),
  dpp_cml = sum(dpp_chd_cml + dpp_stroke_cml + dpp_nonmodelled_cml + dpp_copd_cml +
      dpp_lung_ca_cml + dpp_breast_ca_cml + dpp_colon_ca_cml)
),
  keyby = .(qimd, scenario, mc)
][, lapply(.SD,
  function(x) {
    signif(quantile(
      x, prob = c(0.5, 0.025, 0.975), na.rm = TRUE
    ), 2)
  }), keyby = .(qimd, scenario), .SDcols = !"mc"]





# Why the first method does not work
lutbl <-
  read_fst("./lifecourse_models/smok_incid_table.fst",
    as.data.table = TRUE)
lutbl[between(age, 3, 21) & between(year, 6, 9), mean(mu), keyby = .(age, year)
  ][, dcast(.SD, age~year)]

lutbl <-
  read_fst("./lifecourse_models/smok_cess_table.fst",
    as.data.table = TRUE)
lutbl[between(age, 16, 21), mean(mu), keyby = age]

lutbl <-
  read_fst("./lifecourse_models/smok_status_table.fst",
    as.data.table = TRUE)
lutbl[between(age, 3, 21) & between(year, 6, 9), mean(mu), keyby = .(age, year)
][, dcast(.SD, age~year)]

