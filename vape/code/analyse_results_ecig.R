# ecig script 3
# use to analyse results and plot 

# vape = dual_total + former_smoker_vaper_total + ecig_exclusive_total,
# non_vape = former_smoker_nonvaper_total + smoke_exclusive_total + neither_total,
# smoke = dual_total + smoke_exclusive_total,
# non_smoke = ecig_exclusive_total + neither_total

library(data.table)
library(fst)
library(ggplot2)

ap_dt = fread("./data/result_ap.csv", data.table = TRUE)
pc_dt = fread("./data/results_pc.csv", data.table = TRUE)
pc_dt[, cohort := cohort + 2000L]

ap_dt[scenario == "sc1", scenario := "baseline"]
ap_dt[scenario == "sc2", scenario := "Combined"]
ap_dt[scenario == "sc3", scenario := "ServicesUP"]
ap_dt[scenario == "sc4", scenario := "MinAge21"]
ap_dt[scenario == "sc5", scenario := "TaxUP"]
ap_dt[scenario == "sc6", scenario := "baseline"]
pc_dt[scenario == "sc1", scenario := "baseline"]
pc_dt[scenario == "sc2", scenario := "Combined"]
pc_dt[scenario == "sc3", scenario := "ServicesUP"]
pc_dt[scenario == "sc4", scenario := "MinAge21"]
pc_dt[scenario == "sc5", scenario := "TaxUP"]
pc_dt[scenario == "sc6", scenario := "baseline"]
# AP ####
qimd_year_ap = ap_dt[,  .(vape_prevalence = weighted.mean(vape_prevalence, N), 
                              smoke_prevalence = weighted.mean(smoke_prevalence, N), 
                              non_smoke_prevalence = weighted.mean(non_smoke_prevalence, N), 
                              non_vape_prevalence = weighted.mean(non_vape_prevalence, N),
                              dual_prevalence = weighted.mean(dual_prevalence, N), 
                              ecig_exclusive_prevalence = weighted.mean(ecig_exclusive_prevalence, N), 
                              smoke_exclusive_prevalence = weighted.mean(smoke_exclusive_prevalence, N), 
                              neither_prevalence = weighted.mean(neither_prevalence, N)), 
                         by = .(scenario, year,  qimd)]
age_year_ap = ap_dt[,  .(vape_prevalence = weighted.mean(vape_prevalence, N), 
                 smoke_prevalence = weighted.mean(smoke_prevalence, N), 
                 non_smoke_prevalence = weighted.mean(non_smoke_prevalence, N), 
                 non_vape_prevalence = weighted.mean(non_vape_prevalence, N),
                 dual_prevalence = weighted.mean(dual_prevalence, N), 
                 ecig_exclusive_prevalence = weighted.mean(ecig_exclusive_prevalence, N), 
                 smoke_exclusive_prevalence = weighted.mean(smoke_exclusive_prevalence, N), 
                 neither_prevalence = weighted.mean(neither_prevalence, N)), 
            by = .(scenario, year, age_group)]

year_ap = ap_dt[,  .(vape_prevalence = weighted.mean(vape_prevalence, N), 
                    smoke_prevalence = weighted.mean(smoke_prevalence, N), 
                    non_smoke_prevalence = weighted.mean(non_smoke_prevalence, N), 
                    non_vape_prevalence = weighted.mean(non_vape_prevalence, N),
                    dual_prevalence = weighted.mean(dual_prevalence, N), 
                    ecig_exclusive_prevalence = weighted.mean(ecig_exclusive_prevalence, N), 
                    smoke_exclusive_prevalence = weighted.mean(smoke_exclusive_prevalence, N), 
                    neither_prevalence = weighted.mean(neither_prevalence, N)), 
               by = .(scenario, year)]

p = ggplot(year_ap, 
           aes(x = year, y = smoke_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/smoke_prevalence_by_policy_ap.jpg", p)

p = ggplot(year_ap, 
           aes(x = year, y = vape_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/ecig_ever_prevalence_by_policy_ap.jpg", p)

p = ggplot(year_ap, 
           aes(x = year, y = dual_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/dual_prevalence_by_policy_ap.jpg", p)

p = ggplot(year_ap, 
           aes(x = year, y = ecig_exclusive_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/ecig_exclusive_prevalence_by_policy_ap.jpg", p)

p = ggplot(year_ap, 
           aes(x = year, y = smoke_exclusive_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/smoke_exclusive_prevalence_by_policy_ap.jpg", p)

p = ggplot(year_ap, 
           aes(x = year, y = neither_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/neither_prevalence_by_policy_ap.jpg", p)

# ggplot(year_ap, 
#        aes(x = year, y = ecig_exclusive_prevalence, color = scenario)) +
#   geom_line(aes(group = scenario))

p=
 ggplot(age_year_ap[scenario == "baseline"], 
           aes(x = year, y = vape_prevalence, color = age_group)) +
  geom_line(aes(group = age_group))
ggsave(filename = "./output/ecig_ever_prevalence_by_age_ap.jpg", p)

p = ggplot(qimd_year_ap[scenario == "baseline" ], 
       aes(x = year, y = ecig_exclusive_prevalence, color = qimd)) +
  geom_line(aes(group = qimd))
ggsave(filename = "./output/ecig_exclusive_prevalence_by_qimd_ap.jpg", p)

p = ggplot(qimd_year_ap[scenario == "baseline" ], 
       aes(x = year, y = dual_prevalence, color = qimd)) +
  geom_line(aes(group = qimd))
ggsave(filename = "./output/dual_prevalence_by_qimd_ap.jpg", p)

# PC ####
qimd_year_pc = pc_dt[, .(vape_prevalence = weighted.mean(vape_prevalence, N),
                           smoke_prevalence = weighted.mean(smoke_prevalence, N), 
                           non_smoke_prevalence = weighted.mean(non_smoke_prevalence, N), 
                           non_vape_prevalence = weighted.mean(non_vape_prevalence, N),
                         dual_prevalence = weighted.mean(dual_prevalence, N), 
                         ecig_exclusive_prevalence = weighted.mean(ecig_exclusive_prevalence, N), 
                         smoke_exclusive_prevalence = weighted.mean(smoke_exclusive_prevalence, N), 
                         neither_prevalence = weighted.mean(neither_prevalence, N)), 
                       by = .(scenario, year, qimd)]
cohort_year_pc = pc_dt[, .(vape_prevalence = weighted.mean(vape_prevalence, N),
        smoke_prevalence = weighted.mean(smoke_prevalence, N), 
        non_smoke_prevalence = weighted.mean(non_smoke_prevalence, N), 
        non_vape_prevalence = weighted.mean(non_vape_prevalence, N),
        dual_prevalence = weighted.mean(dual_prevalence, N), 
        ecig_exclusive_prevalence = weighted.mean(ecig_exclusive_prevalence, N), 
        smoke_exclusive_prevalence = weighted.mean(smoke_exclusive_prevalence, N), 
        neither_prevalence = weighted.mean(neither_prevalence, N)), 
    by = .(scenario, year, age_group)]

year_pc = pc_dt[, .(vape_prevalence = weighted.mean(vape_prevalence, N),
                           smoke_prevalence = weighted.mean(smoke_prevalence, N), 
                           non_smoke_prevalence = weighted.mean(non_smoke_prevalence, N), 
                           non_vape_prevalence = weighted.mean(non_vape_prevalence, N),
                    dual_prevalence = weighted.mean(dual_prevalence, N), 
                    ecig_exclusive_prevalence = weighted.mean(ecig_exclusive_prevalence, N), 
                    smoke_exclusive_prevalence = weighted.mean(smoke_exclusive_prevalence, N), 
                    neither_prevalence = weighted.mean(neither_prevalence, N)), 
                       by = .(scenario, year)]

p = ggplot(year_pc[, .(scenario, year, vape_prevalence, smoke_prevalence, non_smoke_prevalence, non_vape_prevalence) ], 
           aes(x = year, y = smoke_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/smoke_prevalence_by_policy_pc.jpg", p)

p = ggplot(year_pc[, .(scenario, year, vape_prevalence, smoke_prevalence, non_smoke_prevalence, non_vape_prevalence) ], 
           aes(x = year, y = vape_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/ecig_ever_prevalence_by_policy_pc.jpg", p)

p = ggplot(year_pc, 
           aes(x = year, y = dual_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/dual_prevalence_by_policy_pc.jpg", p)

p = ggplot(year_pc, 
           aes(x = year, y = ecig_exclusive_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/ecig_exclusive_prevalence_by_policy_pc.jpg", p)

p = ggplot(year_pc, 
           aes(x = year, y = neither_prevalence, color = scenario)) +
  geom_line(aes(group = scenario))
ggsave(filename = "./output/neither_prevalence_by_policy_pc.jpg", p)

p=
  ggplot(cohort_year_pc[scenario == "baseline"], 
         aes(x = year, y = vape_prevalence, color = age_group)) +
  geom_line(aes(group = age_group))
ggsave(filename = "./output/ecig_ever_prevalence_by_age_pc.jpg", p)

p = ggplot(qimd_year_pc[scenario == "baseline" ], 
           aes(x = year, y = vape_prevalence, color = qimd)) +
  geom_line(aes(group = qimd))
ggsave(filename = "./output/ecig_exclusive_prevalence_by_qimd_pc.jpg", p)

p = ggplot(qimd_year_pc[scenario == "baseline" ], 
           aes(x = year, y = dual_prevalence, color = qimd)) +
  geom_line(aes(group = qimd))
ggsave(filename = "./output/dual_prevalence_by_qimd_pc.jpg", p)

p = ggplot(cohort_year_pc[scenario == "baseline"], 
           aes(x = year, y = vape_prevalence, color = age_group)) +
  geom_line(aes(group = cohort))