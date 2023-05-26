# ecig script 2
## use to fit ecig AP and ecig PC model


# decompose age-period-cohort (age-period; period-cohort) - addictive so have cohort effect (difficult to quit)
# APC - projection is acceptable so we used it as the starting point for smoking
# Bayesian APC - start with prior, 

# by year, sex, agegroup, imd - vaping trend with uncertainty (by weight)


# TODO: whichfrst get the gateway effect
# TODO: ecigstp get willingess to quit

library(data.table)
library(fst)
library(ggplot2)

df_13_19 = read_fst("./data/HSE_ts_13_19.fst", as.data.table = TRUE)
df_ecig = read_fst("./data/HSE_ts_13_19_e_cig.fst", as.data.table =  TRUE)
smoke_dt = fread("./data/smk_output.csv", data.table = TRUE) # add the .N for this file

setnames(smoke_dt,  old = c("agegrp10"),"age_group")
smoke_dt[, year := year - 2000L]
smoke_dt = smoke_dt[! (age_group == "All" | sex == "All"| qimd == "All") ]
smoke_dt[, aux := as.numeric(substr(age_group, 1,2)) + 5 ]
smoke_dt[, `:=`(cohort = year - aux,
                aux = NULL)]

hse = df_13_19[df_ecig, on = .(id, year)]
hse = hse[age >= 16 ] # Correct with this step?
hse[ecigevr == "no", ecignw := "no"] # TODO: add documents ( ecignw is only asked when ecigever is yes)
hse[, ecig := as.integer(as.character(factor(ecignw, c("yes", "no"), 1:0)))]

hse = hse[year %in% c(16:19), .(id, year, sex, age, qimd, smok_status, ecig, wt_int)]
#write.csv(hse, "./output/all.csv") # uncomment when needed

#hse = fread("./output/all.csv", data.table = TRUE)
hse[, age_group := paste0((age %/% 10) * 10, "-", ((age %/% 10) * 10) + 9)] #create age group: 10 - 29; 30 -49 
hse[, cohort := year - age] # year - middle of the age group

hse[, summary(ecig)] # Only valid for years 16-19
hse[smok_status == "4", mean(ecig, na.rm = T), keyby = age][, plot(age, V1)]
hse[smok_status == "3", mean(ecig, na.rm = T), keyby = age][, plot(age, V1)]
hse[smok_status == "1", mean(ecig, na.rm = T), keyby = age][, plot(age, V1)]

# AP model ####
# current smoker
m1.3 <- glm(ecig ~ age_group + year, hse[smok_status == "4"], family = binomial, weights = wt_int)

m1 <- glm(ecig ~ age_group + year + sex + qimd, hse[smok_status == "4"], family = binomial, weights = wt_int)
m1.2 <- glm(ecig ~ age_group + year + sex , hse[smok_status == "4"], family = binomial, weights = wt_int)
m1.2.1 <- glm(ecig ~ age_group + year + qimd , hse[smok_status == "4"], family = binomial, weights = wt_int)
m1.3.1 <- glm(ecig ~ age_group + qimd, hse[smok_status == "4"], family = binomial, weights = wt_int)
m1.3.2 <- glm(ecig ~ year + sex, hse[smok_status == "4"], family = binomial, weights = wt_int)
m1.3.3 <- glm(ecig ~ year + qimd, hse[smok_status == "4"], family = binomial, weights = wt_int)
m <- glm(ecig ~ age_group * year, hse[smok_status == "4"], family = binomial, weights = wt_int)
AIC(m1, m1.2, m1.3, m)
AIC(m1.3.2, m1.3.3, m1.3.1, m1.2.1)

ap_m_current = m1.3 # ecig ~ age_group + year
summary(ap_m_current)


m2.2 <- glm(ecig ~ age_group + year + sex , hse[smok_status == "1"], family = binomial, weights = wt_int)

m2 <- glm(ecig ~ age_group + year + sex + qimd, hse[smok_status == "1"], family = binomial, weights = wt_int)
m2.2.1 <- glm(ecig ~ age_group + year + qimd , hse[smok_status == "1"], family = binomial, weights = wt_int)
m2.3 <- glm(ecig ~ age_group + year, hse[smok_status == "1"], family = binomial, weights = wt_int)
m2.3.1 <- glm(ecig ~ age_group + qimd, hse[smok_status == "1"], family = binomial, weights = wt_int)
m2.3.2 <- glm(ecig ~ year + sex, hse[smok_status == "1"], family = binomial, weights = wt_int)
m2.3.3 <- glm(ecig ~ year + qimd, hse[smok_status == "1"], family = binomial, weights = wt_int)
AIC(m2.2, m2.2.1, m2.3, m2, m2.3.1, m2.3.2, m2.3.3)

ap_m_never = m2.2 # ecig ~ age_group + year + sex
summary(ap_m_never)

m3.1 <- glm(ecig ~ age_group + year + sex + qimd, hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)

m3.2 <- glm(ecig ~ age_group + year + sex , hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
m3.2.1 <- glm(ecig ~ age_group + year + qimd , hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
m3.3 <- glm(ecig ~ age_group + year, hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
m3.3.1 <- glm(ecig ~ age_group + qimd, hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
m3.3.2 <- glm(ecig ~ year + sex, hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
m3.3.3 <- glm(ecig ~ year + qimd, hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
AIC(m3.1, m3.2, m3.3, m3.2.1, m3.3.1, m3.3.2, m3.3.3)

ap_m_former = m3.1 # ecig ~ age_group + period + sex + qimd
summary(ap_m_former)
round(exp(coefficients(ap_m_former)),2)

# 
fit_ap = CJ(year = c(16:73), age_group = unique(smoke_dt$age_group), sex = unique(smoke_dt$sex), qimd = unique(smoke_dt$qimd))
fit_ap_current = copy(fit_ap) 
fit_ap_current = fit_ap_current[, c("ecig", "se", "scale") := predict.glm(ap_m_current, newdata = fit_ap_current[, .(age_group, year)], type = "response", se.fit = TRUE)]

fit_ap_current[, mean(ecig), keyby = year][, plot(year, V1)]
fit_ap_current[, smok_status := "4"]
fwrite(fit_ap_current, file = "./vape/data/model_current_ap.csv")

fit_ap_never = copy(fit_ap) 
fit_ap_never = fit_ap_never[, c("ecig", "se", "scale") := predict.glm(ap_m_never, newdata = fit_ap_never[, .(age_group, year, sex)], type = "response", se.fit = TRUE) ]
fit_ap_never[, mean(ecig), keyby = year][, plot(year, V1)]
fit_ap_never[, smok_status := "1"]
fwrite(fit_ap_never, file = "./vape/data/model_never_ap.csv")

fit_ap_former = copy(fit_ap) 
fit_ap_former = fit_ap_former[, c("ecig", "se", "scale") := predict.glm(ap_m_former, newdata = fit_ap_former, type = "response", se.fit = TRUE) ]
fit_ap_former[, smok_status := "3"]
fwrite(fit_ap_former, file = "./vape/data/model_former_ap.csv")

# PC model ####
# current smoker
# scatter plot 
hse[smok_status == "4", weighted.mean(ecig, wt_int), keyby = cohort][, scatter.smooth(cohort, V1)]
pc_m = glm(ecig ~ poly(cohort, 3), hse[smok_status == "4"], family = binomial, weights = wt_int)
lines(x = -80:10, y = predict(pc_m, newdata = data.frame('cohort' = -80:10), type = "response"), col = 'red')


hse[smok_status == "3", weighted.mean(ecig, wt_int), keyby = cohort][, scatter.smooth(cohort, V1)]
pc_m_former = glm(ecig ~ cohort, hse[smok_status == "3"], family = binomial, weights = wt_int)
lines(x = -80:10, y = predict(pc_m_former, newdata = data.frame('cohort' = -80:10), type = "response"), col = 'red')


hse[smok_status == "1", weighted.mean(ecig, wt_int), keyby = cohort][, scatter.smooth(cohort, V1)]
pc_m_never = glm(ecig ~ poly(cohort,2), hse[smok_status == "1"], family = binomial, weights = wt_int)
lines(x = -80:10, y = predict(pc_m_never, newdata = data.frame('cohort' = -80:10), type = "response"), col = 'green')

pc_m1.3 <- glm(ecig ~ poly(cohort,3) + year, hse[smok_status == "4"], family = binomial, weights = wt_int)

pc_m1 <- glm(ecig ~ poly(cohort,3) + year + sex + qimd, hse[smok_status == "4"], family = binomial, weights = wt_int)
pc_m1.2 <- glm(ecig ~ poly(cohort,3) + year + sex , hse[smok_status == "4"], family = binomial, weights = wt_int)
pc_m1.2.1 <- glm(ecig ~ poly(cohort,3) + year + qimd , hse[smok_status == "4"], family = binomial, weights = wt_int)
pc_m1.2.2 <- glm(ecig ~ poly(cohort,3) + sex + qimd , hse[smok_status == "4"], family = binomial, weights = wt_int)
pc_m1.3.1 <- glm(ecig ~ poly(cohort,3) + qimd, hse[smok_status == "4"], family = binomial, weights = wt_int)
pc_m1.3.2 <- glm(ecig ~ poly(cohort,3) + sex, hse[smok_status == "4"], family = binomial, weights = wt_int)
AIC(pc_m1, pc_m1.3, pc_m1.3.2, pc_m1.3.1, pc_m1.2, pc_m1.2.1, pc_m1.2.2)

pc_m1.3_1 = glm(ecig ~ poly(cohort,2) + year, hse[smok_status == "4"], family = binomial, weights = wt_int)
AIC(pc_m1.3_1, pc_m1.3)
pc_m_current = pc_m1.3 # ecig ~ poly(cohort,3) + year 
summary(pc_m_current)

pc_m3.1 <- glm(ecig ~ cohort + year + sex + qimd, hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)

pc_m3.2 <- glm(ecig ~ cohort + year + sex , hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
pc_m3.2.1 <- glm(ecig ~ cohort + year + qimd , hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
pc_m3.2.2 <- glm(ecig ~ cohort + sex + qimd , hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
pc_m3.3 <- glm(ecig ~ cohort + year, hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
pc_m3.3.1 <- glm(ecig ~ cohort + qimd, hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
pc_m3.3.2 <- glm(ecig ~ year + sex, hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
pc_m3.3.3 <- glm(ecig ~ year + qimd, hse[smok_status %in% c("3", "2")], family = binomial, weights = wt_int)
AIC(pc_m3.1, pc_m3.2, pc_m3.3, pc_m3.2.1, pc_m3.2.2, pc_m3.3.3, pc_m3.3.2, pc_m3.3.1)

pc_m_former = pc_m3.1 # ecig ~  cohort + year  + sex + qimd
summary(pc_m_former)

pc_m2.3.2_1 = glm(ecig ~ cohort + sex, hse[smok_status == "1"], family = binomial, weights = wt_int)

pc_m2 <- glm(ecig ~ poly(cohort, 2) + year + sex + qimd, hse[smok_status == "1"], family = binomial, weights = wt_int)
pc_m2.2 <- glm(ecig ~ poly(cohort, 2) + year + sex , hse[smok_status == "1"], family = binomial, weights = wt_int)
pc_m2.2.1 <- glm(ecig ~ poly(cohort, 2) + year + qimd , hse[smok_status == "1"], family = binomial, weights = wt_int)
pc_m2.2.2 <- glm(ecig ~ poly(cohort, 2) + sex + qimd , hse[smok_status == "1"], family = binomial, weights = wt_int)
pc_m2.3 <- glm(ecig ~ poly(cohort, 2) + year, hse[smok_status == "1"], family = binomial, weights = wt_int)
pc_m2.3.1 <- glm(ecig ~ poly(cohort, 2) + qimd, hse[smok_status == "1"], family = binomial, weights = wt_int)
pc_m2.3.2 <- glm(ecig ~ poly(cohort, 2) + sex, hse[smok_status == "1"], family = binomial, weights = wt_int)
pc_m2.3.3 <- glm(ecig ~ year + qimd, hse[smok_status == "1"], family = binomial, weights = wt_int)
AIC(pc_m2, pc_m2.2, pc_m2.3, pc_m2.2.1, pc_m2.2.2, pc_m2.3.1, pc_m2.3.2, pc_m2.3.3)
summary(pc_m2.2)
summary(pc_m2.3.2)
pc_m2.2_1 = glm(ecig ~ cohort + year + sex , hse[smok_status == "1"], family = binomial, weights = wt_int)
AIC(pc_m2.2, pc_m2.3.2, pc_m2.2_1, pc_m2.3.2_1)

pc_m_never = pc_m2.3.2_1 # ecig ~ poly(cohort, 2) + sex
summary(pc_m_never)

fit_pc = CJ(year = c(16:73), sex = unique(smoke_dt$sex), qimd = unique(smoke_dt$qimd), cohort = -100: 100)
fit_pc_current = copy(fit_pc)
fit_pc_current = fit_pc_current[, c("ecig", "se", "scale") := predict.glm(pc_m_current, newdata = fit_pc_current[, .(cohort, year)], type = "response", se.fit = TRUE)]
fit_pc_current[, smok_status := "4"]
# fit_pc_current[, mean(ecig), keyby = year][, plot(year, V1)]
fwrite(fit_pc_current, file = "./vape/data/model_current_pc.csv")


fit_pc_former = copy(fit_pc)
fit_pc_former = fit_pc_former[, c("ecig", "se", "scale") := predict.glm(pc_m_former, newdata = fit_pc_former[, .(cohort, year, sex, qimd)], type = "response", se.fit = TRUE)]
fit_pc_former[, smok_status := "3"]
# fit_pc_former[, mean(ecig), keyby = year][, plot(year, V1)]
fwrite(fit_pc_former, file = "./vape/data/model_former_pc.csv")


fit_pc_never = copy(fit_pc)
fit_pc_never = fit_pc_never[, c("ecig", "se", "scale") := predict.glm(pc_m_never, newdata = fit_pc_never[, .(cohort, sex)], type = "response", se.fit = TRUE)]
fit_pc_never[, smok_status := "1"]
# fit_pc_never[, mean(ecig), keyby = year][, plot(year, V1)]
fwrite(fit_pc_never, file = "./vape/data/model_never_pc.csv")


