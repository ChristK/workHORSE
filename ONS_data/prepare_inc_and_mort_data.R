## workHORSE is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - workHORSE: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## workHORSE is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

library(data.table)
library(CKutils)
library(fst)
library(ungroup)
library(ggplot2)
source("./preparatory_work/aux_functions.R") # for absorb_dt

# Create a file of mortality by agegroup ----
pop <- fread("./ONS_data/population_ONS.csv")
pop[, age := as.integer(gsub("\\+", "", age))]
to_agegrp(pop, 5, 90)
pop <- pop[, lapply(.SD, sum), by = .(sex, qimd, year, agegrp)][, age := NULL]
deaths <- fread("./ONS_data/all_deaths_ONS.csv")
deaths[, agegrp := factor(agegroup, agegrp_name(0, 90), agegrp_name(0, 90))]
ttt <- fread("./ONS_data/cause_specific_deaths_ONS.csv")
ttt[, agegroup := factor(agegroup, agegrp_name(0, 90), agegrp_name(0, 90))]
ttt[grepl("^C18", cause), cause := "C18"]
ttt[grepl("^C19", cause), cause := "C19"]
ttt[grepl("^C34", cause), cause := "C34"]
ttt[grepl("^C50", cause), cause := "C50"]
ttt[grepl("^E11", cause), cause := "E11"]
ttt[grepl("^F0", cause), cause := "F00_03"]
ttt[grepl("^I2", cause), cause := "I20_25"]
ttt[grepl("^I48", cause), cause := "I48"]
ttt[grepl("^I6", cause), cause := "I60_69"]
ttt[grepl("^J4", cause), cause := "J40_47"]
ttt <- ttt[, .(deaths = sum(deaths)), by = .(year, sex, cause, qimd, agegrp = agegroup)]
ttt <- dcast(ttt, year+sex+qimd+agegrp~cause, value.var = "deaths", drop = FALSE, fill = 0L)
absorb_dt(deaths, ttt)
absorb_dt(deaths, pop)
deaths[, sex := factor(sex, 1:2, c("men", "women"))]
deaths[, qimd := factor(
  qimd,
  levels = 1:5,
  labels = c("1 most deprived", "2", "3", "4", "5 least deprived")
)]
fwrite(deaths, "./ONS_data/mortality_by_agegroup.csv") # messes up agegrp level order
write_fst(deaths, "./ONS_data/mortality_by_agegroup.fst")
rm(deaths, pop, ttt)





# Ungroup agegroup bins to get 1-year estimates ----
# by default estimates up to age 100 are calculated
max_age <- 100 # collapses ages older than max_age to max_age

pop <- fread("./ONS_data/population_ONS.csv")
pop[, age := as.integer(gsub("\\+", "", age))]

setkey(pop, year, age)
ttt <- data.table(NULL)
for (i in unique(pop$sex)) {
  for (j in unique(pop$qimd)) {
    tt <- pop[sex == i & qimd == j, .(age, pops, year)]
    tt <- dcast(tt, age~year, value.var = "pops")
    P1 <- pclm2D(tt[, age], setDF(tt[, -1]), 11, control = list(lambda = c(NA, NA), opt.method = "BIC"))
    tt <- P1$fit
    # tt[] <- as.integer(round(tt))
    tt <- data.table(tt)
    tt[, `:=` (age = 0:100, sex = i, qimd = j)]
    tt <- melt(tt, c("age", "sex", "qimd"), variable.name = "year", value.name = "pops")
    ttt <- rbind(ttt, tt)
  }
}
ttt[, year := as.integer(as.character(year))]
fwrite(ttt, "./ONS_data/pop_smoothed.csv")

ttt[pop[age < 90, ], pops := as.numeric(i.pops), on = c("age", "sex", "qimd", "year")] # replace ages < 90 with observed values

ttt[age > max_age, age := max_age]

pop <- copy(ttt)
rm(ttt)
pop[, sex := factor(sex, 1:2, c("men", "women"))]
pop[, qimd := factor(
  qimd,
  levels = 1:5,
  labels = c("1 most deprived", "2", "3", "4", "5 least deprived")
)]
fwrite(pop, "./ONS_data/pop.csv")

# ggplot(pop[
#   (age %% 10) == 5 & between(age, 20, 100) & year > 2000, ],
#   aes(
#     y = pops,
#     x = year
#   )) +
#   geom_point() +
#   facet_grid(age~sex + qimd, scales = "free") # , scales = "free"


# Mortality ----
deaths <- fread("./ONS_data/all_deaths_ONS.csv")
deaths[, agegroup := factor(agegroup, agegrp_name(0, 90), agegrp_name(0, 90))]

lkp_age <- data.table(age = 0:90)
to_agegrp(lkp_age, 5, 90, "age", "agegroup")

x1 <- fread("~/pCloudDrive/My Models/workHORSE/ONS_data/population_ONS.csv")
x1[, age := as.integer(gsub("\\+", "", age))]
to_agegrp(x1, 5L, 90, "age", "agegroup")
x1 <- x1[, .(pops = sum(pops)), keyby = .(year, agegroup, sex, qimd)]


setkey(deaths, year, agegroup)
ttt <- data.table(NULL)
for (i in unique(deaths$sex)) {
  for (j in unique(deaths$qimd)) {
    tt <- deaths[sex == i & qimd == j, .(agegroup, deaths, year)]
    tt[lkp_age[, first(age), keyby = agegroup], age := i.V1, on = "agegroup"]
    tt <- dcast(tt, age~year, value.var = "deaths")
    # x2 <- x1[sex == i & qimd == j, .(agegroup, pops, year)]
    # x2[lkp_age[, first(age), keyby = agegroup], age := i.V1, on = "agegroup"]
    # x2 <- dcast(x2, age~year, value.var = "pops")[, age := NULL]

    P1 <- pclm2D(tt[, age], setDF(tt[, -1]), 11, control = list(lambda = c(NA, NA), opt.method = "AIC"))
    tt <- P1$fit
    # tt[] <- as.integer(round(tt))
    tt <- data.table(tt)
    tt[, `:=` (age = 0:100, sex = i, qimd = j)]
    tt <- melt(tt, c("age", "sex", "qimd"), variable.name = "year", value.name = "all_deaths")
    ttt <- rbind(ttt, tt)
  }
}
ttt[age > max_age, age := max_age]
ttt[, year := as.integer(as.character(year))]

deaths <- ttt[, .("deaths" = sum(all_deaths)), keyby = .(age, sex, qimd, year)]
deaths[, sex := factor(sex, 1:2, c("men", "women"))]
deaths[, qimd := factor(
  qimd,
  levels = 1:5,
  labels = c("1 most deprived", "2", "3", "4", "5 least deprived")
)]
deaths[, year := as.integer(as.character(year))]

deaths[pop, pops := i.pops, on = c("age", "sex", "qimd", "year")]

ttt <- fread("~/pCloudDrive/My Models/workHORSE/ONS_data/cause_specific_deaths_ONS.csv")
ttt[, agegroup := factor(agegroup, agegrp_name(0, 90), agegrp_name(0, 90))]
ttt[grepl("^C18", cause), cause := "C18"]
ttt[grepl("^C19", cause), cause := "C19"]
ttt[grepl("^C34", cause), cause := "C34"]
ttt[grepl("^C50", cause), cause := "C50"]
ttt[grepl("^E11", cause), cause := "E11"]
ttt[grepl("^F0", cause), cause := "F00-3"]
ttt[grepl("^I2", cause), cause := "I20-5"]
ttt[grepl("^I48", cause), cause := "I48"]
ttt[grepl("^I6", cause), cause := "I60-9"]
ttt[grepl("^J4", cause), cause := "J40-7"]
ttt <- ttt[, .(deaths = sum(deaths)), by = .(year, sex, cause, qimd, agegroup)]

tt <- CJ(agegroup = unique(ttt$agegroup),
          sex = unique(ttt$sex),
          qimd = unique(ttt$qimd),
          year = unique(ttt$year),
          cause = unique(ttt$cause),
          deaths = 0L)
tt[ttt, deaths := i.deaths, on = c("agegroup", "sex", "qimd", "year", "cause")]

setkey(tt, year, agegroup)

tttt <- data.table(NULL)
for (i in unique(tt$sex)) {
  for (j in unique(tt$qimd)) {
    ttt <- CJ(age = 0:100,
              sex = unique(tt$sex),
              qimd = unique(tt$qimd),
              year = unique(tt$year))
    for (k in unique(tt$cause)) {
      print(paste0(i, "_", j, "_",k))
      tt1 <- tt[sex == i & qimd == j & cause == k, .(agegroup, deaths, year)]
      tt1[lkp_age[, first(age), keyby = agegroup], age := i.V1, on = "agegroup"]
      tt1 <- dcast(tt1, age~year, value.var = "deaths")
      P1 <- pclm2D(tt1[, age], setDF(tt1[, -1]), 11, control = list(lambda = c(NA, NA), opt.method = "AIC"))
      tt1 <- P1$fit
      # tt1[] <- as.integer(round(tt1))
      tt1 <- data.table(tt1)
      tt1[, `:=` (age = 0:100, sex = i, qimd = j)]
      tt1 <- melt(tt1, c("age", "sex", "qimd"), variable.name = "year",
                  value.name = paste0(k, "_deaths"))
      tt1[, year := as.integer(as.character(year))]
      tt1[age > max_age, age := max_age]
      tt1 <- tt1[, .(sum(get(paste0(k, "_deaths")))),
                 keyby = .(age, sex, qimd, year)]
      setnames(tt1, "V1", paste0(k, "_deaths"))
      ttt <- ttt[tt1, on = c("age", "sex", "qimd", "year")]
    }
    tttt <- rbind(tttt, ttt)
  }
}

tttt[, sex := factor(sex, 1:2, c("men", "women"))]
tttt[, qimd := factor(
  qimd,
  levels = 1:5,
  labels = c("1 most deprived", "2", "3", "4", "5 least deprived")
)]

deaths[tttt, on = c("age", "sex", "qimd", "year")]
absorb_dt(deaths, tttt, c("age", "sex", "qimd", "year"))
fwrite(deaths, "./ONS_data/deaths.csv")


# Cancer incidence ----
ttt <- fread("./ONS_data/cancer_incidence_ONS.csv")
ttt[, agegroup := factor(agegroup, c("0-44", agegrp_name(45, 90)),
                         c("00-44", agegrp_name(45, 90)))]
ttt[grepl("^C18", cause), cause := "C18"]
ttt[grepl("^C19", cause), cause := "C19"]
ttt[grepl("^C34", cause), cause := "C34"]
ttt[grepl("^C50", cause), cause := "C50"]

ttt <- ttt[, .(cases = sum(incidence)), by = .(year, sex, cause, qimd, agegroup)]

tt <- CJ(agegroup = unique(ttt$agegroup),
         sex = unique(ttt$sex),
         qimd = unique(ttt$qimd),
         year = unique(ttt$year),
         cause = unique(ttt$cause),
         cases = 0L)
tt[ttt, cases := i.cases, on = c("agegroup", "sex", "qimd", "year", "cause")]

setkey(tt, year, agegroup)

lkp_age <- data.table(agegroup = unique(tt$agegroup), age = c(0, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90))

tttt <- data.table(NULL)
for (i in unique(tt$sex)) {
  for (j in unique(tt$qimd)) {
    ttt <- CJ(age = 0:100,
              sex = unique(tt$sex),
              qimd = unique(tt$qimd),
              year = unique(tt$year))
    for (k in unique(tt$cause)) {
      print(paste0(i, "_", j, "_",k))
      tt1 <- tt[sex == i & qimd == j & cause == k, .(agegroup, cases, year)]
      tt1[lkp_age, age := i.age, on = "agegroup"]
      tt1 <- dcast(tt1, age~year, value.var = "cases")
      P1 <- pclm2D(tt1[, age], setDF(tt1[, -1]), 11, control = list(lambda = c(NA, NA), opt.method = "AIC"))
      tt1 <- P1$fit
      # tt1[] <- as.integer(round(tt1))
      tt1 <- data.table(tt1)
      tt1[, `:=` (age = 0:100, sex = i, qimd = j)]
      tt1 <- melt(tt1, c("age", "sex", "qimd"), variable.name = "year",
                  value.name = paste0(k, "_cases"))
      tt1[, year := as.integer(as.character(year))]
      tt1[age > max_age, age := max_age]
      tt1 <- tt1[, .(sum(get(paste0(k, "_cases")))),
                 keyby = .(age, sex, qimd, year)]
      setnames(tt1, "V1", paste0(k, "_cases"))
      ttt <- ttt[tt1, on = c("age", "sex", "qimd", "year")]
    }
    tttt <- rbind(tttt, ttt)
  }
}

tttt[, sex := factor(sex, 1:2, c("men", "women"))]
tttt[, qimd := factor(
  qimd,
  levels = 1:5,
  labels = c("1 most deprived", "2", "3", "4", "5 least deprived")
)]

absorb_dt(tttt, pop, c("age", "sex", "qimd", "year"))
fwrite(tttt, "./ONS_data/cancer_inc.csv")

# Cancer Incidence for validation 5 yr agegroups
ttt <- fread("./ONS_data/cancer_incidence_ONS.csv")
ttt[, agegroup := factor(agegroup, c("0-44", agegrp_name(45, 90)),
  c("00-44", agegrp_name(45, 90)))]
ttt[grepl("^C18", cause), cause := "C18"]
ttt[grepl("^C34", cause), cause := "C34"]
ttt[grepl("^C50", cause), cause := "C50"]

ttt <- ttt[, .(cases = sum(incidence)), by = .(year, sex, cause, qimd, agegroup)]
ttt[, sex := factor(sex, 1:2, c("men", "women"))]
ttt[, qimd := factor(
  qimd,
  levels = 1:5,
  labels = c("1 most deprived", "2", "3", "4", "5 least deprived")
)]
ttt[cause == "C18", cause := "colon_ca"]
ttt[cause == "C34", cause := "lung_ca"]
ttt[cause == "C50", cause := "breast_ca"]
ttt <- dcast(ttt, agegroup + sex + qimd + year ~ cause)
ttt[, c("C19", "C20", "C33") := NULL]
setnames(ttt, c("colon_ca", "lung_ca", "breast_ca"),
  c("cases_colon_ca", "cases_lung_ca", "cases_breast_ca"))
setnafill(ttt, "c", 0L, cols = c("cases_colon_ca", "cases_lung_ca", "cases_breast_ca"))
ttt <- ttt[!agegroup == "90+"]
ttt[, unique(agegroup)]


pop <- fread("./ONS_data/pop.csv")[between(age, 0, 89)]
to_agegrp(pop, 5L, 89, "age", "agegroup")
pop[, unique(agegroup)]
pop <- pop[, .(pops = sum(pops)), by = .(year, sex, qimd, agegroup)]
CKutils::absorb_dt(ttt, pop)
ttt <- na.omit(ttt)
ttt[, `:=` (
  rate_breast_ca = cases_breast_ca/pops,
  rate_colon_ca = cases_colon_ca/pops,
  rate_lung_ca = cases_lung_ca/pops
)]
setnames(ttt, "agegroup", "agegrp5")
write_fst(ttt, "./ONS_data/cancer_incd_by_agegrp5_for_validation.fst")

# Prepare data for validation ----
# Prepare mortality file for validation
deaths <- read_fst("./ONS_data/mortality_by_agegroup.fst", as.data.table = TRUE)
deaths[, c("C19", "C20", "C33", "E11", "F00_03", "I48", "N390", "agegrp") := NULL]
deaths <- deaths[agegroup %in% c("30-34", "35-39", "40-44",
  "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89")]

deaths <- groupingsets(
  deaths,
  j = lapply(.SD, sum),
  by = c("year", "sex", "agegroup", "qimd"),
  # .SDcols = xps,
  sets = list(c("year", "sex", "agegroup", "qimd"),
    c("year", "sex"), c("year", "agegroup"), c("year", "qimd"))
)
for (j in seq_len(ncol(deaths)))
  set(deaths, which(is.na(deaths[[j]])), j, "All")

deaths[, Mx_all       := deaths / pops]
deaths[, Mx_chd       := I20_25 / pops]
deaths[, Mx_stroke    := I60_69 / pops]
deaths[, Mx_colon_ca  := C18 / pops]
deaths[, Mx_lung_ca   := C34 / pops]
deaths[, Mx_breast_ca := C50 / pops]
deaths[, Mx_copd      := J40_47 / pops]
deaths[, Mx_nonmodelled := (deaths - I20_25 - I60_69 - C18 -
    C34 - C50 - J40_47) / pops]
deaths[, deaths_nonmodelled := deaths - I20_25 - I60_69 - C18 -
    C34 - C50 - J40_47]
setnames(deaths,
  c("deaths", "C18", "C34", "C50", "I20_25", "I60_69", "J40_47", "agegroup"),
  c("deaths_all", "deaths_colon_ca", "deaths_lung_ca", "deaths_breast_ca", "deaths_chd", "deaths_stroke", "deaths_copd", "agegrp5"))
deaths[, agegrp5 := factor(agegrp5)]
write_fst(deaths, "./ONS_data/mortality_by_agegrp5_for_validation.fst")

# 20-yr agegrpups
deaths <- read_fst("./ONS_data/mortality_by_agegroup.fst", as.data.table = TRUE)
deaths[, c("C19", "C20", "C33", "E11", "F00_03", "I48", "N390") := NULL]
deaths[agegroup %in% c("30-34", "35-39", "40-44", "45-49"), agegrp20 := "30-49"]
deaths[agegroup %in% c("50-54", "55-59", "60-64", "65-69"), agegrp20 := "50-69"]
deaths[agegroup %in% c("70-74", "75-79", "80-84", "85-89"), agegrp20 := "70-89"]
deaths[, agegrp20 := factor(agegrp20)]
deaths <- deaths[!is.na(agegrp20)]
deaths[, c("agegroup", "agegrp") := NULL]
deaths <- deaths[, lapply(.SD, sum), by = .(agegrp20, year, sex, qimd)]
deaths <- groupingsets(
  deaths,
  j = lapply(.SD, sum),
  by = c("year", "sex", "agegrp20", "qimd"),
  # .SDcols = xps,
  sets = list(c("year", "sex", "agegrp20", "qimd"),
    c("year", "sex"), c("year", "agegrp20"), c("year", "qimd"))
)
for (j in seq_len(ncol(deaths)))
  set(deaths, which(is.na(deaths[[j]])), j, "All")
deaths[, Mx_all       := deaths / pops]
deaths[, Mx_chd       := I20_25 / pops]
deaths[, Mx_stroke    := I60_69 / pops]
deaths[, Mx_colon_ca  := C18 / pops]
deaths[, Mx_lung_ca   := C34 / pops]
deaths[, Mx_breast_ca := C50 / pops]
deaths[, Mx_copd      := J40_47 / pops]
deaths[, Mx_nonmodelled := (deaths - I20_25 - I60_69 - C18 -
    C34 - C50 - J40_47) / pops]
deaths[, deaths_nonmodelled := deaths - I20_25 - I60_69 - C18 -
    C34 - C50 - J40_47]
setnames(deaths,
  c("deaths", "C18", "C34", "C50", "I20_25", "I60_69", "J40_47"),
  c("deaths_all", "deaths_colon_ca", "deaths_lung_ca", "deaths_breast_ca", "deaths_chd", "deaths_stroke", "deaths_copd"))
write_fst(deaths, "./ONS_data/mortality_by_agegrp20_for_validation.fst")

rm(deaths)


# Cancer Incidence for validation 20 yr agegroups
ttt <- fread("./ONS_data/cancer_incidence_ONS.csv")
ttt[, agegroup := factor(agegroup, c("0-44", agegrp_name(45, 90)),
  c("00-44", agegrp_name(45, 90)))]
ttt[grepl("^C18", cause), cause := "C18"]
ttt[grepl("^C34", cause), cause := "C34"]
ttt[grepl("^C50", cause), cause := "C50"]

ttt <- ttt[, .(cases = sum(incidence)), by = .(year, sex, cause, qimd, agegroup)]
ttt[, sex := factor(sex, 1:2, c("men", "women"))]
ttt[, qimd := factor(
  qimd,
  levels = 1:5,
  labels = c("1 most deprived", "2", "3", "4", "5 least deprived")
)]
ttt[cause == "C18", cause := "colon_ca"]
ttt[cause == "C34", cause := "lung_ca"]
ttt[cause == "C50", cause := "breast_ca"]
ttt <- dcast(ttt, agegroup + sex + qimd + year ~ cause)
ttt[, c("C19", "C20", "C33") := NULL]
setnames(ttt, c("colon_ca", "lung_ca", "breast_ca"),
  c("cases_colon_ca", "cases_lung_ca", "cases_breast_ca"))
setnafill(ttt, "c", 0L, cols = c("cases_colon_ca", "cases_lung_ca", "cases_breast_ca"))
ttt <- ttt[!agegroup == "90+"]
ttt[agegroup %in% c("00-44", "45-49"), agegrp20 := "30-49"]
ttt[agegroup %in% c("50-54", "55-59", "60-64", "65-69"), agegrp20 := "50-69"]
ttt[agegroup %in% c("70-74", "75-79", "80-84", "85-89"), agegrp20 := "70-89"]

ttt[, unique(agegrp20)]
ttt[, agegroup := NULL]
ttt <- ttt[, lapply(.SD, sum), by = .(agegrp20, year, sex, qimd)]

pop <- fread("./ONS_data/pop.csv")[between(age, 30, 89)]
to_agegrp(pop, 20L, 89, "age", "agegrp20")
pop[, unique(agegrp20)]
pop <- pop[, .(pops = sum(pops)), by = .(year, sex, qimd, agegrp20)]
CKutils::absorb_dt(ttt, pop)
ttt <- na.omit(ttt)

ttt <- groupingsets(
  ttt,
  j = lapply(.SD, sum),
  by = c("year", "sex", "agegrp20", "qimd"),
  # .SDcols = xps,
  sets = list(c("year", "sex", "agegrp20", "qimd"),
    c("year", "sex"), c("year", "agegrp20"), c("year", "qimd"))
)
for (j in seq_len(ncol(ttt)))
  set(ttt, which(is.na(ttt[[j]])), j, "All")

ttt[, `:=` (
  rate_breast_ca = cases_breast_ca/pops,
  rate_colon_ca = cases_colon_ca/pops,
  rate_lung_ca = cases_lung_ca/pops
)]
write_fst(ttt, "./ONS_data/cancer_incd_by_agegrp20_for_validation.fst")

# get_ons_incd("lung_ca", "rate", 20)

# ggplot(deaths[
#   (age %% 10) == 5 & between(age, 20, 100) & year > 2000, ],
#   aes(
#     y = deaths/pops,
#     x = year
#   )) +
#   geom_point() +
#   facet_grid(age~sex + qimd, scales = "free") # , scales = "free"
