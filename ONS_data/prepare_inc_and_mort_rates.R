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
source("./preparatory_work/aux_functions.R") # for absorb_dt

# by default estimates up to age 100 are calculated
max_age <- 100 # collapses ages older than max_age to max_age

if (file.exists("./ONS_data/pop.csv")) {
  pop <- fread("./ONS_data/pop.csv", stringsAsFactors = TRUE)
} else {
  pop <-
    fread("./ONS_data/population_ONS.csv")
  pop[, age := as.integer(gsub("\\+", "", age))]

  setkey(pop, year, age)
  ttt <- data.table(NULL)
  for (i in unique(pop$sex)) {
    for (j in unique(pop$qimd)) {
      tt <- pop[sex == i & qimd == j, .(age, pops, year)]
      tt <- dcast(tt, age ~ year, value.var = "pops")
      P1 <-
        pclm2D(tt[, age], setDF(tt[, -1]), 11, control = list(lambda = c(NA, NA)))
      tt <- P1$fit
      # tt[] <- as.integer(round(tt))
      tt <- data.table(tt)
      tt[, `:=` (age = 0:100,
                 sex = i,
                 qimd = j)]
      tt <-
        melt(tt,
             c("age", "sex", "qimd"),
             variable.name = "year",
             value.name = "pops")
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
}

ggplot(pop[
  (age %% 10) == 5 & between(age, 20, 100) & year > 2000, ],
  aes(
    y = pops,
    x = year
  )) +
  geom_point() +
  facet_grid(age~sex + qimd, scales = "free") # , scales = "free"


# Mortality ----
deaths <- fread("~/pCloudDrive/My Models/workHORSE/ONS_data/all_deaths_ONS.csv")
deaths[, agegroup := factor(agegroup, agegrp_name(0, 90), agegrp_name(0, 90))]

lkp_age <- data.table(age = 0:90)
to_agegrp(lkp_age, 5, 90, "age", "agegroup")

x1 <- fread("~/pCloudDrive/My Models/workHORSE/ONS_data/population_ONS.csv")
x1[, age := as.integer(gsub("\\+", "", age))]
to_agegrp(x1, 5L, 90, "age", "agegroup")
x1 <- x1[, .(pops = sum(pops)), keyby = .(year, agegroup, sex, qimd)]

pop2 <- copy(pop)
pop2[, `:=` (sex = as.integer(sex), qimd = as.integer(qimd))]

setkey(deaths, year, agegroup)
ttt <- data.table(NULL)
for (i in unique(deaths$sex)) {
  for (j in unique(deaths$qimd)) {
    print(paste0(i, "_", j))
    tt <- deaths[sex == i & qimd == j, .(agegroup, deaths, year)]
    tt[lkp_age[, first(age), keyby = agegroup], age := i.V1, on = "agegroup"]
    tt <- dcast(tt, age~year, value.var = "deaths")
    x2 <- x1[sex == i & qimd == j, .(agegroup, pops, year)]
    x2[lkp_age[, first(age), keyby = agegroup], age := i.V1, on = "agegroup"]
    x2 <- setDF(dcast(x2, age~year, value.var = "pops")[, age := NULL])

    P1 <- pclm2D(tt[, age], setDF(tt[, -1]), 11, x2,
                 control = list(lambda = c(NA, NA), opt.method = "BIC"))
    tt <- P1$fit
    # tt[] <- as.integer(round(tt))
    tt <- data.table(tt)
    tt[, `:=` (age = 0:100, sex = i, qimd = j)]
    tt <- melt(tt, c("age", "sex", "qimd"), variable.name = "year", value.name = "mx_all")
    tt[, year := as.integer(as.character(year))]
    tt[pop2[sex == i &
              qimd == j, ], `:=` (deaths = mx_all * i.pops, pops = i.pops), on = .(age, sex, qimd, year)]
    ttt <- rbind(ttt, tt)
  }
}
ttt[age > max_age, age := max_age]

deaths <- ttt[, .("deaths" = sum(deaths), "pops" = sum(pops),
                  "mx_all" = weighted.mean(mx_all, pops)),
              keyby = .(age, sex, qimd, year)]
deaths[, sex := factor(sex, 1:2, c("men", "women"))]
deaths[, qimd := factor(
  qimd,
  levels = 1:5,
  labels = c("1 most deprived", "2", "3", "4", "5 least deprived")
)]


ggplot(deaths[
  between(age, 20, 100), ],
  aes(
    y = mx_all,
    x = age,
    col = factor(year)
  )) +
  geom_line() +
  facet_grid(sex~qimd, scales = "free") # , scales = "free"

ggplot(deaths[
  (age %% 10) == 0 & between(age, 20, 120) & year > 2000, ],
  aes(
    y = mx_all,
    x = year
  )) +
  geom_point() +
  facet_grid(age~sex + qimd, scales = "free") # , scales = "free"

ggplot(deaths[
  (age %% 10) == 5 & between(age, 20, 100) & year > 2000, ],
  aes(
    y = deaths,
    x = year
  )) +
  geom_point() +
  facet_grid(age~sex + qimd, scales = "free") # , scales = "free"

ttt <- fread("~/pCloudDrive/My Models/workHORSE/ONS_data/cause_specific_deaths_ONS.csv")
ttt[, agegroup := factor(agegroup, agegrp_name(0, 90), agegrp_name(0, 90))]
ttt[grepl("^C18", cause), cause := "colon_ca"]
ttt[grepl("^C19", cause), cause := "C19"]
ttt[grepl("^C34", cause), cause := "lung_ca"]
ttt[grepl("^C50", cause), cause := "breast_ca"]
ttt[grepl("^E11", cause), cause := "E11"]
ttt[grepl("^F0", cause), cause := "F00_03"]
ttt[grepl("^I2", cause), cause := "chd"]
ttt[grepl("^I48", cause), cause := "I48"]
ttt[grepl("^I6", cause), cause := "stroke"]
ttt[grepl("^J4", cause), cause := "copd"]
ttt <- ttt[cause %in% c("colon_ca", "lung_ca", "breast_ca", "chd", "stroke", "copd"),
           .(deaths = sum(deaths)), by = .(year, sex, cause, qimd, agegroup)]

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
      x2 <- x1[sex == i & qimd == j, .(agegroup, pops, year)]
      x2[lkp_age[, first(age), keyby = agegroup], age := i.V1, on = "agegroup"]
      x2 <- setDF(dcast(x2, age~year, value.var = "pops")[, age := NULL])


      P1 <- pclm2D(tt1[, age], setDF(tt1[, -1]), 11, x2,
                   control = list(lambda = c(NA, NA), opt.method = "BIC"))
      tt1 <- P1$fit
      # tt1[] <- as.integer(round(tt1))
      tt1 <- data.table(tt1)
      tt1[, `:=` (age = 0:100, sex = i, qimd = j)]
      tt1 <- melt(tt1, c("age", "sex", "qimd"), variable.name = "year",
                  value.name = paste0("mx_", k))
      tt1[, year := as.integer(as.character(year))]
      tt1[get(paste0("mx_", k)) < 0, (paste0("mx_", k)) := 0]
      tt1[pop2[sex == i &
                 qimd == j, ],
          c(paste0("deaths_", k), "pops") := .(get(paste0("mx_", k)) * i.pops, i.pops),
          on = .(age, sex, qimd, year)]

      ggplot(tt1[
        (age %% 10) == 5 & between(age, 20, 100) & year > 2000, ],
        aes(
          y = mx_colon_ca,
          x = year
        )) +
        geom_point() +
        facet_grid(age~sex + qimd, scales = "free") # , scales = "free"


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
ttt <- fread("~/pCloudDrive/My Models/workHORSE/ONS_data/cancer_incidence_ONS.csv")
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
    ttt <- CJ(age = 0:90,
              sex = unique(tt$sex),
              qimd = unique(tt$qimd),
              year = unique(tt$year))
    for (k in unique(tt$cause)) {
      print(paste0(i, "_", j, "_",k))
      tt1 <- tt[sex == i & qimd == j & cause == k, .(agegroup, cases, year)]
      tt1[lkp_age, age := i.age, on = "agegroup"]
      tt1 <- dcast(tt1, age~year, value.var = "cases")
      P1 <- pclm2D(tt1[, age], setDF(tt1[, -1]), 11, control = list(lambda = c(NA, NA)))
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
