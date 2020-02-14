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

f <- fread("./ONS_data/pop_size/2016 SNPP Population females.csv", header = T)
setnames(f, tolower(names(f)))
f[, component := NULL]
f[, sex := "women"]
f[age_group == "90 and over", age_group := 90]
f[, age := as.integer(age_group)] # All_ages ic converted to NA. That's desirable
f[, age_group := NULL]
f <- na.omit(f)

m <- fread("./ONS_data/pop_size/2016 SNPP Population males.csv", header = T)
setnames(m, tolower(names(m)))
m[, component := NULL]
m[, sex := "men"]
m[age_group == "90 and over", age_group := 90]
m[, age := as.integer(age_group)]
m[, age_group := NULL]
m <- na.omit(m)

pops <- rbind(m, f)
pops <- melt(pops, c("area_code", "area_name", "age", "sex"), variable.name = "year", value.name = "pops")
setnames(pops, c("area_code", "area_name"), c("LAD17CD", "LAD17NM"))
pops[, `:=` (LAD17CD = factor(LAD17CD),
             LAD17NM = factor(LAD17NM),
             sex = factor(sex),
             year = as.integer(as.character(year)),
             age = as.integer(as.character(age))
)]

# add years from 2013 to 2017
pops <- pops[year > 2017L]
localities_indx <- read_fst("./synthpop/lsoa_to_locality_indx.fst", as.data.table = TRUE)
dt <- read_fst("./synthpop/lsoa_mid_year_population_estimates.fst",
               as.data.table = TRUE)
dt[localities_indx, on = "LSOA11CD", `:=`(LAD17CD = i.LAD17CD, LAD17NM = i.LAD17NM)]
dt[, c("LSOA11CD", "LAD11CD", "LAD11NM") := NULL]
dt <-
  melt(
    dt,
    grep("^[0-9]", names(dt), value = TRUE, invert = TRUE),
    variable.name = "age",
    value.name = "pops",
    variable.factor = FALSE
  )[year > 2002L]
dt <- dt[, .(pops = sum(pops)), keyby = .(year, age, sex, LAD17CD, LAD17NM)]
dt[, `:=` (LAD17CD = factor(LAD17CD),
             LAD17NM = factor(LAD17NM),
             sex = factor(sex),
             year = as.integer(as.character(year)),
             age = as.integer(as.character(age))
)]

pops <- rbind(pops, dt)



setkey(pops, year, age, sex, LAD17CD)

write_fst(pops, "./ONS_data/pop_size/pop_proj.fst", 100)
