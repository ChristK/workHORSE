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

# Builds the ONS population size files used to weight the synthetic population
# - pop_proj.fst: year, age, sex, LAD23CD, LAD23NM, pops (296 LAs, 2003-2047)
# - pop_proj_england.fst: year, age, sex, pops (England, 2003-2047)
# 2003-2010: LSOA11 mid-year estimates summed to LAD23
# 2011-2025: ONS mid-year estimates (mid-2011 to mid-2025 time series)
# 2026-2047: LAs from the 2022-based SNPP (migration category variant) and
#            England from the 2024-based NPP (principal projection), so from
#            2026 England is not the sum of its LAs
# Run from the project root after download_ons_pop.R and
# update_locality_indx_lad23.R. See ./ONS_data/pop_size/README.md

library(data.table)
library(fst)
library(ggplot2)

first_year <- 2003L # init_year_long - maxlag in sim_design.yaml
mye_from   <- 2011L # first year of the ONS mid-year estimates time series
proj_from  <- 2026L # projections start after the latest estimates (mid-2025)
last_year  <- 2047L # last year of the 2022-based SNPP

# LSOA11 mid-year estimates summed to LAD23 ----
indx <- read_fst("./synthpop/lsoa_to_locality_indx.fst",
                 columns = c("LSOA11CD", "LAD23CD"), as.data.table = TRUE)
lsoa <- read_fst("./synthpop/lsoa_mid_year_population_estimates.fst",
                 as.data.table = TRUE)[between(year, first_year, mye_from - 1L)]
lsoa[indx, on = "LSOA11CD", LAD23CD := i.LAD23CD]
lsoa <- melt(lsoa, grep("^[0-9]", names(lsoa), value = TRUE, invert = TRUE),
             variable.name = "age", value.name = "pops",
             variable.factor = FALSE)
lsoa <- lsoa[, .(pops = sum(pops)),
             keyby = .(year = as.integer(year), age = as.integer(age),
                       sex = as.character(sex),
                       LAD23CD = as.character(LAD23CD))]

# ONS extracts written by download_ons_pop.R ----
mye <- read_fst("./ONS_data/pop_size/ons_mye_2011_2025_lad23.fst",
                as.data.table = TRUE)
snpp <- read_fst("./ONS_data/pop_size/ons_snpp2022_migcat_lad23.fst",
                 as.data.table = TRUE)
npp <- read_fst("./ONS_data/pop_size/ons_npp2024_ppp_england.fst",
                as.data.table = TRUE)

# Checks on the sources against published England totals (ONS bulletins and
# Nomis NM_2002_1 for the mid-year estimates)
stopifnot(
  abs(lsoa[year == 2003L, sum(pops)] - 49925517) < 1,
  abs(lsoa[year == 2010L, sum(pops)] - 52642452) < 1,
  abs(mye[year == 2011L, sum(pops)] - 53107169) < 1,
  abs(mye[year == 2025L, sum(pops)] - 58834812) < 1,
  abs(snpp[year == 2047L, sum(pops)] / 64.4e6 - 1) < 0.005, # SNPP bulletin
  npp[year == 2047L, between(sum(pops), 58.6e6, 62.1e6)],   # NPP base and peak
  abs(mye[year == 2024L, sum(pops)] / npp[year == 2024L, sum(pops)] - 1) < 0.01
)

# Stitch ----
nam <- c("year", "age", "sex", "LAD23CD", "pops")
lad <- rbind(lsoa[, ..nam],
             mye[between(year, mye_from, proj_from - 1L), ..nam],
             snpp[between(year, proj_from, last_year), ..nam])
lad[unique(mye[, .(LAD23CD, LAD23NM)]), on = "LAD23CD", LAD23NM := i.LAD23NM]

eng <- rbind(lad[year < proj_from, .(pops = sum(pops)), keyby = .(year, age, sex)],
             npp[between(year, proj_from, last_year), .(year, age, sex, pops)])

# Every LA x year x age x sex cell exactly once, no NA, no negatives
grid <- CJ(year = first_year:last_year, age = 0:90, sex = c("men", "women"),
           LAD23CD = unique(mye$LAD23CD))
stopifnot(
  nrow(grid) == 296L * 91L * 2L * 45L,
  nrow(lad) == nrow(grid), !anyNA(lad), all(lad$pops >= 0),
  !anyDuplicated(lad, by = c("year", "age", "sex", "LAD23CD")),
  nrow(grid[!lad, on = c("year", "age", "sex", "LAD23CD")]) == 0L,
  nrow(eng) == 91L * 2L * 45L, !anyNA(eng), all(eng$pops >= 0),
  !anyDuplicated(eng, by = c("year", "age", "sex"))
)

lad[, `:=` (LAD23CD = factor(LAD23CD),
            LAD23NM = factor(LAD23NM),
            sex = factor(sex))]
eng[, sex := factor(sex)]
setcolorder(lad, c("year", "age", "sex", "LAD23CD", "LAD23NM", "pops"))
setkey(lad, year, age, sex, LAD23CD)
setkey(eng, year, age, sex)

write_fst(lad, "./ONS_data/pop_size/pop_proj.fst", 100)
write_fst(eng, "./ONS_data/pop_size/pop_proj_england.fst", 100)

# Junction diagnostics ----
# Year-on-year % change at each change of source (ages 20-89, the ages the
# model uses), against the mean annual change over the 3 preceding years.
# 'excess' isolates the discontinuity from the underlying trend.
junction_steps <- function(dt, by) {
  tt <- dt[between(age, 20L, 89L), .(pops = sum(pops)), keyby = c(by, "year")]
  tt[, chg := 100 * (pops / shift(pops) - 1), by = by]
  rbindlist(lapply(c(mye_from, proj_from), function(y) {
    tt[, .(junction = paste0(y - 1L, "->", y),
           pops_before = pops[year == y - 1L],
           pops_after = pops[year == y],
           step_pct = chg[year == y],
           trend_pct = mean(chg[between(year, y - 3L, y - 1L)])),
       by = by]
  }))[, excess_pct := step_pct - trend_pct]
}
diag <- rbind(
  junction_steps(lad[, .(year, age, LAD23CD, LAD23NM, pops)],
                 c("LAD23CD", "LAD23NM")),
  junction_steps(eng[, .(year, age, LAD23CD = "E92000001", LAD23NM = "England",
                         pops)], c("LAD23CD", "LAD23NM"))
)
diag[, flag := abs(excess_pct) > 5]

dir.create("./ONS_data/pop_size/diagnostics", showWarnings = FALSE)
fwrite(diag[order(junction, -abs(excess_pct))],
       "./ONS_data/pop_size/diagnostics/junction_steps.csv")
p <- ggplot(diag[LAD23CD != "E92000001"], aes(excess_pct)) +
  geom_histogram(binwidth = 0.25) +
  geom_vline(aes(xintercept = excess_pct), diag[LAD23CD == "E92000001"],
             colour = "red") +
  facet_wrap(~ junction, ncol = 1) +
  labs(x = "Step minus mean annual change of previous 3 years (%), ages 20-89",
       y = "Local authorities",
       title = "Population source junctions (red line: England)")
ggsave("./ONS_data/pop_size/diagnostics/junction_steps_plot.png", p,
       width = 16, height = 12, units = "cm", dpi = 150)

print(diag[LAD23CD == "E92000001"])
print(diag[, .(n_flagged = sum(flag), median_excess = median(excess_pct),
               max_abs_excess = max(abs(excess_pct))), keyby = junction])
print(head(diag[order(-abs(excess_pct))], 10))
