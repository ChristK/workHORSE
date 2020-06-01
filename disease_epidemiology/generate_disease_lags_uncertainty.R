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

cat("Initialising workHORSE model...\n\n")
if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
  library(CKutils)
}
dependencies(c("data.table", "fst"))

set.seed(2506120L)

mc_max  <- 1e3
mc_max2 <- 10L # To allow 2 or more diseases have the same lag in years but be independent

dt <- CJ(mc = 1:mc_max, disease_enum = 1:mc_max2)
for (i in 3:9) {
  set(dt, NULL, paste0("lag_", i), 2L + rbinom(nrow(dt), 8L, (i - 2) / 8L))
}
write_fst(dt, "./disease_epidemiology/disease_lags_l.fst", 100L)

# create a table with row numbers for each mc
dt[, rn := .I]
tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
write_fst(tt, "./disease_epidemiology/disease_lags_indx.fst", 100L)

# tt <- 2L + rbinom(1e5, 8L, (9L - 2L) / 8L)
# mean(tt)
# hist(tt)
