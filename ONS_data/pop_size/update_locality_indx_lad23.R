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

# Adds April 2023 local authority codes and names (LAD23CD, LAD23NM) to
# ./synthpop/lsoa_to_locality_indx.fst and writes the LAD17 -> LAD23 lookup.
# The script that originally generated the index is not part of this repo, so
# the index is updated in place; all existing columns are kept.
# Every LAD17 lies wholly within one LAD23 (all reorganisations since 2017
# merged whole districts), so LSOA11 -> LAD17 -> LAD23 is exact.
# Run from the project root after download_ons_pop.R.

library(data.table)
library(fst)

# LAD23 codes that replaced groups of LAD17 districts (2019, 2020, 2021 and
# 2023 reorganisations). All other LAD17 codes are unchanged in LAD23.
merged <- rbindlist(list(
  # Bournemouth, Christchurch and Poole (2019)
  data.table(LAD23CD = "E06000058",
             LAD17CD = c("E06000028", "E06000029", "E07000048")),
  # Dorset (2019)
  data.table(LAD23CD = "E06000059",
             LAD17CD = c("E07000049", "E07000050", "E07000051", "E07000052",
                         "E07000053")),
  # Buckinghamshire (2020)
  data.table(LAD23CD = "E06000060",
             LAD17CD = c("E07000004", "E07000005", "E07000006", "E07000007")),
  # North Northamptonshire (2021)
  data.table(LAD23CD = "E06000061",
             LAD17CD = c("E07000150", "E07000152", "E07000153", "E07000156")),
  # West Northamptonshire (2021)
  data.table(LAD23CD = "E06000062",
             LAD17CD = c("E07000151", "E07000154", "E07000155")),
  # Cumberland (2023)
  data.table(LAD23CD = "E06000063",
             LAD17CD = c("E07000026", "E07000028", "E07000029")),
  # Westmorland and Furness (2023)
  data.table(LAD23CD = "E06000064",
             LAD17CD = c("E07000027", "E07000030", "E07000031")),
  # North Yorkshire (2023)
  data.table(LAD23CD = "E06000065",
             LAD17CD = sprintf("E07%06d", 163:169)),
  # Somerset (2023; Taunton Deane and West Somerset merged in 2019 first)
  data.table(LAD23CD = "E06000066",
             LAD17CD = sprintf("E07%06d", 187:191)),
  # East Suffolk (2019)
  data.table(LAD23CD = "E07000244", LAD17CD = c("E07000205", "E07000206")),
  # West Suffolk (2019)
  data.table(LAD23CD = "E07000245", LAD17CD = c("E07000201", "E07000204"))
))
stopifnot(nrow(merged) == 41L, !anyDuplicated(merged$LAD17CD))

indx <- read_fst("./synthpop/lsoa_to_locality_indx.fst", as.data.table = TRUE)

lookup <- unique(indx[, .(LAD17CD = as.character(LAD17CD),
                          LAD17NM = as.character(LAD17NM),
                          RGN11NM = as.character(RGN11NM))])
lookup[, LAD23CD := LAD17CD]
lookup[merged, on = "LAD17CD", LAD23CD := i.LAD23CD]

# Official LAD23 names (e.g. Shepway is now Folkestone and Hythe)
lad23_nm <- unique(read_fst("./ONS_data/pop_size/ons_mye_2011_2025_lad23.fst",
                            columns = c("LAD23CD", "LAD23NM"),
                            as.data.table = TRUE))
snpp_cd <- unique(read_fst("./ONS_data/pop_size/ons_snpp2022_migcat_lad23.fst",
                           columns = "LAD23CD")$LAD23CD)
lookup[lad23_nm, on = "LAD23CD", LAD23NM := i.LAD23NM]

stopifnot(
  nrow(lookup) == 326L,
  !anyNA(lookup$LAD23NM),                         # every LAD17 has a LAD23
  lookup[, uniqueN(LAD23CD)] == 296L,
  setequal(lookup$LAD23CD, lad23_nm$LAD23CD),     # same 296 LAs as the MYE
  setequal(lookup$LAD23CD, snpp_cd),              # and as the SNPP
  lookup[, uniqueN(RGN11NM), by = LAD23CD][, all(V1 == 1L)], # one region each
  !any(lookup$LAD23NM %in% lookup$RGN11NM),       # areas are matched by name
  lookup[, uniqueN(LAD23CD), by = LAD23NM][, all(V1 == 1L)]
)

indx[lookup, on = "LAD17CD", `:=`(LAD23CD = i.LAD23CD, LAD23NM = i.LAD23NM)]
indx[, `:=`(LAD23CD = factor(LAD23CD), LAD23NM = factor(LAD23NM))]
setcolorder(indx, c("LAD23CD", "LAD23NM"), after = "LAD17NM")
stopifnot(!anyNA(indx$LAD23CD), nrow(indx) == 32844L)
setkey(indx, LSOA11CD)

write_fst(indx, "./synthpop/lsoa_to_locality_indx.fst", 100)
fwrite(lookup[order(LAD23CD, LAD17CD),
              .(LAD17CD, LAD17NM, LAD23CD, LAD23NM)],
       "./ONS_data/pop_size/lad17_to_lad23_lookup.csv")
