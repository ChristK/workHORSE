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

# Downloads the ONS population sources used by transform_pops.R and writes
# compact extracts (year, age, sex, [LAD23CD, LAD23NM], pops) that are
# committed. Raw downloads go to ./ONS_data/pop_size/raw/ (gitignored).
# Run from the project root. See ./ONS_data/pop_size/README.md for details.

library(data.table)
library(readxl)
library(fst)

raw_dir <- "./ONS_data/pop_size/raw"
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

src <- list(
  # 2022-based SNPP, migration category variant (ONS's headline projection for
  # this release, in place of a principal projection), April 2023 LA boundaries
  snpp = c(
    url = "https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/populationandmigration/populationprojections/datasets/localauthoritiesinenglandz1/2022basedmigrationcategoryvarianton2023localauthoritygeographies/2022snpppopulationsyoamigcat23.zip",
    file = "2022snpppopulationsyoamigcat23.zip"
  ),
  # Mid-year estimates, mid-2011 to mid-2025 detailed time series, April 2023
  # LA boundaries
  mye = c(
    url = "https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/estimatesofthepopulationforenglandandwales/mid2011tomid2025detailedtimeseries/myebtablesenglandwales20112025.xlsx",
    file = "myebtablesenglandwales20112025.xlsx"
  ),
  # 2024-based national population projections, England (zip of all variants)
  npp = c(
    url = "https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/populationandmigration/populationprojections/datasets/z3zippedpopulationprojectionsdatafilesengland/2024based/en1.zip",
    file = "npp2024_en1.zip"
  )
)

for (s in src) {
  f <- file.path(raw_dir, s[["file"]])
  if (!file.exists(f)) {
    download.file(s[["url"]], f, mode = "wb", method = "libcurl",
                  headers = c("User-Agent" = "Mozilla/5.0 (workHORSE data update)"))
  }
}

# SNPP 2022-based, migration category variant, LAD23 ----
unzip(file.path(raw_dir, src$snpp[["file"]]), exdir = file.path(raw_dir, "snpp"))
snpp <- rbindlist(lapply(c("males", "females"), function(x)
  fread(file.path(raw_dir, "snpp", paste0("2022 SNPP Population ", x, ".csv")),
        header = TRUE)))
stopifnot(snpp[, all(COMPONENT == "Population")])
snpp <- snpp[AGE_GROUP != "All ages"]
snpp[AGE_GROUP == "90 and over", AGE_GROUP := "90"]
snpp <- melt(snpp, c("AREA_CODE", "AREA_NAME", "COMPONENT", "SEX", "AGE_GROUP"),
             variable.name = "year", value.name = "pops",
             variable.factor = FALSE)
snpp <- snpp[, .(year    = as.integer(year),
                 age     = as.integer(AGE_GROUP),
                 sex     = fifelse(SEX == "males", "men", "women"),
                 LAD23CD = AREA_CODE,
                 LAD23NM = AREA_NAME,
                 pops)]
stopifnot(!anyNA(snpp), snpp[, uniqueN(LAD23CD)] == 296L,
          snpp[, setequal(unique(year), 2022:2047)],
          snpp[, setequal(unique(age), 0:90)])
setkey(snpp, year, age, sex, LAD23CD)
write_fst(snpp, "./ONS_data/pop_size/ons_snpp2022_migcat_lad23.fst", 100)

# MYE mid-2011 to mid-2025, LAD23 ----
mye <- as.data.table(read_excel(file.path(raw_dir, src$mye[["file"]]),
                                sheet = "MYEB1", skip = 1))
mye <- mye[country == "E"] # drop Wales
mye <- melt(mye, c("ladcode23", "laname23", "country", "sex", "age"),
            measure.vars = patterns("^population_"),
            variable.name = "year", value.name = "pops",
            variable.factor = FALSE)
mye <- mye[, .(year    = as.integer(sub("^population_", "", year)),
               age     = as.integer(age), # 90 is 90+
               sex     = fifelse(sex == "m", "men", "women"),
               LAD23CD = ladcode23,
               LAD23NM = laname23,
               pops)]
stopifnot(!anyNA(mye), mye[, uniqueN(LAD23CD)] == 296L,
          mye[, setequal(unique(year), 2011:2025)],
          mye[, setequal(unique(age), 0:90)])
setkey(mye, year, age, sex, LAD23CD)
write_fst(mye, "./ONS_data/pop_size/ons_mye_2011_2025_lad23.fst", 100)

# NPP 2024-based principal projection, England ----
unzip(file.path(raw_dir, src$npp[["file"]]),
      files = "en_ppp_machine_readable.xlsx", exdir = file.path(raw_dir, "npp"))
npp <- as.data.table(read_excel(
  file.path(raw_dir, "npp", "en_ppp_machine_readable.xlsx"),
  sheet = "Population"))
npp <- melt(npp, c("Sex", "Age"), variable.name = "year", value.name = "pops",
            variable.factor = FALSE)
# Ages are 0-104, "105 - 109" and "110 and over". Collapse 90+ into 90
npp[, age := pmin(as.integer(sub("^([0-9]+).*$", "\\1", Age)), 90L)]
npp <- npp[, .(pops = sum(pops)),
           keyby = .(year = as.integer(year), age,
                     sex = fifelse(Sex == "Males", "men", "women"))]
stopifnot(!anyNA(npp), npp[, setequal(unique(age), 0:90)],
          npp[, min(year) == 2024L && max(year) >= 2047L])
write_fst(npp, "./ONS_data/pop_size/ons_npp2024_ppp_england.fst", 100)
