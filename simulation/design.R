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

# Simulation design
if (!exists("design")) {
  design <- list()
  design$iteration_n    <- 20L # override saved file. If changed in file all outputs will be deleted
  design$clusternumber  <- parallel::detectCores()/2L # Change to your number of CPU cores (explicit parallelisation)
  design$n_cpus         <- 1L  # Change to your number of CPU cores (implicit parallelisation)
  design$logs           <- FALSE
  design$process_output <- TRUE
  design$scenarios      <- ""
}



design$cols_for_output <- c("pid", "year", "pid_mrk", "age", "sex", "qimd",  "wt", # pid NEEDS TO be first
                               "lqimd", "income", "education", "ncc", "ethnicity") # variables to be returned in the output file
design$strata_for_output     <- c("scenario", "year", "agegrp", "sex", "qimd", "ethnicity")
design$exposures             <- c("age", "sex", "qimd", "active_days", "fruit",
                                  "veg", "smok_status", "smok_cig", "ets", "alcohol",
                                  "bmi", "sbp", "tchol")
design$n                     <- 2e5L # Define the sample size
design$init_year             <- 2013L # I subtract 2000 below. DO NOT CHANGE. GUI needs to hook to a secondary parameter
design$init_year_long        <- design$init_year #This remains 2013 throughout
design$sim_horizon_max       <- 2041L - design$init_year_long
design$ageL                  <- 30L  # Define lower age limit to diseases-model simulation (min = 30)
design$ageH                  <- 89L  # Define lower age limit to diseases-model simulation (max = 84)
design$cvd_lag               <- 4L
design$copd_lag              <- 5L
design$cancer_lag            <- 9L # 9 is the max acceptable
design$nonmodelled_lag       <- 5L # 9 is the max acceptable
design$maxlag                <- 10L # 9 is the max acceptable
design$smoking_relapse_limit <- 3L
design$stochastic            <- TRUE
design$kismet                <- TRUE
design$jumpiness             <- 1.0 # increase for more erratic jumps in trajectories
design$export_xps            <- FALSE
design$simsmok_calibration   <- TRUE
design$output_dir            <- "./output/" #tempdir()
design$cancer_cure           <- 10L # need to be >2 & <= 10.
                                    # Otherwise, in dqsample in init_prevalence_fn.r dqsample(2:2, 4, T)
                                    # and the healthcare cost function
design$export_xps            <- FALSE
design$validation            <- FALSE
design$max_prvl_for_outputs  <- 2L # Need to be > 0. If 1 incidence is
# not recorded, only prevalence. If 2, duration 1 is incidence and duration
# 1 + 2 is prevalence

# design$diseases    <- c("chd", "stroke", "lung_ca")
# design$friendly_disease_names <- c("CHD", "Stroke", "Lung cancer")
# design$fatal_diseases <- c("nonmodelled", "chd", "stroke", "lung_ca")
# design$friendly_fatal_diseases_names <- c("Other COD", "CHD", "Stroke", "Lung cancer")
# design$scenarios <- gsub(".R", "", sort(list.files("./Scenarios", pattern = "^sc.*\\.R$")))


# Calculate disease_enum logic for lagged exposures
# This ensures independency of age with same median for different diseases
tt <- design[grepl("_lag$", names(design))]
l <- lapply(names(tt), function(x) {
  disease_enum <- (cumsum(unlist(tt) == tt[[x]]) * duplicated(tt))[[x]]
  if (disease_enum == 0L) disease_enum <- 1L
  names(disease_enum) <- NULL
  return(disease_enum)
})
names(l) <- paste0(names(tt), "_enum")
design <- c(design, l)

design$init_year <- design$init_year - 2000L

rm(tt, l)

update_design_fromGUI <- function(design, x = parameters) {
  design$national_qimd       <<- x$national_qimd_checkbox # T = use national qimd, F = use local qimd
  design$init_year_fromGUI   <<- fromGUI_timeframe(x)["init year"] - 2000L
  design$sim_horizon_fromGUI <<- fromGUI_timeframe(x)["horizon"]
  design$locality            <<- fromGUI_location(x)
}
