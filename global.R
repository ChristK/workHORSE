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

# file.remove(list.files("./output/", full.names = TRUE, recursive = TRUE))

cat("Initialising workHORSE model...\n\n")
if (!require(CKutils)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
  library(CKutils)
}
if (!require(workHORSEmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  roxygen2::roxygenise("./Rpackage/workHORSE_model_pkg/", clean = TRUE)
  # TODO remove before deployment
  remotes::install_local("./Rpackage/workHORSE_model_pkg/", force = TRUE)
  library(workHORSEmisc)
}

options(rgl.useNULL = TRUE)  # suppress error by demography in rstudio server
dependencies(yaml::read_yaml("./dependencies.yaml"))

design <- Design$new("./simulation/sim_design.yaml")

options(future.fork.enable = TRUE) # TODO remove for production
registerDoFuture()



options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

strata_for_gui <- c("mc", "friendly_name", design$sim_prm$strata_for_output)

def_col <- viridis(16, option = "D")
def_col_small <- def_col[c(1, 9, 5, 3, 7, 2, 8, 4, 6)]

def_sym <- c("circle-dot", "square", "diamond", "cross", "triangle-up",
  "pentagon", " star", "hexagon-open-dot", "triangle-down")

# Produce scenario tabs 2-9 using scenario 1 as a template if not exist ------
for (i in 2:9) {
  if (!file.exists(paste0("ui/scenario", i, "_tab.R")) ||
      !identical(file.size("ui/scenario1_tab.R"), file.size(paste0("ui/scenario", i, "_tab.R")))) {
    tt <- readLines(file.path("ui", "scenario1_tab.R"))
    tt <- gsub("Scenario 1",
               paste0("Scenario ", i),
               gsub(
                 "input.level == 1",
                 paste0("input.level == ", i),
                 gsub("_sc1",  paste0("_sc", i),
                      gsub("def_col_small\\[\\[1\\]\\]",
                           paste0("def_col_small\\[\\[", i, "\\]\\]"),
                      gsub("def_sym\\[\\[1\\]\\]",
                           paste0("def_sym\\[\\[", i, "\\]\\]"), tt)
                      ))))
    writeLines(tt, paste0("ui/scenario", i, "_tab.R"))
  }
}



# load localities structure
localities_indx <- read_fst("./synthpop/lsoa_to_locality_indx.fst", as.data.table = TRUE)
localitities_list <- list(
  "Country" = list("England"),
  "Region" = sort(as.character((localities_indx[, unique(RGN11NM)]))),
  "Local Authority" = sort(as.character((localities_indx[, unique(LAD17NM)])))
)
rm(localities_indx)



# sum(gtools::rdirichlet(1, rep(1, 180))) # to simulate a 180 length vector that sums to 1
# extract_uptake_table(setDT(readRDS("./output/input_parameters.rds")), 2)[]
# To write when running from shiny server or docker see https://groups.google.com/forum/#!topic/shiny-discuss/srWETT6uL-I


# Table related ----
# default global search value
if (!exists("default_search_tbl_summary")) default_search_tbl_summary <- ""

# default column search values
if (!exists("default_search_columns_tbl_summary")) default_search_columns_tbl_summary <- NULL

# dockerfile_object <- containerit::dockerfile()
# dockerfile_object
# containerit::write(dockerfile_object, file = "./dockerfile")
