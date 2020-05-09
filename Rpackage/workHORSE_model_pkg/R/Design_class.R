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


#' R6 Class representing a simulation design
#'
#' A design has a sim_prm list that holds the simulation parameters.

Design <-
  R6::R6Class(classname = "Design",
          public = list(

            #' @field sim_prm The simulation parameters.
            sim_prm = NA,

            #' @description
            #' Create a new design object.
            #' @param l Either a path to a yaml file or a list with appropriate format.
            #' @return A new `Design` object.
            initialize = function(l) {

              data_type <- typeof(l)
              if (data_type == "character") {
                l <- read_yaml(base::normalizePath("test", mustWork = TRUE))

              } else if (data_type != "list") {

                stop("You can initialise the object only with an R object of type `list` or a path to a YAML configuration file")
              }

              # Validation
              stopifnot(
                c(
                  "iteration_n"           ,
                  "clusternumber"         ,
                  "n_cpus"                ,
                  "logs"                  ,
                  "process_output"        ,
                  "scenarios"             ,
                  "cols_for_output"       ,
                  "strata_for_output"     ,
                  "exposures"             ,
                  "n"                     ,
                  "init_year"             ,
                  "init_year_long"        ,
                  "sim_horizon_max"       ,
                  "ageL"                  ,
                  "ageH"                  ,
                  "cvd_lag"               ,
                  "copd_lag"              ,
                  "cancer_lag"            ,
                  "nonmodelled_lag"       ,
                  "maxlag"                ,
                  "smoking_relapse_limit" ,
                  "stochastic"            ,
                  "kismet"                ,
                  "jumpiness"             ,
                  "export_xps"            ,
                  "simsmok_calibration"   ,
                  "output_dir"            ,
                  "cancer_cure"           ,
                  "validation"            ,
                  "max_prvl_for_outputs"  ,
                  "cvd_lag_enum"          ,
                  "copd_lag_enum"         ,
                  "cancer_lag_enum"       ,
                  "nonmodelled_lag_enum"  ,
                  "locality"              ,
                  "init_year_fromGUI"     ,
                  "sim_horizon_fromGUI"
                ) %in% names(l)
              )

              self$sim_prm = l

              invisible(self)
            },

            #' @description
            #' Create a new design object.
            #' @param path Path including filename and extension to save a yaml file with the simulation parameters.
            #' @return The `Design` object.
            save_to_disk = function(path) {
              write_yaml(self$sim_prm, base::normalizePath(path, mustWork = FALSE))

              invisible(self)
            },

            #' @description
            #' Print the simulation parameters.
            #' @return The `Design` object.
            print = function() {
              print(self$sim_prm)

              invisible(self)
            }
          ))

# design <- Design$new("./validation/design_for_validation.yaml") # or Design$new(design)
