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
#' @description
#' A design has a sim_prm list that holds the simulation parameters.
#'
#' @details
#' To be completed...
#'
#' @export
Design <-
  R6::R6Class(classname = "Design",
          public = list(

            #' @field sim_prm The simulation parameters.
            sim_prm = NA,

            #' @description Create a new design object.
            #' @param sim_prm Either a path to a yaml file or a list with appropriate
            #'   format.
            #' @return A new `Design` object.
            #' @examples
            #' design <- Design$new("./validation/design_for_trends_validation.yaml")
            initialize = function(sim_prm) {

              data_type <- typeof(sim_prm)
              if (data_type == "character") {
                sim_prm <- read_yaml(base::normalizePath(sim_prm, mustWork = TRUE))

              } else if (data_type != "list") {

                stop("You can initialise the object only with an R object of
                     type `list` or a path to a YAML configuration file")
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
                  "n_primers"             ,
                  "n_synthpop_aggregation"
                ) %in% names(sim_prm),

                any(sim_prm$clusternumber == 1L, sim_prm$n_cpus == 1L),
                sapply(sim_prm, function(x) if (is.numeric(x)) x >= 1 else TRUE)
              )

              tt <- sim_prm[grepl("_lag$", names(sim_prm))]
              l <- lapply(names(tt), function(x) {
                disease_enum <- (cumsum(unlist(tt) == tt[[x]]) * duplicated(tt))[[x]]
                if (disease_enum == 0L) disease_enum <- 1L
                names(disease_enum) <- NULL
                return(disease_enum)
              })
              names(l) <- paste0(names(tt), "_enum")
              sim_prm <- c(sim_prm, l)

              sim_prm$init_year <- sim_prm$init_year_long - 2000L

              # place holders to be updated from self$update_fromGUI(parameters)

              sim_prm$national_qimd       <- TRUE
              sim_prm$init_year_fromGUI   <- sim_prm$init_year
              sim_prm$sim_horizon_fromGUI <- sim_prm$sim_horizon_max
              sim_prm$locality            <- "England"

              # Create synthpop_dir_ if it doesn't exists
              sim_prm$output_dir <-
                base::normalizePath(sim_prm$output_dir, mustWork = FALSE)
              if (!dir.exists(sim_prm$output_dir)) {
                dir.create(sim_prm$output_dir, recursive = TRUE)
                message(paste0("Directory ", sim_prm$output_dir, " was created"))
              }


              self$sim_prm = sim_prm

              invisible(self)
            },

            #' @description Create a new design object.
            #' @param path Path including filename and extension to save a yaml
            #'   file with the simulation parameters.
            #' @return The `Design` object.
            save_to_disk = function(path) {
              write_yaml(self$sim_prm, base::normalizePath(path, mustWork = FALSE))

              invisible(self)
            },

            #' @description Updates the design object from GUI.
            #' @param GUI_prm A GUI parameter object.
            #' @return The `Design` object.
            update_fromGUI = function(GUI_prm = parameters) {
              self$sim_prm$national_qimd       <- GUI_prm$national_qimd_checkbox
              # T = use national qimd, F = use local qimd
              self$sim_prm$init_year_fromGUI   <- fromGUI_timeframe(GUI_prm)["init year"] - 2000L
              self$sim_prm$sim_horizon_fromGUI <- fromGUI_timeframe(GUI_prm)["horizon"]
              self$sim_prm$locality            <- fromGUI_location(GUI_prm)
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

