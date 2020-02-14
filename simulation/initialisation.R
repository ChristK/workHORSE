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
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
  library(CKutils)
}
if (!require(workHORSEmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  roxygen2::roxygenise("./Rpackage/workHORSE_model_pkg/") # TODO remove before deployment
  remotes::install_local("./Rpackage/workHORSE_model_pkg/", force = TRUE)
  library(workHORSEmisc)
}

output_dir <- function(x = character(0)) {
  normalizePath(paste0(design$output_dir, x),  mustWork = FALSE)
}

synthpop_dir <- "/mnt/storage_slow/synthpop/"

ramdisk_dir <- function(x = character(0)) {
  normalizePath(paste0("/mnt/RAM_disk/", x),  mustWork = FALSE)
}

options(rgl.useNULL = TRUE)  # suppress error by demography in rstudio server
# Then try/install packages...

dependencies(
  c(
    # "gamlss", # only necessary when fitting the models
    # "gamlss.tr", # only necessary when fitting the models
    # "mc2d", # only necessary for generating fixed_mc
    "foreach",
    "doParallel",
    "doRNG",
    # "mc2d", # for rpert()
    "gamlss.dist", # For distr in prevalence.R
    "dqrng",
    "qs",
    "fst",
    "wrswoR",
    "ggplot2",
    "shiny",
    "shinydashboard",
    # "shinydashboardPlus",
    "shinyBS",
    "bsplus", # from https://ijlyttle.github.io/bsplus/index.html
    "shinythemes",
    "shinyjs",
    "shinyWidgets", # http://shinyapps.dreamrs.fr/shinyWidgets/
    "DT",
    "plotly",
    "viridis",
    "colourpicker",
    "dichromat",
    "htmltools",
    "promises", # needed for %...>%
    "future",
    "future.apply",
    "data.table"
  ), TRUE, FALSE, FALSE, FALSE
)

options(future.fork.enable = TRUE) # enable fork in Rstudio TODO remove for production
# plan(multiprocess, workers = design$clusternumber)
plan(list(tweak(multiprocess, workers = 2L), tweak(multiprocess, workers = design$clusternumber/2L))) # for future apply to work within future


# library(promises)
# library(gamlss)
# library(parallel)
options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

# Check that every synthpop file has a metafile and
# that no orphan exists. Then function deletes orphans
# check_synthpop_has_metafile(synthpop_dir)


run_simulation <- function(parameters, iteration_n) {
  # NOTE iteration should be a vector, i.e. 1:20

  # ******************************************************************
  # IMPACT NCD workHORSE
  # ******************************************************************

  # parameters <- qread("./parameters_test.qs") # TODO remove before release
  # iteration_n = 1:20L
  parameters_dt <- workHORSEmisc::fromGUI_to_dt(parameters)

  update_design_fromGUI(design, parameters)


  # This needs to be here to ensure that a synthpop will be generated
  # if necessary even when the user does not press the "confirm area
  # selection" button in GUI
  if (!all(sapply(iteration_n, check_synthpop_exists,
                  locality = parameters$locality_select,
                  n = design$n,
                  synthpop_dir,
                  sim_horizon = design$sim_horizon_max,
                  init_year = design$init_year_long,
                  max_lag = design$maxlag,
                  smoking_relapse_limit = design$smoking_relapse_limit,
                  ageL = design$ageL,
                  ageH = design$ageH,
                  jumpiness = design$jumpiness,
                  simsmok_calibration = design$simsmok_calibration,
                  n_cpu =  design$n_cpus))) {

    X <- iteration_n # Number of iterations
    maxjobs <- design$clusternumber/2 # parallel::detectCores()/2L
    N <- length(X)
    i.list <- splitIndices(N, N/maxjobs)
    result.list <- vector("list", N)

    for (i in seq_along(i.list)) {
      i.vec <- i.list[[i]]
      result.list[i.vec] <- {
        if (Sys.info()[1] == "Windows") {
          cl <- makeCluster(maxjobs) # used for clustering. Windows compatible
          registerDoParallel(cl)
        } else {
          # cl <- makeCluster(maxjobs, type = "FORK")
          registerDoParallel(maxjobs)  # used for forking. Only Linux/OSX compatible
        }

        foreach(
          mc_iter = X[i.vec],
          # .combine = rbind,
          .inorder = TRUE,
          .verbose = TRUE,
          .packages = c(    "gamlss.dist", # For distr in prevalence.R
                            "dqrng",
                            "qs",
                            "fst",
                            "CKutils",
                            "workHORSEmisc",
                            "data.table"),
          .export = NULL, # ls(),
          .noexport = NULL # c("time_mark")
        ) %dorng%
          {
            generate_synthpop(mc = mc_iter,
                              locality = parameters$locality_select,
                              n = design$n,
                              synthpop_dir,
                              sim_horizon = design$sim_horizon_max,
                              init_year = design$init_year_long,
                              max_lag = design$maxlag,
                              smoking_relapse_limit = design$smoking_relapse_limit,
                              ageL = design$ageL,
                              ageH = design$ageH,
                              jumpiness = design$jumpiness,
                              simsmok_calibration = design$simsmok_calibration,
                              include_diseases = TRUE,
                              design = design,
                              n_cpu =  design$n_cpus)

            NULL
          }
        if (exists("cl")) stopCluster(cl)
      }
    }
  }



  if (design$logs) time_mark("start parallelisation")

  # Parallelisation ----
  # To reduce memory requirements I split the pop files in chunks
  chunks_per_file <- 3L
  par_grid <- CJ(iteration_n, chunk = 1:chunks_per_file) # grid to avoid nested parallelisation

  X <- (1:nrow(par_grid)) # Number of iterations
  maxjobs <- design$clusternumber/2L # parallel::detectCores()/2L
  N <- length(X)
  i.list <- splitIndices(N, N/maxjobs)
  result.list <- vector("list", N)

  for (i in seq_along(i.list)) {
    i.vec <- i.list[[i]]
    result.list[i.vec] <- {
      if (Sys.info()[1] == "Windows") {
        cl <- makeCluster(maxjobs) # used for clustering. Windows compatible
        registerDoParallel(cl)
      } else {
        # cl <- makeCluster(maxjobs, type = "FORK")
        registerDoParallel(maxjobs)  # used for forking. Only Linux/OSX compatible
      }

      foreach(
        mc_iter = X[i.vec],
        .inorder = FALSE,
        .verbose = TRUE,
        .packages = c("data.table", "workHORSEmisc", "gamlss.dist", # For distr in prevalence.R
                      "dqrng", "qs", "fst", "wrswoR", "CKutils")
      ) %dopar% {

        if (design$logs) {
          if (!dir.exists(normalizePath(output_dir("logs/")))) dir.create(normalizePath(output_dir("logs/")), FALSE, TRUE)
          sink(file = normalizePath(output_dir(paste0("logs/log", mc_iter, ".txt")
          ), mustWork = FALSE),
          append = TRUE,
          type = "output",
          split = FALSE)
        }
        # mc_iter = 1
        sp_iter <- par_grid$iteration_n[[mc_iter]]
        chunk_iter <- par_grid$chunk[[mc_iter]]

        lags_mc <- get_lag_mc(sp_iter, design) # TODO lag of 10 crashes shift_byID
        max_lag_mc <- max(unlist(lags_mc))
        POP <- get_synthpop(sp_iter, design, max_lag_mc, synthpop_dir, chunk_iter,
                            chunks_per_file)
        fromGUI_lqimd(parameters, POP, design)


        # Run scenarios
        if (design$logs) print("scenario outputs")
        output_chunk <- list()
        output_chunk <- sapply(
          fromGUI_scenario_names(parameters_dt)$true_scenario,
          run_scenario,
          POP,
          parameters_dt,
          lags_mc,
          sp_iter,
          design,
          output_chunk,
          USE.NAMES = FALSE
        )


        rm(POP)
        # invisible(gc(verbose = FALSE, full = TRUE))

        invisible(lapply(output_chunk, setDT)) # does not work if I merge with the one below
        invisible(lapply(output_chunk, function(x) {
          oldnam <- grep("_sc.*$", names(x), value = TRUE)
          newnam <- gsub("_sc.*$", "", oldnam)
          setnames(x, oldnam, newnam)
        }))
        output_chunk <- rbindlist(output_chunk, idcol = "scenario")

        # POP[year == 30, sum(t2dm_prvl>0)]
        # output_chunk[year == 30 & scenario == "sc1", sum(t2dm_prvl>0)]

        setkey(output_chunk, scenario, pid, year) # for identify_longdeads
        output_chunk <- output_chunk[output_chunk$year >= design$init_year_fromGUI &
                                       between(output_chunk$age, design$ageL, design$ageH) &
                                       !identify_longdeads(all_cause_mrtl, pid_mrk),
        ]
        output_chunk[, pid_mrk  := mk_new_simulant_markers(pid)] # Necessary because of pruning above
        output_chunk[, scenario := factor(scenario)]
        generate_health_econ(output_chunk)
        output_chunk[, c("ncc", "pid", "income", "education") := NULL] # Consider deleting income/education also

        # gen incd
        for (nam in grep("_prvl$", names(output_chunk), value = TRUE)) {
          newnam <- gsub("_prvl$", "_incd", nam)
          set(output_chunk, NULL, newnam, fifelse(output_chunk[[nam]] == 1L, 1L, 0L))
        }

        if (design$logs) print("transform prvl/dgn/mrtl to be summed" )
        invisible(output_chunk[, lapply(.SD, clamp_int, 0L, 1L, TRUE),
                               .SDcols = patterns("_prvl$|_dgn$|_mrtl$")])
        if ("lqimd" %in% names(output_chunk)) output_chunk[, ("lqimd") := NULL] else output_chunk[, c("nqimd") := NULL]

        output_chunk[, year := year + 2000L]
        to_agegrp(output_chunk, 20L, 89L, "age", "agegrp", TRUE, 30L)

        # Scale-up to ONS population projections
        for (nam in grep("_prvl$|_dgn$|_incd|_mrtl$|_cost$|^eq5d", names(output_chunk), value = TRUE)) {
          set(output_chunk, NULL, nam, output_chunk[[nam]] * output_chunk$wt)
        }

        # summarise by strata
        output_chunk[, c("age", "pid_mrk") := NULL]
        output_chunk <- output_chunk[, lapply(.SD, sum), keyby = eval(design$strata_for_output)]
        setnames(output_chunk, "wt", "pops")
        output_chunk[, mc := sp_iter]
        tt <- fromGUI_scenario_names(parameters_dt)
        setnames(tt, "true_scenario", "scenario")
        absorb_dt(output_chunk, tt)


        write_fst(output_chunk, ramdisk_dir(paste0("chunk_",sp_iter, "_", chunk_iter, ".fst")), 90)

        rm(output_chunk)
        # invisible(gc(verbose = FALSE, full = TRUE))
        if (design$logs) sink()
        NULL
      }}
    stopImplicitCluster()
    if (exists("cl")) stopCluster(cl)
  }
  if (design$logs) time_mark("End of parallelisation")

  while (sink.number() > 0L) sink()

  filenam <- ramdisk_dir(paste0("chunk_", par_grid$iteration_n, "_", par_grid$chunk, ".fst"))
  filenam[!file.exists(filenam)]
  output <- rbindlist(lapply(filenam, read_fst, as.data.table = TRUE))
  file.remove(filenam)
  strata <- c("mc", design$strata_for_output, "friendly_name")

  output <- output[, lapply(.SD, sum), keyby = eval(strata)]

  # spread overhead policy costs to individuals
  for (sc_nam in levels(output$scenario)) {
    l <- fromGUI_scenario_parms(sc_nam, parameters_dt)
    for (nam in grep("_ovrhd$", names(l), value = TRUE)) {
      newnam <- gsub("^sc_ls_|_cost_ovrhd$", "", nam)
      newnam <- paste0(newnam, "_ovrhd_cost")
      output[, (newnam) := pops * l[[nam]] / sum(pops), by = year]
    }
  }

  # calculate useful cost indices
  output[, policy_cost := invitation_cost + attendees_cost +
           active_days_cost + bmi_cost + alcohol_cost + smoking_cost +
           alcoholreduc_ovrhd_cost + smkcess_ovrhd_cost +
           wghtloss_ovrhd_cost + pa_ovrhd_cost]

  setkey(output, year, friendly_name)

  fwrite_safe(output, output_dir("output.csvy"))

  # calculate net cost/effect
  baseline <- output[scenario == fromGUI_baseline_scenario(parameters_dt)]

  strata <- setdiff(strata, c("scenario", "friendly_name"))


  output[baseline, on = strata, `:=` (
    net_utility = eq5d - i.eq5d,
    net_policy_cost = policy_cost - i.policy_cost,
    net_healthcare_cost = healthcare_cost - i.healthcare_cost,
    net_socialcare_cost = socialcare_cost - i.socialcare_cost,
    net_informal_care_cost = informal_care_cost - i.informal_care_cost,
    net_productivity_cost = productivity_cost - i.productivity_cost, # (higher is better)
    cpp_cvd = i.cvd_incd - cvd_incd,
    cpp_chd = i.chd_incd - chd_incd,
    cpp_stroke = i.stroke_incd - stroke_incd,
    cpp_poststroke_dementia = i.poststroke_dementia_incd - poststroke_dementia_incd,
    cpp_copd = i.copd_incd - copd_incd,
    cpp_lung_ca = i.lung_ca_incd - lung_ca_incd,
    cpp_colon_ca = i.colon_ca_incd - colon_ca_incd,
    cpp_breast_ca = i.breast_ca_incd - breast_ca_incd,
    cpp_htn = i.htn_incd - htn_incd,
    cpp_af = i.af_incd - af_incd,
    cpp_t2dm = i.t2dm_incd - t2dm_incd,

    cypp_cvd = i.cvd_prvl - cvd_prvl,
    cypp_chd = i.chd_prvl - chd_prvl,
    cypp_stroke = i.stroke_prvl - stroke_prvl,
    cypp_poststroke_dementia = i.poststroke_dementia_prvl - poststroke_dementia_prvl,
    cypp_copd = i.copd_prvl - copd_prvl,
    cypp_lung_ca = i.lung_ca_prvl - lung_ca_prvl,
    cypp_colon_ca = i.colon_ca_prvl - colon_ca_prvl,
    cypp_breast_ca = i.breast_ca_prvl - breast_ca_prvl,
    cypp_htn = i.htn_prvl - htn_prvl,
    cypp_af = i.af_prvl - af_prvl,
    cypp_t2dm = i.t2dm_prvl - t2dm_prvl,

    dpp_nonmodelled = i.nonmodelled_mrtl - nonmodelled_mrtl,
    dpp_chd = i.chd_mrtl - chd_mrtl,
    dpp_stroke = i.stroke_mrtl - stroke_mrtl,
    dpp_copd = i.copd_mrtl - copd_mrtl,
    dpp_lung_ca = i.lung_ca_mrtl - lung_ca_mrtl,
    dpp_colon_ca = i.colon_ca_mrtl - colon_ca_mrtl,
    dpp_breast_ca = i.breast_ca_mrtl - breast_ca_mrtl,
    dpp_all_cause = i.all_cause_mrtl - all_cause_mrtl
  )]

  output[, `:=` (
    total_hcp_cost = net_policy_cost + net_healthcare_cost, # healthcare perspective
    total_hscp_cost = net_policy_cost + net_healthcare_cost +  net_socialcare_cost,
    societal_cost = net_policy_cost + net_healthcare_cost +  net_socialcare_cost +
      net_informal_care_cost - net_productivity_cost
  )]

  output[, mc := factor(mc)] # helpful later when I calculate cumsums
  # write_fst(output, "./output/test_results.fst") # TODO remove before release

  return(output)
}
