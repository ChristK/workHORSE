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

#' @export
fromGUI_to_dt <- function(parameters) {
  parameters_dt <- lapply(parameters, unclass)
  parameters_dt <- data.table(
    "input_names" = names(parameters_dt),
    "value" = parameters_dt,
    key = "input_names"
  )
  parameters_dt[, c("input_names", "scenario") := tstrsplit(input_names, "_sc", fixed = TRUE, keep = 1:2)]
  parameters_dt[!is.na(scenario), scenario := paste0("sc", scenario)]
  parameters_dt[is.na(scenario), scenario := "global"]
  setcolorder(parameters_dt)
  parameters_dt <-
    parameters_dt[lapply(lapply(value, `[`), length) > 0, ]
  setkey(parameters_dt, input_names)
  tt <- parameters_dt["scenarios_number_slider", value[[1]]]
  parameters_dt <-
    parameters_dt[scenario %in% c("global", paste0("sc", seq_len(tt))), ]
  parameters_dt <-
    parameters_dt[!input_names %like% "^shinyjs-delay-|^btn-"]
  # write_fst(parameters_dt, "./parameters_dt.fst")

  # parameters_dt[, transpose(lapply(lapply(value, `[`), length))][, table(V1)]
  parameters_dt[, value1 := sapply(value, `[`, 1)]
  parameters_dt[, value2 := sapply(value, `[`, 2)]
  fromGUI_scenario_names(parameters_dt)
  return(parameters_dt)
}
# parameters_dt <- fromGUI_to_dt(parameters)

#' @export
fromGUI_location <- function(parameters) {
  par <- parameters$locality_select
  # if (!exists("localities_indx")) localities_indx <- read_fst("./synthpop/lsoa_to_locality_indx.fst", as.data.table = TRUE)
  # localities_indx <- unique(localities_indx[, .(LAD17CD, LAD17NM, SHA11CD, SHA11NM)])
  # if (length(par) == 1L && par %in% localities_indx$LAD17NM) {
  #   return(localities_indx[LAD17NM == par, as.character(LAD17CD)])
  # } else if (length(par) == 1L && par == "England") {
  #   return("national")
  # }
}
# fromGUI_location(parameters)

#' @export
fromGUI_timeframe <- function(parameters) {
  out <- c(parameters$simulation_period_slider[1],
           parameters$simulation_period_slider[2] - parameters$simulation_period_slider[1])
  names(out) <- c("init year", "horizon")
  out
}
# fromGUI_timeframe(parameters)["init year"]
# fromGUI_timeframe(parameters)["horizon"]

#' @export
fromGUI_lqimd <- function(parameters, POP, design) {
  if (!parameters$national_qimd_checkbox) {
    setnames(POP, c("qimd", "lqimd"), c("nqimd", "qimd"))
    design$cols_for_output <<- c(setdiff(design$cols_for_output, "lqimd"), "nqimd")
  }
}

#' @export
fromGUI_scenario_count <- function(parameters_dt) {
  if (!is.data.table(parameters_dt)) parameters_dt <- fromGUI_to_dt(parameters_dt)
  # TODO make it work for more than one parallel and one serial ensemble scenarios
  ns <-
    parameters_dt[input_names %flike% "scenarios_number_slider",
                  as.integer(value1)]
  # get parallel ensemble scenarios
  nsp <-
    parameters_dt[input_names %flike% "parallel_ensemble_checkbox",
                  sum(as.logical(value1))]
  nss <-
    parameters_dt[input_names %flike% "serial_ensemble_checkbox",
                  sum(as.logical(value1))]
  return(ns - nsp - nss + 2L)
}
# fromGUI_scenario_count(parameters_dt)
# system.time({fromGUI_scenario_count(parameters_dt)})
# system.time({fromGUI_scenario_count(parameters)})

#' @export
fromGUI_scenario_names <- function(parameters_dt) {
  if (!"friendly_name" %in% names(parameters_dt)) {
    tt = parameters_dt[input_names %flike% "friendly", .(scenario, friendly_name = unlist(value))]
    absorb_dt(parameters_dt, tt)
  }
  # taking into account ensembles
  if (!"true_scenario" %in% names(parameters_dt)) {
    parameters_dt[, true_scenario := scenario]
    parameters_dt[scenario == "global", true_scenario := NA]
    # handle ensemble scenarios
    if (parameters_dt[input_names %flike% "parallel_ensemble_checkbox" &
                      as.logical(value1), .N > 0L]) {
      tt <- parameters_dt[input_names %flike% "parallel_ensemble_checkbox" &
                            as.logical(value1),
                          .(scenario, true_scenario = paste(scenario, collapse = ""),
                            friendly_name = first(sort(unique(friendly_name))))]
      parameters_dt[tt, on = "scenario", `:=` (true_scenario = i.true_scenario,
                                               friendly_name = i.friendly_name)]
    }
    if (parameters_dt[input_names %flike% "serial_ensemble_checkbox" &
                      as.logical(value1), .N > 0L]) {
      tt <- parameters_dt[input_names %flike% "serial_ensemble_checkbox" &
                            as.logical(value1),
                          .(scenario, true_scenario = paste(scenario, collapse = ""),
                            friendly_name = first(sort(unique(friendly_name))))]
      parameters_dt[tt, on = "scenario", `:=` (true_scenario = i.true_scenario,
                                               friendly_name = i.friendly_name)]
    }
  }
  setorder(unique(na.omit(parameters_dt[, .(true_scenario, friendly_name)])))[]
}
# fromGUI_scenario_names(parameters_dt)$true_scenario
# parameters_dt[, c("true_scenario", "friendly_name") := NULL]

#' @export
fromGUI_baseline_scenario <- function(parameters_dt) {
  parameters_dt[input_names == "baseline" & as.logical(value1), unique(true_scenario)]
}
# fromGUI_baseline_scenario(parameters_dt)

#' @export
fromGUI_uptake_table <- function(scenario_nam, parameters_dt) {
  # recreate age groups
  if (parameters_dt[input_names == "age_eligibility_slider" & scenario == scenario_nam, value[[1]][[2]]] - parameters_dt[input_names == "age_eligibility_slider" & scenario == scenario_nam, value[[1]][[1]]] < 11) {
    agegrp_nam <-
      paste0(parameters_dt[input_names == "age_eligibility_slider" & scenario == scenario_nam, value[[1]][[1]]], "-", parameters_dt[input_names == "age_eligibility_slider" & scenario == scenario_nam, value[[1]][[2]]])
  } else {
    agegrp_nam <-
      agegrp_name(as.integer(parameters_dt[input_names == "age_eligibility_slider" & scenario == scenario_nam, value[[1]][[1]]] /
                               10) * 10L,
                  as.integer(parameters_dt[input_names == "age_eligibility_slider" & scenario == scenario_nam, value[[1]][[2]]] /
                               10) * 10L,
                  10L)
  }

  tt <- CJ(qimd = factor(c("1 most deprived", "2", "3", "4", "5 least deprived")),
           sex = factor(c("men", "women")),
           agegrp10 = factor(agegrp_nam))

  # out <- parameters_dt[input_names %like% "uptake_qrisk_", .(input_names, "val" = unlist(get(sc_nam)))]
  # if (!"val" %in% names(out)) out[, val := NA_real_]
  out <- parameters_dt[input_names %like% "uptake_qrisk_" & scenario == scenario_nam, .(input_names, "val" = unlist(value))]
  if (nrow(out) == 0L) stop("Parameters don't exist.")

  out[, c("qrisk_cat", "ord") := tstrsplit(input_names, "_", fixed = TRUE, keep = 3:4)]
  out <- setcolorder(dcast(out, as.integer(ord)~qrisk_cat, value.var = "val"), c("ord", "low", "mid", "high"))
  if (is.unsorted(out$ord)) setkey(out, ord)
  out[, ord := NULL]
  out <- cbind(tt, out)
  out <- melt(out, 1:3, variable.name = "Qrisk2_cat",  value.nam = "val")
  return(out)
}

#' @export
fromGUI_uptake_table_agegrps <-
  function(scenario_nam, parameters_dt) {
    if (parameters_dt[input_names == "age_eligibility_slider" &
                      scenario == scenario_nam, value[[1]][[2]]] -
        parameters_dt[input_names == "age_eligibility_slider" &
                      scenario == scenario_nam, value[[1]][[1]]] < 11) {
      agegrp_nam <-
        paste0(parameters_dt[input_names == "age_eligibility_slider" &
                               scenario == scenario_nam, value[[1]][[1]]], "-",
               parameters_dt[input_names == "age_eligibility_slider" &
                               scenario == scenario_nam, value[[1]][[2]]])
    } else {
      agegrp_nam <-
        agegrp_name(
          floor(parameters_dt[input_names == "age_eligibility_slider" &
                                scenario == scenario_nam, value[[1]][[1]]] /
                  10) * 10L,
          floor(parameters_dt[input_names == "age_eligibility_slider" &
                                scenario == scenario_nam, value[[1]][[2]]] /
                  10) * 10L,
          10L,
          match_input = TRUE,
          match_input_max_age = parameters_dt[input_names == "age_eligibility_slider" &
                                                scenario == scenario_nam, value[[1]][[2]]]
        )
    }
    out <-
      data.table(age = parameters_dt[input_names == "age_eligibility_slider" &
                                       scenario ==
                                       scenario_nam,
                                     value[[1]][[1]]]:parameters_dt[input_names ==
                                                                      "age_eligibility_slider" &
                                                                      scenario == scenario_nam, value[[1]][[2]]],
                 agegrp10 = factor(agegrp_nam))
    out
  }
# fromGUI_uptake_table_agegrps(scenario_nam, parameters_dt)

#' @export
fromGUI_statins_table <- function(scenario_nam, parameters_dt) {
  out <- parameters_dt[input_names %like% "statin_px_qrisk" & scenario == scenario_nam, .(input_names, "val" = as.numeric(unlist(value)))]
  if (nrow(out) == 0L) stop("Parameters don't exist.")
  out[, c("Qrisk2_cat", "qimd") := tstrsplit(input_names, "_", fixed = TRUE, keep = 4:5)]
  replace_from_table(out, "qimd", as.character(1:5), c("1 most deprived", "2", "3", "4", "5 least deprived"))
  out[, `:=`(qimd = factor(qimd),
              Qrisk2_cat = factor(Qrisk2_cat, c("low","mid","high")),
              input_names = NULL)]
  return(out)
}

#' @export
fromGUI_antihtn_table <- function(scenario_nam, parameters_dt) {
  out <- parameters_dt[input_names %like% "antihtn_px_qrisk" & scenario == scenario_nam, .(input_names, "val" = as.numeric(unlist(value)))]
  if (nrow(out) == 0L) stop("Parameters don't exist.")
  out[, c("Qrisk2_cat", "qimd") := tstrsplit(input_names, "_", fixed = TRUE, keep = 4:5)]
  replace_from_table(out, "qimd", as.character(1:5), c("1 most deprived", "2", "3", "4", "5 least deprived"))
  out[, `:=` (qimd = factor(qimd),
              Qrisk2_cat = factor(Qrisk2_cat, c("low","mid","high")),
              input_names = NULL)]
  return(out)
}


#' @export
fromGUI_scenario_parms <- function(scenario_nam, parameters_dt) {
  l <- list()
  l$sc_is_ensemble <- parameters_dt[true_scenario == scenario_nam &
                                      input_names %like% "^parallel_ensemble_checkbox|^serial_ensemble_checkbox",
                                    any(as.logical(value1))]
  if (l$sc_is_ensemble) {
    # TODO logic for ensembles
    l$sc_is_baseline <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names == "baseline", any(as.logical(value1))]

  } else {

    l$sc_is_health_econ_perspective <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names == "health_econ_perspective_checkbox", unlist(value)]
    l$sc_is_baseline <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names == "baseline", as.logical(value1)]
    l$sc_init_year <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "ininit_year_slider", as.integer(value1)]
    l$sc_last_year <-
      parameters_dt[input_names %flike% "simulation_period_slider", as.integer(value2)]
    l$sc_eligib_age <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "age_eligibility_slider", unlist(value)]
    l$sc_eligib_freq <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "frequency_eligibility_slider", as.integer(value1)]
    l$sc_eligib_htn <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "invite_known_hypertensives_checkbox", as.logical(value1)]
    l$sc_eligib_diab <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "invite_known_diabetics_checkbox", as.logical(value1)]
    l$sc_eligib_noone <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "cancel_program_checkbox", as.logical(value1)]
    l$sc_invit_detailed <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "coverage_detailed_checkbox", as.logical(value1)]
    if (l$sc_invit_detailed) { # if by qimd
      l$sc_invit_qimd1 <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_qimd1_slider", as.numeric(value1)/100]
      l$sc_invit_qimd2 <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_qimd2_slider", as.numeric(value1)/100]
      l$sc_invit_qimd3 <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_qimd3_slider", as.numeric(value1)/100]
      l$sc_invit_qimd4 <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_qimd4_slider", as.numeric(value1)/100]
      l$sc_invit_qimd5 <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_qimd5_slider", as.numeric(value1)/100]

      l$sc_invit_qimd1_cost <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_cost_qimd1", as.numeric(value1)]
      l$sc_invit_qimd2_cost <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_cost_qimd2", as.numeric(value1)]
      l$sc_invit_qimd3_cost <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_cost_qimd3", as.numeric(value1)]
      l$sc_invit_qimd4_cost <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_cost_qimd4", as.numeric(value1)]
      l$sc_invit_qimd5_cost <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_cost_qimd5", as.numeric(value1)]
    } else {
      l$sc_invit_qimdall <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_qimd0_slider", as.numeric(value1)/100]
      l$sc_invit_qimdall_cost <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "coverage_cost_qimd0", as.numeric(value1)]
    }

    l$sc_qrisk_ignore_tchol <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "ignore_cholesterol_checkbox", as.logical(value1)]
    l$sc_qrisk_ignore_sbp <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "ignore_sbp_checkbox", as.logical(value1)]
    l$sc_qrisk_ignore_bmi <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "ignore_bmi_checkbox", as.logical(value1)]
    l$sc_uptake_detailed <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "uptake_detailed_checkbox", as.logical(value1)]
    if (l$sc_uptake_detailed) { # if detailed uptake
      l$sc_uptake <- fromGUI_uptake_table(scenario_nam = scenario_nam, parameters_dt = parameters_dt)
      setnafill(l$sc_uptake, "c", 0, cols = "val")
      l$sc_uptake[, val := val/sum(val)]
      setnames(l$sc_uptake, "val", "uptake_wt")

      l$sc_uptake_structural0s <-
        parameters_dt[true_scenario == scenario_nam &
                        input_names %flike% "uptake_structural_zeroes", as.logical(value1)]
    }
    l$sc_uptake_all <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "uptake_slider", as.numeric(value1)/100]
    l$sc_uptake_all_cost <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "uptake_cost", as.numeric(value1)]
    l$sc_px_statins <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "statin_px_slider", as.numeric(value1)/100]
    l$sc_px_antihtn <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "antihtn_px_slider", as.numeric(value1)/100]
    l$sc_px_detailed <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "px_detailed_checkbox", as.logical(value1)]
    if (l$sc_px_detailed) {
      l$sc_px_statins_wt <- fromGUI_statins_table(scenario_nam = scenario_nam,
                                                  parameters_dt = parameters_dt)
      setnafill(l$sc_px_statins_wt, "c", 0L, cols = "val")
      l$sc_px_statins_wt[val == 0, val := 1e-6]
      l$sc_px_statins_wt[, val := val/sum(val)] # no reaso to * l$sc_px_statins
      setnames(l$sc_px_statins_wt, "val", "px_statins_wt")

      l$sc_px_antihtn_wt <- fromGUI_antihtn_table(scenario_nam = scenario_nam,
                                                  parameters_dt = parameters_dt)
      setnafill(l$sc_px_antihtn_wt, "c", 0, cols = "val")
      l$sc_px_antihtn_wt[val == 0, val := 1e-6]

      l$sc_px_antihtn_wt[, val := val/sum(val)] # no reason to l$sc_px_antihtn
      setnames(l$sc_px_antihtn_wt, "val", "px_antihtn_wt")
    }

    l$sc_ls_smkcess <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "smkcess_slider", as.numeric(value1)/100]
    l$sc_ls_smkcess_cost_ind <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %like% "smkcess_cost$", as.numeric(value1)]
    l$sc_ls_smkcess_cost_ovrhd <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "smkcess_cost_ovrhd", as.numeric(value1)]
    l$sc_ls_wghtpct <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "wghtpct_slider", as.numeric(value1)/100]
    l$sc_ls_wghtreduc <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "wghtreduc_slider", as.numeric(value1)/100]
    l$sc_ls_wghtloss_cost_ind <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %like% "wghtloss_cost$", as.numeric(value1)]
    l$sc_ls_wghtloss_cost_ovrhd <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "wghtloss_cost_ovrhd", as.numeric(value1)]
    l$sc_ls_papct <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "papct_slider", as.numeric(value1)/100]
    l$sc_ls_papincr <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "papincr_slider", as.integer(value1)]
    l$sc_ls_pa_cost_ind <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %like% "pa_cost$", as.numeric(value1)]
    l$sc_ls_pa_cost_ovrhd <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "pa_cost_ovrhd", as.numeric(value1)]
    l$sc_ls_alcoholpct <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "alcoholpct_slider", as.numeric(value1)/100]
    l$sc_ls_alcoholreduc <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "alcoholreduc_slider", as.integer(value1)/100]
    l$sc_ls_alcoholreduc_cost_ind <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %like% "alcoholreduc_cost$", as.numeric(value1)]
    l$sc_ls_alcoholreduc_cost_ovrhd <-
      parameters_dt[true_scenario == scenario_nam &
                      input_names %flike% "alcoholreduc_cost_ovrhd", as.numeric(value1)]


  }
  return(l)
}
