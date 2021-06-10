source("./global.R")
parameters <- qread("./output/parameters.qs")
input <- qread("./output/input.qs")
parameters <- fromGUI_prune(parameters) # TODO delete for production
design$update_fromGUI(parameters)
parameters_dt <- fromGUI_to_dt(parameters)
scenario_nam <- "sc2"
mc_iter <- 1L

design$get_lags(mc_iter)
dt <- SynthPop$new(mc_iter, design) #run until ncc , with problems on 5(frt, veg, smok all, hdl, alcohol)
dt0 <- SynthPop$new(0, design)

# basic_sc_nam <-
#   parameters_dt[true_scenario == scenario_nam, unique(true_scenario), keyby = scenario]$scenario
# scenario_parms <-
#   lapply(basic_sc_nam, fromGUI_scenario_parms, parameters_dt)
# names(scenario_parms) <- basic_sc_nam
# # Sort scenarios in chronological order (important for serial ensembles)
# basic_sc_nam <-
#   names(sort(unlist(
#     lapply(scenario_parms, `[[`, "sc_init_year")
#   )))
# scenario_parms <- scenario_parms[basic_sc_nam]
# for (sc in basic_sc_nam) print(sc)
# 
# scenario_parms[[sc]]$sc_soc_qimd1_change <- 5L
# scenario_parms[[sc]]$sc_soc_qimd_fatality_change <- TRUE
# scenario_parms[[sc]]$sc_soc_qimd_rf_change <-
#   c("tob", "ets", "fv", "alc", "pa", "sbp", "bmi", "tchol")
# # scenario_parms <- scenario_parms[[sc]]
# # dt <- dt$pop

output_chunk  <- list()
output_chunk <- sapply(
  sort(fromGUI_scenario_names(parameters_dt)$true_scenario),
  # parameters_dt[grepl("^sc[1-9]", scenario), sort(unique(scenario))],
  run_scenario,
  mc_iter,
  dt,
  parameters_dt,
  design,
  output_chunk,
  USE.NAMES = FALSE
)

invisible(lapply(output_chunk, setDT))
invisible(lapply(output_chunk, function(x) {
  oldnam <- grep("_sc.*$", names(x), value = TRUE)
  newnam <- gsub("_sc.*$", "", oldnam)
  setnames(x, oldnam, newnam)
}))
output_chunk <- rbindlist(output_chunk, idcol = "scenario")

setkey(output_chunk, scenario, pid, year) # for identify_longdead
output_chunk <-
  output_chunk[output_chunk$year >= design$sim_prm$init_year_fromGUI &
      between(output_chunk$age,
        design$sim_prm$ageL,
        design$sim_prm$ageH) &
      !identify_longdead(all_cause_mrtl, pid_mrk),]
output_chunk[, pid_mrk  := mk_new_simulant_markers(pid)]
output_chunk[, scenario := factor(scenario)]

export_smok(mc_iter, output_chunk)

generate_health_econ(output_chunk,
  ceiling(mc_iter /
      design$sim_prm$n_synthpop_aggregation))

output_chunk[, c("ncc", "pid", "income", "education") := NULL]

output_chunk2 <- copy(output_chunk)
library(collapse)
output_chunk <- copy(output_chunk2)

# gen incd
for (nam in grep("_prvl$", names(output_chunk), value = TRUE)) {
  newnam <- gsub("_prvl$", "_incd", nam)
  set(output_chunk,
    NULL,
    newnam,
    fifelse(output_chunk[[nam]] == 1L, 1L, 0L)) # OR 0L^abs((x - 1L)*2L)
}

# No speed benefit
# vars <- grep("_prvl$", names(output_chunk), value = TRUE)
# newnam <- gsub("_prvl$", "_incd", vars)
# for (i in seq_along(vars)) {
#   set(output_chunk,
#     NULL,
#     newnam[[i]],
#     0L^abs((output_chunk[[vars[[i]]]] - 1L) * 2L)
#   )
# }


if (design$sim_prm$logs)
  print("transform prvl/dgn/mrtl to be summed")
invisible(output_chunk[, lapply(.SD, fclamp_int, 0L, 1L, TRUE),
  .SDcols = patterns("_prvl$|_dgn$|_mrtl$")])
if ("lqimd" %in% names(output_chunk)) {
  output_chunk[, ("lqimd") := NULL]
} else {
  output_chunk[, ("nqimd") := NULL]
}

output_chunk[, year := year + 2000L]
to_agegrp(output_chunk, 20L, 89L, "age", "agegrp", TRUE, 30L)
output_chunk20 <- copy(output_chunk)

# Scale-up to ONS population projections
# output_chunk[, sum(wt), keyby = .(year, scenario)]
w <- output_chunk$wt
for (nam in grep("_prvl$|_dgn$|_incd|_mrtl$|_cost$|^eq5d",
  names(output_chunk),
  value = TRUE)) {
  set(output_chunk, NULL, nam, output_chunk[[nam]] * w)
}

# summarise by strata
output_chunk[, c("age", "pid_mrk") := NULL]

output_chunk <-
  output_chunk[, lapply(.SD, sum),
    keyby = eval(design$sim_prm$strata_for_output)]

output_chunk4 <- # TODO
  collapv(output_chunk20,
    by = design$sim_prm$strata_for_output,
    fsum,
    w = "wt")
output_chunk4[, c("age", "pid_mrk") := NULL]
setkeyv(output_chunk4, design$sim_prm$strata_for_output)

setcolorder(output_chunk)
setcolorder(output_chunk4)
all.equal(output_chunk, output_chunk4)

setnames(output_chunk, "wt", "pops")
output_chunk[, mc := mc_aggr]


fwrite_safe(output_chunk, output_dir("intermediate_out.csv"))
