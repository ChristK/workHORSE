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

# TODO adjust RR for incidence (not mortality)

mc_max <- 1e3
ageL <- 30L
ageH <- 89L


cat("Initialising workHORSE model...\n\n")
if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
  library(CKutils)
}
dependencies(c("data.table", "future", "fst"))

options(
  future.fork.enable = TRUE, # enable fork in Rstudio
  future.globals.maxSize = +Inf,
  future.globals.onReference = "ignore")
plan(multiprocess, workers = 20L)

# dqRNGkind("pcg64") # dqRNGkind("Xoroshiro128+") ~10% faster
SEED <- 4925886L # sample(1e7, 1)
set.seed(SEED) # Current is to ensure reproducibility.
# dqset.seed(SEED)# Ensure that seed hasn't been used elsewher in the model

stochRRtabl <- # need to run by id
  function(m, ci, stochastic = TRUE) {
    # lognormal
    if (stochastic)
      kk <- stats::runif(1)
    else
      kk <- 0.5
    rr <- exp(stats::qnorm(kk, log(m), abs(log(m) - log(ci)) / 1.96))
    rr[!is.finite(rr)] <- 1 # fix for rtruncnorm above
    return(rr)
  }

generate_scalar_rr_l <- function(suffix, mean_rr, ci_rr, mc_max,
                                 distr = c("lognorm", "norm")) {
  distr <- match.arg(distr)
  if (distr == "lognorm") {
    dt <- data.table(
      1:mc_max,
      exp(stats::qnorm(runif(mc_max), log(mean_rr), abs(log(mean_rr) - log(ci_rr)) / 1.96)),
      mean_rr)
  } else {
    dt <- data.table(
      1:mc_max,
      stats::qnorm(runif(mc_max), mean_rr, abs(mean_rr - ci_rr) / 1.96),
      mean_rr)
  }

  colnam <- paste0("rr_", suffix)
  setnames(dt, c("mc", colnam, paste0("mean_rr_", suffix)))
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")

  setkey(dt, mc)
  if (mean_rr >= 1) dt[get(colnam) < 1, (colnam) := 1]
  if (mean_rr <= 1) dt[get(colnam) > 1, (colnam) := 1]

  write_fst(dt, filenam, 100)
}


suffix <- "af_stroke"
tt <- setDT(expand.grid(age = ageL:ageH, sex = c("men", "women")))
to_agegrp(tt, 5L, 80L, "age", "agegroup")
dt <-
  fread(
    "./simulation/rr/input/af.rrstroke.csv",
    stringsAsFactors = TRUE,
    key = c("agegroup")
  )

generate_rr_l <- function(dt, tt, suffix, mc_max) {
  # uses no interpolation or smoothing
  setnames(dt, gsub("\\.", "_", names(dt)))
  if (dt[, all(mean_rr >= 1)]) {
    constrain <- "above_1"
  } else if (dt[, all(mean_rr <= 1)]) {
    constrain <- "below_1"
  } else constrain <- "mix"
  nam <- paste0("mean_rr_", suffix)
  colnam <- paste0("rr_", suffix)
  colby <- "mc"
  setnames(dt, "mean_rr", nam)
  if ("sex" %in% names(dt)) {
    dt[, sex := factor(sex, 1:2, c("men", "women"))]
    colby <- c(colby, "sex")
  }
  dt[agegroup == "80-84", agegroup := "80+"]
  dt <- clone_dt(dt, mc_max, "mc")
  dt[, (colnam) := stochRRtabl(get(nam), ci_rr, TRUE), by = mc]
  dt[is.na(get(colnam)), (colnam) := 1]
  # dt[, .(unique(get(nam)), median(get(colnam))), by = .(agegroup, sex)]
  dt <- na.omit(tt[dt, on = .NATURAL, allow.cartesian = TRUE])
  # dt[mc == 1, plot(age, get(colnam), col = sex, pch = ".", main = (colnam))]
  dt[, agegroup := NULL]
  setkey(dt, mc, age)
  setcolorder(dt)
  if (constrain == "above_1") dt[get(colnam) < 1, (colnam) := 1]
  if (constrain == "below_1") dt[get(colnam) > 1, (colnam) := 1]
  dt[ci_rr == get(nam), (colnam) := ci_rr] # case when rr == 1 and ci == 1
  dt[, ci_rr := NULL]
  if (dt[, anyNA(get(colnam))]) {
    warning("RR has NAs")
  } else if (dt[, all(get(colnam) >= 1)]) {
    message("All RR above 1")
  } else if (dt[, all(get(colnam) <= 1)]) {
    message("All RR below 1")
  } else warning("RR crosses 1")
  dt
}

generate_rr_l_lin_interp <- function(dt, tt, suffix, mc_max) {
  # uses linear interpolation for interpolation and smoothing
  setnames(dt, gsub("\\.", "_", names(dt)))
  dt[, c("mean_rr_lag", "ci_rr_lag") := shift(.(mean_rr, ci_rr)),
    keyby = eval(ifelse("sex" %in% names(dt), "sex", NULL))]
  dt[mean_rr == mean_rr_lag & ci_rr == ci_rr_lag,
    c("mean_rr", "ci_rr") := NA]
  dt[, c("mean_rr_lag", "ci_rr_lag") := NULL]

  if (dt[, all(mean_rr >= 1, na.rm = TRUE)]) {
    constrain <- "above_1"
  } else if (dt[, all(mean_rr <= 1, na.rm = TRUE)]) {
    constrain <- "below_1"
  } else constrain <- "mix"
  nam <- paste0("mean_rr_", suffix)
  colnam <- paste0("rr_", suffix)
  colby <- "mc"
  setnames(dt, "mean_rr", nam)
  if ("sex" %in% names(dt)) {
    dt[, sex := factor(sex, 1:2, c("men", "women"))]
    colby <- c(colby, "sex")
  }
  dt[agegroup == "80-84", agegroup := "80+"]
  dt <- clone_dt(dt, mc_max, "mc")
  dt[, (colnam) := stochRRtabl(get(nam), ci_rr, TRUE), by = mc]
  dt[is.na(get(nam)), (colnam) := NA]
  # dt[, .(unique(get(nam)), median(get(colnam))), by = .(agegroup, sex)]
  dt <- tt[dt, on = .NATURAL, allow.cartesian = TRUE]
  dt[age != substr(agegroup, start = 1, stop = 2), (colnam) := NA]
  dt[, (colnam) := approx(age, get(colnam), xout = age, rule = 2)$y, by = eval(colby)]
  # dt[mc == 1, plot(age, get(paste0(colnam, 2)), col = sex, pch = "*", main = (colnam))]
  # dt[mc == 1, points(age, get(colnam), col = sex, main = (colnam))]

  dt[, agegroup := NULL]
  setkey(dt, mc, age)
  setcolorder(dt)
  if (constrain == "above_1") dt[get(colnam) < 1, (colnam) := 1]
  if (constrain == "below_1") dt[get(colnam) > 1, (colnam) := 1]
  dt[ci_rr == get(nam), (colnam) := ci_rr] # case when rr == 1 and ci == 1
  dt[, ci_rr := NULL]
  if (dt[, anyNA(get(colnam))]) {
    warning("RR has NAs")
  } else if (dt[, all(get(colnam) >= 1)]) {
    message("All RR above 1")
  } else if (dt[, all(get(colnam) <= 1)]) {
    message("All RR below 1")
  } else warning("RR crosses 1")
  dt
}

generate_rr_l_loess <- function(dt, tt, suffix, mc_max) {
  # uses loess for interpolation and smoothing
  setnames(dt, gsub("\\.", "_", names(dt)))
  if (dt[, all(mean_rr >= 1)]) {
    constrain <- "above_1"
  } else if (dt[, all(mean_rr <= 1)]) {
    constrain <- "below_1"
  } else constrain <- "mix"
  nam <- paste0("mean_rr_", suffix)
  colnam <- paste0("rr_", suffix)
  colby <- "mc"
  setnames(dt, "mean_rr", nam)
  if ("sex" %in% names(dt)) {
    dt[, sex := factor(sex, 1:2, c("men", "women"))]
    colby <- c(colby, "sex")
    }
  dt[agegroup == "80-84", agegroup := "80+"]
  dt <- clone_dt(dt, mc_max, "mc")
  dt[, (colnam) := stochRRtabl(get(nam), ci_rr, TRUE), by = mc]
  dt[is.na(get(colnam)), (colnam) := 1]
  # dt[, .(unique(get(nam)), median(get(colnam))), by = .(agegroup, sex)]
  dt <- na.omit(tt[dt, on = .NATURAL, allow.cartesian = TRUE])
  dt[, (colnam) := predict(loess(as.formula(paste0(colnam, " ~ age")), .SD,
                                 span = 0.75)), by = eval(colby)]
  # print(dt[, plot(age, get(colnam), col = sex, pch = ".", main = (colnam))])
  dt[, agegroup := NULL]
  setkey(dt, mc, age)
  setcolorder(dt)
  if (constrain == "above_1") dt[get(colnam) < 1, (colnam) := 1]
  if (constrain == "below_1") dt[get(colnam) > 1, (colnam) := 1]
  dt[ci_rr == get(nam), (colnam) := ci_rr] # case when rr == 1 and ci == 1
  dt[, ci_rr := NULL]
  if (dt[, anyNA(get(colnam))]) {
    warning("RR has NAs")
  } else if (dt[, all(get(colnam) >= 1)]) {
    message("All RR above 1")
  } else if (dt[, all(get(colnam) <= 1)]) {
    message("All RR below 1")
  } else warning("RR crosses 1")
  dt
}

# calculate 95% CI of a ratio from a p value
ci_from_p_forratio <- function(mean_rr, p) {
  z <- -0.862+ sqrt(0.743 - 2.404*log(p))
  est <- log(mean_rr)
  se <- abs(est/z)
  ci <- list()
  ci$lower_ci <- exp(est - 1.96 * se)
  ci$upper_ci <- exp(est + 1.96 * se)
  ci
}

# generate_rr_l <- generate_rr_l_lin_interp

# AF risk (Christiansen et al. 2016). For no risk factor and many limitations
future({
  suffix <- "af_stroke"
  tt <- setDT(expand.grid(age = ageL:ageH, sex = c("men", "women")))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/af.rrstroke.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "sbp_chd"
  tt <- setDT(expand.grid(age = ageL:ageH, sex = c("men", "women")))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/sbp.rrchd.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "sbp_stroke"
  tt <- setDT(expand.grid(age = ageL:ageH, sex = c("men", "women")))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/sbp.rrstroke.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "bmi_chd"
  tt <- setDT(expand.grid(age = ageL:ageH))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/bmi.rrchd.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "bmi_stroke"
  tt <- setDT(expand.grid(age = ageL:ageH))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/bmi.rrstroke.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "tchol_chd"
  tt <- setDT(expand.grid(age = ageL:ageH))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/tchol.rrchd.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "tchol_stroke"
  tt <- setDT(expand.grid(age = ageL:ageH))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/tchol.rrstroke.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "t2dm_chd"
  tt <- setDT(expand.grid(age = ageL:ageH))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/t2dm.rrchd.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "t2dm_stroke"
  tt <- setDT(expand.grid(age = ageL:ageH))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/t2dm.rrstroke.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "pa_chd"
  tt <- setDT(expand.grid(age = ageL:ageH, active_days = 0:4))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/pa.rrchd.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "pa_stroke"
  tt <- setDT(expand.grid(age = ageL:ageH, active_days = 0))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/pa.rrstroke.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "pa_colon_ca"
  tt <- setDT(expand.grid(age = ageL:ageH, active_days = 0:4))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/pa.rrcolon_ca.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "pa_breast_ca"
  tt <- setDT(expand.grid(age = ageL:ageH, active_days = 0:4))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/pa.rrbreast_ca.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "tobacco_chd"
  tt <- setDT(expand.grid(age = ageL:ageH, sex = c("men", "women"), smok_status = factor(3:4)))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/tobacco.rrchd.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )[, smok_status := factor(smok_status)]
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  suffix <- "tobacco_stroke"
  tt <- setDT(expand.grid(age = ageL:ageH, sex = c("men", "women"), smok_status = factor(4)))
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/tobacco.rrstroke.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )[, smok_status := factor(smok_status)]
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})
future({
  # This is not a RR but a probability to develop dementia within a year of
  # having stroke
  suffix <- "stroke_dementia"
  tt <- data.table("age" = ageL:ageH)
  to_agegrp(tt, 5L, 80L, "age", "agegroup")
  dt <-
    fread(
      "./simulation/rr/input/stroke.rrdementia.csv",
      stringsAsFactors = TRUE,
      key = c("agegroup")
    )
  dt <- generate_rr_l(dt, tt, suffix, mc_max)
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})

future({
  generate_scalar_rr_l("ets_chd", 1.26, 1.38, mc_max)
  generate_scalar_rr_l("ets_stroke", 1.25, 1.38, mc_max)
  generate_scalar_rr_l("ets_copd", 1.66, 2.00, mc_max)
  generate_scalar_rr_l("ets_lung_ca", 1.33, 1.54, mc_max)
  generate_scalar_rr_l("ets_breast_ca", 1.072, 1.126, mc_max)
  generate_scalar_rr_l("fv_chd", 0.96, 0.99, mc_max)
  generate_scalar_rr_l("fv_stroke", 0.95, 0.97, mc_max)
  generate_scalar_rr_l("fv_stroke", 0.95, 0.97, mc_max)
  generate_scalar_rr_l("fruit_lung_ca", 0.96^0.8, 0.98^0.8, mc_max) # ^0.8 to adjust risk of reported 100gr servings to 80gr UK servings
  generate_scalar_rr_l("t2dm_colon_ca", 1.527, 2.304, mc_max)
  generate_scalar_rr_l("t2dm_breast_ca", 1.513, 2.206, mc_max)
  generate_scalar_rr_l("tobacco_nonmodelled", 1.99, 2.14, mc_max)
  generate_scalar_rr_l("t2dm_nonmodelled", 2.08, 2.26, mc_max)
  generate_scalar_rr_l("pa_nonmodelled", 1.34, 1.47, mc_max) # physical inactivity
  generate_scalar_rr_l("alcohol_nonmodelled", 1.46, 1.65, mc_max)
  generate_scalar_rr_l("sbp_nonmodelled", 1.29, 1.38, mc_max)
  generate_scalar_rr_l("sbp_nonmodelled", 1.29, 1.38, mc_max)
  generate_scalar_rr_l("statins_t2dm", 1.09, 1.17, mc_max)
  generate_scalar_rr_l("statins_tchol", 0.27 * 0.43 / 0.36,
    0.27 * 0.3958 / 0.36, mc_max, "norm")


  #	Inflate all_cause (proxy to non_cvd) mortality for hypertensives
  # RR (= 1.3) roughly based on stringhini_socioeconomic_2017 figure 4

  # theoretical minimum distribution.
  # level of sbp below no risk exist. From Singh et al
  suffix <- "sbp_tmred_chd"
  dt <- data.table(
    1:mc_max,
    rnorm(
      mc_max,
      mean = runif(mc_max, 110, 115),
      sd = runif(mc_max, 4, 6)),
    112.5)
  colnam <- paste0("rr_", suffix)
  setnames(dt, c("mc", colnam, paste0("mean_rr_", suffix)))
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  setkey(dt, mc)
  write_fst(dt, filenam, 100)

  # theoretical minimum distribution.
  # level of sbp below no risk exist. From Singh et al
  suffix <- "sbp_tmred_stroke"
  dt <- data.table(
    1:mc_max,
    rnorm(
      mc_max,
      mean = runif(mc_max, 110, 115),
      sd = runif(mc_max, 4, 6)),
    112.5)
  colnam <- paste0("rr_", suffix)
  setnames(dt, c("mc", colnam, paste0("mean_rr_", suffix)))
  filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
  setkey(dt, mc)
  write_fst(dt, filenam, 100)
  # get_rr_mc(1, "chd", "sbp_tmred", TRUE)[]
})



# Ethnicity RR for CVD
# NOTE: the paper below is for 12 condition that together comprises CVD and studies
# time to 1st event. So the RR of progression from CHD to stroke is not considered.
# Nevertheless under the assumption that CHD and stoke incidence are independent,
# (but depend on common risk factors) I can use it here.
# George J, Mathur R, Shah AD, Pujades-Rodriguez M, Denaxas S, Smeeth L, et al.
# Ethnicity and the first diagnosis of a wide range of cardiovascular diseases:
# Associations in a linked electronic health record cohort of 1 million patients.
# PLOS ONE. 2017 Jun 9;12(6):e0178945. fig(3, fully adjusted for risk factors
# and medication)
# The paper has separate risks for (un)stable angina, CHD and MI. Below I will
# combine them using MC to get a RR for the combined CHD, weighted for each
# separate condition. Results for stroke were not statistically significant

future({
  # CHD
  # South Asian
  mf <- 10L
  # Stable angina (17793+477 = 18270 events)
  weight_sa <- 18270L * mf
  rr_sa <- exp(rnorm(weight_sa, log(1.61), abs(log(1.46) - log(1.78))/3.92))
  # summary(rr_sa)
  # quantile(rr_sa, c(0.025, 0.5, 0.975))

  # Unstable angina (5739 + 193 =  5932 events)
  weight_ua <- 5932L * mf
  rr_ua <- exp(rnorm(weight_sa, log(1.80), abs(log(1.54) - log(2.11))/3.92))

  # CHD NOS (8140 + 244  = 8384  events)
  weight_nos <- 8384L * mf
  rr_nos <- exp(rnorm(weight_sa, log(1.63), abs(log(1.42) - log(1.87))/3.92))

  # MI (14535 + 350 = 14885  events)
  weight_mi <- 14885L * mf
  rr_mi <- exp(rnorm(weight_sa, log(1.61), abs(log(1.44) - log(1.81))/3.92))

  # UCD (3420 + 51 = 3471  events) was non significant
  weight_ucd <- 3471L * mf
  rr_ucd <- exp(rnorm(weight_sa, log(1.09), abs(log(0.81) - log(1.46))/3.92))

  rr_chd <- sample(c(rr_sa, rr_ua, rr_nos, rr_mi), mc_max, FALSE)
  quantile(rr_chd, c(0.025, 0.5, 0.975))
  dt1 <- data.table(
    mc = seq_along(rr_chd),
    ethnicity = "pakistani",
    rr_ethnicity_chd = rr_chd,
    mean_rr_ethnicity_chd = mean(rr_chd)
  )

  dt2 <- copy(dt1)
  dt2[, ethnicity := "indian"]
  dt3 <- copy(dt1)
  dt3[, ethnicity := "bangladeshi"]
  dt4 <- copy(dt1)
  dt4[, ethnicity := "other asian"]


  # CHD
  # Black
  # Stable angina (17793 + 171 = 18441 events)
  weight_sa <- 17964L * mf
  rr_sa <- exp(rnorm(weight_sa, log(0.71), abs(log(0.60) - log(0.83))/3.92))

  # Unstable angina (5739 + 72 =  5811 events)
  weight_ua <- 5811L * mf
  rr_ua <- exp(rnorm(weight_sa, log(0.70), abs(log(0.55) - log(0.90))/3.92))

  # CHD NOS (8140 + 82 = 8222  events)
  weight_nos <- 8222L * mf
  rr_nos <- exp(rnorm(weight_sa, log(0.60), abs(log(0.48) - log(0.76))/3.92))

  # MI (14535 + 84 = 14619  events)
  weight_mi <- 14619L * mf
  rr_mi <- exp(rnorm(weight_sa, log(0.47), abs(log(0.37) - log(0.58))/3.92))

  # UCD (3420 + 28 = 3448  events) was non significant
  weight_ucd <- 3448L * mf
  rr_ucd <- exp(rnorm(weight_sa, log(0.70), abs(log(0.48) - log(1.04))/3.92))

  rr_chd <- sample(c(rr_sa, rr_ua, rr_nos, rr_mi), mc_max, FALSE)
  quantile(rr_chd, c(0.025, 0.5, 0.975))

  dt5 <- data.table(
    mc = seq_along(rr_chd),
    ethnicity = "black caribbean",
    rr_ethnicity_chd = rr_chd,
    mean_rr_ethnicity_chd = mean(rr_chd)
  )

  dt6 <- copy(dt5)
  dt6[, ethnicity := "black african"]
  dt <- rbind(dt1, dt2, dt3, dt4, dt5, dt6)

  dt[, ethnicity := factor(ethnicity, c("white", "indian", "pakistani",
    "bangladeshi", "other asian", "black caribbean", "black african", "chinese",
    "other"   ))]
  setkey(dt, mc)

  filenam <- "./simulation/rr/ethnicity_chd_rr_l.fst"
  filenam_indx <- "./simulation/rr/ethnicity_chd_rr_indx.fst"
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})

# Alcohol from GBD 2017
# formulas from above cannot be used because of J-shape relationship
# CHD
future({
  tt <- fread("./simulation/rr/input/alcohol.rrchd.csv")
  tt <- melt(tt, 1:2, variable.name = "agegroup")
  tt[, c("mean_rr", "ci_rr") := tstrsplit(value, " ")[c(1, 4)]]
  tt[, c("ci_rr") := sub(")", "", ci_rr)]
  tt[, value := NULL]
  tt[, mean_rr := as.numeric(mean_rr)]
  tt[, ci_rr := as.numeric(ci_rr)]

  tt[, uniqueN(mean_rr), keyby = .(alcohol, sex)] # No variation by age
  tt <- tt[, .("mean_rr" = unique(mean_rr), "ci_rr" = unique(ci_rr)), keyby = .(alcohol, sex)]

  ttt <- data.table(alcohol = 0:72)
  replace_from_table(ttt, "alcohol", 0:72, c(rep(c(0, 12, 24, 36, 48, 60), each = 12),  72), "alcohol_grp")
  dt <- ttt[tt, on = "alcohol_grp==alcohol", allow.cartesian = TRUE]
  dt[, alcohol_grp := NULL]
  dt[, sex := factor(sex, 1:2, c("men", "women"))]

  dt <- clone_dt(dt, mc_max, "mc")
  dt[, rr_alcohol_chd := stochRRtabl(mean_rr, ci_rr, TRUE), by = mc]
  dt[, ci_rr := NULL]
  setnames(dt, "mean_rr", "mean_rr_alcohol_chd")
  dt[is.na(rr_alcohol_chd), .N]
  dt[!alcohol %in% c(12 * (0:6)), rr_alcohol_chd := NA]
  # dt[, rr_alcohol_chd := predict(loess(rr_alcohol_chd ~ alcohol, .SD,
  #   span = 0.7), alcohol), by = .(mc, sex)]
  dt[, rr_alcohol_chd := approx(alcohol, rr_alcohol_chd,
    alcohol)$y, by = .(mc, sex)]

  dt[rr_alcohol_chd < 1, rr_alcohol_chd := 1]

  # dt[mc == 1 & sex == "men", {plot(alcohol, rr_alcohol_chd)
  # lines(alcohol, rr_alcohol_chd_sm)}]
  setkey(dt, mc, alcohol)

  filenam <- "./simulation/rr/alcohol_chd_rr_l.fst"
  filenam_indx <- "./simulation/rr/alcohol_chd_rr_indx.fst"
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})

# Stroke
future({
  tt <- fread("./simulation/rr/input/alcohol.rrstroke.csv")
  tt <- melt(tt, 1:2, variable.name = "agegroup")
  tt[, c("mean_rr", "ci_rr") := tstrsplit(value, " ")[c(1, 4)]]
  tt[, c("ci_rr") := sub(")", "", ci_rr)]
  tt[, value := NULL]
  tt[, mean_rr := as.numeric(mean_rr)]
  tt[, ci_rr := as.numeric(ci_rr)]

  tt[, uniqueN(mean_rr), keyby = .(alcohol, sex)][, table(V1)] # No variation by age
  tt <- tt[, .("mean_rr" = unique(mean_rr), "ci_rr" = unique(ci_rr)), keyby = .(alcohol, sex)]

  ttt <- data.table(alcohol = 0:72)
  replace_from_table(ttt, "alcohol", 0:72, c(rep(c(0, 12, 24, 36, 48, 60), each = 12),  72), "alcohol_grp")
  dt <- ttt[tt, on = "alcohol_grp==alcohol", allow.cartesian = TRUE]
  dt[, alcohol_grp := NULL]
  dt[, sex := factor(sex, 1:2, c("men", "women"))]

  dt <- clone_dt(dt, mc_max, "mc")
  dt[, rr_alcohol_stroke := stochRRtabl(mean_rr, ci_rr, TRUE), by = mc]
  dt[, ci_rr := NULL]
  setnames(dt, "mean_rr", "mean_rr_alcohol_stroke")
  dt[is.na(rr_alcohol_stroke), .N]
  dt[!alcohol %in% c(12 * (0:6)), rr_alcohol_stroke := NA]
  # dt[, rr_alcohol_stroke := predict(loess(rr_alcohol_stroke ~ alcohol, .SD,
  #   span = 0.7), alcohol), by = .(mc, sex)]
  dt[, rr_alcohol_stroke := approx(alcohol, rr_alcohol_stroke,
    alcohol)$y, by = .(mc, sex)]

  dt[rr_alcohol_stroke < 1, rr_alcohol_stroke := 1]

  # dt[mc == 10 & sex == "men", {plot(alcohol, rr_alcohol_stroke)
  # lines(alcohol, rr_alcohol_stroke_sm)}]
  setkey(dt, mc, alcohol)

  filenam <- "./simulation/rr/alcohol_stroke_rr_l.fst"
  filenam_indx <- "./simulation/rr/alcohol_stroke_rr_indx.fst"
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})


# Colon_ca
# Alcohol
future({
  tt <- fread("./simulation/rr/input/alcohol.rrcolon_ca.csv")
  tt <- melt(tt, 1, variable.name = "agegroup")
  tt[, c("mean_rr", "ci_rr") := tstrsplit(value, " ")[c(1, 4)]]
  tt[, c("ci_rr") := sub(")", "", ci_rr)]
  tt[, value := NULL]
  tt[, mean_rr := as.numeric(mean_rr)]
  tt[, ci_rr := as.numeric(ci_rr)]
  tt[, uniqueN(mean_rr), keyby = .(alcohol)][, table(V1)] # No variation by age
  tt <- tt[, .("mean_rr" = unique(mean_rr), "ci_rr" = unique(ci_rr)), keyby = .(alcohol)]

  ttt <- data.table(alcohol = 0:72)
  replace_from_table(ttt, "alcohol", 0:72, c(rep(c(0, 12, 24, 36, 48, 60), each = 12),  72), "alcohol_grp")
  dt <- ttt[tt, on = "alcohol_grp==alcohol", allow.cartesian = TRUE]
  dt[, alcohol_grp := NULL]

  dt <- clone_dt(dt, mc_max, "mc")
  dt[, rr_alcohol_colon_ca := stochRRtabl(mean_rr, ci_rr, TRUE), by = mc]
  dt[, ci_rr := NULL]
  setnames(dt, "mean_rr", "mean_rr_alcohol_colon_ca")
  dt[is.na(rr_alcohol_colon_ca), .N]
  dt[!alcohol %in% c(12 * (0:6)), rr_alcohol_colon_ca := NA]
  # dt[, rr_alcohol_colon_ca := predict(loess(rr_alcohol_colon_ca ~ alcohol, .SD,
  #   span = 0.7), alcohol), by = .(mc)]
  dt[, rr_alcohol_colon_ca := approx(alcohol, rr_alcohol_colon_ca,
    alcohol)$y, by = .(mc)]

  dt[rr_alcohol_colon_ca < 1, rr_alcohol_colon_ca := 1]

  # dt[mc == 10 & sex == "men", {plot(alcohol, rr_alcohol_colon_ca)
  # lines(alcohol, rr_alcohol_colon_ca_sm)}]
  setkey(dt, mc, alcohol)

  filenam <- "./simulation/rr/alcohol_colon_ca_rr_l.fst"
  filenam_indx <- "./simulation/rr/alcohol_colon_ca_rr_indx.fst"
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})

# Breast_ca
future({
  tt <- fread("./simulation/rr/input/alcohol.rrbreast_ca.csv")
  tt <- melt(tt, 1, variable.name = "agegroup")
  tt[, c("mean_rr", "ci_rr") := tstrsplit(value, " ")[c(1, 4)]]
  tt[, c("ci_rr") := sub(")", "", ci_rr)]
  tt[, value := NULL]
  tt[, mean_rr := as.numeric(mean_rr)]
  tt[, ci_rr := as.numeric(ci_rr)]
  tt[, uniqueN(mean_rr), keyby = .(alcohol)][, table(V1)] # No variation by age
  tt <- tt[, .("mean_rr" = unique(mean_rr), "ci_rr" = unique(ci_rr)), keyby = .(alcohol)]

  ttt <- data.table(alcohol = 0:72)
  replace_from_table(ttt, "alcohol", 0:72, c(rep(c(0, 12, 24, 36, 48, 60), each = 12),  72), "alcohol_grp")
  dt <- ttt[tt, on = "alcohol_grp==alcohol", allow.cartesian = TRUE]
  dt[, alcohol_grp := NULL]

  dt <- clone_dt(dt, mc_max, "mc")
  dt[, rr_alcohol_breast_ca := stochRRtabl(mean_rr, ci_rr, TRUE), by = mc]
  dt[, ci_rr := NULL]
  setnames(dt, "mean_rr", "mean_rr_alcohol_breast_ca")
  dt[is.na(rr_alcohol_breast_ca), .N]
  dt[!alcohol %in% c(12 * (0:6)), rr_alcohol_breast_ca := NA]
  # dt[, rr_alcohol_breast_ca := predict(loess(rr_alcohol_breast_ca ~ alcohol, .SD,
  #   span = 0.7), alcohol), by = .(mc)]
  dt[, rr_alcohol_breast_ca := approx(alcohol, rr_alcohol_breast_ca,
    alcohol)$y, by = .(mc)]
  dt[rr_alcohol_breast_ca < 1, rr_alcohol_breast_ca := 1]

  # dt[mc == 10, {plot(alcohol, rr_alcohol_breast_ca)
  # # lines(alcohol, rr_alcohol_breast_ca_sm)
  # }]
  dt[, sex := factor("women", c("men", "women"))]
  setkey(dt, mc, alcohol)

  filenam <- "./simulation/rr/alcohol_breast_ca_rr_l.fst"
  filenam_indx <- "./simulation/rr/alcohol_breast_ca_rr_indx.fst"
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})


# Tobacco (packyears)
# COPD
future({
  tt <- fread("./simulation/rr/input/packyears.rrcopd.csv")

  ttt <- data.table(packyears = 0:60)
  # replace_from_table(ttt, "alcohol", 0:60, c(rep(c(0, 12, 24, 36, 48, 60), each = 12),  72), "alcohol_grp")
  dt <- tt[ttt, on = "packyears", nomatch = NA]
  dt[, lci_rr := NULL]
  setnames(dt, "uci_rr", "ci_rr")

  dt <- clone_dt(dt, mc_max, "mc")
  dt[!is.na(mean_rr), rr_packyears_copd := stochRRtabl(mean_rr, ci_rr, TRUE), by = mc]
  dt[, ci_rr := NULL]
  setnames(dt, "mean_rr", "mean_rr_packyears_copd")
  setnafill(dt, "locf", cols = "mean_rr_packyears_copd")
  # dt[, rr_packyears_copd := predict(loess(rr_packyears_copd ~ packyears, .SD,
  #   span = 0.7), packyears), by = .(mc)]
  dt[, rr_packyears_copd := approx(packyears, rr_packyears_copd,
    packyears)$y, by = .(mc)]
  dt[packyears == 0, rr_packyears_copd := 1]
  dt[rr_packyears_copd < 1, rr_packyears_copd := 1]

  # dt[mc == 536, {
  #   plot(packyears, rr_packyears_copd)
  #   lines(packyears, rr_packyears_copd)
  # }]
  setkey(dt, mc, packyears)

  filenam <- "./simulation/rr/packyears_copd_rr_l.fst"
  filenam_indx <- "./simulation/rr/packyears_copd_rr_indx.fst"
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})

# Tobacco (packyears)
# Colon ca
future({
  tt <- data.table("packyears" = c(0, 35, 60), "mean_rr" = c(1, 1.269, 1.505),
    "ci_rr" = c(1, ci_from_p_forratio(1.269, 0.0001)[[1]],
      ci_from_p_forratio(1.505, 0.0001)[[1]]))
  ttt <- data.table(packyears = 0:60)
  # replace_from_table(ttt, "alcohol", 0:60, c(rep(c(0, 12, 24, 36, 48, 60), each = 12),  72), "alcohol_grp")
  dt <- tt[ttt, on = "packyears", nomatch = NA]

  dt <- clone_dt(dt, mc_max, "mc")
  dt[!is.na(mean_rr), rr_packyears_colon_ca := stochRRtabl(mean_rr, ci_rr, TRUE), by = mc]
  dt[, ci_rr := NULL]
  setnames(dt, "mean_rr", "mean_rr_packyears_colon_ca")
  setnafill(dt, "locf", cols = "mean_rr_packyears_colon_ca")
  dt[, rr_packyears_colon_ca := approx(packyears, rr_packyears_colon_ca,
    packyears)$y, by = .(mc)]
  dt[packyears == 0, rr_packyears_colon_ca := 1]
  dt[rr_packyears_colon_ca < 1, rr_packyears_colon_ca := 1]

  # dt[mc == 15, {
  #   plot(packyears, rr_packyears_colon_ca)
  #   lines(packyears, rr_packyears_colon_ca)
  # }]
  setkey(dt, mc, packyears)

  filenam <- "./simulation/rr/packyears_colon_ca_rr_l.fst"
  filenam_indx <- "./simulation/rr/packyears_colon_ca_rr_indx.fst"
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})

# BMI
# Colon_ca
future({
  tt <- fread("./simulation/rr/input/bmi.rrcolon_ca.csv")
  tt <- melt(tt, 1, variable.name = "agegroup")
  tt[, c("mean_rr", "ci_rr") := tstrsplit(value, " ")[c(1, 4)]]
  tt[, c("ci_rr") := sub(")", "", ci_rr)]
  tt[, value := NULL]
  tt[, mean_rr := as.numeric(mean_rr)]
  tt[, ci_rr := as.numeric(ci_rr)]
  tt[, uniqueN(mean_rr), keyby = .(sex)][, table(V1)] # No variation by age
  tt <- tt[, .("mean_rr" = unique(mean_rr), "ci_rr" = unique(ci_rr)), keyby = .(sex)]

  dt <- clone_dt(tt, mc_max, "mc")
  dt[!is.na(mean_rr), rr_bmi_colon_ca := stochRRtabl(mean_rr, ci_rr, TRUE), by = mc]
  dt[, ci_rr := NULL]
  setnames(dt, "mean_rr", "mean_rr_bmi_colon_ca")

  setkey(dt, mc, sex)

  filenam <- "./simulation/rr/bmi_colon_ca_rr_l.fst"
  filenam_indx <- "./simulation/rr/bmi_colon_ca_rr_indx.fst"
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})

# BMI
# Breast_ca
future({
  tt <- fread("./simulation/rr/input/bmi.rrbreast_ca.csv")
  tt <- melt(tt, 1, variable.name = "agegroup")
  tt[, c("mean_rr", "ci_rr") := tstrsplit(value, " ")[c(1, 4)]]
  tt[, c("ci_rr") := sub(")", "", ci_rr)]
  tt[, value := NULL]
  tt[, mean_rr := as.numeric(mean_rr)]
  tt[, ci_rr := as.numeric(ci_rr)]
  tt[, uniqueN(mean_rr), keyby = .(sex)][, table(V1)] # varies by age
  tt[agegroup == "80-84", agegroup := "80+"]

  ttt <- setDT(expand.grid(age = ageL:ageH, sex = c("men", "women")))
  to_agegrp(ttt, 5L, 80L, "age", "agegroup")
  absorb_dt(ttt, tt)

  dt <- clone_dt(ttt, mc_max, "mc")
  dt[!is.na(mean_rr), rr_bmi_breast_ca := stochRRtabl(mean_rr, ci_rr, TRUE), by = mc]
  dt[, ci_rr := NULL]
  setnames(dt, "mean_rr", "mean_rr_bmi_breast_ca")
  # No smoothing as the RR are for pre and post menopausal women
  dt[, agegroup := NULL]
  setkey(dt, mc, sex)

  filenam <- "./simulation/rr/bmi_breast_ca_rr_l.fst"
  filenam_indx <- "./simulation/rr/bmi_breast_ca_rr_indx.fst"
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})

# Tobacco
# Breast_ca
future({
  dt <- data.table("smok_status" = 2:4,
    mean_rr = c(1.09, 1.09, 1.11),
    ci_rr = c(1.12, 1.12, 1.16))

  dt <- clone_dt(dt, mc_max, "mc")
  dt[!is.na(mean_rr), rr_tobacco_breast_ca := stochRRtabl(mean_rr, ci_rr, TRUE), by = mc]
  setnames(dt, "mean_rr", "mean_rr_tobacco_breast_ca")
  dt[, smok_status := factor(smok_status, 1:4)]

  setkey(dt, mc, smok_status)

  filenam <- "./simulation/rr/tobacco_breast_ca_rr_l.fst"
  filenam_indx <- "./simulation/rr/tobacco_breast_ca_rr_indx.fst"
  write_fst(dt, filenam, 100)
  # create a table with row numbers for each mc
  dt[, rn := .I]
  tt <- dt[, .(from = min(rn), to = max(rn)), keyby = mc]
  write_fst(tt, filenam_indx, 100L)
})

# xx <- read_fst( "./simulation/rr/sbp_chd_rr_l.fst", as.data.table = T)
