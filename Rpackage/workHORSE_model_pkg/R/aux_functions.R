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


`:=` = function(...)
  NULL # due to NSE notes in R CMD check

.onUnload <- function(libpath) {
  library.dynam.unload("workHORSEmisc", libpath)
}

# Ensures that when fwrites appent file colnames of file to be written, match those already in the file
#' @export
fwrite_safe <- function(x,
                        file = "",
                        append = TRUE,
                        ...) {
  if (append) {
    if (file.exists(file)) {
      col_names_disk <- names(fread(file, nrows = 0))
      col_names_file <- names(x)
      col_names <- outersect(col_names_disk, col_names_file)
      if (length(col_names) > 0)
        x[, (col_names) := NA]
      setcolorder(x, col_names_disk)
    }
  }
  fwrite(x, file, append, ...)
}


#' @export
inflate <- function(x, percentage_rate, year, baseline_year) {
  x * (1 + percentage_rate/100)^(year - baseline_year)
}
# inflate(1000, 3, 2011:2020, 2013)

#' @export
deflate <- function(x, percentage_rate, year, baseline_year) {
  x * (1 - percentage_rate/100)^(year - baseline_year)
}

# Necessary aux functions
# plots mean x by y stratified by z
#' @export
plot_cor <- function(x, y, z, wt, dt) {
  xlab  <- toupper(x)
  ylab  <- toupper(y)
  title <- paste0(ylab, "~", xlab, "|", toupper(z))

  dt[, .(y = weighted.mean(get(y), get(wt), na.rm = TRUE)), keyby = .(x = round(get(x)), z = get(z))
  ][, qplot(
    x,
    y,
    data = na.omit(.SD),
    facets = as.formula("~z"),
    geom = "smooth",
    main = title,
    xlab = xlab,
    ylab = ylab
  )]
}

# helper func for gamlss::fitDistr
#' @export
distr_best_fit <- function(dt, var, wt, distr_family, distr_extra = NULL, pred = FALSE, seed = NULL, trace = TRUE) {
  if (pred) {
    print("Selection based on minimum prediction global deviance")
    if (!is.null(seed)) set.seed(seed)
    lns <- sample(nrow(dt), round(nrow(dt) * 0.8))
    dt_trn   <- dt[lns, ] # train dataset
    dt_crv   <- dt[!lns, ]  # cross-validation dataset
    marg_distr <- gamlss::fitDistPred(dt_trn[[var]], type = distr_family, weights = dt_trn[[wt]],
                                      extra = distr_extra,
                                      try.gamlss = TRUE, trace = trace, newdata = dt_crv[[var]])
  } else {
    print("Selection based on BIC")
    marg_distr <- gamlss::fitDist(dt[[var]], log(nrow(dt)), type = distr_family, weights = dt[[wt]],
                                  extra = distr_extra,
                                  try.gamlss = TRUE, trace = trace)
  }
  marg_distr
}

#' @export
distr_validation <- function(marg_distr, dt, title, discrete = FALSE, smooth = 0.35) {
  params <- vector("list", length(marg_distr$parameters) + 1L)
  names(params) <- c(marg_distr$parameters, "p")
  for (i in marg_distr$parameters)
    params[[i]] <- get(i, marg_distr)

  params$p <- seq(0.001, 0.999, length.out = nrow(dt))
  distr_nam <- marg_distr$family[[1]]
  y <- do.call(paste0("q", distr_nam), params)
  y <- y[between(y, dt[, min(var)], dt[, max(var)])]
  y_wt <- rep(1 / length(y), length(y))
  # validation plots
  # see http://www.csss.washington.edu/files/working-papers/2002/wp27.pdf
  out <- reldist_diagnostics(dt$var, y,
                             dt[, wt / sum(wt)], y_wt,
                             main = title,
                             discrete = discrete,
                             smooth = smooth)
  out
}

#' @export
centile_predictAll <- function(gamlss_model, orig_data, newdata, cent = 0.5) {
  stopifnot("gamlss" %in% class(gamlss_model))
  tt <- predictAll(gamlss_model, newdata, data = orig_data)
  tt$p <- cent
  distr_nam <- attr(tt, "family")[[1]]
  out <- do.call(paste0("q", distr_nam), tt)
  return(out)
}

#' @export
mean_predictAll <- function(gamlss_model, orig_data, newdata) {
  stopifnot("gamlss" %in% class(gamlss_model))
  tt <- predictAll(gamlss_model, newdata, data = orig_data)
  setDT(tt)
  tt[, n := 1e5]
  setcolorder(tt, "n")
  distr_nam <- attr(tt, "family")[[1]]
  out <- apply(tt, 1, function(x) {
    mean(do.call(paste0("r", distr_nam), as.list(tt)))
  })
  return(out)
}

#' @export
plot_synthpop_val <- function(dt, x, grp, wt, title, x_label, standardised_to_grp = c(FALSE, TRUE), print_to_screen = c(FALSE, TRUE)) {
  if (standardised_to_grp) {
    dt[, weight := get(wt) / sum(get(wt)), by = c("type", grp)]
  } else {
    dt[, weight := get(wt) / sum(get(wt)), by = c("type")]
  }

  if ("year" %in% grp) dt[, year := year + 2000L]
  x <- enquo(x)
  p <- ggplot(dt, aes(
    !!x,
    colour = type,
    weight = weight,
    linetype = type
  )) +
    geom_density() +
    facet_grid(grp) +
    xlab(x_label) +
    ggtitle(title)
  if (print_to_screen) print(p)
  suppressWarnings(
    cowplot::ggsave2(paste0(gsub(" ", "_", title), "_density.png"), p,
                     width = 16,
                     height = 9,
                     units = "cm",
                     scale = 2,
                     dpi = 300,
                     # compression = "lzw",
                     path = "./preparatory_work/plots"
    )
  )

  p <- ggplot(dt, aes(
    !!x,
    colour = type,
    weight = weight,
    linetype = type
  )) +
    stat_ecdf() +
    facet_grid(grp) +
    xlab(x_label) +
    ggtitle(title)
  if (print_to_screen) print(p)

  if ("year" %in% grp) dt[, year := year - 2000L]

  suppressWarnings(
    cowplot::ggsave2(paste0(gsub(" ", "_", title), "_cdf.png"), p,
                     width = 16,
                     height = 9,
                     units = "cm",
                     scale = 2,
                     dpi = 300,
                     # compression = "lzw",
                     path = "./preparatory_work/plots"
    )
  )
}

#' @export
validate_gamlss_tbl <- function(dt, gamlss_tbl, mc = 10L, colname, distr_nam = distr_nam) {
  stopifnot(is.data.table(dt), is.data.table(gamlss_tbl), mc >= 1)
  nam_var <- intersect(names(dt), names(gamlss_tbl))
  nam_param <- setdiff(names(gamlss_tbl), nam_var)
  dt[, age := as.integer(age)]
  x <- copy(dt)
  x[, `:=`(type, "Observed")]
  z <- copy(dt)
  z[gamlss_tbl, (nam_param) := mget(nam_param), on = nam_var]
  if (z[is.na(mu), .N] > 0) stop("NAs produced in join")
  z[, `:=`(type, "Modelled")]
  z <- rbindlist(rep(list(z), mc))
  z[, p := dqrunif(.N, 0, 0.999)]
  z[, (colname) := do.call((distr_nam), .SD), .SDcols = c("p", nam_param)]
  z[, c("p", (nam_param)) := NULL]
  out <- rbind(x, z, use.names = TRUE, fill = TRUE)
}

# dt <- dt[age >= 20L]
# gamlss_tbl <- copy(smok_incid_model_tbl)
# colname <- "smok_incid"
# z[, .SD, .SDcols = c("p", nam_param)]
# z[, do.call("qBI", .SD), .SDcols = c("p", nam_param)]

#' @export
shift_bypid <-
  function(x, lag, id, replace = NA) {
    if (typeof(x) == "integer") {
      return(shift_bypidInt(x, lag, replace, id))
    } else if (typeof(x) == "logical") {
      return(shift_bypidBool(x, lag, replace, id))
    } else if (typeof(x) == "double") {
      return(shift_bypidNum(x, lag, replace, id))
    } else
      stop("class of x not supported")
  }

#' @export
lung_ca_rr <- function(
  smok_status = 3,
  num_cig = 3,
  smok_dur = 3,
  quit_yrs = 40,
  age = 62,
  education = 4L,
  ethnicity = 1L,
  bmi = 27,
  copd = FALSE,
  personal_history_of_cancer = FALSE,
  family_history_of_lung_cancer = FALSE) {

  ethn_vec <- c(
    0,         # "white"
    -0.5241286, # "indian"
    -0.5241286, # "pakistani"
    -0.5241286, # "bangladeshi"
    -0.5241286, # "other asian"
    0.3211605, # "black caribbean"
    0.3211605, # "black african"
    -0.5241286, # "chinese"
    0          # "other"

  )
  # Original education in PLCO
  # 1=less than high school completed; 2=high school graduate;
  # 3=post high school training; 4=some college; 5=college graduate;
  # 6=postgraduate or professional degree
  edu_vec <- c(5L, 5L, 2L, 2L, 1L, 1L, 1L)
  smok_vec <- c(0, 2.542472, 2.542472, 2.799727)

  # prb_a is the probability of lung cancer if not a smoker (assumes also no COPD)
  a_nocopd <- -7.02198 + 0.079597 * (age - 62) -
    0.0879289 * (edu_vec[education] - 4) -
    0.028948 * (bmi - 27) + ethn_vec[ethnicity] +
    0.4845208 * personal_history_of_cancer +
    0.5856777 * family_history_of_lung_cancer

  # a_copd       <- a_nocopd + 0.3457265 # needs special treatment as it depends on another disease
  # a            <- a_nocopd + 0.3457265 * copd
  exp_a_nocopd <- exp(a_nocopd)
  # exp_a_copd   <- exp(a_copd) # TODO: simplification as if not a smoker the likelyhood of COPD is decreasing but there is residual
  # exp_a        <- exp(a)

  prb_a_nocopd <- exp_a_nocopd / (1 + exp_a_nocopd)
  # prb_a_copd <- exp_a_copd / (1 + exp_a_copd)
  # prb_a <- exp_a / (1 + exp_a)

  # prb_b is the probability of lung cancer for observed smoking status
  x <- fifelse(smok_status != 1L, {
    smok_vec[smok_status] +
      (-0.1815486 * (((num_cig / 100) ^ -1) - 4.021541613)) +
      0.0305566 * (smok_dur - 27)
  }, 0)

  x <- x - fifelse(smok_status %in% (2:3), 0.0321362 * (quit_yrs - 8.593417626), 0)

  b        <- a_nocopd + x + 0.3457265 * copd
  b_nocopd <- a_nocopd + x
  b_copd   <- a_nocopd + x + 0.3457265 # COPD coef

  exp_b_nocopd <- exp(b_nocopd)
  exp_b_copd   <- exp(b_copd)
  exp_b        <- exp(b)
  prb_b_nocopd <- exp_b_nocopd / (1 + exp_b_nocopd)
  prb_b_copd   <- exp_b_copd / (1 + exp_b_copd)
  prb_b        <- exp_b / (1 + exp_b)

  # Then RR is:
  rr <- list(
    "rr_for_parf"  = clamp(prb_b/prb_a_nocopd, 1, 20), # otherwise light ex smokers had RR < 1
    "rr_no_copd"   = clamp(prb_b_nocopd/prb_a_nocopd, 1, 20),
    "rr_with_copd" = clamp(prb_b_copd/prb_a_nocopd, 1, 20)
  )

  return(rr)
}

# get observed mortality from ONS data
#' @export
get_ons_mrtl <- function(disease, type = c("rate", "absolute"), agegrp_width = c(5L, 20L)) {
  if (type == "rate") prefix <- "Mx_"
  if (type == "absolute") prefix <- "deaths_"

  colnams <- c("year", "sex", "qimd", paste0("agegrp", agegrp_width), paste0(prefix, disease), "pops")
  if (agegrp_width == 5) tt <- read_fst("./ONS_data/mortality_by_agegrp5_for_validation.fst", columns = colnams, as.data.table = TRUE)
  if (agegrp_width == 20) tt <- read_fst("./ONS_data/mortality_by_agegrp20_for_validation.fst", columns = colnams, as.data.table = TRUE)
  tt
}
# get_ons_mrtl("chd", "rate", 20)

# get observed cancer incidence from ONS data
#' @export
get_ons_incd <- function(disease, type = c("rate", "absolute"), agegrp_width = 20L) {
  if (type == "rate") prefix <- "rate_"
  if (type == "absolute") prefix <- "cases_"

  colnams <- c("year", "sex", "qimd", paste0("agegrp", agegrp_width), paste0(prefix, disease), "pops")
  if (agegrp_width == 5) tt <- read_fst("./ONS_data/cancer_incd_by_agegrp5_for_validation.fst", columns = colnams, as.data.table = TRUE)
  if (agegrp_width == 20) tt <- read_fst("./ONS_data/cancer_incd_by_agegrp20_for_validation.fst", columns = colnams, as.data.table = TRUE)
  tt
}

# get disease epi parameters for mc
#' @export
get_disease_epi_mc <-
  function(mc,
           disease,
           epi_par = c("incidence", "prevalence", "fatality", "duration"),
           what = c("value", "median_value", "p"), # p is the percentile of the uncertainty
           stochastic = TRUE) {
    # stochastic = F returns median_value even if what == "value"
    stopifnot(between(mc, 1, 2e3)) # TODO make 2e3 derived from disease_epi_indx.fst
    epi_par2 <- epi_par <- match.arg(epi_par)
    what <- match.arg(what)
    if (epi_par == "fatality") {epi_par <- paste0("case_", epi_par)}
    if (epi_par == "duration") {epi_par <- paste0(epi_par, "_years")} else {epi_par <- paste0(epi_par, "_rates")}
    if (!stochastic && what == "value") {what <- "median_value"}
    indx <- read_fst("./disease_epidemiology/disease_epi_indx.fst", from = mc, to = mc, as.data.table = TRUE)

    if (epi_par2 != "duration") {

      colnam <- c("age", "sex", "qimd", paste0(what, "_", disease, "_", epi_par))
      out <- read_fst("./disease_epidemiology/disease_epi_l.fst", columns = colnam,
                      from = indx$from, to = indx$to, as.data.table = TRUE)
      setnames(out, tail(colnam, 1), epi_par2)

    } else {# if duration

      dur_colnam <- paste0(what, "_", disease, "_", epi_par)
      prev_colnam <- paste0(what, "_", disease, "_", "prevalence_rates")
      colnam <- c("age", "sex", "qimd", dur_colnam, prev_colnam)

      out <- read_fst("./disease_epidemiology/disease_epi_l.fst", columns = colnam,
                      from = indx$from, to = indx$to, as.data.table = TRUE)
      age_range <- out[, range(age)]
      setnames(out, prev_colnam, "prevalence_rates")
      setnames(out, dur_colnam, "duration")

      out <- out[prevalence_rates > 0] # To exclude men in breast ca. May have

      # get the prevalent cases. Each row is now a case
      n <- out[, 10/min(prevalence_rates)]# logic to get at least 10 cases per strata and prevent 0
      out <- out[rep(seq_len(.N), prevalence_rates * n)]
      out[, `:=` (pid = .I - 1L, prevalence_rates = NULL)]
      # Now project the pids that each case will have copd
      out <- out[rep(seq_len(.N), duration)]

      out[, age := age + seq_len(.N) - 1L, by = pid]
      out[, disease_years := seq_len(.N), by = pid]
      out <- out[between(age, age_range[1], age_range[2]), .("duration" = mean(disease_years)), keyby = .(age, sex, qimd)]
    }

    invisible(out)
  }
# get_disease_epi_mc(1, "chd", "i", "v")

# get disease disease epi parameters for mc


#' @export
get_rr_mc <-
  function(mc,
           disease,
           risk_factor,
           stochastic = TRUE) {
    stopifnot(between(mc, 1, 2e3)) # TODO make 2e3 derived from filenam_indx
    suffix <- paste0(risk_factor, "_", disease)
    filenam <- paste0("./simulation/rr/", suffix, "_rr_l.fst")
    filenam_indx <- paste0("./simulation/rr/", suffix, "_rr_indx.fst")
    if (file.exists(filenam_indx)) {
      indx <- read_fst(filenam_indx, from = mc, to = mc, as.data.table = TRUE)
      out <- read_fst(filenam, from = indx$from, to = indx$to, as.data.table = TRUE)
    } else {
      out <- read_fst(filenam, from = mc, to = mc, as.data.table = TRUE)
    }
    out[, mc := NULL]
    colnam_rr <- paste0("rr_", suffix)
    colnam_mean_rr <- paste0("mean_rr_", suffix)
    if (stochastic) {
      out[, (colnam_mean_rr) := NULL]
      setnames(out, colnam_rr, paste0(risk_factor, "_rr"))
    } else {
      out[, (colnam_rr) := NULL]
      setnames(out, colnam_mean_rr, paste0(risk_factor, "_rr"))
    }
    if (identical(dim(out), c(1L, 1L))) out <- as.numeric(out)
    out[]
  }
# get_rr_mc(1, "chd", "sbp", TRUE)

get_lag_mc_hlp <-
  function(mc,
           disease_enum, # 1:10
           lag) { # 2:9
    stopifnot(between(disease_enum, 1, 10), between(lag, 2, 9), between(mc, 1, 2e3))
    colnam <- paste0("lag_", lag)
    filenam <- "./disease_epidemiology/disease_lags_l.fst"
    filenam_indx <- "./disease_epidemiology/disease_lags_indx.fst"
    indx <- read_fst(filenam_indx, from = mc, to = mc, as.data.table = TRUE)
    out  <- read_fst(filenam, colnam, from = indx$from, to = indx$to,
                     as.data.table = TRUE)[disease_enum, get(colnam)]
    out
  }

#' @export
get_lag_mc <- function(mc, design) {
  # arguement checks in the helper function get_lag_mc_hlp
  nam <- grep("_lag$", names(design), value = TRUE)

  if (design$stochastic) {
    out <- vector("list", length(nam))
    names(out) <- nam
    invisible(lapply(nam, function(x) {
      out[x] <<- get_lag_mc_hlp(mc, design[[paste0(x, "_enum")]], design[[x]])
    }))
  } else {
    out <- design[nam]
  }
  out$plco_lag <- 6L # for lung ca plco formula
  out
}


#' @export
get_lifetable_all <-
  function(mc, disease, design, type = c("qx", "mx")) {
    if (disease %in% c("allcause", "nonmodelled")) {
      disease2 <- "chd"
    } else {
      disease2 <- disease
    }
    prb <- get_disease_epi_mc(mc, disease = disease2, "fatality", "p", stochastic = design$stochastic)
    # NOTE fatality column is poorely named. It is a probability which is correlated with incidence & prevalence
    colnam <- c("year", "age", "sex", "qimd", "disease",
                paste0(type, c("_total", "_total_1", "_total_99", "_total_10",
                               "_total_20", "_total_30", "_total_40", "_total_60",
                               "_total_70", "_total_80", "_total_90")))
    disease_ <- disease
    indx <- read_fst( "./lifecourse_models/mortality_projections_indx.fst", as.data.table = TRUE)[disease == disease_]

    lifetable_all <- read_fst("./lifecourse_models/mortality_projections.fst", colnam, from = indx$from, to = indx$to, as.data.table = TRUE)[
      between(year, design$init_year, design$init_year + design$sim_horizon_max) & between(age, design$ageL, design$ageH)
    ]
    if (design$stochastic) {
      absorb_dt(lifetable_all, prb)
      setnames(lifetable_all, colnam, gsub(paste0("^", type, "_"), "", colnam))

      lifetable_all[fatality < 0.1, qx_mc :=
                      qunif(normalise(c(fatality, 0, 0.1))[1],
                            total_1, total_10)]
      lifetable_all[between(fatality, 0.1, 0.2), qx_mc :=
                      qunif(normalise(c(fatality, 0.1, 0.2))[1],
                            total_10, total_20)]
      lifetable_all[between(fatality, 0.2, 0.3), qx_mc :=
                      qunif(normalise(c(fatality, 0.2, 0.3))[1],
                            total_20, total_30)]
      lifetable_all[between(fatality, 0.3, 0.4), qx_mc :=
                      qunif(normalise(c(fatality, 0.3, 0.4))[1],
                            total_30, total_40)]
      lifetable_all[between(fatality, 0.4, 0.5), qx_mc :=
                      qunif(normalise(c(fatality, 0.4, 0.5))[1],
                            total_40, total)]
      lifetable_all[between(fatality, 0.5, 0.6), qx_mc :=
                      qunif(normalise(c(fatality, 0.5, 0.6))[1],
                            total, total_60)]
      lifetable_all[between(fatality, 0.6, 0.7), qx_mc :=
                      qunif(normalise(c(fatality, 0.6, 0.7))[1],
                            total_60, total_70)]
      lifetable_all[between(fatality, 0.7, 0.8), qx_mc :=
                      qunif(normalise(c(fatality, 0.7, 0.8))[1],
                            total_70, total_80)]
      lifetable_all[between(fatality, 0.8, 0.9), qx_mc :=
                      qunif(normalise(c(fatality, 0.8, 0.9))[1],
                            total_80, total_90)]
      lifetable_all[between(fatality, 0.9, 1), qx_mc :=
                      qunif(normalise(c(fatality, 0.9, 1))[1],
                            total_90, total_99)]
      lifetable_all[fatality == 0.5, qx_mc := total]
      return(lifetable_all[, c("year", "age", "sex", "qimd", "qx_mc"), with = FALSE])
    } else {
      setnames(lifetable_all, paste0(type, "_total"), "qx_mc")
      return(lifetable_all[, c("year", "age", "sex", "qimd", "qx_mc"), with = FALSE])
    }
  }

# get_lifetable_all(4, "chd", design, "qx")[age == 60 & qimd == "3" & sex == "men", plot(year, qx_mc)]


# Get population estimates adjusted for mortality
#' @export
generate_pop_adj_for_mrtl <- function(mc, dt, design, update_dt = FALSE) {
  tt <- get_lifetable_all(mc = mc, "allcause", design = design, "qx")
  orig_pops <- pops <- dt[, .(N = as.numeric(.N)), keyby = .(year, age, sex, qimd)]
  absorb_dt(pops, tt)

  out <- data.table()

  ttt <- pops[year == design$init_year & between(age, design$ageL, design$ageH)][, N_adj := as.numeric(N)]
  out <- rbind(out, ttt)
  for (i in seq_along(unique(pops$year))) {
    ttt1 <- copy(ttt)[, `:=`(year = year + 1L, age = age + 1L)] # + 1 here correct, not + i.
    ttt <- pops[year == design$init_year + i & between(age, design$ageL, design$ageH)]
    ttt[ttt1, on = c("age", "sex", "qimd"), N := i.N_adj]
    ttt[, N_adj := N * (1 - qx_mc)]
    out <- rbind(out, ttt)
  }
  out[orig_pops, N := i.N, on = c("year", "age", "sex", "qimd")]
  if (update_dt) {
    dt[out, on = c("year", "age", "sex", "qimd"), pops_adj := N_adj]
    message("dt was updated.")
  }
  out[, qx_mc := NULL]
  out[]
}
# generate_pop_adj_for_mrtl(mc, dt, design, TRUE)
#
# dt[is.na(pops_adj) & between(year, design$init_year, design$init_year + design$sim_horizon) & between(age, design$ageL, design$ageH), .N]

#' @export
simulate_fatality <-
  function(mc, dt, design, disease, prvl_constr = 2) {
    # To calculate fatality based on mortality projections I need to have an estimate of
    # disease prevalence by year. Then fatality rate = qx rate / prev rate
    # To estimate prevalence I need to run a macrosimulation. This is what I do below
    tt <- generate_pop_adj_for_mrtl(mc = mc, dt, design)
    absorb_dt(tt, get_lifetable_all(mc = mc, disease, design, "qx")) # disease mortality
    absorb_dt(tt,
              get_disease_epi_mc(mc = mc, disease, "p", "v", design$stochastic)[, year := design$init_year])
    colnam <- paste0("prb_", disease, "_inc")
    setnames(dt, colnam, "V1___")
    ttt <-
      dt[, .(incd = sum(V1___) / .N), # No need for N_adj for mortality
         keyby = .(year, age, sex, qimd)]
    setnames(dt, "V1___", colnam)
    absorb_dt(tt, ttt)

    # Estimate N_adj as pop adjustment for all non-disease
    # mortl rather than all mrtl
    tt[, N_adj := N_adj + N * qx_mc]

    # With the method below very low or even -ve prevalence is produced.
    # I will constrain prevalence to 1/2 of the observed in init_year.

    ttt <- tt[year == design$init_year]
    # With the method below very low or even -ve prevalence is produced.
    # I will constrain prevalence between 1/2 and *2 of the observed in init_year. If this
    # limit is reached then I will adjust qx_mc
    constrain <-
      ttt[, .(
        age,
        sex,
        qimd,
        prvl_limL = prevalence * 1 / prvl_constr,
        prvl_limU = clamp(prevalence * prvl_constr)
      )]

    # this method does not work for age == 30. Hence I will calculate the
    # ratio of prevalence between age 30 and 31 for init_year and I will
    # apply it in all concecutive years
    ttt30 <-
      ttt[age == design$ageL, .(age, sex, qimd, prvl30 = prevalence)]
    ttt31 <-
      ttt[age == design$ageL + 1L, .(sex, qimd, prvl31 = prevalence)]
    absorb_dt(ttt30, ttt31)
    ttt30[, `:=` (
      prvl30_ratio = prvl30 / prvl31,
      age = NULL,
      prvl30 = NULL,
      prvl31 = NULL
    )]

    for (i in (seq_along(unique(tt$year)) - 1L)) {
      ttt <- tt[year == design$init_year + i]
      ttt[, cases := N_adj * (prevalence + incd - qx_mc)]
      ttt[, `:=`(year = year + 1L, age = age + 1L)] # move to next year
      ttt31 <-
        ttt[age == design$ageL + 1L, .(sex, qimd, prevalence)]
      ttt31[ttt30, on = c("sex", "qimd"), `:=` (prevalence = prevalence * i.prvl30_ratio,
                                                age = design$ageL)]
      ttt[age == design$ageH + 1L, age := design$ageL]
      suppressMessages(absorb_dt(ttt, ttt31, on = c("age", "sex", "qimd")))
      tt[ttt, on = c("year", "age", "sex", "qimd"), prevalence := cases / N_adj]
      # check constrain
      ttt <-
        tt[year == design$init_year + i + 1L, .(year, age, sex, qimd, prevalence)]
      if (ttt[constrain, on = .NATURAL,][prevalence < prvl_limL, .N] > 0) {
        ttt[constrain, on = .NATURAL, prevalence := fifelse(prevalence < i.prvl_limL,
                                                            i.prvl_limL, prevalence)]
        suppressMessages(absorb_dt(tt, ttt, c("year", "age", "sex", "qimd")))
      }
      if (ttt[constrain, on = .NATURAL,][prevalence > prvl_limU, .N] > 0) {
        ttt[constrain, on = .NATURAL, prevalence := fifelse(prevalence > i.prvl_limU,
                                                            i.prvl_limU, prevalence)]
        suppressMessages(absorb_dt(tt, ttt, c("year", "age", "sex", "qimd")))
      }
    }
    tt[, fatality := clamp(qx_mc / clamp(prevalence))]
    tt[, .(year, age, sex, qimd, fatality)]
  }

#' @export
generate_rns <- function(mc, dt, colnams) {
  dqRNGkind("pcg64") # dqRNGkind("Xoroshiro128+") ~10% faster
  SEED <- 4719349L # sample(1e7, 1)
  set.seed(SEED + mc)
  dqset.seed(SEED, mc)

  nrows <- nrow(dt)
  for (nam in colnams)
    set(dt, NULL, nam, dqrunif(nrows))
  invisible(dt)
}


# Given a correlation matrix (Pearson), produces a matrix of
# correlated uniforms

#' @export
generate_corr_unifs <- function(n, M) {
  # generate normals, check correlations
  # from http://comisef.wikidot.com/tutorial:correlateduniformvariates
  stopifnot(is.matrix(M))
  # Check that matrix is semi-positive definite
  stopifnot(min(eigen(M)$values) >= 0)

  M_original <- M


  X <- matrix(dqrnorm(n * dim(M)[[2]]), n)
  colnames(X) <- colnames(M)

  # adjust correlations for uniforms
  for (i in 1:dim(M)[[2]]){
    for (j in 1:dim(M)[[2]]){
      if (i != j){
        M[i, j] <- 2 * sin(pi * M[i, j] / 6)
        M[j, i] <- 2 * sin(pi * M[j, i] / 6)
      }
    }
  }

  # induce correlation, check correlations
  C <- chol(M)
  Y <- X %*% C
  cor(Y)

  Y <- pnorm(Y)
  message(paste0("Mean square error is: ", signif(sum((cor(Y) - M_original) ^ 2), 3)))
  return(Y)
}


# Estimate health utility decreaments
#' @export
generate_eq5d_decr <- function(dt) {
  # From Sullivan et al. 2011
  # TODO add copd and cancers
  utility_pop_norms <-
    c(0.922, 0.922, 0.922, 0.922, 0.922, 0.922, 0.922, 0.914, 0.914, 0.914, 0.914, 0.914,
      0.914, 0.914, 0.914, 0.914, 0.914, 0.888, 0.888, 0.888, 0.888, 0.888, 0.888, 0.888,
      0.888, 0.888, 0.888, 0.854, 0.854, 0.854, 0.854, 0.854, 0.854, 0.854, 0.854, 0.854,
      0.854, 0.814, 0.814, 0.814, 0.814, 0.814, 0.814, 0.814, 0.814, 0.814, 0.814, 0.775,
      0.775, 0.775, 0.775, 0.775, 0.775, 0.775, 0.775, 0.775, 0.775, 0.706, 0.706, 0.706,
      0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706,
      0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706, 0.706)
  utility_income    <- c(0.012572, 0.038783, 0.0396568, 0.0396568, 0.0408501)
  utility_education <- c(0.0060444, 0.0056836, 0.0028418, 0.0028418, 0.0028418, 0.0056836, 0)
  utility_ncc       <- c(0, 0, -0.0528, -0.0415, -0.0203, 0.0083, 0.04087, 0.06687, 0.11589, 0.13444, 0.18361)

  new_ncc <- clamp_int(dt$ncc + (dt$htn_prvl > 0L) + (dt$af_prvl > 0L) +
                         (dt$t2dm_prvl > 0L) +  (dt$chd_prvl > 0L) +
                         (dt$stroke_prvl > 0L) + (dt$poststroke_dementia_prvl > 0L) +
                         (dt$copd_prvl > 0L) + (dt$lung_ca_prvl > 0L) +
                         (dt$colon_ca_prvl > 0L) + (dt$breast_ca_prvl > 0L),
                       0L, 10L)


  out <- utility_pop_norms[dt$age - 17L] + utility_income[dt$income] + utility_education[dt$education] +
    0.0010046 * (dt$sex == "men") +
    utility_ncc[new_ncc + 1L] - 0.0460 * (dt$htn_prvl > 0L) -
    0.0384 * (dt$af_prvl > 0L) - 0.0714 * (dt$t2dm_prvl > 0L) -
    0.0679 * (dt$chd_prvl > 0L) - 0.0578 * (dt$stroke_prvl > 0L) -
    0.0957 * (dt$copd_prvl > 0L) - 0.1192427 * (dt$lung_ca_prvl > 0L) -
    0.0673908 * (dt$colon_ca_prvl > 0L) - 0.0194279 * (dt$breast_ca_prvl > 0L) -
    0.2165659 * (dt$poststroke_dementia_prvl > 0L)

  # half eq5d for year of death
  out <- out/((dt$all_cause_mrtl > 0) + 1L) # out/1 for alive and out/2 for deads

  clamp(out, 0, 1, TRUE)
  dt[, eq5d := out]
}
# generate_eq5d_decr(output)
# output[, summary(eq5d)]


#' @export
generate_healthcare_costs <- function(dt) {
  # cancers cost by year since diagnosis age < 65 and >=65
  # (first element is for 0)
  lung_ca_cost_young <- c(0, 16356.18,	5276.19,	4650.95,	3104.12,	2964.66,
                          2964.66,	2964.66,	2964.66,	2964.66,	2964.66)
  lung_ca_cost_old <- c(0, 14092.30,	5020.51,	4584.70,	3910.65,	3536.44,
                        3536.44,	3536.44,	3536.44,	3536.44,	3536.44)
  colon_ca_cost_young <- c(0, 21763.69, 5827.05,	4284.87,	3401.63,	2775.23,
                           2118.61,	2277.83,	1961.72,	1592.15,	1592.15)
  colon_ca_cost_old <- c(0, 20270.32,	4917.08,	3954.82,	3278.44,	3218.01,
                         3185.47,	2720.61,	3056.47,	2598.58,	2598.58)
  breast_ca_cost_young <- c(0, 13877.30,	4272.08,	2528.85,	2070.96,	1984.96,
                            1912.91,	1695.59,	1664.21,	1529.40,	1529.40)
  breast_ca_cost_old <- c(0, 11332.18,	3108.77,	2638.09,	2653.20,	2540.47,
                          2582.31,	2464.93,	2491.66,	2646.23,	2646.23)

  out <- 1237 + 72 * (dt$htn_prvl > 0L) + 1102 * (dt$af_prvl > 0L) +
    586 * (dt$t2dm_prvl > 0L) + 1667 * (dt$chd_prvl > 0L) +
    (10603 + 9040) * (dt$stroke_prvl == 1L) +
    5569 * (dt$stroke_prvl > 1L) + 2404 * (dt$copd_prvl > 0L) +
    lung_ca_cost_young[1L + dt$lung_ca_prvl] * (dt$age < 65) +
    lung_ca_cost_old[1L + dt$lung_ca_prvl] * (dt$age >= 65) +
    colon_ca_cost_young[1L + dt$colon_ca_prvl] * (dt$age < 65) +
    colon_ca_cost_old[1L + dt$colon_ca_prvl] * (dt$age >= 65) +
    breast_ca_cost_young[1L + dt$breast_ca_prvl] * (dt$age < 65) +
    breast_ca_cost_old[1L + dt$breast_ca_prvl] * (dt$age >= 65) +
    2289 * (dt$poststroke_dementia_prvl > 0L)

  # half cost for year of death
  out <-
    out / ((dt$all_cause_mrtl > 0) + 1L) # out/1 for alive and out/2 for dead
  dt[, healthcare_cost := out]
}

#' @export
generate_productivity_costs <- function(dt) {
  tt <- fread("./simulation/health_econ/costs_productivity_use.csv",
              stringsAsFactors = TRUE)
  tt[, eq5d_r := as.integer(100L * eq5d_r)]
  dt[, eq5d_r := as.integer(100L * round(eq5d/0.05)*0.05)]
  absorb_dt(dt, tt)
  # half cost for year of death
  dt[all_cause_mrtl > 0, productivity_cost := productivity_cost/2]
  dt[, eq5d_r := NULL]
  dt
}

#' @export
generate_socialcare_costs <- function(dt) {
  socialcare_costs <- c(34, 67, 101, 134, 168, 168, 168, 168, 168, 168,
                        168, 168, 168, 168, 168, 168, 168, 168, 168, 168,
                        168, 168, 168, 168, 168, 168, 168, 168, 199, 229,
                        260, 291, 322, 322, 322, 322, 322, 322, 322, 322,
                        322, 322, 322, 322, 322, 322, 322, 307, 292, 277,
                        263, 248, 287, 327, 366, 405, 445, 500, 556, 611,
                        667, 722, 909, 1096,1282, 1469, 1656, 2115, 2574,
                        3034, 3493, 3953, 3953, 3953, 3953, 3953, 3953,
                        3953, 3953, 3953, 3953, 3953, 3953)
  # for ages 18 to 100
  out <- socialcare_costs[dt$age - 17L] + 699 * (dt$stroke_prvl > 0L) +
    4862 * (dt$poststroke_dementia_prvl > 0L)

  # half cost for year of death
  out <-
    out / ((dt$all_cause_mrtl > 0) + 1L) # out/1 for alive and out/2 for deads
  dt[, socialcare_cost := out]
  dt
}

#' @export
generate_informal_care_costs <- function(dt) {
  tt <- fread("./simulation/health_econ/costs_informal_care.csv",
              stringsAsFactors = TRUE)
  tt[, eq5d_r := as.integer(100L * eq5d_r)]
  dt[, eq5d_r := as.integer(100*round(eq5d/0.05)*0.05)]
  absorb_dt(dt, tt)
  # half cost for year of death
  dt[all_cause_mrtl > 0, informal_care_cost := informal_care_cost/2]
  dt[, eq5d_r := NULL]
  dt
}

#' @export
generate_health_econ <- function(dt) {
  generate_eq5d_decr(dt)
  generate_healthcare_costs(dt)
  generate_socialcare_costs(dt)
  generate_productivity_costs(dt)
  generate_informal_care_costs(dt)
  invisible(dt)
}

#' @export
set_eligible <- function(scenario_nam, dt, parameters_dt) {
  l <- fromGUI_scenario_parms(scenario_nam, parameters_dt)
  colnam <- "eligible_sc"
  set(dt, NULL, colnam, 0L)
  if (!l$sc_eligib_noone) {
    dt[between(year + 2000L, l$sc_init_year, l$sc_last_year) &
         between(age, l$sc_eligib_age[[1]], l$sc_eligib_age[[2]]) &
         dead == FALSE &
         # statin_px_curr_xps == 0L &
         # af_dgn == 0L &
         htn_dgn < fifelse(l$sc_eligib_htn, Inf, 1) &
         t2dm_dgn < fifelse(l$sc_eligib_diab, Inf, 1L) &
         # ckd_prvl_curr_xps <= 3L &
         # ra_prvl == 0L &
         chd_dgn == 0L &
         stroke_dgn == 0L,
       (colnam) := 1L]
  }
  invisible(dt)
  # dt[year == 20 &
  #      between(age, l$sc_eligib_age[[1]], l$sc_eligib_age[[2]]), prop_if(eligible_sc1 == 1)]
}
# set_eligible("sc2", POP, parameters_dt)
# POP[between(year, 18, 35) & between(age, 40, 74) & dead == FALSE, prop_if(eligible_sc1 == 1)]
# Around 70% of the population should be eligible. From https://www.healthcheck.nhs.uk/commissioners-and-providers/data/total-eligible-population/
# POP[between(year, 18, 35) & between(age, 40, 74) & dead == FALSE, 1 - prop_if(eligible_sc1 == 1), keyby = age]

#' @export
set_invitees <- function(scenario_nam, dt, parameters_dt) {
  l <- fromGUI_scenario_parms(scenario_nam, parameters_dt)
  colnam <- "invitees_sc"
  colnam_cost <- "invitation_cost_sc"
  elig_colnam <- "eligible_sc"
  set(dt, NULL, colnam, 0L)

  # TODO find better approach

  if (l$sc_invit_detailed) {
    tt1 <- data.table(
      year = (l$sc_init_year:l$sc_last_year) - 2000L,
      qimd = "1 most deprived",
      mu = l$sc_invit_qimd1,
      freq = l$sc_eligib_freq
    )
    elig <- vector("numeric", nrow(tt1))
    for (i in seq_len(nrow(tt1))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt1[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt1[i, mu])
      if (i > tt1[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt1[i - 1, mu] + tt1[i - tt1[i, freq], mu])
    }
    tt1[, mu := clamp(mu / elig)]
    tt1[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    tt2 <-
      data.table(
        year = (l$sc_init_year:l$sc_last_year) - 2000L,
        qimd = "2",
        mu = l$sc_invit_qimd2,
        freq = l$sc_eligib_freq
      )
    elig <- vector("numeric", nrow(tt2))
    for (i in seq_len(nrow(tt2))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt2[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt2[i, mu])
      if (i > tt2[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt2[i - 1, mu] + tt2[i - tt2[i, freq], mu])
    }
    tt2[, mu := clamp(mu / elig)]
    tt2[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    tt3 <-
      data.table(
        year = (l$sc_init_year:l$sc_last_year) - 2000L,
        qimd = "3",
        mu = l$sc_invit_qimd3,
        freq = l$sc_eligib_freq
      )
    elig <- vector("numeric", nrow(tt3))
    for (i in seq_len(nrow(tt3))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt3[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt3[i, mu])
      if (i > tt3[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt3[i - 1, mu] + tt3[i - tt3[i, freq], mu])
    }
    tt3[, mu := clamp(mu / elig)]
    tt3[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    tt4 <-
      data.table(
        year = (l$sc_init_year:l$sc_last_year) - 2000L,
        qimd = "4",
        mu = l$sc_invit_qimd4,
        freq = l$sc_eligib_freq
      )
    elig <- vector("numeric", nrow(tt4))
    for (i in seq_len(nrow(tt4))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt4[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt4[i, mu])
      if (i > tt4[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt4[i - 1, mu] + tt4[i - tt4[i, freq], mu])
    }
    tt4[, mu := clamp(mu / elig)]
    tt4[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    tt5 <-
      data.table(
        year = (l$sc_init_year:l$sc_last_year) - 2000L,
        qimd = "5 least deprived",
        mu = l$sc_invit_qimd5,
        freq = l$sc_eligib_freq
      )
    elig <- vector("numeric", nrow(tt5))
    for (i in seq_len(nrow(tt5))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt5[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt5[i, mu])
      if (i > tt5[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt5[i - 1, mu] + tt5[i - tt5[i, freq], mu])
    }
    tt5[, mu := clamp(mu / elig)]
    tt5[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    tt <- rbind(tt1, tt2, tt3, tt4, tt5)

    ttcost <-
      data.table(
        V1 = 1L,
        V2 = c("1 most deprived", "2", "3", "4", "5 least deprived"),
        V3 = c(
          l$sc_invit_qimd1_cost,
          l$sc_invit_qimd2_cost,
          l$sc_invit_qimd3_cost,
          l$sc_invit_qimd4_cost,
          l$sc_invit_qimd5_cost
        )
      )
    setnames(ttcost, c(colnam, "qimd", colnam_cost))

  } else {
    # if input not by qimd
    # The method below is robust in changes of mu and/or freq to
    # define probability of being invited given that eligible
    # population is reduced the more you invite for a check
    tt <- data.table(
      year = (l$sc_init_year:l$sc_last_year) - 2000L,
      mu = l$sc_invit_qimdall,
      freq = l$sc_eligib_freq
    )
    elig <- vector("numeric", nrow(tt))
    for (i in seq_len(nrow(tt))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt[i, mu])
      if (i > tt[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt[i - 1, mu] + tt[i - tt[i, freq], mu])
    }
    tt[, mu := clamp(mu / elig)]
    tt[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    ttcost <- data.table(1L, l$sc_invit_qimdall_cost)
    setnames(ttcost, c(colnam, colnam_cost))
  }

  absorb_dt(dt, tt)
  setnafill(dt, "c", 0, cols = c("mu", "freq"))
  dt[, (colnam) := identify_invitees(eligible_sc, mu, freq, pid_mrk)]
  dt[, c("mu", "freq") := NULL]
  absorb_dt(dt, ttcost)
  setnafill(dt, "c", 0, cols = colnam_cost)
  invisible(dt)
  # dt[ eligible_sc1 == 1L &
  #      between(age, l$sc_eligib_age[[1]], l$sc_eligib_age[[2]]), prop_if(invitees_sc1== 1), keyby = year]

}
# set_invitees("sc2", POP, parameters_dt)

#' @export
set_attendees <- function(scenario_nam, dt, parameters_dt) {
  l <- fromGUI_scenario_parms(scenario_nam, parameters_dt)
  colnam       <- "attendees_sc"
  colnam_cost  <- "attendees_cost_sc"
  invit_colnam <- "invitees_sc"
  set(dt, NULL, colnam_cost, 0L)


  if (l$sc_uptake_detailed) {
    agegrp <- fromGUI_uptake_table_agegrps(scenario_nam = scenario_nam,
                                           parameters_dt = parameters_dt)
    absorb_dt(dt, agegrp)
    dt[, "Qrisk2_cat" := Qrisk2(.SD, l$sc_qrisk_ignore_bmi, l$sc_qrisk_ignore_sbp,
                                l$sc_qrisk_ignore_tchol)$Qrisk2_cat]
    if (l$sc_uptake_structural0s) {
      absorb_dt(dt, l$sc_uptake)
      setnafill(dt, "c", 0, cols = "uptake_wt")
      dt[, c("Qrisk2_cat", "agegrp10") := NULL]
      tt <- sort(dt[invitees_sc == 1L & uptake_wt > 0,
                    sample_int_expj(.N, as.integer(round(l$sc_uptake_all * .N)),
                                    uptake_wt)])
      tt <-
        dt[invitees_sc == 1L &
             uptake_wt > 0, .(year, pid)][tt,][, (colnam) := 1L]
      absorb_dt(dt, tt, on = c("pid","year"))
    } else { # if no sructural 0s
      absorb_dt(dt, l$sc_uptake[uptake_wt == 0, uptake_wt := 1e-6])
      setnafill(dt, "c", 0, cols = "uptake_wt")
      dt[, c("Qrisk2_cat", "agegrp10") := NULL]

      tt <- sort(dt[invitees_sc == 1L,
                    sample_int_expj(.N, as.integer(round(l$sc_uptake_all * .N)),
                                    uptake_wt)])
      tt <-
        dt[invitees_sc == 1L, .(year, pid)][tt,][, (colnam) := 1L]
      absorb_dt(dt, tt, on = c("pid","year"))
    }
  } else {
    set(dt, NULL, "uptake_wt", l$sc_uptake_all)
    dt[invitees_sc == 1L, (colnam) := rbinom(.N, 1, uptake_wt)]
  }

  setnafill(dt, "c", 0, cols = colnam)
  dt[attendees_sc == 1L, (colnam_cost) := l$sc_uptake_all_cost]
  dt[, uptake_wt := NULL]
  invisible(dt)

  # dt[ invitees_sc1 == 1L &
  #      between(age, l$sc_eligib_age[[1]], l$sc_eligib_age[[2]]), prop_if(attendees_sc1== 1), keyby = year]
}

#' @export
set_px <- function(scenario_nam, dt, parameters_dt) {
  l <- fromGUI_scenario_parms(scenario_nam, parameters_dt)
  dt[, "Qrisk2_cat" := Qrisk2(.SD, FALSE,  FALSE, FALSE)$Qrisk2_cat]
  atte_colnam <- "attendees_sc"

  # for statins
  colnam     <- "statin_px_sc"
  colnam_bio <- "tchol_sc"

  if (l$sc_px_detailed) {
    absorb_dt(dt, l$sc_px_statins_wt)

    # Below assumes people on statin_px_curr_xps but undertreated will titrate
    # statin treatment.

    # Adjusted prb where the denominator changed from all attendees to those eligible for statins
    tt <- dt[attendees_sc == 1L,
             clamp(l$sc_px_statins/prop_if(tchol_curr_xps >= 5))]


    tt <- sort(dt[attendees_sc == 1L & tchol_curr_xps >= 5,
                  sample_int_expj(.N, as.integer(round(tt * .N)),
                                  px_statins_wt)]) # rows that will have statin px
    tt <-
      dt[attendees_sc == 1L & tchol_curr_xps >= 5,
         .(year, pid)][tt,][, (colnam) := 1L]
    absorb_dt(dt, tt, on = c("pid","year"))
  } else {
    tt <- dt[attendees_sc == 1L,
             clamp(l$sc_px_statins/prop_if(tchol_curr_xps >= 5 & Qrisk2_cat != "low"))]

    set(dt, NULL, "px_statins_wt", tt)
    dt[attendees_sc == 1L & tchol_curr_xps >= 5 & Qrisk2_cat != "low",
       (colnam) := rbinom(.N, 1L, px_statins_wt)]
  }
  setnafill(dt, "c", 0, cols = colnam)
  dt[, (colnam) := hc_effect(statin_px_sc, 0.9749866, pid_mrk)]
  # 0.9749866 comes from the following study in Wales.
  # King W, Lacey A, White J, Farewell D, Dunstan F, Fone D. Socioeconomic
  # inequality in medication persistence in primary and secondary prevention of
  # coronary heart disease – A population-wide electronic cohort study. PLOS ONE
  # 2018;13:e0194081.
  # from 33228 individuals px a statin in Wales for primary prevention, 5378 had
  # discontinued it within 7 years without socioeconomic gradient statin <-
  # c(0.96, (33228 - 5378)/33228) year <- c(0.1, 7) m1 <- glm(statin~ -1 + year,
  # family = gaussian(link = "log")) exp(m1$coefficients) x <- predict(m1,
  # newdata = data.table(year = 0:7), type = "re") shift(x, -1)/x plot(0:70,
  # predict(m1, newdata = data.table(year = 0:70), type = "re"), ylim = c(0, 1))
  # Every year 0.9749866 of those taking statin continue next year


  # estimate tchol change
  # atorvastatin effect
  # from Law MR, et al. Quantifying effect of statins on low density lipoprotein cholesterol,
  # ischaemic heart disease, and stroke: systematic review and meta-analysis. BMJ 2003;326:1423.
  # table 2. 43% (0.3958 - 0.46875) reduction of ldl. to convert to tc, tc/ldl = 0.27/0.36 from
  # Edwards JE, et al. Statins in hypercholesterolaemia: A dose-specific meta-analysis of lipid
  # changes in randomised, double blind trials. BMC Family Practice 2003;4:18.

  atorv_eff <- 0.27*0.43/0.36 # TODO need to be in the MC parameters
  # atorv_eff <- rnorm(
  #       0.27*0.43/0.36,
  #       (0.27*0.46875/0.36 - 0.27*0.3958/0.36)/(2*1.96))
  # adherence <- rpert(1e6, 0.5, 0.8, 1, 8) # proportion of prescribed dose taken
  # or to avoid dependency for rpert
  # adherence <- rBE(1e6, 0.8, 0.1) # proportion of prescribed dose taken

  dt[, (colnam_bio) := tchol_curr_xps]
  dt[statin_px_sc == 1L & statin_px_curr_xps == 0L, (colnam_bio) := tchol_curr_xps * (1 - atorv_eff * rBE(.N, 0.8, 0.2))]

  # for bpmed
  colnam     <- "bpmed_px_sc"
  colnam_bio <- "sbp_sc"

  if (l$sc_px_detailed) {
    absorb_dt(dt, l$sc_px_antihtn_wt)

    # Adjusted prb where the denominator changed from all attendees to those eligible for bpmed
    tt <- dt[attendees_sc == 1L,
             clamp(l$sc_px_antihtn/prop_if(sbp_curr_xps >= 135))]

    tt <- sort(dt[attendees_sc == 1L & sbp_curr_xps >= 135,
                  sample_int_expj(.N, as.integer(round(tt * .N)),
                                  px_antihtn_wt)])
    tt <-
      dt[attendees_sc == 1L & sbp_curr_xps >= 135,
         .(year, pid)][tt,][, (colnam) := 1L]
    absorb_dt(dt, tt, on = c("pid","year"))
  } else {
    tt <- dt[attendees_sc == 1L,
             clamp(l$sc_px_antihtn/prop_if(sbp_curr_xps >= 135))]

    set(dt, NULL, "px_antihtn_wt", tt)
    dt[attendees_sc == 1L & sbp_curr_xps >= 135,
       (colnam) := rbinom(.N, 1L, px_antihtn_wt)]
  }
  setnafill(dt, "c", 0, cols = colnam)
  dt[, (colnam) := hc_effect(bpmed_px_sc, 0.9749866, pid_mrk)] # assume same prb as statins

  # Estimate sbp change
  dt[, (colnam_bio) := sbp_curr_xps]
  dt[bpmed_px_sc == 1L & bpmed_curr_xps == 0L, (colnam_bio) := sbp_curr_xps - clamp((sbp_curr_xps - 135) * rBE(.N, 0.8, 0.2), 0, 1e3)] # Assume that antihtn medication can potentially achieve sbp 135 for all. Not 110 to account for residual risk

  dt[, c("Qrisk2_cat", "px_statins_wt", "px_antihtn_wt") := NULL]

  dt
}

#' @export
set_lifestyle <- function(scenario_nam, dt, parameters_dt, design = design) {
  l <- fromGUI_scenario_parms(scenario_nam, parameters_dt)
  atte_colnam <- "attendees_sc"

  # PA
  set(dt, NULL, "hc_eff", 0L)
  colnam      <- "active_days_sc"
  colnam_cost <- "active_days_cost_sc"
  dt[, (colnam) := active_days_curr_xps]
  set(dt, NULL, colnam_cost, 0)
  dt[attendees_sc == 1L, hc_eff := rbinom(.N, 1, l$sc_ls_papct)]
  dt[hc_eff == 1L, (colnam_cost) := l$sc_ls_pa_cost_ind] # Cost only the year of referral
  dt[, hc_eff := hc_effect(hc_eff, 0.8, pid_mrk)] # TODO 0.8 to advanced settings
  dt[hc_eff == 1L, (colnam) := clamp(active_days_sc + l$sc_ls_papincr, 0, 7)]

  # Weight management
  set(dt, NULL, "hc_eff", 0L)
  colnam      <- "bmi_sc"
  colnam_cost <- "bmi_cost_sc"
  dt[, (colnam) := bmi_curr_xps]
  set(dt, NULL, colnam_cost, 0)
  dt[attendees_sc == 1L & bmi_curr_xps > 30, hc_eff := rbinom(.N, 1, l$sc_ls_wghtpct)]
  dt[hc_eff == 1L, (colnam_cost) := l$sc_ls_wghtloss_cost_ind] # Cost only the year of referral
  dt[, hc_eff := hc_effect(hc_eff, 0.8, pid_mrk)] # TODO 0.8 to advanced settings
  dt[hc_eff == 1L, (colnam) := bmi_sc * (1 - l$sc_ls_wghtreduc)]

  # Alcohol
  set(dt, NULL, "hc_eff", 0L)
  colnam      <- "alcohol_sc"
  colnam_cost <- "alcohol_cost_sc"
  dt[, (colnam) := alcohol_curr_xps]
  set(dt, NULL, colnam_cost, 0)
  dt[attendees_sc == 1L & alcohol_curr_xps >= 16, hc_eff := rbinom(.N, 1, l$sc_ls_alcoholpct)]
  dt[hc_eff == 1L, (colnam_cost) := l$sc_ls_alcoholreduc_cost_ind] # Cost only the year of referral
  dt[, hc_eff := hc_effect(hc_eff, 0.8, pid_mrk)] # TODO 0.8 to advanced settings
  dt[hc_eff == 1L, (colnam) := alcohol_sc * (1 - l$sc_ls_alcoholreduc)]


  # Smoking cessation
  set(dt, NULL, "hc_eff", 0L)
  colnam_status   <- "smok_status_sc"
  colnam_quit_yrs <- "smok_quit_yrs_sc"
  colnam_dur      <- "smok_dur_sc"
  colnam_cig      <- "smok_cig_sc"
  colnam_cost     <- "smoking_cost_sc"
  dt[, (c(colnam_status, colnam_quit_yrs, colnam_dur, colnam_cig)) :=
       .(smok_status_curr_xps, smok_quit_yrs_curr_xps, smok_dur_curr_xps, smok_cig_curr_xps)]
  set(dt, NULL, colnam_cost, 0)
  dt[attendees_sc == 1L & smok_status_curr_xps == "4",
     hc_eff := rbinom(.N, 1, l$sc_ls_smkcess)]
  dt[hc_eff == 1L, (colnam_cost) := l$sc_ls_smkcess_cost_ind] # Cost only the year of referral

  # Handle smok_relapse probabilities
  tbl <-
    read_fst("./lifecourse_models/smok_relapse_table.fst",
             as.data.table = TRUE)
  tbl <- dcast(tbl, sex + qimd ~ smok_quit_yrs, value.var = "pr")
  nam <- tbl[, paste0(sex, " ", qimd)]
  tbl <- as.matrix(tbl[, mget(paste0(1:15))], rownames = nam)

  dt[, (c(colnam_status, colnam_quit_yrs, colnam_dur)) :=
       simsmok_cessation(smok_status_sc, smok_quit_yrs_sc,
                         smok_dur_sc, sex, qimd, pid_mrk, hc_eff,
                         dqrunif(.N), tbl, design$smoking_relapse_limit)]

  dt[, smok_status_sc := factor(smok_status_sc)]
  # needed for QRisk and QDrisk
  dt[, smoke_cat_sc := 0L]
  dt[smok_status_sc == "3", smoke_cat_sc := 1L]
  dt[smok_status_sc == "4", smoke_cat_sc := 3L]
  dt[smok_status_sc == "4" & smok_cig_sc < 10L, smoke_cat_sc := 2L]
  dt[smok_status_sc == "4" & smok_cig_sc > 19L, smoke_cat_sc := 4L]

  dt[, hc_eff := NULL]
  invisible(dt)
}


#' @export
run_scenario <- function(scenario_nam, dt, parameters_dt, lags_mc, mc, design,
                         output, timing = c(TRUE, FALSE)) {
  if (timing[[1]]) ptm <- proc.time()
  # The order is important
  set_eligible(scenario_nam, dt, parameters_dt)
  set_invitees(scenario_nam, dt, parameters_dt)
  set_attendees(scenario_nam, dt, parameters_dt)
  set_px(scenario_nam, dt, parameters_dt) # Too slow
  set_lifestyle(scenario_nam, dt, parameters_dt, design)

  # NOTE commented code below is wrong
  # TODO I can calculate the effect of xps change to disease prb for efficiency
  # No need to recalculate disease probability for everyone
  # only apply disease impact on attendees (works only with kismet == TRUE)
  # if (design$kismet) {
  #   nam <- grep("^prb_|_prvl$|_dgn$", names(dt), value = TRUE)
  #   nam <- grep("^rn_|^ckd|^ra_|^cst", nam, value = TRUE, invert = TRUE)
  #   nam2 <- paste0(nam, "_sc")
  #   for (colnam in nam) {
  #     set(dt, NULL, paste0(colnam, "_sc"), dt[[colnam]])
  #   }
  #
  #
  #   tt <- dt[attendees_sc > 0L, unique(pid)]
  #   if (length(tt) > 0) {
  #     dtf <- dt[pid %in% tt]
  #     setkey(dtf, pid, year) # to be sure
  #     nam <- grep("_prvl_sc$|_dgn_sc$", names(dtf), value = TRUE)
  #     for (i in nam) set(dtf, NULL, i, NULL)
  #   }  else dtf <- dt # inefficient but can't break from *apply
  # } else dtf <- dt

  af_model(scenario_nam, mc, dt, design, timing = timing[[2]])
  htn_model(scenario_nam, mc, dt, design, timing = timing[[2]])
  t2dm_model(scenario_nam, mc, dt, design, lags_mc, timing = timing[[2]])
  chd_model(scenario_nam, mc, dt, design, lags_mc, timing = timing[[2]])
  stroke_model(scenario_nam, mc, dt, design, lags_mc, timing = timing[[2]])
  poststroke_dementia_model(scenario_nam, mc, dt, design, lags_mc, timing = timing[[2]])
  copd_model(scenario_nam, mc, dt, design, lags_mc, timing = timing[[2]])
  lung_ca_model(scenario_nam, mc, dt, design, lags_mc, timing = timing[[2]])
  colon_ca_model(scenario_nam, mc, dt, design, lags_mc, timing = timing[[2]])
  breast_ca_model(scenario_nam, mc, dt, design, lags_mc, timing = timing[[2]])
  nonmodelled_model(scenario_nam, mc, dt, design, lags_mc, timing = timing[[2]])

  # if (design$kismet) absorb_dt(dt, dtf[, .SD, .SDcols = patterns("_sc$|pid|year")],
  #                              on = c("pid", "year"))
  # rm(dtf)

  # TODO export scenario xps for additional validation here
  # if (design$export_xps) export_xps(mc_iter, POP, TRUE, "xps_output.csv")
  # if (design$validation) {
  #   export_all_incd(mc_iter, POP, TRUE)
  #   export_all_prvl(mc_iter, POP, TRUE)
  #   export_mrtl(mc_iter, POP, TRUE)
  # }

  output <- gen_output(scenario_nam, design, lags_mc, dt, output)

  dt[, (grep("_sc$", names(dt), value = TRUE)) := NULL]

  if (timing[[1]]) print(proc.time() - ptm)
  invisible(output)
}

#' @export
finalise_synthpop <- function(mc, dt, design, lags_mc, timing = c(TRUE, FALSE)) {
  message("Finalising synthpop...")
  if (timing[[1]]) ptm <- proc.time()
  init_prevalence(mc, dt, design, timing = timing[[2]])
  af_model("", mc, dt, design, timing = timing[[2]])
  htn_model("", mc, dt, design, timing = timing[[2]])
  t2dm_model("", mc, dt, design, lags_mc, timing = timing[[2]])
  nonmodelled_model("", mc, dt, design, lags_mc, timing = timing[[2]])
  chd_model("", mc, dt, design, lags_mc, timing = timing[[2]])
  stroke_model("", mc, dt, design, lags_mc, timing = timing[[2]])
  poststroke_dementia_model("", mc, dt, design, lags_mc, timing = timing[[2]])
  copd_model("", mc, dt, design, lags_mc, timing = timing[[2]])
  lung_ca_model("", mc, dt, design, lags_mc, timing = timing[[2]])
  colon_ca_model("", mc, dt, design, lags_mc, timing = timing[[2]])
  breast_ca_model("", mc, dt, design, lags_mc, timing = timing[[2]])
  if (timing[[1]]) print(proc.time() - ptm)
}

# Function for timing log
#' @export
time_mark <- function(x, file_nam = output_dir("times.txt")) {
  sink(
    file = file_nam,
    append = TRUE,
    type = "output",
    split = FALSE
  )
  cat(paste0(x, " at: ", Sys.time(), "\n"))
  sink()
}

#' @export
export_xps <- function(mc_, dt,
                       write_to_disk = TRUE, filenam = "val_xps_output.csv",
                       reweighted_to_hse = FALSE) {
  to_agegrp(dt, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
  setnames(dt, "t2dm_prvl_curr_xps", "t2dm_prvl_original")
  if ("t2dm_prvl" %in% names(dt)) {
    dt[, t2dm_prvl_curr_xps := fifelse(t2dm_prvl == 0L, 0L, 1L)]
  } else {
    dt[, t2dm_prvl_curr_xps := fifelse(
      t2dm_prvl_original == 0L, 0L, 1L)]
  }
  dt[, smok_never_curr_xps := fifelse(smok_status_curr_xps == "1", 1L, 0L)]
  dt[, smok_active_curr_xps := fifelse(smok_status_curr_xps == "4", 1L, 0L)]
  if (reweighted_to_hse) {
    wt <- read_fst("./synthpop/hse_sociodemographics.fst",
                   as.data.table = TRUE)
    absorb_dt(dt, wt)
  } else {
    dt[, hse_wt := 1]
  }

  xps <- grep("_curr_xps$", names(dt), value = TRUE)
  xps <- xps[-which(xps == "smok_status_curr_xps")]
  out_xps <- groupingsets(
    dt,
    j = lapply(.SD, weighted.mean, hse_wt),
    by = c("year", "sex", "agegrp20", "qimd"),
    .SDcols = xps,
    sets = list(c("year", "sex", "agegrp20", "qimd"),
                c("year", "sex"), c("year", "agegrp20"), c("year", "qimd"))
  )[, `:=` (year = year + 2000L, mc = mc_)]
  for (j in seq_len(ncol(out_xps)))
    set(out_xps, which(is.na(out_xps[[j]])), j, "All")
  dt[, c("agegrp20", "smok_never_curr_xps", "smok_active_curr_xps",
         "t2dm_prvl_curr_xps", "hse_wt") := NULL]
  setnames(dt, "t2dm_prvl_original", "t2dm_prvl_curr_xps")

  setkey(out_xps, year)
  if (write_to_disk) fwrite_safe(out_xps, output_dir(filenam))
  invisible(out_xps)
}

#' @export
export_mrtl <- function(mc_, dt, write_to_disk = TRUE, filenam = "val_mrtl_output.csv") {
  to_agegrp(dt, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
  out <- groupingsets(
    dt,
    j = .N,
    by = c("year", "sex", "agegrp20", "qimd", "all_cause_mrtl"),
    sets = list(c("year", "sex", "agegrp20", "qimd", "all_cause_mrtl"),
                c("year", "sex", "all_cause_mrtl"), c("year", "agegrp20", "all_cause_mrtl"), c("year", "qimd", "all_cause_mrtl"))
  )[, `:=` (year = year + 2000L, mc = mc_)]

  out <- dcast(out, year + agegrp20 + sex + qimd + mc ~ all_cause_mrtl, value.var = "N")
  setnafill(out, "c", 0L, cols = paste0(0:7))
  out[, pops := do_cols_dt(out, paste0(0:7), "+")]
  for (j in seq_len(ncol(out)))
    set(out, which(is.na(out[[j]])), j, "All")
  out[, `0` := NULL]
  setnames(out,
           paste0(1:7),
           paste0("deaths_", c("nonmodelled", "chd", "stroke", "copd", "lung_ca", "colon_ca", "breast_ca")))

  dt[, c("agegrp20") := NULL]
  setkey(out, year)
  if (write_to_disk) fwrite_safe(out, output_dir(filenam))
  invisible(out)
}

#' @export
export_incd <- function(mc_, dt, write_to_disk = TRUE, filenam = "val_incd_output.csv") {
  to_agegrp(dt, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
  xps <- grep("_ca_prvl$", names(dt), value = TRUE)

  out <- groupingsets(
    dt,
    j = lapply(.SD, function(x) sum(x == 1L)),
    by = c("year", "sex", "agegrp20", "qimd"),
    .SDcols = xps,
    sets = list(c("year", "sex", "agegrp20", "qimd"),
                c("year", "sex"), c("year", "agegrp20"), c("year", "qimd"))
  )[, `:=` (year = year + 2000L, mc = mc_)]
  for (j in seq_len(ncol(out)))
    set(out, which(is.na(out[[j]])), j, "All")

  pop <- groupingsets(
    dt,
    j = .N,
    by = c("year", "sex", "agegrp20", "qimd"),
    sets = list(c("year", "sex", "agegrp20", "qimd"),
                c("year", "sex"), c("year", "agegrp20"), c("year", "qimd"))
  )[, `:=` (year = year + 2000L)]
  for (j in seq_len(ncol(pop)))
    set(pop, which(is.na(pop[[j]])), j, "All")
  setnames(pop, "N", "pops")
  absorb_dt(out, pop)
  xps2 <-  gsub("_prvl$", "", xps)
  xps2 <- paste0("cases_", xps2)
  setnames(out, xps, xps2)

  dt[, c("agegrp20") := NULL]
  setkey(out, year)
  if (write_to_disk) fwrite_safe(out, output_dir(filenam))
  invisible(out)
}


#' @export
export_all_incd <- function(mc_, dt, write_to_disk = TRUE, filenam = "val_all_incd_output.csv") {
  to_agegrp(dt, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
  xps <- grep("_prvl$", names(dt), value = TRUE)

  out <- groupingsets(
    dt,
    j = lapply(.SD, function(x) sum(x == 1L)),
    by = c("year", "sex", "agegrp20", "qimd"),
    .SDcols = xps,
    sets = list(c("year", "sex", "agegrp20", "qimd"),
                c("year", "sex"), c("year", "agegrp20"), c("year", "qimd"))
  )[, `:=` (year = year + 2000L, mc = mc_)]
  for (j in seq_len(ncol(out)))
    set(out, which(is.na(out[[j]])), j, "All")

  pop <- groupingsets(
    dt,
    j = .N,
    by = c("year", "sex", "agegrp20", "qimd"),
    sets = list(c("year", "sex", "agegrp20", "qimd"),
                c("year", "sex"), c("year", "agegrp20"), c("year", "qimd"))
  )[, `:=` (year = year + 2000L)]
  for (j in seq_len(ncol(pop)))
    set(pop, which(is.na(pop[[j]])), j, "All")
  setnames(pop, "N", "pops")
  absorb_dt(out, pop)
  xps2 <-  gsub("_prvl$", "", xps)
  xps2 <- paste0("cases_", xps2)
  setnames(out, xps, xps2)

  dt[, c("agegrp20") := NULL]
  setkey(out, year)
  if (write_to_disk) fwrite_safe(out, output_dir(filenam))
  invisible(out)
}

#' @export
export_all_prvl <- function(mc_, dt, write_to_disk = TRUE, filenam = "val_all_prvl_output.csv") {
  to_agegrp(dt, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
  xps <- grep("_prvl$", names(dt), value = TRUE)

  out <- groupingsets(
    dt,
    j = lapply(.SD, function(x) sum(x > 0L)),
    by = c("year", "sex", "agegrp20", "qimd"),
    .SDcols = xps,
    sets = list(c("year", "sex", "agegrp20", "qimd"),
                c("year", "sex"), c("year", "agegrp20"), c("year", "qimd"))
  )[, `:=` (year = year + 2000L, mc = mc_)]
  for (j in seq_len(ncol(out)))
    set(out, which(is.na(out[[j]])), j, "All")

  pop <- groupingsets(
    dt,
    j = .N,
    by = c("year", "sex", "agegrp20", "qimd"),
    sets = list(c("year", "sex", "agegrp20", "qimd"),
                c("year", "sex"), c("year", "agegrp20"), c("year", "qimd"))
  )[, `:=` (year = year + 2000L)]
  for (j in seq_len(ncol(pop)))
    set(pop, which(is.na(pop[[j]])), j, "All")
  setnames(pop, "N", "pops")
  absorb_dt(out, pop)
  xps2 <-  gsub("_prvl$", "", xps)
  xps2 <- paste0("cases_", xps2)
  setnames(out, xps, xps2)

  dt[, c("agegrp20") := NULL]
  setkey(out, year)
  if (write_to_disk) fwrite_safe(out, output_dir(filenam))
  invisible(out)
}

#' @export
sim_init_prvl <- function(mc, disease, design = design, dt = POP) {
  tbl <- get_disease_epi_mc(mc, disease, "p", "v", design$stochastic)
  incd <-  get_disease_epi_mc(mc, disease, "i", "v", design$stochastic)
  absorb_dt(tbl, incd)
  tbl[, prevalence := prevalence - incidence]
  tbl[, year := design$init_year] # given init year is 13 and close enough to 11 which I have epi for
  col_nam <- setdiff(names(tbl), names(dt))
  absorb_dt(dt, tbl)
  nam <- paste0(disease, "_prvl")
  dt[year == design$init_year, (nam) := as.integer(dqrunif(.N) < prevalence)] # faster than rbinom
  setnafill(dt, "c", 0L, cols = nam)
  dt[, (col_nam) := NULL]
}

# Calibrate mortality
# The first run only uses nonmodelled mortality and all other mortalities are set to 0
# Then calculate expected disease mortality based on population size by age sex qimd
#' @export
calibrate_to_mrtl <- function(mc_, dt, disease_, design, calibration_factor = 0) {
  # for disease mortality
  disease_prvl <- paste0(disease_, "_prvl")
  mrtl <- get_lifetable_all(mc_, disease_, design, "qx")

  # FIX for prvl cancer >
  # if (grepl("_ca$", disease_)) {
  #
  #
  #
  #   # dt[get(disease_prvl) > design$cancer_cure, (disease_prvl) := 0L]
  all_mrtl <- get_lifetable_all(mc_, "nonmodelled", design, "qx")
  mrtl[all_mrtl, on = .(year, age, sex, qimd), qx_mc := clamp(qx_mc + i.qx_mc * calibration_factor, 0, Inf)]
  # }

  absorb_dt(dt, mrtl)

  # Only keep prevalent cases
  dt_prvl <- dt[get(disease_prvl) > 0L, ]
  # clone it to minimise rounding error from low prevalence
  mltp_factor <- 10L
  dt_prvl <- clone_dt(dt_prvl, mltp_factor )
  # update paid to eradicate duplicates
  setkey(dt_prvl, .id, pid, year)

  dt_prvl[, pid := rleidv(.SD, cols = c(".id", "pid"))]

  # start loop over years
  for (yr in design$init_year:max(dt$year)) {
    # sample simulants to kill
    pid_to_kill <- dt_prvl[year == yr & (disease_prvl) > 0L,
                           .(pid = pid * rbinom(.N, 1L, qx_mc)),
                           keyby = .(age, sex, qimd)][pid > 0, pid]
    # I use rbinom and not sample because sample only uses integer
    # remove them from the prevalence of following years, so they can't be selected
    dt_prvl[year > yr & pid %in% pid_to_kill, (disease_prvl) := 0L]
  }

  prvl <- dt_prvl[get(disease_prvl) > 0L, .(prvl = .N / mltp_factor), keyby = .(year, age, sex, qimd)]
  absorb_dt(prvl, mrtl)
  # tt_pop_size <- generate_pop_adj_for_mrtl(mc_, dt, design)
  tt_pop_size <- dt[, .N, by = .(year, age, sex, qimd)]
  absorb_dt(prvl, tt_pop_size)
  prvl[, ftlt := clamp(qx_mc * N / prvl)]
  prvl[is.na(ftlt), ftlt := 0] # When 0/0
  prvl[, c("qx_mc", "prvl", "N") := NULL]
  dt[, c("qx_mc") := NULL]
  setnames(prvl, "ftlt", paste0("prb_", disease_, "_mrtl"))
  prvl
}


