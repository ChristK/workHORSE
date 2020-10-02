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

# Generate uncertainty (1e3 mc iterations) assuming 97.5 percentile is up to 20%
# higher than the median for costs
mc_max <- 1e3
uncertainty_space <- 1.2


library(data.table)
library(fst)
library(future.apply)
library(dqrng)
library(extraDistr)
if (!require(CKutils)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_github("ChristK/CKutils", force = TRUE)
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

options(
  future.fork.enable = TRUE, # enable fork in Rstudio
  future.globals.maxSize = +Inf,
  future.globals.onReference = "ignore"
)
plan(multiprocess, workers = 20L)


dqRNGkind("pcg64")
SEED <- 8950040 # sample(1e7, 1)
set.seed(SEED) # Current is to ensure reproducibility.
dqset.seed(SEED) # Ensure that seed hasn't been used elsewhere in the model

fit.betapr <- # NOT VECTORISED
  function(q, p = c(0.025, 0.5, 0.975), fnscale = 1, uncertainty_space = 1.2) {
    if (length(q) == 1L)
      q <- c(q * (2 - uncertainty_space), q, q * uncertainty_space)
    ofn <- function(x)
      sum((q - qbetapr(p, x[1], x[2], x[3])) ^ 2)
    osol <-
      stats::optim(
        c(1, 1, 1),
        ofn,
        method = "L-BFGS-B",
        control = list(
          "fnscale" = fnscale,
          "maxit"   = 1e7,
          "ndeps"   = rep(1e-3, 3)
        ),
        lower = c(1, 1, 1),
        upper = c(Inf, Inf, Inf)
      )
    out <- as.list(osol$par)
    out <- stats::setNames(out, c("shape1", "shape2", "scale"))
    return(out)
  }

# Informal care ----
informal_costs <- fread("./simulation/health_econ/input/costs_informal_care.csv")
informal_costs[, mltp_factor := 100]
informal_costs[informal_care_cost * uncertainty_space > 1e3, mltp_factor := 5e3]
informal_costs[informal_care_cost * uncertainty_space > 1e4, mltp_factor := 2e4]
informal_costs[, summary(informal_care_cost/mltp_factor)]

informal_costs[, c("shape1", "shape2", "scale") := {
  l <- future_lapply(informal_care_cost/mltp_factor, fit.betapr, future.seed = TRUE)
  list(sapply(l, `[[`, 1), sapply(l, `[[`, 2), sapply(l, `[[`, 3))
}]

# calibration
informal_costs[, mltp_factor := informal_care_cost/qbetapr(0.5, shape1, shape2, scale)]

# internal validation
# informal_costs[, informal_care_cost_predl := mltp_factor * qbetapr(0.025, shape1, shape2, scale)]
# informal_costs[, informal_care_cost_pred  := mltp_factor * qbetapr(0.5, shape1, shape2, scale)]
# informal_costs[, informal_care_cost_predu := mltp_factor * qbetapr(0.975, shape1, shape2, scale)]
#
# informal_costs[, plot(informal_care_cost, informal_care_cost_pred - informal_care_cost)]
# informal_costs[, summary(informal_care_cost_pred/informal_care_cost)]
# informal_costs[, summary(informal_care_cost_predl/informal_care_cost)]
# informal_costs[, summary(informal_care_cost_predu/informal_care_cost)]

informal_costs_l <- clone_dt(informal_costs, mc_max, "mc")

informal_costs_l[, median_informal_care_cost := informal_care_cost]
informal_costs_l[, p := dqrunif(1), by = mc]
informal_costs_l[,informal_care_cost := mltp_factor * qbetapr(p, shape1, shape2, scale)]
informal_costs_l[, c("p", "shape1", "shape2", "scale", "mltp_factor") := NULL]
informal_costs_l[, sex := factor(sex)]

setkey(informal_costs_l, mc)
write_fst(informal_costs_l, "./simulation/health_econ/informal_care_costs_l.fst", 100L)
# create a table with row numbers for each mc
informal_costs_l[, rn := .I]
tt <- informal_costs_l[, .(from = min(rn), to = max(rn)), keyby = mc]
write_fst(tt, "./simulation/health_econ/informal_care_costs_indx.fst", 100L)
rm(informal_costs_l, informal_costs, tt)

# Productivity ----
productivity_costs <- fread("./simulation/health_econ/input/costs_productivity_use.csv")
productivity_costs[, mltp_factor := 1]
productivity_costs[productivity_cost * uncertainty_space > 1e1, mltp_factor := 5]
productivity_costs[productivity_cost * uncertainty_space > 1e2, mltp_factor := 5e1]
productivity_costs[productivity_cost * uncertainty_space > 1e3, mltp_factor := 5e2]
productivity_costs[productivity_cost * uncertainty_space > 1e4, mltp_factor := 5e3]
productivity_costs[, summary(productivity_cost/mltp_factor)]

productivity_costs[, c("shape1", "shape2", "scale") := {
  l <- future_lapply(productivity_cost/mltp_factor, fit.betapr, future.seed = TRUE)
  list(sapply(l, `[[`, 1), sapply(l, `[[`, 2), sapply(l, `[[`, 3))
}]

# calibration
productivity_costs[, mltp_factor := productivity_cost/qbetapr(0.5, shape1, shape2, scale)]

# internal validation
# productivity_costs[, productivity_cost_predl := mltp_factor * qbetapr(0.025, shape1, shape2, scale)]
# productivity_costs[, productivity_cost_pred  := mltp_factor * qbetapr(0.5, shape1, shape2, scale)]
# productivity_costs[, productivity_cost_predu := mltp_factor * qbetapr(0.975, shape1, shape2, scale)]
# productivity_costs[, plot(productivity_cost, productivity_cost_pred - productivity_cost)]
# productivity_costs[, summary(productivity_cost_pred/productivity_cost)]
# productivity_costs[, summary(productivity_cost_predl/productivity_cost)]
# productivity_costs[, summary(productivity_cost_predu/productivity_cost)]

productivity_costs_l <- clone_dt(productivity_costs, mc_max, "mc")

productivity_costs_l[, median_productivity_cost := productivity_cost]
productivity_costs_l[, p := dqrunif(1), by = mc]
productivity_costs_l[, productivity_cost := mltp_factor * qbetapr(p, shape1, shape2, scale)]
productivity_costs_l[, c("p", "shape1", "shape2", "scale", "mltp_factor") := NULL]
productivity_costs_l[, sex := factor(sex)]

setkey(productivity_costs_l, mc)
write_fst(productivity_costs_l, "./simulation/health_econ/productivity_costs_l.fst", 100L)
# create a table with row numbers for each mc
productivity_costs_l[, rn := .I]
tt <- productivity_costs_l[, .(from = min(rn), to = max(rn)), keyby = mc]
write_fst(tt, "./simulation/health_econ/productivity_costs_indx.fst", 100L)
rm(productivity_costs_l, productivity_costs, tt)

# Social care ----
socialcare_costs <- fread("./simulation/health_econ/input/costs_socialcare.csv")
socialcare_costs[, mltp_factor := 1]
socialcare_costs[socialcare_cost * uncertainty_space > 1e1, mltp_factor := 5]
socialcare_costs[socialcare_cost * uncertainty_space > 1e2, mltp_factor := 5e1]
socialcare_costs[socialcare_cost * uncertainty_space > 1e3, mltp_factor := 5e2]
socialcare_costs[socialcare_cost * uncertainty_space > 1e4, mltp_factor := 5e3]
socialcare_costs[, summary(socialcare_cost/mltp_factor)]

socialcare_costs[, c("shape1", "shape2", "scale") := {
  l <- future_lapply(socialcare_cost/mltp_factor, fit.betapr, future.seed = TRUE)
  list(sapply(l, `[[`, 1), sapply(l, `[[`, 2), sapply(l, `[[`, 3))
}]

# calibration
socialcare_costs[, mltp_factor := socialcare_cost/qbetapr(0.5, shape1, shape2, scale)]

# internal validation
# socialcare_costs[, socialcare_cost_predl := mltp_factor * qbetapr(0.025, shape1, shape2, scale)]
# socialcare_costs[, socialcare_cost_pred  := mltp_factor * qbetapr(0.5, shape1, shape2, scale)]
# socialcare_costs[, socialcare_cost_predu := mltp_factor * qbetapr(0.975, shape1, shape2, scale)]
# socialcare_costs[, plot(socialcare_cost, socialcare_cost_pred - socialcare_cost)]
# socialcare_costs[, summary(socialcare_cost_pred/socialcare_cost)]
# socialcare_costs[, summary(socialcare_cost_predl/socialcare_cost)]
# socialcare_costs[, summary(socialcare_cost_predu/socialcare_cost)]

socialcare_costs_l <- clone_dt(socialcare_costs, mc_max, "mc")

socialcare_costs_l[, median_socialcare_cost := socialcare_cost]
socialcare_costs_l[, p := dqrunif(1), by = mc]
socialcare_costs_l[, socialcare_cost := mltp_factor * qbetapr(p, shape1, shape2, scale)]
socialcare_costs_l[, c("p", "shape1", "shape2", "scale", "mltp_factor") := NULL]

setkey(socialcare_costs_l, mc, age)
write_fst(socialcare_costs_l, "./simulation/health_econ/socialcare_costs_l.fst", 100L)
# create a table with row numbers for each mc
socialcare_costs_l[, rn := .I]
tt <- socialcare_costs_l[, .(from = min(rn), to = max(rn)), keyby = mc]
write_fst(tt, "./simulation/health_econ/socialcare_costs_indx.fst", 100L)
rm(socialcare_costs_l, socialcare_costs, tt)

# added diseases
socialcare_costs_added_diseases <-
  fread("./simulation/health_econ/input/costs_socialcare_added_diseases.csv")
socialcare_costs_added_diseases[, mltp_factor := 1]
socialcare_costs_added_diseases[socialcare_cost * uncertainty_space > 1e1, mltp_factor := 5]
socialcare_costs_added_diseases[socialcare_cost * uncertainty_space > 1e2, mltp_factor := 5e1]
socialcare_costs_added_diseases[socialcare_cost * uncertainty_space > 1e3, mltp_factor := 5e2]
socialcare_costs_added_diseases[socialcare_cost * uncertainty_space > 1e4, mltp_factor := 5e3]
socialcare_costs_added_diseases[, summary(socialcare_cost/mltp_factor)]

socialcare_costs_added_diseases[, c("shape1", "shape2", "scale") := {
  l <- future_lapply(socialcare_cost/mltp_factor, fit.betapr, future.seed = TRUE)
  list(sapply(l, `[[`, 1), sapply(l, `[[`, 2), sapply(l, `[[`, 3))
}]

# calibration
socialcare_costs_added_diseases[, mltp_factor := socialcare_cost/qbetapr(0.5, shape1, shape2, scale)]

# internal validation
# socialcare_costs_added_diseases[, socialcare_cost_predl := mltp_factor * qbetapr(0.025, shape1, shape2, scale)]
# socialcare_costs_added_diseases[, socialcare_cost_pred  := mltp_factor * qbetapr(0.5, shape1, shape2, scale)]
# socialcare_costs_added_diseases[, socialcare_cost_predu := mltp_factor * qbetapr(0.975, shape1, shape2, scale)]
# socialcare_costs_added_diseases[, plot(socialcare_cost, socialcare_cost_pred - socialcare_cost)]
# socialcare_costs_added_diseases[, summary(socialcare_cost_pred/socialcare_cost)]
# socialcare_costs_added_diseases[, summary(socialcare_cost_predl/socialcare_cost)]
# socialcare_costs_added_diseases[, summary(socialcare_cost_predu/socialcare_cost)]

socialcare_costs_added_diseases_l <- clone_dt(socialcare_costs_added_diseases, mc_max, "mc")

socialcare_costs_added_diseases_l[, median_socialcare_cost := socialcare_cost]
socialcare_costs_added_diseases_l[, p := dqrunif(1), by = mc]
socialcare_costs_added_diseases_l[, socialcare_cost := mltp_factor * qbetapr(p, shape1, shape2, scale)]
socialcare_costs_added_diseases_l[, c("p", "shape1", "shape2", "scale", "mltp_factor") := NULL]

setkey(socialcare_costs_added_diseases_l, mc)
write_fst(socialcare_costs_added_diseases_l, "./simulation/health_econ/socialcare_costs_added_diseases_l.fst", 100L)
# create a table with row numbers for each mc
socialcare_costs_added_diseases_l[, rn := .I]
tt <- socialcare_costs_added_diseases_l[, .(from = min(rn), to = max(rn)), keyby = mc]
write_fst(tt, "./simulation/health_econ/socialcare_costs_added_diseases_indx.fst", 100L)
rm(socialcare_costs_added_diseases_l, socialcare_costs_added_diseases, tt)

# Healthcare ----
healthcare_costs <- fread("./simulation/health_econ/input/costs_healthcare.csv")
healthcare_costs[, mltp_factor := 1]
healthcare_costs[healthcare_cost * uncertainty_space > 1e1, mltp_factor := 5]
healthcare_costs[healthcare_cost * uncertainty_space > 1e2, mltp_factor := 5e1]
healthcare_costs[healthcare_cost * uncertainty_space > 1e3, mltp_factor := 5e2]
healthcare_costs[healthcare_cost * uncertainty_space > 1e4, mltp_factor := 5e3]
healthcare_costs[, summary(healthcare_cost/mltp_factor)]

healthcare_costs[, c("shape1", "shape2", "scale") := {
  l <- future_lapply(healthcare_cost/mltp_factor, fit.betapr, future.seed = TRUE)
  list(sapply(l, `[[`, 1), sapply(l, `[[`, 2), sapply(l, `[[`, 3))
}]

# calibration
healthcare_costs[, mltp_factor := healthcare_cost/qbetapr(0.5, shape1, shape2, scale)]

# internal validation
# healthcare_costs[, healthcare_cost_predl := mltp_factor * qbetapr(0.025, shape1, shape2, scale)]
# healthcare_costs[, healthcare_cost_pred  := mltp_factor * qbetapr(0.5, shape1, shape2, scale)]
# healthcare_costs[, healthcare_cost_predu := mltp_factor * qbetapr(0.975, shape1, shape2, scale)]
# healthcare_costs[, plot(healthcare_cost, healthcare_cost_pred - healthcare_cost)]
# healthcare_costs[, summary(healthcare_cost_pred/healthcare_cost)]
# healthcare_costs[, summary(healthcare_cost_predl/healthcare_cost)]
# healthcare_costs[, summary(healthcare_cost_predu/healthcare_cost)]

healthcare_costs_l <- clone_dt(healthcare_costs, mc_max, "mc")

healthcare_costs_l[, median_healthcare_cost := healthcare_cost]
healthcare_costs_l[, p := dqrunif(1), by = mc]
healthcare_costs_l[, healthcare_cost := mltp_factor * qbetapr(p, shape1, shape2, scale)]
healthcare_costs_l[, c("p", "shape1", "shape2", "scale", "mltp_factor") := NULL]

setkey(healthcare_costs_l, mc, disease, years_since_diagnosis)
write_fst(healthcare_costs_l, "./simulation/health_econ/healthcare_costs_l.fst", 100L)
# create a table with row numbers for each mc
healthcare_costs_l[, rn := .I]
tt <- healthcare_costs_l[, .(from = min(rn), to = max(rn)), keyby = mc]
write_fst(tt, "./simulation/health_econ/healthcare_costs_indx.fst", 100L)
rm(healthcare_costs_l, healthcare_costs, tt)



