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

options(future.fork.enable = TRUE) # enable fork in Rstudio
# TODO remove for production
# plan(multiprocess, workers = design$clusternumber)
plan(list(
  tweak(multiprocess, workers = 2L),
  tweak(multiprocess, workers = design$sim_prm$clusternumber / 2L)
))
# for future apply to work within future


options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

strata_for_gui <- c("mc", "friendly_name", design$sim_prm$strata_for_output)



def_col <- viridis(16, option = "D")
def_col_small <- def_col[c(1, 9, 5, 3, 7, 2, 8, 4, 6)]

def_sym <- c("circle-dot", "square", "diamond", "cross", "triangle-up", "pentagon", " star", "hexagon-open-dot", "triangle-down")

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

most_cost_effective <- function(dt) {
 names(dt[year == max(year), .(nmb_cml, mc, friendly_name)
     ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
      ][, friendly_name[which.max(nmb_cml)], by = mc
       ][, first(sort(counts(V1), decreasing = TRUE))])
}

rank_cost_effective <- function(dt) {
  names(dt[year == max(year), .(nmb_cml, mc, friendly_name)][,
   sum_dt(.SD, c("mc", "friendly_name"), character(0))][,
    friendly_name[which.max(nmb_cml)], by = mc][,
     sort(counts(V1), decreasing = TRUE)])
}

most_effective <- function(dt) {
  names(dt[year == max(year),
            .(
              cpp_cml = sum(net_utility_cml)

            ), by = .(mc, friendly_name)
            ][, friendly_name[which.max(cpp_cml)], by = mc
              ][, first(sort(counts(V1), decreasing = TRUE))])
}

most_equitable <- function(dt) {
  # TODO combine rel and abs ineq
  tt <- dt[year == max(year), .(net_utility_cml, nmb_cml, eq5d_cml, pops_cml, mc, friendly_name, qimd)
                   ][, sum_dt(.SD, c("mc", "friendly_name", "qimd"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name"))
  calc_sei(tt, c("mc", "friendly_name"))
  names(tt[, .(nmb_cml = sum(nmb_cml), sei = mean(sei)), by = c("mc", "friendly_name")][, friendly_name[which.max(sei)], by = mc][, head(sort(counts(V1), decreasing = TRUE), 1L)])
}

most_equitable_rel <- function(dt) {
  tt <- dt[year == max(year), .(net_utility_cml, nmb_cml, eq5d_cml, pops_cml, mc, friendly_name, qimd)
           ][, sum_dt(.SD, c("mc", "friendly_name", "qimd"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name"))
  calc_rei(tt, c("mc", "friendly_name"))
  names(tt[, .(nmb_cml = sum(nmb_cml), rei = mean(rei)), by = c("mc", "friendly_name")][, friendly_name[which.max(rei)], by = mc][, head(sort(counts(V1), decreasing = TRUE), 1L)])

}

most_equitable_abs <- function(dt) {
  tt <- dt[year == max(year), .(net_utility_cml, nmb_cml, eq5d_cml, pops_cml, mc, friendly_name, qimd)
           ][, sum_dt(.SD, c("mc", "friendly_name", "qimd"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name"))
  calc_sei(tt, c("mc", "friendly_name"))
  names(tt[, .(nmb_cml = sum(nmb_cml), sei = mean(sei)), by = c("mc", "friendly_name")][, friendly_name[which.max(sei)], by = mc][, head(sort(counts(V1), decreasing = TRUE), 1L)])
}

benefit_cost_ratio <- function(dt) {
  dt[year == max(year),
     .(
       cbr_cml = sum(nmb_cml)/sum(intervention_cost_cml)
     ), by = .(mc, friendly_name)
     ][, mean(cbr_cml), by = friendly_name][, signif(first(sort(V1, decreasing = TRUE)), 2)]
}

incr_cost_effect_ratio <- function() " '£X, fill me!' "

tot_net_monet_benef <- function(dt) {
  dt[year == max(year),
     .(
       nmb_cml = sum(nmb_cml)
     ), by = .(mc, friendly_name)
     ][, mean(nmb_cml), by = friendly_name][, signif(first(sort(V1, decreasing = TRUE)), 2)]
}
reduce_rel_index_ineq <- function() " 'X, fill me!' "   # TODO
increase_abs_index_ineq <- function() " 'X, fill me!' " # TODO
social_care_cost_sav <- function() " '£X, fill me!' "   # TODO
prod_benef <- function() " '£X, fill me!' "             # TODO
inform_care_cost_sav <- function() " '£X, fill me!' "   # TODO



benefit_cost_ratio_cml <-
  function(dt,
           perspective = c("Societal perspective",
                           "Health and social care perspective",
                           "Healthcare perspective"),
           wtp = input$out_wtp_box) {
    # perspective input$health_econ_perspective_checkbox
    # wtp input$out_wtp_box
    if (perspective == "Societal perspective") {
      dt[, .(bcr_cml = (
        -(
          net_healthcare_cost_cml + net_socialcare_cost_cml + net_informal_care_cost_cml - net_productivity_cost_cml
        ) + (net_utility_cml * wtp)
      ) / net_policy_cost_cml)]
    }
  }


most_benefit_cost_ratio <- function(dt,
                                    perspective = c("Societal perspective",
                                                        "Health and social care perspective",
                                                        "Healthcare perspective"),
                                    wtp = input$out_wtp_box, how_many = 1L) {

  names(dt[year == max(year), .(net_healthcare_cost_cml, net_socialcare_cost_cml,
                                net_informal_care_cost_cml, net_productivity_cost_cml,
                                net_utility_cml, net_policy_cost_cml,
                                mc, friendly_name)
           ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
             ][, bcr_cml := benefit_cost_ratio_cml(.SD, perspective = perspective, wtp = wtp)
               ][, friendly_name[which.max(bcr_cml)], by = mc
                 ][, head(sort(counts(V1), decreasing = TRUE), how_many)])
}

most_benefit_cost_ratio_value <- function(dt,
                                    perspective = c("Societal perspective",
                                                    "Health and social care perspective",
                                                    "Healthcare perspective"),
                                    wtp = input$out_wtp_box) {

  dt[year == max(year), .(net_healthcare_cost_cml, net_socialcare_cost_cml,
                                net_informal_care_cost_cml, net_productivity_cost_cml,
                                net_utility_cml, net_policy_cost_cml,
                                mc, friendly_name)
           ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
             ][, bcr_cml := benefit_cost_ratio_cml(.SD, perspective = perspective, wtp = wtp)
               ][, friendly_name[which.max(bcr_cml)], by = mc
                 ][, head(sort(counts(V1), decreasing = TRUE), 1)]

}


rank_benefit_cost_ratio <- function(dt,
                                    perspective = c("Societal perspective",
                                                    "Health and social care perspective",
                                                    "Healthcare perspective"),
                                    wtp = input$out_wtp_box) {

  names(dt[year == max(year), .(net_healthcare_cost_cml, net_socialcare_cost_cml,
                                net_informal_care_cost_cml, net_productivity_cost_cml,
                                net_utility_cml, net_policy_cost_cml,
                                mc, friendly_name)
           ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
             ][, bcr_cml := benefit_cost_ratio_cml(.SD, perspective = perspective, wtp = wtp)
               ]
        [, friendly_name[which.max(bcr_cml)], by = mc][, sort(counts(V1), decreasing = TRUE)])
}

icer_cml <-
  function(dt,
           perspective = c("Societal perspective",
                           "Health and social care perspective",
                           "Healthcare perspective")) {
    # perspective input$health_econ_perspective_checkbox
    # wtp input$out_wtp_box
    if (perspective == "Societal perspective") {
      dt[, .(icer_cml = round(net_policy_cost_cml + net_healthcare_cost_cml +
                           net_socialcare_cost_cml + net_informal_care_cost_cml -
                           net_productivity_cost_cml) / net_utility_cml)]
    }
  }

scn_icer_cml <-
  function(dt,
           friendly_name_,
           perspective = c("Societal perspective",
                           "Health and social care perspective",
                           "Healthcare perspective")) {
    dt[year == max(year) & friendly_name == friendly_name_,
       .(
         net_healthcare_cost_cml,
         net_socialcare_cost_cml,
         net_informal_care_cost_cml,
         net_productivity_cost_cml,
         net_utility_cml,
         net_policy_cost_cml,
         mc
       )][, sum_dt(.SD, c("mc"), character(0))
          ][, icer_cml := icer_cml(.SD, perspective = perspective)
            ][, mean(icer_cml, na.rm = TRUE)]
  }




net_monetary_benefit_cml <-
  function(dt,
           perspective = c("Societal perspective",
                           "Health and social care perspective",
                           "Healthcare perspective"),
           wtp = input$out_wtp_box) {
    # perspective input$health_econ_perspective_checkbox <- "Societal perspective"
    # wtp input$out_wtp_box
    if (perspective == "Societal perspective") {
      dt[, .(nmb_cml = (-(net_healthcare_cost_cml + net_socialcare_cost_cml + net_informal_care_cost_cml - net_productivity_cost_cml) + (net_utility_cml * wtp)) - net_policy_cost_cml)]
    }
  }

most_benefit_cost_ratio_value <- function(dt,
                                          perspective = c("Societal perspective",
                                                          "Health and social care perspective",
                                                          "Healthcare perspective"),
                                          wtp = input$out_wtp_box) {

  dt[year == max(year), .(net_healthcare_cost_cml, net_socialcare_cost_cml,
                          net_informal_care_cost_cml, net_productivity_cost_cml,
                          net_utility_cml, net_policy_cost_cml,
                          mc, friendly_name)
     ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
       ][, nmb_cml := net_monetary_benefit_cml(.SD, perspective = perspective, wtp = wtp)
         ][, friendly_name[which.max(nmb_cml)], by = mc
           ][, head(sort(counts(V1), decreasing = TRUE), 1)]

}


sum_dt <- function(dt, by = c("year", "friendly_name", "mc"), cols_to_exclude = c("scenario", "agegrp", "sex" , "qimd", "ethnicity")) {
  cols_to_exclude <- intersect(cols_to_exclude, names(dt))

  dt[, lapply(.SD, function(x) {
    if (is.numeric(x)) sum(x) else x
  }), keyby = by,
  .SDcols = -cols_to_exclude]
}

median_dt <- function(dt, by = c("year", "friendly_name"), cols_to_exclude = c("scenario", "agegrp", "sex" , "qimd", "ethnicity", "mc"), digits = 0) {
  cols_to_exclude <- intersect(cols_to_exclude, names(dt))

  dt[, lapply(.SD, function(x) {
    if (is.numeric(x)) round(median(x), digits) else x
  }), keyby = by,
  .SDcols = -cols_to_exclude]
}

mean_dt <- function(dt, by = c("year", "friendly_name"), cols_to_exclude = c("scenario", "agegrp", "sex" , "qimd", "ethnicity", "mc"), digits = 0) {
  cols_to_exclude <- intersect(cols_to_exclude, names(dt))

  dt[, lapply(.SD, function(x) {
    if (is.numeric(x)) round(mean(x), digits) else x
  }), keyby = by,
  .SDcols = -cols_to_exclude]
}

# rigit scores
calc_rigit_scores <- function(dt, by = c("mc", "year", "friendly_name")) {
  keydt <- key(dt)
  keydt <- intersect(keydt, names(dt))
  setkey(dt, mc, qimd, friendly_name)
  dt[, qimd2 := cumsum(pops_cml)/sum(pops_cml),
     keyby = by]
  dt[, qimd2 := (qimd2 + c(0, qimd2[seq_len(.N - 1)]))/2,
     keyby = by]
  # dt[, qimd2_median := median(qimd2, na.rm = T), keyby = .(year, scenario, qimd)]
  if (length(keydt)) setkeyv(dt, keydt)
  dt
}

# Absolute equity slope index
calc_sei <- function(dt, by = c("mc", "year", "friendly_name")) {
  dt[, sei := lm(net_utility_cml~qimd2)$coefficients["qimd2"], by = by]
}

# Relative equity slope index (exclude qimd = 1)
calc_rei <- function(dt, by = c("mc", "year", "friendly_name")) {
  dt[, proport_utility_cml := net_utility_cml/eq5d_cml]
  dt[, rei := lm(proport_utility_cml~qimd2)$coefficients["qimd2"], by = by]
}

# source(file.path("metamodel", "metamodel_predict.R"), local = TRUE)

# sum(gtools::rdirichlet(1, rep(1, 180))) # to simulate a 180 length vector that sums to 1
# extract_uptake_table(setDT(readRDS("./output/input_parameters.rds")), 2)[]
# To write when running from shiny server or docker see https://groups.google.com/forum/#!topic/shiny-discuss/srWETT6uL-I


# Table related ----
# default global search value
if (!exists("default_search_tbl_summary")) default_search_tbl_summary <- ""

# default column search values
if (!exists("default_search_columns_tbl_summary")) default_search_columns_tbl_summary <- NULL

