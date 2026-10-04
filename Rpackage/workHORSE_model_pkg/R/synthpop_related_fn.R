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


# get all unique LADs (April 2023 boundaries) included in locality vector.
# KEEP!!!
#' @export
get_unique_LADs <- function(locality) {
  indx_hlp <-
    read_fst("./synthpop/lsoa_to_locality_indx.fst",
             columns = c("LAD23CD", "LAD23NM", "RGN11NM"),
             as.data.table = TRUE)

  unknown <- setdiff(locality, c("England", levels(indx_hlp$LAD23NM),
                                 levels(indx_hlp$RGN11NM)))
  if (length(unknown) > 0L)
    stop("Unknown locality: ", paste(unknown, collapse = ", "),
         ". Local authorities use April 2023 boundaries.")

  if ("England" %in% locality) {
    lads <- indx_hlp[, unique(LAD23CD)] # national
  } else {
    lads <-
      indx_hlp[LAD23NM %in% locality |
                 RGN11NM %in% locality, unique(LAD23CD)]
  }
  return(as.character(lads))
}

# Get ONS population (estimates and projections) by year, age and sex for the
# input localities. England uses the 2024-based NPP from 2026, so it is not the
# sum of its LAs (2022-based SNPP). KEEP!!!
#' @export
get_pop_proj <- function(locality) {
  if ("England" %in% locality) { # national
    tt <- read_fst("./ONS_data/pop_size/pop_proj_england.fst",
                   as.data.table = TRUE)
  } else {
    lads <- get_unique_LADs(locality)
    tt <- read_fst("./ONS_data/pop_size/pop_proj.fst",
                   columns = c("year", "age", "sex", "LAD23CD", "pops"),
                   as.data.table = TRUE)[LAD23CD %in% lads]
  }
  tt[, .(pops = sum(pops)), keyby = .(year, age, sex)]
}

# Get dt projections for the input localities. KEEP!!!
#' @export
get_pop_size <- function(design, parameters) {
  tt <- get_pop_proj(parameters$locality_select)
  tt <- tt[between(age, design$ageL, design$ageH) &
             between(year, parameters$ininit_year_slider_sc1,
                     parameters$inout_year_slider)]

  return(tt)
}

