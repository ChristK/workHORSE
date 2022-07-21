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

library(demography)
library(data.table)
library(CKutils)
library(Rcpp)
library(fst)
library(future.apply)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
if (!require(workHORSEmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_local("./Rpackage/workHORSE_model_pkg/")
  library(workHORSEmisc)
}
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multicore)

hor <- 35L # maximum simulation horizon

# load data ------------------------
deaths <- read_fst("./ONS_data/mortality_by_agegroup.fst", as.data.table = TRUE)
deaths[, Mx_allcause  := deaths / pops]
deaths[, Mx_chd       := I20_25 / pops]
deaths[, Mx_stroke    := I60_69 / pops]
deaths[, Mx_colon_ca  := C18 / pops]
deaths[, Mx_lung_ca   := C34 / pops]
deaths[, Mx_breast_ca := C50 / pops]
deaths[, Mx_copd      := J40_47 / pops]
deaths[sex == "women", Mx_nonmodelled := (deaths - I20_25 - I60_69 - C18 -
                                            C34 - C50 - J40_47) / pops]
deaths[sex == "men", Mx_nonmodelled := (deaths - I20_25 - I60_69 - C18 -
                                          C34 - J40_47) / pops] # men breast cancer is nonmodelled



# mortality projections ----
nam <- sort(grep("^Mx_", names(deaths), value = TRUE))
lifetable_all <- data.table(NULL)
for (disease_ in nam) {
  for (sex_ in unique(deaths$sex)) {
    for (qimd_ in unique(deaths$qimd)) {
      message(paste0(gsub("Mx_", "", disease_),
             "_",
             sex_,
             "_",
             gsub(" ", "_", qimd_)))
      if (!(disease_ == "Mx_breast_ca" && sex_ == "men")) {
        power_ <- 0.4
        if (disease_ == "Mx_copd") power_ <- 0.1
        tt <-
          mortality_projection(
            deaths[qimd == qimd_ & sex == sex_, c("year", "agegrp", "sex", "qimd", ..disease_, "pops")],
            disease_,
            gsub("Mx_", "", disease_), hor, "rapca", 0.4, power = power_)
        p <- plot_projection(tt)
        print(p)
        cowplot::ggsave2(
          filename = paste0(gsub("Mx_", "", disease_),
                            "_",
                            sex_,
                            "_",
                            gsub(" ", "_", qimd_),
                            ".png"),
          p,
          height = 9,
          width = 16,
          units = "cm",
          dpi = 300,
          scale = 2.5,
          path = "./validation/mrtl_projections"
        )
        lifetable_all <- rbind(lifetable_all, tt)
      }
    }
  }
}


# Code like the one below would produce coherent projections, however the fir is not good
# tt <- future_lapply(levels(deaths$qimd), function(x) {
#   mortality_projection(deaths[qimd == x & sex == "men"], "Mx_chd", "chd", hor, "M", 0.7)
# })
# tt <- rbindlist(tt)
# plot_projection(tt)
# lifetable_all <- rbind(lifetable_all, tt)


# calculate qx from mx ----
# lifetable_all[, sex := tolower(sex)] # hardcoded to demography:::lt()
setkey(lifetable_all, age)
tt <- lifetable_all[, lapply(.SD, function(x) {
  if (max(x) <= 1L) { # to exclude age col
    return(demography:::lt(x , min(age), 1L, sex)$qx)
  } else return(x)
}),
.SDcols = patterns("^mx_|age"),
by = .(year, sex, qimd, type)]
setnames(tt, gsub("mx_", "qx_", names(tt)))
absorb_dt(lifetable_all, tt)

setcolorder(lifetable_all,
            c("year", "age", "sex", "qimd", "type", "qx_total", "mx_total"))

lifetable_all[, year := year - 2000L]
lifetable_all[, qimd := gsub("_", " ", qimd)]
setnames(lifetable_all, "type", "disease")
setkey(lifetable_all, disease, year, age, sex, qimd)

write_fst(lifetable_all, "./lifecourse_models/mortality_projections.fst")

# create a table with row numbers for each mc
lifetable_all[, rn := .I]
tt <- lifetable_all[, .(from = min(rn), to = max(rn)), keyby = .(disease)]
write_fst(tt, "./lifecourse_models/mortality_projections_indx.fst", 100L)
