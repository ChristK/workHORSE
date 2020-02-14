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

setwd("~/My Models/workHORSE_WS4/")

# file.remove(list.files("/mnt/storage_slow/synthpop/", "^synthpop_", full.names = TRUE))
# roxygen2::roxygenise("./Rpackage/workHORSE_model_pkg/") # TODO remove before deployment
# remotes::install_local("./Rpackage/workHORSE_model_pkg/", force = TRUE)
# library(workHORSEmisc)

source("./global.R")
localities_indx <- read_fst("./synthpop/lsoa_to_locality_indx.fst", as.data.table = TRUE)
localitities_list <- list(
  "Country" = list("England"),
  "Region" = sort(as.character((localities_indx[, unique(RGN11NM)]))),
  "Local Authority" = sort(as.character((localities_indx[, unique(LAD17NM)])))
)
rm(localities_indx)
names(localitities_list$`Local Authority`) <- localitities_list$`Local Authority`

# design$locality <- localitities_list$Country[[1]] # England
# design$locality <- localitities_list$`Local Authority`[["Liverpool"]]
design$locality <-
  localitities_list$`Local Authority`[c("Corby", "Daventry",
                                         "East Northamptonshire", "South Northamptonshire",
                                         "Kettering", "Northampton", "Wellingborough"
                                         )] # Northamptonshire
# grep("Kettering", localitities_list$`Local Authority`, value = T)



# check_synthpop_has_metafile(synthpop_dir)

### Run foreach inside of a for loop, ensuring that it never receives
### a first argument with a length more than maxjobs. This avoids some
### memory problems (swapping, or getting jobs killed on the cluster)
### when using mclapply(1:N, FUN) where N is large.
# From https://r.789695.n4.nabble.com/mclapply-memory-leak-td4711759.html

X <- 1:100 # Number of iterations
maxjobs <- design$clusternumber/2 # parallel::detectCores()/2L
N <- length(X)
i.list <- splitIndices(N, N/maxjobs)
result.list <- vector("list", N)



for(i in seq_along(i.list)) {
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
                          locality = design$locality,
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

check_synthpop_has_metafile(synthpop_dir)



# plan(multiprocess, workers = design$clusternumber/2)

### Run mclapply inside of a for loop, ensuring that it never receives
### a first argument with a length more than maxjobs. This avoids some
### memory problems (swapping, or getting jobs killed on the cluster)
### when using mclapply(1:N, FUN) where N is large.
# maxjobs_foreach <- function(X, FUN, maxjobs = parallel::detectCores()/2L){
#   N <- length(X)
#   i.list <- splitIndices(N, N/maxjobs)
#   result.list <- vector("list", X)
#   for(i in seq_along(i.list)){
#     i.vec <- i.list[[i]]
#     result.list[i.vec] <- mclapply(X[i.vec], FUN)
#   }
#   result.list
# }


# tt <- read_fst("/mnt/storage_slow/synthpop/synthpop_96bdd875788e8df12e19c709804bd243_4.fst", as.data.table = T)
# nam <- grep("^rn_", names(tt), value = T)
# tt[, (nam) := NULL]
# write_fst(tt, "/mnt/storage_slow/synthpop/000.fst", 100)
# ll <- list.files("/mnt/storage_slow/synthpop", "^synthpop_.*fst$", full.names = TRUE)
# ll <- grep("_indx.fst$", ll, invert = TRUE, value = TRUE)
# for (i in 1:100) {
#   tt <- read_fst(ll[i], as.data.table = T)
#   if (key(tt)[[1]] == "pid") {
#     tt[, rn := .I]
#     tt <- tt[, .(from = min(rn), to = max(rn)), keyby = .(pid)]
#     nam <- ll[i]
#     nam <- gsub(".fst$", "_indx.fst", nam)
#     write_fst(tt, nam, 100L)
#   } else stop()
# }

# Scale-up to ONS population projections
# NOTE scaling-up assumes that the relative distribution of qimd remains stable
# because the population projections are at LAD level
# Also scaling cannot be performed in chunks without introducing signifficant
# error. This is because chunks are not random. One may have all aged 40-50
# and another none. We first need to merge chunks before we scale up

# Calculate weights so that their sum is the population of the area based on ONS
# ll <- list.files("/mnt/storage_slow/synthpop", "^synthpop_.*fst$", full.names = TRUE)
# ll <- grep("_indx.fst$", ll, invert = TRUE, value = TRUE)
# for (i in 1:100) {
#   dt <- read_fst(ll[i], as.data.table = T)
#   calc_pop_weights(dt, design$locality)
#   write_fst(dt, ll[i], 100L)
# }
