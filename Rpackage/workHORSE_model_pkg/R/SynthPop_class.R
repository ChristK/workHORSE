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



# From
# https://stackoverflow.com/questions/33424233/how-do-i-tell-an-r6-class-what-to-do-with-square-brackets
# Allows data.table syntax to the R6class object directly. Assumes it has a
# field 'pop' that is a data.table
#' @export
`[.SynthPop` = function(x, ...)
  x$pop[...]


#' R6 Class representing a synthetic population
#'
#' @description
#' A synthpop has a `pop` field that contains the life course of simulants in a
#' `data.table`.
#'
#' @details
#' To be completed...
#'
#' @export
SynthPop <-
  R6::R6Class(
    classname = "SynthPop",

    # public ------------------------------------------------------------------
    public = list(
      #' @field mc The Monte Carlo iteration of the synthetic population. Every
      #'   integer generates a unique synthetic population.
      mc = NA,

      #' @field metadata Metadata of the synthpop.
      metadata = NA,

      #' @field pop The data.table that contains the life-course of simulants.
      #'   If the file exists, it is loaded from disk. If it doesn't, it is
      #'   first generated, then saved to disk, and then loaded from disk.
      pop = NA,


      #' @description Create a new SynthPop object.
      #' @param mc_ The Monte Carlo iteration of the synthetic population. Each
      #'   integer generates a unique synthetic population. If `mc = 0` an
      #'   object with an empty synthpop is initiated.
      #' @param design_ A \code{\link[workHORSEmisc]{Design}} object.
      #' @param synthpop_dir_ The directory where 'SynthPop' objects are stored.
      #'   The synthpop file in \code{\link[fst]{fst-package}} format. If
      #'   `filename` already exists, then the synthpop is loaded from there.
      #'   Otherwise it is generated from scratch and then saved as `filename`
      #'   in \code{\link[fst]{fst-package}} format. Two additional files are
      #'   saved for each 'synthpop'. A metadata file, and an index file.
      #' @return A new `SynthPop` object.
      #' @examples
      #' design <- Design$new("./validation/design_for_trends_validation.yaml")
      #' POP$write_synthpop(1:6, design, "./test")
      #' POP <- SynthPop$new(4L, design, "./test")
      #' POP$print()
      #' POP$count_synthpop()
      #'
      #' POP$delete_synthpop(1L)
      #' POP$delete_synthpop(5:6)
      #' POP$get_filename()
      initialize = function(mc_, design_, synthpop_dir_) {
        stopifnot(length(mc_) == 1L, is.numeric(mc_), ceiling(mc_) >= 0L)
        stopifnot("Design" %in% class(design_))
        stopifnot(is.character(synthpop_dir_))

        mc_ <- ceiling(mc_)
        # Create synthpop_dir_ if it doesn't exists
        synthpop_dir_ <-
          base::normalizePath(synthpop_dir_, mustWork = FALSE)
        if (!dir.exists(synthpop_dir_)) {
          dir.create(synthpop_dir_, recursive = TRUE)
          message(paste0("Directory ", synthpop_dir_, " was created"))
        }

        # get unique lsoas
        lsoas <- private$get_unique_LSOAs(design_)

        private$checksum <- private$gen_checksum(design_, lsoas)

        self$mc <- mc_
        private$design <- design_
        private$synthpop_dir <- synthpop_dir_

        if (mc_ > 0) {
          private$filename <- private$gen_synthpop_filename(mc_,
                                                            synthpop_dir_,
                                                            private$checksum)

          # logic for the synthpop load
          files_exist <- sapply(private$filename, file.exists)
          if (all(!files_exist)) {
            # No files exist. Create the synthpop and store the file on disk (no
            # parallelism)
            private$gen_synthpop(mc_,
                                 private$filename,
                                 design_,
                                 lsoas)

          } else if (file.exists(private$filename$metafile) &&
                     !all(files_exist)) {
            # Metafile exists but not all three files. It means that most likely
            # a generate_synthpop() is still running. So the function waits
            # until the file is created before it proceeds to load it. Note that
            # if this is not the case then the loop is infinite!!!
            while (!all(sapply(private$filename, file.exists)))
              Sys.sleep(5)

            # Ensure the file write is complete (size stable)
            sz1 <- file.size(private$filename$synthpop)
            Sys.sleep(3)
            sz2 <- file.size(private$filename$synthpop)
            while (sz1 != sz2) {
              sz1 <- file.size(private$filename$synthpop)
              Sys.sleep(3)
              sz2 <- file.size(private$filename$synthpop)
            }

          } else if (!file.exists(private$filename$metafile) &&
                     !all(files_exist)) {
            # Metafile doesn't exist but some other files exist. In this case
            # delete everything and start from cratch
            self$delete_incomplete_synthpop()
            private$gen_synthpop(mc_,
                                 private$filename,
                                 design_,
                                 lsoas)
          }
          # No need to provision for case when all file present. The following
          # lines handle this case anyway

          self$pop <-
            private$get_synthpop(mc_, private$filename, design_)
          self$metadata <- yaml::read_yaml(private$filename$metafile)

          if (design_$sim_prm$logs)
            self$print()
        }
        invisible(self)
      },

      #' @description
      #' Delete (all) synthpop files in the synthpop directory.
      #' @param mc_ If `mc_ = NULL`, delete all files in the synthpop directory.
      #'   If `mc_` is an integer vector delete the specific synthpop files
      #'   including the metadata and index files.
      #' @return The invisible `SynthPop` object.
      delete_synthpop = function(mc_) {
        if (is.null(mc_)) {
          file.remove(list.files(
            private$synthpop_dir,
            full.names = TRUE,
            recursive = TRUE
          ))
        } else if (length(mc_) == 1L &&
                   is.numeric(mc_) && ceiling(mc_) > 0L) {
          file.remove(unlist(
            private$gen_synthpop_filename(mc_,
                                          private$synthpop_dir,
                                          private$checksum)
          ))

        } else if (length(mc_) > 1L &&
                   all(is.numeric(mc_)) && all(ceiling(mc_) > 0L)) {
          filnam <-
            lapply(mc_,
                   private$gen_synthpop_filename,
                   private$synthpop_dir,
                   private$checksum)
          file.remove(unlist(filnam))

        } else
          message("mc_ need to be NULL or numeric. Nothing was deleted.")

        return(invisible(self))
      },

      #' @description
      #' Check that every synthpop file has a metafile and an index file. Delete
      #' any orphan files.
      #' @return The invisible `SynthPop` object.
      delete_incomplete_synthpop =
        function() {
          files <- list.files(private$synthpop_dir, "^synthpop_.*\\.fst$")
          files <-
            grep("_indx.fst$", files, invert = TRUE, value = TRUE) # remove indx files
          files <- sub("\\.fst$", "", files)
          metafiles <-
            list.files(private$synthpop_dir, "^synthpop_.*_meta\\.yaml$")
          metafiles <- sub("_meta\\.yaml$", "", metafiles)

          to_remove <- setdiff(metafiles, files)
          if (length(to_remove) > 0) {
            to_remove <- paste0(to_remove, "_meta.yaml")
            file.remove(file.path(private$synthpop_dir, to_remove))
          }

          to_remove <- setdiff(files, metafiles)
          if (length(to_remove) > 0) {
            to_remove2 <- paste0(to_remove, ".fst")
            file.remove(file.path(private$synthpop_dir, to_remove2))
            to_remove3 <- paste0(to_remove, "_indx.fst")
            file.remove(file.path(private$synthpop_dir, to_remove3))
          }

          return(invisible(self))
        },

      #' @description
      #' Count the synthpop files in a directory. It includes files without
      #' metafiles and index files.
      #' @return The invisible `SynthPop` object.
      count_synthpop =
        function() {
          files <- list.files(private$synthpop_dir, "^synthpop_.*\\.fst$")
          files <- # remove indx files
            grep("_indx.fst$", files, invert = TRUE, value = TRUE)
          print(length(files))
          return(invisible(self))
        },

      #' @description
      #' Get the synthpop file paths.
      #' @param x One of "all", "synthpop", "metafile", or "indxfile". Can be
      #'   abbreviated.
      #' @return The invisible `SynthPop` object.
      get_filename = function(x = c("all", "synthpop", "metafile", "indxfile")) {
        x <- match.arg(x)
        switch(
          x,
          all      = print(private$filename),
          synthpop = print(private$filename[["synthpop"]]),
          metafile = print(private$filename[["metafile"]]),
          indxfile = print(private$filename[["indxfile"]])
        )
        invisible(self)
      },

      #' @description
      #' Get the synthpop design.
      #' @return The invisible `SynthPop` object.
      get_design = function() {
        # print(private$design)
        # invisible(self)
        private$design
      },

      #' @description
      #' Get the synthpop dir.
      #' @return The invisible `SynthPop` object.
      get_dir = function() {
        print(private$synthpop_dir)
        invisible(self)
      },

      #' @description
      #' Generate synthpop files in parallel and writes them to disk. It skips
      #' files that are already on disk.
      #' @param mc_ An integer vector for the Monte-Carlo iteration of the
      #'   synthetic population. Each integer generates a unique synthetic
      #'   population.
      #' @return The invisible `SynthPop` object.
      write_synthpop = function(mc_) {
        stopifnot(all(is.numeric(mc_)), all(ceiling(mc_) > 0L))

        mc_ <- ceiling(mc_)

        # get unique lsoas
        lsoas <- private$get_unique_LSOAs(private$design)

        checksum <- private$gen_checksum(private$design, lsoas)

        if (Sys.info()["sysname"] == "Windows") {
          cl <- makeCluster(private$design$sim_prm$clusternumber)
          registerDoParallel(cl)
        } else {
          registerDoParallel(private$design$sim_prm$clusternumber)
        }

        on.exit((if (exists("cl")) stopCluster(cl)), add = TRUE)
        on.exit(self$delete_incomplete_synthpop(), add = TRUE)

        foreach(
          mc_iter = mc_,
          .inorder = FALSE,
          .verbose = TRUE,
          .packages = c(
            "gamlss.dist",
            # For distr in prevalence.R
            "dqrng",
            "qs",
            "fst",
            "CKutils",
            "workHORSEmisc",
            "data.table"
          ),
          .export = NULL,
          .noexport = NULL # c("time_mark")
        ) %dopar% {
          data.table::setDTthreads(private$design$sim_prm$n_cpus)
          fst::threads_fst(private$design$sim_prm$n_cpus)
          filename <-
            private$gen_synthpop_filename(mc_iter, private$synthpop_dir, checksum)

          # logic for the synthpop load
          files_exist <- sapply(filename, file.exists)
          if (all(!files_exist)) {
            # No files exist. Create the synthpop and store the file on disk
            private$gen_synthpop(mc_iter,
                                 filename,
                                 private$design,
                                 lsoas)

          } else if (file.exists(filename$metafile) &&
                     !all(files_exist)) {
            # Metafile exists but not all three files. It means that most likely
            # a generate_synthpop() is still running. So the function waits
            # until the file is created before it proceeds to load it. Note that
            # if this is not the case then the loop is infinite!!!
            while (!all(sapply(filename, file.exists)))
              Sys.sleep(5)

            # Ensure the file write is complete (size stable)
            sz1 <- file.size(filename$synthpop)
            Sys.sleep(3)
            sz2 <- file.size(filename$synthpop)
            while (sz1 != sz2) {
              sz1 <- file.size(filename$synthpop)
              Sys.sleep(3)
              sz2 <- file.size(filename$synthpop)
            }

          } else if (!file.exists(filename$metafile) &&
                     !all(files_exist)) {
            # Metafile doesn't exist but some other files exist. In this case
            # delete everything and start from cratch
            self$delete_incomplete_synthpop()
            private$gen_synthpop(mc_iter,
                                 filename,
                                 private$design,
                                 lsoas)
          }
          # No need to provision for case when all files present.

          return(NULL)
        }
        invisible(self)
      },

      #' @description
      #' Prints the synthpop object metadata.
      #' @return The invisible `SynthPop` object.
      print = function() {
        print(c(
          "path" = private$filename$synthpop[[]],
          "mc" = self$mc,
          self$metadata
        ))
        invisible(self)
      }
    ),



    # private -----------------------------------------------------------------
    private = list(
      filename = NA,
      checksum = NA,
      # The design object with the simulation parameters.
      design = NA,
      synthpop_dir = NA,

      # Special deep copy for data.table. Use POP$clone(deep = TRUE) to
      # dispatch. Otherwise a reference is created
      deep_clone = function(name, value) {
        if ("data.table" %in% class(value)) {
          data.table::copy(value)
        } else if ("R6" %in% class(value)) {
          value$clone()
        } else {
          # For everything else, just return it. This results in a shallow
          # copy of s3.
          value
        }
      },


      # get all unique LSOAs included in locality vector
      get_unique_LSOAs = function(design_) {
        indx_hlp <-
          read_fst("./synthpop/lsoa_to_locality_indx.fst",
                   as.data.table = TRUE)

        if ("England" %in% design_$sim_prm$locality) {
          lsoas <- indx_hlp[, unique(LSOA11CD)] # national
        } else {
          lsoas <-
            indx_hlp[LAD17NM %in% design_$sim_prm$locality |
                       RGN11NM %in% design_$sim_prm$locality, unique(LSOA11CD)]
        }
        return(sort(lsoas))
      },

      # get a smaller design list only with characteristics that are important
      # for synthpop creation and define the uniqueness of the object. I.e. if
      # these parameters are different the synthpop has to have different
      # filename and vice-versa
      get_unique_characteristics = function(design_) {
        design_$sim_prm[c(
          "n",
          "sim_horizon_max",
          "init_year_long",
          "maxlag",
          "smoking_relapse_limit",
          "ageL",
          "ageH",
          "jumpiness",
          "simsmok_calibration"
        )]
      },

      # gen synthpop unique checksum for the given set of inputs
      gen_checksum =
        function(design_,
                 lsoas_) {
          # get a md5 checksum based on function arguments
          # First get function call arguments
          fcall <- private$get_unique_characteristics(design_)

          locality_years_age_id <-
            digest::digest(paste(lsoas_, fcall, sep = ",", collapse = ","),
                           serialize = FALSE)
          return(locality_years_age_id)
        },

      # gen synthpop filename for the given set of inputs
      gen_synthpop_filename =
        function(mc_,
                 synthpop_dir_,
                 checksum_) {
          return(
            list(
              "synthpop" = normalizePath(
                paste0(synthpop_dir_,
                       "/synthpop_",
                       checksum_,
                       "_",
                       mc_,
                       ".fst"),
                mustWork = FALSE
              ),
              "metafile" = normalizePath(
                paste0(
                  synthpop_dir_,
                  "/synthpop_",
                  checksum_,
                  "_",
                  mc_,
                  "_meta.yaml"
                ),
                mustWork = FALSE
              ),
              "indxfile" = normalizePath(
                paste0(synthpop_dir_,
                       "/synthpop_",
                       checksum_,
                       "_",
                       mc_,
                       "_indx.fst"),
                mustWork = FALSE
              )
            )
          )
        },

      # generate synthpop demographics
      gen_synthpop_demog =
        function(design_, lsoas_) {
          # locality <- c("East Midlands", "Adur", "Allerdale")

          # load dt
          dt <-
            read_fst("./synthpop/lsoa_mid_year_population_estimates.fst",
                     as.data.table = TRUE)[LSOA11CD %in% lsoas_ &
                                             year %in% design_$sim_prm$init_year_long]
          # delete unwanted ages
          dt[, c(paste0(0:(design_$sim_prm$ageL - 1L)), c(paste0((
            design_$sim_prm$ageH + 1L
          ):90))) := NULL]

          dt <-
            melt(
              dt,
              grep("^[0-9]", names(dt), value = TRUE, invert = TRUE),
              variable.name = "age",
              value.name = "population_size",
              variable.factor = FALSE
            )
          dt[, age := as.integer(age)]
          dt[, year := as.integer(year)]
          dt[, population_size := population_size / sum(population_size)]

          # load ethnicity proportions by lsoa
          rows <-
            read_fst("./synthpop/ethn2011_pct_indx.fst", as.data.table = TRUE)
          rows <- rows[LSOA11CD %in% lsoas_, c(min(from), max(to))]
          ethn <- read_fst(
            "./synthpop/ethn2011_pct.fst",
            from = rows[1],
            to = rows[2],
            as.data.table = TRUE
          )[LSOA11CD %in% lsoas_]

          absorb_dt(dt, ethn)
          .ethn_nam <- c(
            "white",
            "indian",
            "pakistani",
            "bangladeshi",
            "other asian",
            "black caribbean",
            "black african",
            "chinese",
            "other"
          )

          for (j in .ethn_nam)
            set(dt, NULL, j,
                dt[, get(j) * population_size])
          dt[, population_size := NULL]
          dt <-
            melt(
              dt,
              measure.vars = .ethn_nam,
              variable.name = "ethnicity",
              variable.factor = TRUE
            )

          # I do not explicitly set.seed because I do so in the gen_synthpop()
          dt <- dt[sample(.N, design_$sim_prm$n, TRUE, value)]
          dt[, value := NULL]

          indx_hlp <-
            read_fst("./synthpop/lsoa_to_locality_indx.fst",
                     as.data.table = TRUE)

          dt[indx_hlp, on = "LSOA11CD", `:=` (
            tds = i.tds,
            tds_quintile = i.tds_quintile,
            imd = i.imd,
            qimd = i.qimd,
            sha = i.SHA11NM,
            CCG17CDH = CCG17CDH
          )]
          return(invisible(dt))
        },



      gen_synthpop = # returns NULL. Writes synthpop on disk
        function(mc_,
                 filename_,
                 design_,
                 lsoas_) {
          # increase design_$sim_prm$jumpiness for more erratic jumps in trajectories

          # In Shiny app this function runs as a future. It is not
          # straightforward to check whether the future has been resolved or
          # not. To circumvent the problem I will save the metafile here (almost
          # function beginning) and the synthpop file at the end. So if both
          # files exist the function has finished. If only metafile exists the
          # function probably still runs.

          # Save synthpop metadata
          if (!file.exists(filename_$metafile)) {
            yaml::write_yaml(private$get_unique_characteristics(design_),
                             filename_$metafile)
          }
          # NOTE In shiny app if 2 users click the  button at the same time, 2
          # functions will run almost concurrently with potential race condition

          # To avoid edge cases when the function stopped prematurely and a metafile
          # was created while the file was not. On.exit ensures that either both
          # files exist or none

          del_incomplete <- function(filename_) {
            if (file.exists(filename_$metafile) &&
                (!file.exists(filename_$synthpop) ||
                 !file.exists(filename_$indxfile))) {
              suppressWarnings(sapply(filename_, file.remove))
              # file.remove(filename_$metafile)
            }
          }

          on.exit(del_incomplete(filename_), add = TRUE)

          dqRNGkind("pcg64")
          SEED <-
            2121870L # sample(1e7, 1) # Hard-coded for reproducibility
          set.seed(SEED + mc_)
          dqset.seed(SEED, mc_)

          dt <- private$gen_synthpop_demog(design_, lsoas_)

          # Calculate local qimd (lqimd) ----
          message("Calculate local qimd (lqimd)")
          if ("England" %in% design_$sim_prm$locality) {
            dt[, lqimd := qimd]
          } else {
            dt[, lqimd := cut(
              imd,
              breaks = quantile(imd, probs = 0:5 / 5, type = 5),
              labels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
              include.lowest = TRUE,
              right = TRUE
            )]
            CCG17CDH_ <- as.character(unique(dt$CCG17CDH))
            # names(sort(table(as.character(dt$CCG17CDH)))[1]) # if more than one ccgs
          }

          # Generate the cohorts of 30 year old to enter the simulation every year ----
          message("Generate the cohorts of 30 year old")

          # new
          tt <- dt[age == design_$sim_prm$ageL, .N]
          design2_ <- design_$clone()
          design2_$sim_prm$n <-
            tt * design_$sim_prm$sim_horizon_max
          design2_$sim_prm$ageH <- design_$sim_prm$ageL
          dt2 <- private$gen_synthpop_demog(design2_, lsoas_)
          dt2[, age := age - rep(1:design_$sim_prm$sim_horizon_max, tt)]

          # as sim progress these will become 30 yo
          # no population growth here as I will calibrate to dt projections
          if ("England" %in% design_$sim_prm$locality) {
            dt2[, lqimd := qimd]

          } else {
            dt2[, lqimd := cut(
              imd,
              breaks = quantile(dt$imd, probs = 0:5 / 5, type = 5),
              labels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
              include.lowest = TRUE,
              right = TRUE
            )]
            # names(sort(table(as.character(dt$CCG17CDH)))[1]) # if more than one ccgs
          }
          dt <- rbind(dt2, dt)
          rm(dt2, design2_)

          # NOTE!! from now on year in the short form i.e. 13 not 2013
          dt[, `:=`(year = year - 2000L,
                    pid  = .I)]
          new_n <- nrow(dt)


          # Generate correlated ranks for the individuals ----
          message("Generate correlated ranks for the individuals")

          cm_mean <- as.matrix(
            read_fst(
              "./lifecourse_models/exposure_corr_mean.fst",
              as.data.table = TRUE
            ),
            rownames = "rn"
          )

          rank_mtx <- generate_corr_unifs(new_n, cm_mean)
          message("generate correlated uniforms")

          # Restrict the range of some RNs to avoid unrealistic exposures
          # This scaling does not affect correlations
          # /0.999 because I multiplied all the columns below
          rank_mtx <- rank_mtx * 0.999
          rank_mtx[, "frtpor_r"] <- rank_mtx[, "frtpor_r"] * 0.99 / 0.999
          rank_mtx[, "vegpor_r"] <- rank_mtx[, "vegpor_r"] * 0.93 / 0.999
          rank_mtx[, "smok_cig_ex_r"] <-
            rank_mtx[, "smok_cig_ex_r"] * 0.99 / 0.999
          rank_mtx[, "totalwu_r"] <- rank_mtx[, "totalwu_r"] * 0.99 / 0.999
          rank_mtx[, "smok_quit_yrs_r"] <-
            rank_mtx[, "smok_quit_yrs_r"] * 0.99 / 0.999
          rank_mtx[, "smok_dur_ex_r"] <-
            rank_mtx[, "smok_dur_ex_r"] * 0.99 / 0.999
          rank_mtx[, "smok_dur_curr_r"] <-
            rank_mtx[, "smok_dur_curr_r"] * 0.88 / 0.999

          # sum((cor(rank_mtx) - cm_mean) ^ 2)
          message("correlated ranks matrix to data.table")

          rank_mtx <- data.table(rank_mtx)

          # NOTE rankstat_* is unaffected by the RW. Stay constant through the lifecourse
          dt[, c(
            "rank_education",
            "rank_income",
            "rank_pa",
            "rank_fruit",
            "rank_veg",
            "rankstat_smok",
            "rankstat_smok_quit_yrs",
            "rankstat_smok_dur_ex",
            "rankstat_smok_dur_curr",
            "rankstat_smok_cig_ex",
            "rankstat_smok_cig_curr",
            "rank_ets",
            "rank_alcohol",
            "rank_bmi",
            "rank_sbp",
            "rank_bpmed",
            "rank_tchol",
            "rank_hdl",
            "rankstat_af_dgn",
            "rankstat_famcvd",
            "rank_t2dm",
            "rankstat_t2dm_dgn",
            "rank_statin_px",
            "rankstat_ckd"
          ) := rank_mtx]

          rm(rank_mtx)

          # add non-correlated RNs
          rank_cols <-
            c(
              "rankstat_ra",
              "rank_cst",
              "rankstat_ncc",
              "rankstat_ca_history",
              "rankstat_famlungca"
            )


          for (nam in rank_cols)
            set(dt, NULL, nam, dqrunif(new_n)) # NOTE do not replace with generate_rns function.

          # Generate education (exception as it remains stable through lifecourse) ----
          message("Generate education")

          tbl <-
            read_fst("./lifecourse_models/education_table.fst",
                     as.data.table = TRUE)
          nam <- intersect(names(dt), names(tbl))
          # logic necessary for new cohorts entering the simulation that currently age < 30
          # These will have the same distribution as if 30 years old
          tt <- tbl[age == min(age)]
          tt <- clone_dt(tt, design_$sim_prm$sim_horizon_max)
          tt[, age := age - .id] # as the sim progress these will become 30 yo
          # increase population by 0.5% every year
          # TODO extract dt increase by LAD from ONS
          tt[, .id := NULL]
          tbl <- rbind(tt, tbl)
          dt[tbl, education := (rank_education > ed1) + (rank_education > ed2) +
               (rank_education > ed3) + (rank_education > ed4) +
               (rank_education > ed5) + (rank_education > ed6) + 1L,
             on = nam]
          dt[, education := factor(
            education,
            levels = 1:7,
            labels = c(
              "NVQ4/NVQ5/Degree or equiv",
              "Higher ed below degree",
              "NVQ3/GCE A Level equiv",
              "NVQ2/GCE O Level equiv",
              "NVQ1/CSE other grade equiv",
              "Foreign/other",
              "No qualification"
            )
          )]
          dt[, rank_education := NULL]

          # Project forward for simulation and back project for lags  ----
          message("Project forward and back project")

          dt <-
            clone_dt(dt,
                     design_$sim_prm$sim_horizon_max +
                       design_$sim_prm$maxlag + 1L)

          dt[.id <= design_$sim_prm$maxlag, `:=` (age  = age  - .id,
                                                  year = year - .id)]
          dt[.id > design_$sim_prm$maxlag, `:=` (
            age  = age  + .id - design_$sim_prm$maxlag - 1L,
            year = year + .id - design_$sim_prm$maxlag - 1L
          )]
          # dt <-
          #   dt[between(age, design_$sim_prm$ageL - design_$sim_prm$maxlag, design_$sim_prm$ageH)]
          # delete unnecessary ages
          del_dt_rows(
            dt,
            !between(
              dt$age,
              design_$sim_prm$ageL - design_$sim_prm$maxlag,
              design_$sim_prm$ageH
            ),
            environment()
          )

          dt[, `:=` (.id = NULL)]

          # to_agegrp(dt, 20L, 85L, "age", "agegrp20", to_factor = TRUE)
          # to_agegrp(dt, 10L, 85L, "age", "agegrp10", to_factor = TRUE)
          # to_agegrp(dt,  5L, 85L, "age", "agegrp5" , to_factor = TRUE)

          # generate population weights
          private$calc_pop_weights(dt, design_$sim_prm$locality)

          # Simulate exposures -----

          # Random walk for ranks ----
          message("Random walk for ranks")

          setkeyv(dt, c("pid", "year"))
          setindexv(dt, c("year", "age", "sex", "sha", "qimd", "ethnicity"))

          dt[, pid_mrk := mk_new_simulant_markers(pid)]

          dt[, lapply(.SD,
                      fscramble_trajectories,
                      pid_mrk,
                      design_$sim_prm$jumpiness),
             .SDcols = patterns("^rank_")]
          # ggplot2::qplot(year, rank_income, data = dt[pid %in% sample(1e5, 1)], ylim = c(0,1))

          # Generate income ----
          message("Generate income")

          tbl <-
            read_fst("./lifecourse_models/income_table.fst", as.data.table = TRUE)
          nam <- intersect(names(dt), names(tbl))
          dt[tbl, income := (rank_income > inc1) + (rank_income > inc2) +
               (rank_income > inc3) + (rank_income > inc4) + 1L,
             on = nam]
          dt[, income := factor(
            income,
            levels = 1:5,
            labels = c("1 Highest", "2", "3", "4", "5 Lowest")
          )]
          dt[, rank_income := NULL]

          # Generate active days ----
          message("Generate active days")

          tbl <-
            read_fst("./lifecourse_models/active_days_table.fst",
                     as.data.table = TRUE)
          nam <- intersect(names(dt), names(tbl))
          dt[tbl, active_days := (rank_pa > pa0) + (rank_pa > pa1) + (rank_pa > pa2) +
               (rank_pa > pa3) + (rank_pa > pa4) + (rank_pa > pa5) + (rank_pa > pa6),
             on = nam]
          dt[, rank_pa := NULL]
          # dt[, active_days := factor(active_days, levels = 0:7)]

          # Generate fruit consumption (ZISICHEL) ----
          message("Generate fruit consumption")

          tbl <-
            read_fst("./lifecourse_models/frtpor_table.fst", as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, fruit :=
               my_qZISICHEL(rank_fruit,
                            mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu) * 80L]  # g/d
          dt[, (col_nam) := NULL]
          dt[, rank_fruit := NULL]

          # Generate veg consumption (DEL) ----
          message("Generate veg consumption")

          tbl <-
            read_fst("./lifecourse_models/vegpor_table.fst", as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, veg :=
               my_qDEL(rank_veg, mu, sigma, nu, n_cpu = design_$sim_prm$n_cpu) * 80L]  # g/d
          dt[, (col_nam) := NULL]
          dt[, rank_veg := NULL]

          # Smoking simulation ----
          message("Smoking simulation")

          # Assign smok_status when pid_mrk == true (the first year an individul enters the simulation (with lags))
          tbl <-
            read_fst("./lifecourse_models/smok_status_table.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, smok_status_ref := my_qMN4(rankstat_smok, mu, sigma, nu)] # for calibration
          dt[(pid_mrk), smok_status := smok_status_ref]
          dt[, (col_nam) := NULL]
          dt[, rankstat_smok := dqrunif(.N)] # this is now used for simsmoke. There shouldn't be correlated any more (colname hardcoded in C++ code).

          # Assign smok_quit_yrs when pid_mrk == true (the first year an individual enters the simulation)
          # I could use these estimates for calibration but I need to calculate mortality first
          tbl <-
            read_fst("./lifecourse_models/smok_quit_yrs_table.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          set(dt, NULL, "smok_quit_yrs", 0L)
          dt[(pid_mrk) &
               smok_status %in% 2:3,
             smok_quit_yrs := my_qDPO(rankstat_smok_quit_yrs, mu, sigma)]
          dt[, rankstat_smok_quit_yrs := NULL]
          dt[, (col_nam) := NULL]

          # Assign smok_dur_ex when pid_mrk == true (the first year an individul enters the simulation)
          tbl <-
            read_fst("./lifecourse_models/smok_dur_ex_table.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          set(dt, NULL, "smok_dur", 0L)
          dt[(pid_mrk) &
               smok_status %in% 2:3, smok_dur := my_qDPO(rankstat_smok_dur_ex, mu, sigma)]
          dt[, rankstat_smok_dur_ex := NULL]
          dt[, (col_nam) := NULL]

          # Assign smok_dur_curr when pid_mrk == true (the first year an individual enters the simulation)
          tbl <-
            read_fst("./lifecourse_models/smok_dur_curr_table.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[(pid_mrk) &
               smok_status == 4, smok_dur := as.integer(round(qNBI(rankstat_smok_dur_curr, mu, sigma)))]
          dt[, rankstat_smok_dur_curr := NULL]
          dt[, (col_nam) := NULL]

          # Ensure smoking histories start from age 12
          dt[age - smok_quit_yrs < 12L, smok_quit_yrs := age - 12L]
          dt[age - smok_dur < 12L, smok_dur := age - 12L]
          dt[age - smok_dur - smok_quit_yrs < 12L ,
             `:=`(smok_dur = as.integer(smok_dur / ((
               smok_dur + smok_quit_yrs
             ) / (age - 12L))),
             smok_quit_yrs = as.integer(smok_quit_yrs / ((
               smok_dur + smok_quit_yrs
             ) / (age - 12L))))]

          # Assign smok_incid probabilities
          tbl <-
            read_fst("./lifecourse_models/smok_incid_table.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          setnames(dt, "mu", "prb_smok_incid")

          # Assign smok_cessation probabilities
          tbl <-
            read_fst("./lifecourse_models/smok_cess_table.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          setnames(dt, "mu", "prb_smok_cess")

          # Handle smok_relapse probabilities
          tbl <-
            read_fst("./lifecourse_models/smok_relapse_table.fst",
                     as.data.table = TRUE)
          tbl <-
            dcast(tbl, sex + qimd ~ smok_quit_yrs, value.var = "pr")
          nam <- tbl[, paste0(sex, " ", qimd)]
          tbl <-
            as.matrix(tbl[, mget(paste0(1:15))], rownames = nam)

          simsmok(dt, tbl, design_$sim_prm$smoking_relapse_limit)
          # dt[!(pid_mrk), table(smok_status)]
          # dt[pid == 1, plot(year, smok_status, ylim = c(0, 4))]
          # dt[pid == 10, .(age, smok_status, smok_quit_yrs, smok_dur)]
          # dt[, sum(smok_status == 4)/.N, keyby = year]

          if (design_$sim_prm$simsmok_calibration) {
            # calculate dif between ref (multinom) and simsmok
            # I will further calibrate to better match HSE
            resample <-
              function(x, ...)
                x[sample.int(length(x), ...)]
            obs <-
              dt[smok_status == 1L, .(nsa = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
            ref <-
              dt[smok_status_ref == 1L, .(nsr = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
            absorb_dt(ref, obs)
            setnafill(ref, "c", 0L, cols = "nsa")
            ref[, `:=`(dif = nsr - nsa,
                       nsr = NULL,
                       nsa = NULL)]
            # Further calibrate to match better with HSE
            ref[sex == "men" &
                  rbinom(.N, 1L, 0.9) == 1L, dif := dif - 2L]
            ref[sex == "men" &
                  rbinom(.N, 1L, 1 * clamp((year - min(year)) / 10, 0, 1)) == 1L, dif := dif + 2L]
            ref[sex == "women" &
                  rbinom(.N, 1L, 0.2) == 1L, dif := dif - 1L]
            ref[age < 49 &
                  rbinom(.N, 1L, 0.5) == 1L, dif := dif + 1L]
            ref[age < 49 &
                  sex == "men" &
                  qimd == "3" &
                  rbinom(.N, 1L, 0.5) == 1L, dif := dif - 1L]
            ref[age < 49 &
                  sex == "women" &
                  qimd %in% c("4", "5 least deprived") &
                  rbinom(.N, 1L, 0.4) == 1L, dif := dif - 1L]
            ref[between(age, 50, 69) &
                  rbinom(.N, 1L, 0.2) == 1L, dif := dif + 1L]
            ref[between(age, 50, 69) &
                  sex == "women" &
                  qimd == "5 least deprived" &
                  rbinom(.N, 1L, 0.4) == 1L, dif := dif + 1L]
            ref[between(age, 70, 89) &
                  rbinom(.N, 1L, 0.4) == 1L, dif := dif + 1L]
            ref[between(age, 70, 89) &
                  sex == "men" &
                  qimd %in% c("1 most deprived", "2") &
                  rbinom(.N, 1L, 0.4) == 1L, dif := dif - 1L]
            ref[between(age, 70, 89) &
                  sex == "women" &
                  qimd %in% c("4d", "2") &
                  rbinom(.N, 1L, 0.4) == 1L, dif := dif + 1L]
            absorb_dt(dt, ref)

            # when not enough never smokers convert those ex smokers with the longer quit years
            tt <-
              dt[smok_status %in% 2:3, .(year, age, sex, qimd, pid, smok_quit_yrs, dif)]
            setnafill(tt, "c", 0, cols = "dif")
            tt[dif < 0, dif := 0L]
            setkey(tt, year, age, sex, qimd, smok_quit_yrs)
            pid_to_conv <-
              tt[dif > 0, .(pid = tail(pid, max(dif))), keyby = .(year, age, sex, qimd)]
            dt[pid_to_conv, on = .(year, pid), `:=`(
              smok_status = 1L,
              smok_quit_yrs = 0L,
              smok_dur = 0L,
              smok_cig = 0L
            )]

            # when too many never smokers convert to smok status 2 (occasional)
            tt <-
              dt[smok_status %in% 1, .(year, age, sex, qimd, pid, dif)]
            setnafill(tt, "c", 0, cols = "dif")
            tt[dif > 0, dif := 0L]
            tt[, dif := -dif]
            setkey(tt, year, age, sex, qimd)
            # Ensure there are enough people to sample from
            ttt <-
              tt[, .(lpid = length(pid), mdif =  max(dif)), by = .(year, age, sex, qimd)][mdif >
                                                                                            lpid, ]
            tt[ttt, on = .NATURAL, dif := i.lpid]
            pid_to_conv <-
              tt[dif > 0, .(pid = resample(pid, max(dif))), keyby = .(year, age, sex, qimd)]
            dt[pid_to_conv, on = .(year, pid), `:=`(
              smok_status = 2L,
              smok_quit_yrs = dt[smok_status == 2L, median(smok_quit_yrs)],
              smok_dur = dt[smok_status == 2L, median(smok_dur)],
              smok_cig = 1L
            )]
            dt[, dif := NULL]

            # Same logic for active smokers
            # I will further calibrate to better match HSE
            # calculate dif between ref (multinom) and simsmok
            obs <-
              dt[smok_status == 4L, .(nsa = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
            ref <-
              dt[smok_status_ref == 4L, .(nsr = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
            absorb_dt(ref, obs)
            setnafill(ref, "c", 0, cols = "nsa")
            ref[, `:=`(dif = nsr - nsa,
                       nsr = NULL,
                       nsa = NULL)]
            # Further calibrate to match better with HSE
            ref[sex == "men" &
                  rbinom(.N, 1L, 0.3) == 1L, dif := dif - 1L]
            ref[sex == "women" &
                  qimd != "1 most deprived" &
                  qimd != "5 least deprived" &
                  rbinom(.N, 1L, 0.8) == 1L, dif := dif - 1L]
            ref[qimd == "1 most deprived" &
                  rbinom(.N, 1L, 0.5) == 1L, dif := dif + 1L]
            ref[qimd == "5 least deprived" &
                  rbinom(.N, 1L, 0.6) == 1L, dif := dif - 1L] # - reduces
            absorb_dt(dt, ref)

            # when not enough active smokers convert those ex smokers with the shortest quit years
            tt <-
              dt[smok_status == 3L, .(year, age, sex, qimd, pid, smok_quit_yrs, dif)]
            setnafill(tt, "c", 0, cols = "dif")
            tt[dif < 0, dif := 0L]
            setkey(tt, year, age, sex, qimd, smok_quit_yrs)
            pid_to_conv <-
              tt[dif > 0, .(pid = head(pid, max(dif))), keyby = .(year, age, sex, qimd)]
            dt[pid_to_conv, on = .(year, pid),
               `:=`(
                 smok_status = 4L,
                 smok_quit_yrs = 0L,
                 smok_dur = smok_dur + 1
               )] # TODO fix smoking duration

            # when too many never smokers convert to smok status 3
            tt <-
              dt[smok_status == 4L, .(year, age, sex, qimd, pid, dif)]
            setnafill(tt, "c", 0, cols = "dif")
            tt[dif > 0, dif := 0L]
            tt[, dif := -dif]
            setkey(tt, year, age, sex, qimd)
            # Ensure there are enough people to sample from
            ttt <-
              tt[, .(lpid = length(pid), mdif =  max(dif)), by = .(year, age, sex, qimd)][mdif >
                                                                                            lpid,]
            tt[ttt, on = .NATURAL, dif := i.lpid]
            pid_to_conv <-
              tt[dif > 0, .(pid = resample(pid, max(dif))), keyby = .(year, age, sex, qimd)]
            dt[pid_to_conv, on = .(year, pid), `:=`(smok_status = 3L, smok_quit_yrs = 1L)]
            dt[, dif := NULL]

            rm(tt, ttt, obs, ref, pid_to_conv)
          }

          # Assign smok_cig_curr when pid_mrk == true (the first year an individual enters the simulation)
          set(dt, NULL, "smok_cig", 0L)

          tbl <-
            read_fst("./lifecourse_models/smok_cig_curr_table.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[smok_status == 4L,
             smok_cig := qZINBI(rankstat_smok_cig_curr, mu, sigma, nu)]
          dt[, (col_nam) := NULL]

          # Assign smok_cig_ex when pid_mrk == true (the first year an individual enters the simulation)
          # dt[smok_status == 2, smok_cig := 1L]


          tbl <-
            read_fst("./lifecourse_models/smok_cig_ex_table.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[(pid_mrk) &
               smok_status == 3L,
             smok_cig := my_qZABNB(rankstat_smok_cig_ex,
                                   mu,
                                   sigma,
                                   nu,
                                   tau,
                                   n_cpu = design_$sim_prm$n_cpu)]
          dt[, (col_nam) := NULL]

          simsmok_cig(dt) # carry forward smok_cig if smok_status == 3
          dt[smok_cig == 0L & smok_status > 1L, smok_cig := 1L]

          if (design_$sim_prm$simsmok_calibration)
            simsmok_postcalibration(dt) # need to be post cig simulation

          dt[, smok_status := factor(smok_status)]
          # needed for QRisk and QDrisk
          dt[, smoke_cat := 0L]
          dt[smok_status == "3", smoke_cat := 1L]
          dt[smok_status == "4", smoke_cat := 3L]
          dt[smok_status == "4" &
               smok_cig < 10, smoke_cat := 2L]
          dt[smok_status == "4" &
               smok_cig > 19, smoke_cat := 4L]

          dt[, c(
            "rankstat_smok",
            "rankstat_smok_cig_curr",
            "rankstat_smok_cig_ex",
            "prb_smok_incid",
            "prb_smok_cess",
            "smok_status_ref"
          ) := NULL]

          # Generate ETS (BI) ----
          message("Generate ETS")

          # Note at the moment this is independent of smoking prevalence
          # TODO calculate how many each smoker pollutes by year, SHA (not qimd) to be used in
          # scenarios. Ideally correct for mortality
          tbl <-
            read_fst("./lifecourse_models/ets_table.fst", as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, ets := as.integer(rank_ets < mu)]
          dt[, rank_ets := NULL]
          dt[, (col_nam) := NULL]
          # View(dt[, prop_if(ets == 1)/prop_if(smok_status == "4"), keyby = .(year, sha)])

          # Generate alcohol (ZINBI) ----
          message("Generate alcohol")

          tbl <-
            read_fst("./lifecourse_models/alcohol_table.fst", as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, alcohol := qZINBI(rank_alcohol, mu, sigma, nu)]
          dt[, rank_alcohol := NULL]
          dt[, (col_nam) := NULL]

          # Generate BMI (BCPEo) ----
          message("Generate BMI")

          tbl <-
            read_fst("./lifecourse_models/bmi_table.fst", as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, bmi := my_qBCPEo(rank_bmi, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
          dt[, rank_bmi := NULL]
          dt[, (col_nam) := NULL]

          # Generate SBP (BCPEo) ----
          message("Generate SBP")

          tbl <-
            read_fst("./lifecourse_models/sbp_table.fst", as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, sbp := my_qBCPEo(rank_sbp, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
          dt[, rank_sbp := NULL]
          dt[, (col_nam) := NULL]

          # Generate BP medication (BI) -----
          message("Generate BP medication")

          dt[, `:=` (sbp_acc = sbp,
                     sbp = round(clamp(sbp, 110, 200), -1))]
          tbl <-
            read_fst("./lifecourse_models/bp_med_table.fst", as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, bpmed := as.integer(rank_bpmed < mu)]
          dt[, rank_bpmed := NULL]
          dt[, (col_nam) := NULL]
          dt[, `:=` (sbp = sbp_acc,
                     sbp_acc = NULL)]

          # TODO calculate probability of dgn HTN

          # Generate tchol (BCT) ----
          message("Generate tchol")

          tbl <-
            read_fst("./lifecourse_models/tchol_table.fst", as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, tchol := my_qBCT(rank_tchol, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
          dt[, rank_tchol := NULL]
          dt[, (col_nam) := NULL]

          # Generate HDL (to tchol ratio) (GB1) ----
          message("Generate HDL (to tchol ratio)")

          # NOTE this very highly correlated with hdl level (~0.76) and
          #  highly to tchol (~-0.47). The latter is captured by the correlated RNs
          tbl <-
            read_fst("./lifecourse_models/hdl_to_tchol_table.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, tchol_hdl_ratio := 1 / qGB1(rank_hdl, mu, sigma, nu, tau)]
          dt[, rank_hdl := NULL]
          dt[, (col_nam) := NULL]

          # Generate statins medication (BI) -----
          message("Generate statins medication")

          dt[, `:=` (tchol_acc = tchol,
                     tchol = round(clamp(tchol, 2, 12), 0))]
          tbl <-
            read_fst("./lifecourse_models/statin_px_table.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, statin_px := as.integer(rank_statin_px < mu)]
          dt[, rank_statin_px := NULL]
          dt[, (col_nam) := NULL]
          dt[, `:=` (tchol = tchol_acc,
                     tchol_acc = NULL)]


          # Generate AF dgn (BI) ----
          message("Generate AF dgn")

          # NOTE AF prevalence ~6.9%. This is higher than the QOF prevalence of 1.6% for 2019 and 1.7% for 2015/6
          # although PHE includes all ages and my estimates only ~30-89
          tbl <-
            read_fst("./lifecourse_models/af_dgn_table.fst", as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, af_dgn := as.integer(rankstat_af_dgn < mu)]

          # calibration to match QOF assuming characteristics similar to HSE
          if ("England" %in% design_$sim_prm$locality) {
            tt <-
              0.0181 / dt[year == 15L, sum(af_dgn) / .N] # to account for difference in age between my population and PHE
            dt[, af_dgn := as.integer(rankstat_af_dgn < (tt * mu))]
          } else {
            ttt <-
              fread("./lifecourse_models/undiagnosed_af_estimates_byCCG.csv")
            dt[ttt, `:=` (QOF_prevalence = i.QOF_prevalence,
                          expected_prevalence = i.expected_prevalence), on = "CCG17CDH"]
            rm(ttt)
            setnafill(dt, "c", 0.0181, cols = "QOF_prevalence") # For missing CCG assume national average
            tt <-
              dt[year == 15L, sum(af_dgn) / .N, by = CCG17CDH]
            dt[tt, af_pct := i.V1, on = "CCG17CDH"]
            # +0.001 to account for difference in age between my population and PHE
            dt[, af_dgn := as.integer(rankstat_af_dgn < (mu * (0.001 + QOF_prevalence / 100) / af_pct))]
            dt[, `:=` (QOF_prevalence = NULL, af_pct = NULL)]
          }

          # Generate AF undgn + diagn (BI) ----
          message("Generate AF undgn + diagn")

          # TODO alt method for Turakhia et al. 2018 Table 3
          # From PHE that estimate .44 undiagnosed AF cases for every diagnosed case
          dt[, af_prvl := as.integer(rankstat_af_dgn < mu)]

          # calibration to match PHE undiagnosed+diagnosed prevalence
          if ("England" %in% design_$sim_prm$locality) {
            tt <-
              0.0255 / dt[year == 15L, sum(af_prvl) / .N] # + 0.001 to account for difference in age between my population and PHE
            dt[, af_prvl := as.integer(rankstat_af_dgn < (tt * mu))]
          } else {
            setnafill(dt, "c", 0.0255, cols = "expected_prevalence") # For missing CCG assume national average
            tt <-
              dt[year == 15L, sum(af_prvl) / .N, by = CCG17CDH]
            dt[tt, af_pct := i.V1, on = "CCG17CDH"]
            # +0.001 to account for difference in age between my population and PHE
            dt[, af_prvl := as.integer(rankstat_af_dgn < (mu * (0.001 + expected_prevalence / 100) / af_pct))]
            dt[, `:=` (expected_prevalence = NULL, af_pct = NULL)]
          }
          dt[, af_prvl := af_prvl * 2L]# to avoid 1 in prevalence and avoid confusion with incidence
          dt[, rankstat_af_dgn := NULL]
          dt[, (col_nam) := NULL]

          # Generate family CVD dgn (BI) ----
          message("Generate family CVD dgn")

          tbl <-
            read_fst("./lifecourse_models/famcvd_table.fst", as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, famcvd := as.integer(rankstat_famcvd < mu)]
          dt[, rankstat_famcvd := NULL]
          dt[, (col_nam) := NULL]


          # Generate CKD 4/5 prevalence -----------
          message("Generate CKD 4/5 prevalence")

          tbl <-
            read_fst("./lifecourse_models/ckd_table.fst", as.data.table = TRUE)
          nam <- intersect(names(dt), names(tbl))
          dt[tbl, ckd_prvl := (rankstat_ckd > ckd0) + (rankstat_ckd > ckd1) +
               (rankstat_ckd > ckd2) + (rankstat_ckd > ckd3) + 1L,
             on = nam]
          dt[, ckd5_prvl := fifelse(ckd_prvl == 5L, 1L, 0L)]
          dt[, rankstat_ckd := NULL]


          # Generate rheumatoid arthritis prevalence -------
          message("Generate rheumatoid arthritis prevalence")

          # FROM Symmons et al 2002 table 1
          RA_prvl <- fread(
            "./lifecourse_models/RA_prvl.csv",
            stringsAsFactors = FALSE,
            colClasses = c("integer", "factor",
                           "numeric", "numeric", "numeric"),
            key = c("sex", "age")
          )
          tt <- dqrunif(1L)
          RA_prvl[, prvl := mc2d::qpert(tt, lci, prvl, uci)]
          RA_prvl[, prvl_smth := predict(loess(prvl ~ age, degree = 1, span = 0.4)), by = sex]
          # RA_prvl[sex == "1", plot(age, prvl_smth)]
          # RA_prvl[sex == "1", lines(age, prvl)]
          # RA_prvl[, .(prvl = MESS::auc(age, prvl, 30, 89, type = "spline"),
          #              prvl_smth = MESS::auc(age, prvl_smth, 30, 89, type = "spline")),
          #          keyby = sex] # AUC almost identical
          RA_prvl <-
            RA_prvl[between(age, min(dt$age), max(dt$age)), .(age, sex = factor(sex,
                                                                                levels = 1:2,
                                                                                labels = c("men", "women")), prvl_smth)]
          dt[RA_prvl, on = c("age", "sex"), ra_prvl := as.integer(rankstat_ra < prvl_smth)]
          dt[, rankstat_ra := NULL]
          rm(RA_prvl)

          # Generate corticosteroids ------
          message("Generate corticosteroids")

          # Assuming cst independent from RA (although both covary with age).
          cst_prvl <- fread(
            "./lifecourse_models/corticosteroids_prvl.csv",
            stringsAsFactors = TRUE,
            colClasses = c("factor", "integer",
                           "numeric"),
            key = c("sex", "age")
          )
          tt <- dqrunif(1L)
          cst_prvl[, prvl := mc2d::qpert(tt, prvl * 0.8, prvl, prvl * 1.2)] # Assume some uncertainty
          dt[cst_prvl, on = c("age", "sex"), cst_prvl := as.integer(rank_cst < prvl)]
          dt[, rank_cst := NULL]
          rm(cst_prvl)

          # Generate T2DM (BI) -----
          message("Generate T2DM")

          # Both dgn and undgn (move to apply on initial year only)
          dt[, `:=` (bmi_acc = bmi,
                     bmi = round(clamp(bmi, 18, 50)))]
          tbl <-
            read_fst("./lifecourse_models/dm_table.fst", as.data.table = TRUE)
          tbl[, year := design_$sim_prm$init_year] # only estimate prevalence for 2013
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, t2dm_prvl := as.integer(rank_t2dm < mu)]
          setnafill(dt, "c", 0L, cols = "t2dm_prvl")
          dt[, `:=` (rank_t2dm = NULL)]
          dt[, (col_nam) := NULL]

          # Generate probability of dgn T2DM (BI) -----
          message("Generate probability of dgn T2DM")

          dt[, `:=` (bmi = round(clamp(bmi, 18, 50), -1))]
          tbl <-
            read_fst("./lifecourse_models/dm_dgn_table.fst", as.data.table = TRUE)
          tbl[, t2dm_prvl := 1L]
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, t2dm_dgn := as.integer(rankstat_t2dm_dgn < mu)]
          dt[, (col_nam) := NULL]


          # t2dm duration (GPO)
          tbl <-
            read_fst("./lifecourse_models/dm_dur_table.fst", as.data.table = TRUE)
          tbl[, t2dm_prvl := 1L]
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[t2dm_prvl == 1L, t2dm_prvl := 2L + rpois(.N, 3L) + qGPO(dqrunif(.N), mu, sigma)] # +2 to avoid confussion with incd
          # rpois(.N, 3L) to assume 3 year mean period from onset till
          # diagnosis because the model was fitted in diagnosed patients
          dt[t2dm_prvl > age, t2dm_prvl := age]
          dt[, (col_nam) := NULL]
          dt[, t2dm_prvl := carry_backward(t2dm_prvl, pid_mrk)]
          # dt[pid == 1, .(year, t2dm_prvl)]

          # use qrisk diabetes for probability of t2dm incidence

          # Generate family history of diabetes ------------
          message("Generate family history of T2DM")

          # If the prob of being diab is ~8%. So the probability of having at least
          # one of 3 family members with diabetes is 1 - (1-0.08)^3.
          # I let family members vary between 2 and 2+rpois(n, 1)
          # This does not account for the future increase of diabetes prevalence.
          # Therefore it underestimates. I assume diabetes prevalence will increase
          # by 2% every year to adjust for that
          tt <-
            dt[year == design_$sim_prm$init_year &
                 age > 40, prop_if(t2dm_prvl > 0L)] # t2dm prevalence
          dt[, fam_t2dm :=
               rbinom(.N, 1, 1 - (1 - ((
                 1.02 ^ (year - design_$sim_prm$init_year)
               ) * tt)) ^ (2 + rpois(.N, 1)))]



          # Estimate number of comorbidities (ncc) calculation ----
          # to be used in QALY
          message("Generate ncc")

          dt[, ncc := clamp(qbinom(rankstat_ncc, ceiling(age / 8L), fifelse(age < 55, 0.25, 0.40)), 0, 10)]
          # calibrated to Sullivan et all 2011 (web table 1)
          # to_agegrp(output, 10L, 89L, "age", "agegrp10", to_factor = TRUE)
          # output[, round(mean(ncc), 1), keyby = agegrp10]
          # target by agegrp 1.1  1.6  2.4  3.1  4.0  4.4 from

          dt[, `:=` (
            bmi = bmi_acc,
            bmi_acc = NULL,
            rankstat_t2dm_dgn = NULL,
            # fam_t2dm = NULL,
            # cst_prvl = NULL,
            pid_mrk = NULL,
            # to be recreated when loading synthpop
            rankstat_ncc = NULL
          )]

          dt[, c("prb_t2dm_incd_nocvd", "t2dm_incd_cvd_mltp") := QDiabetes(.SD)]
          dt[, t2dm_incd_cvd_mltp := t2dm_incd_cvd_mltp / prb_t2dm_incd_nocvd]


          # Estimate history of cancer
          tbl <-
            read_fst("./lifecourse_models/history_of_cancer.fst",
                     as.data.table = TRUE)
          col_nam <-
            setdiff(names(tbl), intersect(names(dt), names(tbl)))
          absorb_dt(dt, tbl)
          dt[, history_of_ca := as.integer(rankstat_ca_history < mu)]
          dt[, (col_nam) := NULL]

          # Estimate family history of lung ca (crude approx)
          # The prob of having lung ca is around 50/1e5.
          # So the probability of having at least
          # one of 3 family members with lung ca is 1 - (1-50/1e5)^3.
          tt <- 1 - (1 - 50 / 1e5) ^ 3
          dt[, fam_lung_ca := as.integer(rankstat_famlungca < tt)]

          exps_tolag <- c(
            "active_days",
            "fruit",
            "veg",
            "smok_status",
            "smok_quit_yrs",
            "smok_dur",
            "smok_cig",
            "ets",
            "alcohol",
            "af_dgn",
            "af_prvl",
            "ckd_prvl",
            "bmi",
            "sbp",
            "bpmed",
            "tchol",
            "statin_px",
            "t2dm_prvl"
          )
          exps_nam <-  paste0(exps_tolag, "_curr_xps")
          setnames(dt, exps_tolag, exps_nam)

          # Ensure pid does not overlap for files from different mc
          if (max(dt$pid + (mc_ - 1L) * new_n) < .Machine$integer.max) {
            dt[, pid := pid + (mc_ - 1L) * new_n]
          }

          # Include diseases ----
          message("Include diseases")
          setkey(dt, pid, year)
          # I need to delete rows here so the RN with the file are
          # reproducible and don't need to save them
          # TODO lag of 10 crashes shift_byID
          lags_mc <-
            get_lag_mc(mc_, design_$sim_prm)
          max_lag_mc <- max(unlist(lags_mc))
          dt <-
            dt[year >= (design_$sim_prm$init_year - max_lag_mc) &
                 between(dt$age, design_$sim_prm$ageL, design_$sim_prm$ageH)]

          setkey(dt, pid, year)
          dt[, pid_mrk := mk_new_simulant_markers(pid)]

          finalise_synthpop(mc_, dt, design_$sim_prm, lags_mc)

          generate_rns(
            mc_,
            dt,
            c(
              "rn_af_dgn"                  ,  "rn_breast_ca_dgn"          ,
              "rn_breast_ca_incd"          ,  "rn_breast_ca_mrtl"         ,
              "rn_chd_dgn"                 ,  "rn_chd_incd"               ,
              "rn_chd_mrtl"                ,  "rn_colon_ca_dgn"           ,
              "rn_colon_ca_incd"           ,  "rn_colon_ca_mrtl"          ,
              "rn_copd_dgn"                ,  "rn_copd_incd"              ,
              "rn_copd_mrtl"               ,  "rn_htn_dgn"                ,
              "rn_lung_ca_dgn"             ,  "rn_lung_ca_incd"           ,
              "rn_lung_ca_mrtl"            ,  "rn_multi_mrtl"             ,
              "rn_nonmodelled_mrtl"        ,  "rn_poststroke_dementia_dgn",
              "rn_poststroke_dementia_incd",  "rn_stroke_dgn"             ,
              "rn_stroke_incd"             ,  "rn_stroke_mrtl"            ,
              "rn_t2dm_dgn"                ,  "rn_t2dm_incd"
            )
          )

          # Run to estimate prevalence ignoring mortality
          output <- list()
          output <-
            gen_output("", design_$sim_prm, lags_mc, dt, output, TRUE)
          output <-
            rbindlist(output)[between(age, design_$sim_prm$ageL,
                                      design_$sim_prm$ageH) &
                                year >= design_$sim_prm$init_year]

          # Estimate fatality rates for all diseases
          .diseases <- c("chd", "stroke", "copd",
                         "lung_ca", "colon_ca", "breast_ca")
          for (disease_nam in .diseases) {
            tt <- simulate_fatality(mc_, output, disease_nam,
                                    design_$sim_prm)
            absorb_dt(dt, tt, on = c("year", "age", "sex", "qimd"))
          }

          # it doesn't matter that I get breast_ca fatality for men
          # there is no breast_ca incidence for men anyway


          # Fatality rate calibration ----
          ftlt_calib <- FALSE
          if (ftlt_calib) {
            # second run with fatalities to calibrate to mortality projections
            output <- list()
            output <-
              gen_output("", design_$sim_prm, lags_mc, dt, output)
            output <- rbindlist(output, idcol = FALSE)

            absorb_dt(dt, output, on = c("year", "pid"))
            rm(output)

            dt[, prb_nonmodelled_mrtl :=
                 fifelse(t2dm_prvl > 0L,
                         prb_nonmodelled_mrtl_not2dm *
                           nonmodelled_mrtl_t2dm_mltp,
                         prb_nonmodelled_mrtl_not2dm)]
            nam <- grep("^prb_.*_mrtl$", names(dt), value = TRUE)
            nam <- c(nam, paste0(.diseases, "_prvl"))

            tt <- dt[between(age, design_$sim_prm$ageL,
                             design_$sim_prm$ageH) &
                       year >= design_$sim_prm$init_year,
                     .SD,
                     .SDcols = c("year", "age", "sex", "qimd", nam)]
            dt[, prb_nonmodelled_mrtl := NULL]

            # delete longdeads
            # setkey(tt, pid, year)
            # tt[, pid_mrk := mk_new_simulant_markers(pid)]
            # tt[, dead := identify_longdeads(all_cause_mrtl, pid_mrk)]
            # tt <- tt[!(dead), .SD, .SDcols = c("year", "age", "sex", "qimd", nam)]


            set(tt, NULL, "pops", 1)
            # prb now represent uncalibrated deaths after sum
            for (disease_nam in .diseases) {
              ftlt <- paste0("prb_", disease_nam, "_mrtl")
              prvl <- paste0(disease_nam, "_prvl")
              set(tt, which(tt[[prvl]] == 0), ftlt, 0)
              set(tt, NULL, prvl, NULL)
            }
            tt <- tt[, lapply(.SD, sum), keyby = .(year, age, sex, qimd)]

            for (disease_nam in c("nonmodelled", .diseases)) {
              expected_deaths <- get_lifetable_all(mc_, disease_nam,
                                                   design_$sim_prm, "mx")
              colnam <- paste0("prb_", disease_nam, "_mrtl")
              setnames(tt, colnam, "col__")
              expected_deaths[tt,
                              calib_mltp := qx_mc * i.pops / i.col__,
                              on = .(year, age, sex, qimd)]
              expected_deaths[!is.finite(calib_mltp), calib_mltp := 1]
              if (disease_nam == "nonmodelled") {
                colnam <- paste0(colnam, "_not2dm")
              }
              setnames(dt, colnam, "col__")

              dt[expected_deaths,
                 (colnam) := clamp(i.calib_mltp * col__),
                 on = .(year, age, sex, qimd)]
              dt[, col__ := NULL]
              tt[, `:=` (col__ = NULL)]
            }
          }

          dt[, ncc := clamp(
            ncc - (chd_prvl > 0) - (stroke_prvl > 0) -
              (poststroke_dementia_prvl > 0) -
              (htn_prvl > 0) - (t2dm_prvl > 0) - (af_prvl > 0) -
              (copd_prvl > 0) - (lung_ca_prvl > 0) -
              (colon_ca_prvl > 0) -
              (breast_ca_prvl > 0),
            0L,
            10L
          )]
          # to be added back in the qaly fn. Otherwise when I prevent disease
          # the ncc does not decrease.

          # Prune & write to disk ----
          # del rn as they are reproducible
          nam <- grep("_mrtl$|^rn_", names(dt), value = TRUE)
          nam <- grep("^prb_", nam, value = TRUE, invert = TRUE)
          nam <- c(nam, c("LSOA11CD",
                          "LAD11CD",
                          "LAD11NM",
                          "tds_quintile",
                          "imd",
                          "sha",
                          "CCG17CDH",
                          "pid_mrk"))
          dt[, (nam) := NULL]

          message("Writing to disk")

          setkey(dt, pid, year) # Just in case
          write_fst(dt,
                    filename_$synthpop,
                    90) # 100 is too slow

          # Write pid based indx
          dt[, rn := .I]
          tt <-
            dt[, .(from = min(rn), to = max(rn)), keyby = .(pid)]
          write_fst(tt, filename_$indxfile, 90L)

          # rm(dt, tt)

          return(invisible(NULL))
        },


      # Load a synthpop file from disk in full or in chunks.
      get_synthpop =
        function(mc_,
                 filename_,
                 design_,
                 exclude_cols = c(),
                 chunk_number = 1L,
                 max_chunk = 1L) {
          stopifnot(chunk_number <=  max_chunk)

          mm <- metadata_fst(filename_$synthpop)
          mm <-
            setdiff(mm$columnNames,
                    exclude_cols)

          # Read the complete file if max_chunk == 1L
          if (max_chunk == 1L) {
            dt <-
              read_fst(filename_$synthpop,
                       columns = mm,
                       as.data.table = TRUE)

            # regenerate RN (need to be before any deletion of rows)
            # TODO Kismet = F
            generate_rns(
              mc_,
              dt,
              c(
                "rn_t2dm_incd",
                "rn_chd_incd",
                "rn_stroke_incd",
                "rn_copd_incd",
                "rn_lung_ca_incd",
                "rn_colon_ca_incd",
                "rn_breast_ca_incd",
                "rn_poststroke_dementia_incd",
                "rn_t2dm_dgn",
                "rn_chd_dgn",
                "rn_stroke_dgn",
                "rn_copd_dgn",
                "rn_htn_dgn",
                "rn_af_dgn",
                "rn_lung_ca_dgn",
                "rn_colon_ca_dgn",
                "rn_breast_ca_dgn",
                "rn_poststroke_dementia_dgn",
                "rn_chd_mrtl",
                "rn_stroke_mrtl",
                "rn_copd_mrtl",
                "rn_lung_ca_mrtl",
                "rn_colon_ca_mrtl",
                "rn_breast_ca_mrtl",
                "rn_nonmodelled_mrtl",
                "rn_multi_mrtl"
              )
            )
          } else {
            # if max_chunk > 1L split the file in chunks based on pid and load
            # only the chunk_number

            indx <-
              read_fst(filename_$indxfile, as.data.table = TRUE)

            # Use ceiling instead of floor to produce max_chunk + 1 groups and
            # discard the last group with the remains when pid is not multiple
            # of max_chunk. Use floor to create exactly max_chunk groups and the
            # last group contains the remains indx$pid is sorted. Additional
            # sorted below added for extra safety
            chunk <-
              split(sort(indx$pid), sort(indx$pid %% max_chunk))[[chunk_number]]
            from_row <- indx[pid == min(chunk), from]
            to_row <- indx[pid == max(chunk), to]

            dt <-
              read_fst(
                filename_$synthpop,
                columns = mm,
                from = from_row,
                to = to_row,
                as.data.table = TRUE
              )

            stopifnot(all.equal(sort(chunk), sort(unique(dt$pid))))

            # generate all RN and pick those that apply to the chunk of data
            dt_rn <- data.table(1:indx[.N, to])
            dt_rn <- generate_rns(
              mc_,
              dt_rn,
              c(
                "rn_t2dm_incd",
                "rn_chd_incd",
                "rn_stroke_incd",
                "rn_copd_incd",
                "rn_lung_ca_incd",
                "rn_colon_ca_incd",
                "rn_breast_ca_incd",
                "rn_poststroke_dementia_incd",
                "rn_t2dm_dgn",
                "rn_chd_dgn",
                "rn_stroke_dgn",
                "rn_copd_dgn",
                "rn_htn_dgn",
                "rn_af_dgn",
                "rn_lung_ca_dgn",
                "rn_colon_ca_dgn",
                "rn_breast_ca_dgn",
                "rn_poststroke_dementia_dgn",
                "rn_chd_mrtl",
                "rn_stroke_mrtl",
                "rn_copd_mrtl",
                "rn_lung_ca_mrtl",
                "rn_colon_ca_mrtl",
                "rn_breast_ca_mrtl",
                "rn_nonmodelled_mrtl",
                "rn_multi_mrtl"
              )
            )[from_row:to_row, ][, V1 := NULL]

            for (nam in names(dt_rn))
              set(dt, NULL, nam, dt_rn[[nam]])
            rm(dt_rn)
          }

          dt <- dt[between(
            year,
            design_$sim_prm$init_year - design_$sim_prm$maxlag,
            design_$sim_prm$init_year + (
              design_$sim_prm$init_year_fromGUI -
                design_$sim_prm$init_year
            ) + design_$sim_prm$sim_horizon_fromGUI
          ) &
            between(age,
                    design_$sim_prm$ageL - design_$sim_prm$maxlag,
                    design_$sim_prm$ageH)]

          # dt <- dt[between(year,
          #                  design$init_year_fromGUI,
          #                  sum(fromGUI_timeframe(parameters)) - 2000L) &
          #            between(age, design$ageL - max_lag, design$ageH)]

          dt[, pid_mrk := mk_new_simulant_markers(pid)]
          # Above necessary because of pruning  and potential merging above



          # simulate disease epidemiology
          lags_mc <-
            get_lag_mc(mc_, design_$sim_prm)
          max_lag_mc <- max(unlist(lags_mc))

          output <- list()
          output <-
            gen_output("", design_$sim_prm, lags_mc, dt, output)
          output <- rbindlist(output, idcol = FALSE)

          absorb_dt(dt, output, on = c("year", "pid"))
          rm(output)
          dt[, dead := identify_longdeads(all_cause_mrtl, pid_mrk)]

          setcolorder(dt, c("pid", "pid_mrk", "year", "age", "sex", "qimd"))
          setkeyv(dt, c("pid", "year"))
          setindexv(dt, c("year", "age", "sex", "qimd", "ethnicity"))
          invisible(dt)

        },

      # Calculate weights so that their sum is the population of the area based
      # on ONS
      calc_pop_weights = function(dt, locality) {
        tt <-
          read_fst("./ONS_data/pop_size/pop_proj.fst", as.data.table = TRUE)
        lads <- get_unique_LADs(locality)
        tt <- tt[LAD17CD %in% lads &
                   between(age, min(dt$age), max(dt$age)) &
                   between(year - 2000L, min(dt$year), max(dt$year)),
                 .(pops = sum(pops)), keyby = .(year, age, sex)]
        tt[, year := year - 2000L]
        dt[, wt := .N, by = .(year, age, sex)]
        absorb_dt(dt, tt)
        dt[, wt := pops / wt]
        dt[, pops := NULL]

        invisible(dt)
      }




    )
  )
