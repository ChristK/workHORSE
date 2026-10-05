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

# The checksum in the file names of the cached synthpops (private$gen_checksum()
# of SynthPop) decides whether a cached synthpop is reused. It has to change
# with everything that changes the generated values, including the versions of
# the code and of the dependencies (synthpop_code_versions()), otherwise stale
# cached synthpops, with their population weights, are silently reused.
# SynthPop$new(0L, design) initialises an empty synthpop and computes the
# checksum, so no synthpop is built here.
#
# The checksum reads ./synthpop/lsoa_to_locality_indx.fst and the population
# files in ./ONS_data/pop_size/ relative to the working directory, i.e. the
# repo root. The tests that need them run from the nearest directory above the
# tests that holds these files and simulation/sim_design.yaml, and are skipped
# when the package is tested outside the repo.

# Nearest directory above the tests that holds the files the checksum reads and
# the simulation design; NULL when the package is tested outside the repo
find_model_root <- function(from = getwd()) {
  needed <- c("simulation/sim_design.yaml",
              "synthpop/lsoa_to_locality_indx.fst",
              "ONS_data/pop_size/pop_proj.fst",
              "ONS_data/pop_size/pop_proj_england.fst")
  d <- normalizePath(from, mustWork = FALSE)
  repeat {
    if (all(file.exists(file.path(d, needed)))) return(d)
    up <- dirname(d)
    if (identical(up, d)) return(NULL)
    d <- up
  }
}

# Evaluate code with the working directory set to dir, restored afterwards
in_dir <- function(dir, code) {
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  force(code)
}

# The design of simulation/sim_design.yaml for a locality. The output and
# synthpop directories are in the temporary directory of the session, so that
# nothing is created in the repo or in the shared synthpop directory.
make_design <- function(root, locality = "England") {
  prm <- yaml::read_yaml(file.path(root, "simulation", "sim_design.yaml"))
  prm$output_dir   <- tempdir()
  prm$synthpop_dir <- tempdir()
  design <- Design$new(prm)
  design$sim_prm$locality <- locality
  design
}

# The private environment of an empty synthpop (mc = 0) with this design. The
# checksum is computed when the object is initialised.
empty_synthpop_private <- function(design) {
  SynthPop$new(0L, design)$.__enclos_env__$private
}
checksum_of <- function(design) empty_synthpop_private(design)$checksum

root <- find_model_root()


test_that("synthpop_code_versions() returns the five versions that determine the generated values", {
  v <- workHORSEmisc:::synthpop_code_versions()
  expect_type(v, "character")
  expect_named(v, c("workHORSEmisc", "gamlss.dist", "dqrng", "CKutils", "R"))
  expect_false(anyNA(v))
  expect_true(all(nzchar(v)))
  for (p in c("workHORSEmisc", "gamlss.dist", "dqrng", "CKutils"))
    expect_identical(v[[p]], as.character(utils::packageVersion(p)), label = p)
  # R is major.minor of the running R, without the patch level
  expect_match(v[["R"]], "^[0-9]+\\.[0-9]+$")
  expect_true(startsWith(format(getRversion()), paste0(v[["R"]], ".")))
})

test_that("the checksum is stable across calls", {
  skip_if_not(!is.null(root), "model files not found above the tests")
  in_dir(root, {
    d <- make_design(root, "Rutland")
    a <- checksum_of(d)
    expect_match(a, "^[0-9a-f]{32}$")
    expect_identical(checksum_of(d), a)
    # a new design object with the same settings
    expect_identical(checksum_of(make_design(root, "Rutland")), a)
    # the checksum that SynthPop prints
    expect_output(SynthPop$new(0L, d)$get_checksum(), a, fixed = TRUE)
  })
})

test_that("the checksum changes when the versions change", {
  skip_if_not(!is.null(root), "model files not found above the tests")
  skip_if_not(utils::packageVersion("testthat") >= "3.2.0",
              "with_mocked_bindings() needs testthat >= 3.2.0")
  in_dir(root, {
    d <- make_design(root, "Rutland")
    v <- workHORSEmisc:::synthpop_code_versions()
    # checksum when synthpop_code_versions() returns `versions`
    checksum_with <- function(versions)
      testthat::with_mocked_bindings(
        checksum_of(d),
        synthpop_code_versions = function() versions,
        .package = "workHORSEmisc"
      )

    a <- checksum_of(d)
    # the helper is the only source of the versions in the checksum
    expect_identical(checksum_with(v), a)
    # a new version of any one of them gives a new checksum
    changed <- vapply(names(v), function(nm) {
      v2 <- v
      v2[[nm]] <- paste0(v2[[nm]], ".1")
      checksum_with(v2)
    }, character(1))
    expect_false(a %in% changed)
    expect_false(anyDuplicated(changed) > 0L)
    # also when an entry is added or removed
    expect_false(checksum_with(c(v, extra = "1.0")) %in% c(a, changed))
    expect_false(checksum_with(v[-1L]) %in% c(a, changed))
  })
})

test_that("England and the nine regions have different checksums", {
  skip_if_not(!is.null(root), "model files not found above the tests")
  in_dir(root, {
    regions <- fst::read_fst("./synthpop/lsoa_to_locality_indx.fst",
                             columns = "RGN11NM", as.data.table = TRUE)
    regions <- sort(unique(as.character(regions$RGN11NM)))
    expect_length(regions, 9L)
    cs <- vapply(c("England", regions),
                 function(loc) checksum_of(make_design(root, loc)),
                 character(1))
    expect_length(unique(cs), 10L)
    # England has the same LSOAs as all the regions together, but uses other
    # population data
    expect_false(checksum_of(make_design(root, regions)) %in% cs)
  })
})

test_that("the synthpop metadata keeps the characteristics and adds the code versions", {
  skip_if_not(!is.null(root), "model files not found above the tests")
  in_dir(root, {
    d    <- make_design(root, "Rutland")
    priv <- empty_synthpop_private(d)
    chr  <- priv$get_unique_characteristics(d)
    md   <- priv$get_metadata(d)
    v    <- as.list(workHORSEmisc:::synthpop_code_versions())
    # the existing fields are unchanged and first, the versions are added
    expect_identical(names(md), c(names(chr), "code_versions"))
    expect_identical(md[names(chr)], chr)
    expect_identical(md$code_versions, v)
    # as written to and read back from the yaml metafile
    f <- tempfile(fileext = "_meta.yaml")
    yaml::write_yaml(md, f)
    back <- yaml::read_yaml(f)
    unlink(f)
    expect_equal(back[names(chr)], chr)
    expect_identical(back$code_versions, v)
  })
})

test_that("gen_synthpop() writes this metadata to the metafile", {
  src <- deparse(body(SynthPop$private_methods$gen_synthpop))
  expect_true(any(grepl("yaml::write_yaml(private$get_metadata(design_)",
                        src, fixed = TRUE)))
})
