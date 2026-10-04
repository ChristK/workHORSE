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

# Writes qGPO_reference_6.1-1.csv, the reference of test-qGPO.R: the quantiles
# of the generalised Poisson distribution computed by gamlss.dist::qGPO().
#
# THIS SCRIPT MUST BE RUN WITH gamlss.dist 6.1-1. That is the last release with
# a correct generalised Poisson: in 6.1-11 (CRAN, September 2026) dGPO()
# returns Poisson densities for sigma > 1e-06 and so qGPO() returns Poisson
# quantiles. The script stops with any other version.
# To install 6.1-1 in a separate library:
#   install.packages("https://cran.r-project.org/src/contrib/Archive/gamlss.dist/gamlss.dist_6.1-1.tar.gz",
#                    repos = NULL, type = "source", lib = "<lib_gd611>")
#
# usage: Rscript make_qGPO_reference.R [<lib_gd611> [<n_cores>]]
# <lib_gd611> is put first in .libPaths(). qGPO() is slow, it calls pGPO() for
# every x: evaluating the points takes about 3 minutes on one core. They can be
# split over <n_cores> forked processes (Unix only; default 1). The result does
# not depend on <n_cores>. The file is written next to this script.

args <- commandArgs(TRUE)
if (length(args) >= 1L) .libPaths(c(args[1L], .libPaths()))
n_cores <- if (length(args) >= 2L) as.integer(args[2L]) else 1L
this_dir <- {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
out_file <- file.path(this_dir, "qGPO_reference_6.1-1.csv")

stopifnot(requireNamespace("data.table", quietly = TRUE),
          requireNamespace("gamlss.dist", quietly = TRUE))
if (as.character(utils::packageVersion("gamlss.dist")) != "6.1.1")
  stop("gamlss.dist 6.1-1 is required, found ", format(utils::packageVersion("gamlss.dist")),
       " in ", find.package("gamlss.dist"))
cat("gamlss.dist", format(utils::packageVersion("gamlss.dist")), "from", find.package("gamlss.dist"), "\n")

# The points ---------------------------------------------------------------
# A grid that covers the ranges of dm_dur_table.fst (mu 6.39-13.33, sigma
# 0.182-0.318) and goes beyond them: mu 0.5-30; sigma from the Poisson limit
# (1e-08, and 1e-05 and 9.99e-05 below the 1e-04 threshold of pGPO) up to 2;
# p from 0 and tiny values up to 1 and values near 1, where qGPO() returns Inf
# (p + 1e-09 >= 1) or reaches its max.value of 10000.
grid <- data.table::CJ(
  p = c(0, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95,
        0.99, 0.999, 0.9999, 0.999999, 0.99999999, 0.999999998, 0.999999999,
        0.9999999995, 1),
  mu = c(0.5, 1, 2, 5, 6.39, 8, 10, 13.33, 20, 30),
  sigma = c(1e-8, 1e-5, 9.99e-5, 1e-4, 1.01e-4, 0.01, 0.05, 0.182, 0.25, 0.318,
            0.5, 1, 2),
  sorted = FALSE
)
# Random draws, rounded so that they are exact decimals (15 digits at most)
set.seed(20260404L)
n_rnd <- 1000L
rnd_wide <- data.table::data.table(
  p = round(runif(n_rnd), 10), mu = round(runif(n_rnd, 0.5, 30), 4),
  sigma = round(runif(n_rnd, 0.01, 2), 4))
rnd_dm <- data.table::data.table( # the ranges of dm_dur_table.fst
  p = round(runif(n_rnd / 2L), 10), mu = round(runif(n_rnd / 2L, 6.39, 13.33), 4),
  sigma = round(runif(n_rnd / 2L, 0.182, 0.318), 4))
pts <- rbind(grid, rnd_wide, rnd_dm)

# The numbers go through text (15 significant digits) and are parsed back with
# the reader of the test (fread), so that the quantiles are computed at exactly
# the values the test reads
txt <- pts[, lapply(.SD, function(x) sprintf("%.15g", x))]
tmp <- tempfile(fileext = ".csv")
data.table::fwrite(txt, tmp, quote = FALSE)
ref <- data.table::fread(tmp)
stopifnot(nrow(ref) == nrow(pts), all(vapply(ref, is.double, NA)))

# qGPO() is slow when the quantile is large: balance the load
ord <- order(ref$mu * ref$sigma * (ref$p + 0.01), decreasing = TRUE)
chunks <- split(ord, ceiling(seq_along(ord) / 10L))
t0 <- Sys.time()
res <- parallel::mclapply(chunks, function(ix)
  gamlss.dist::qGPO(ref$p[ix], ref$mu[ix], ref$sigma[ix]),
  mc.cores = n_cores, mc.preschedule = FALSE)
stopifnot(!vapply(res, inherits, NA, "try-error"))
ref[unlist(chunks, use.names = FALSE), q := unlist(res, use.names = FALSE)]
cat(sprintf("qGPO() evaluated at %d points in %.1f min\n", nrow(ref),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
stopifnot(!anyNA(ref$q), all(ref$q == round(ref$q)), all(ref$q >= 0))

# Write q as text as well and check the file reads back as computed
txt[, q := ifelse(is.infinite(ref$q), "Inf", sprintf("%.0f", ref$q))]
data.table::fwrite(txt, out_file, quote = FALSE)
chk <- data.table::fread(out_file)
stopifnot(identical(chk$p, ref$p), identical(chk$mu, ref$mu),
          identical(chk$sigma, ref$sigma), identical(as.numeric(chk$q), ref$q))
cat(sprintf("%d rows written to %s (q = Inf in %d, q = 10000 in %d)\n", nrow(chk),
            out_file, sum(is.infinite(ref$q)), sum(ref$q == 10000)))
