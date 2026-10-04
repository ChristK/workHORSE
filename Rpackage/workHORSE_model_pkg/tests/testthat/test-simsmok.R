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

# Direct tests of the smoking trajectory functions in src/simsmok.cpp with small
# synthetic inputs:
# * simsmok(): baseline model, edits the data.frame in place.
# * simsmok_cessation(): health check (HC) cessation scenario, returns the
#   scenario's smok_status, smok_quit_yrs and smok_dur.
#
# pr_relapse is a 10 x 15 matrix (rows: men qimd 1-5, then women qimd 1-5;
# columns: years since quitting 1-15). pr_distinct() gives every cell a
# different value, so the outcome of a relapse draw tells which cell was read.
# With relapse_rn just below the value of the cell that should be read the
# simulant relapses, with relapse_rn just above it the simulant does not.
# Reading a cell with a smaller value fails the first check, and reading a cell
# with a larger value fails the second.

pr_distinct <- function() outer(1:10, 1:15, function(i, j) 0.001 * i + 0.01 * j)
pr_row      <- function(sex, qimd) ifelse(sex == 1L, qimd, qimd + 5L) # 1-based row of pr_relapse
eps         <- 1e-7


# ---- simsmok(): baseline model ----

# data.frame for simsmok() with two rows per simulant: the first year (pid_mrk
# TRUE) holds the inputs, the second year is simulated. The vectors are fresh
# because simsmok() modifies them in place.
sm_df <- function(status, quit_yrs, dur, rn, sex = 1L, qimd = 1L, incid = 0, cess = 0) {
  n   <- length(status)
  two <- function(first, second)
    c(rbind(rep(first, length.out = n), rep(second, length.out = n)))
  data.frame(
    smok_status    = two(status, 0L),
    prb_smok_incid = two(0, incid),
    prb_smok_cess  = two(0, cess),
    rankstat_smok  = two(0.5, rn),
    pid_mrk        = rep(c(TRUE, FALSE), n),
    sex            = two(sex, sex),
    qimd           = two(qimd, qimd),
    smok_quit_yrs  = two(quit_yrs, 0L),
    smok_dur       = two(dur, 0L)
  )
}

# ex smokers (status 2 and 3) with 0 to 3 years since quitting in the first
# year, men and women of every qimd, relapse draw just below and just above the
# expected cell
sm_cases <- function(pr) {
  cs     <- expand.grid(sex = 1:2, qimd = 1:5, q = 0:3, side = c(-1, 1), status = 2:3)
  cs$col <- pmax(cs$q, 1L) # 0 years since quitting counts as the first year after quitting
  cs$rn  <- pr[cbind(pr_row(cs$sex, cs$qimd), cs$col)] + cs$side * eps
  cs
}
sm_run_cases <- function(cs, pr) {
  df <- sm_df(cs$status, cs$q, 10L, cs$rn, cs$sex, cs$qimd)
  simsmok(df, pr, 3L)
  df
}

test_that("simsmok(): the relapse draw reads the column of the years since quitting, 0 years counting as the first column", {
  pr  <- pr_distinct()
  cs  <- sm_cases(pr)
  df  <- sm_run_cases(cs, pr)
  yr1 <- seq(1L, 2L * nrow(cs), by = 2L)
  yr2 <- yr1 + 1L
  relapse <- cs$side < 0
  expect_identical(df$smok_status[yr2],   ifelse(relapse, 4L, cs$status))
  expect_identical(df$smok_quit_yrs[yr2], ifelse(relapse, 0L, cs$q + 1L))
  expect_identical(df$smok_dur[yr2],      ifelse(relapse, 11L, 10L))
  # the first year of each simulant is not simulated
  expect_identical(df$smok_status[yr1],   cs$status)
  expect_identical(df$smok_quit_yrs[yr1], cs$q)
  expect_identical(df$smok_dur[yr1],      rep(10L, nrow(cs)))
})

test_that("simsmok(): beyond the relapse cut-off an ex smoker never relapses", {
  q  <- c(4L, 5L, 10L, 15L)
  df <- sm_df(rep(3L, 4L), q, 10L, rn = 0, sex = c(1L, 2L, 1L, 2L), qimd = c(1L, 2L, 5L, 5L))
  simsmok(df, pr_distinct(), 3L) # all cells are > 0 = relapse_rn
  expect_identical(df$smok_status[c(2L, 4L, 6L, 8L)],   rep(3L, 4L))
  expect_identical(df$smok_quit_yrs[c(2L, 4L, 6L, 8L)], q + 1L)
  expect_identical(df$smok_dur[c(2L, 4L, 6L, 8L)],      rep(10L, 4L))
})

test_that("simsmok(): never smokers start and smokers quit with their incidence and cessation probabilities", {
  df <- sm_df(status = c(1L, 1L, 4L, 4L), quit_yrs = 0L, dur = c(0L, 0L, 7L, 7L),
              rn = c(0.01, 0.9, 0.01, 0.9), incid = 0.5, cess = 0.5)
  simsmok(df, pr_distinct(), 3L)
  sim <- c(2L, 4L, 6L, 8L)
  expect_identical(df$smok_status[sim],   c(4L, 1L, 3L, 4L))
  expect_identical(df$smok_quit_yrs[sim], c(0L, 0L, 1L, 0L))
  expect_identical(df$smok_dur[sim],      c(1L, 0L, 7L, 8L))
})

test_that("simsmok(): repeated calls give identical results (no reads outside pr_relapse)", {
  # includes ex smokers with 0 years since quitting. pr_relapse is rebuilt every
  # time so that a read outside the matrix would not see the same memory
  run   <- function() sm_run_cases(sm_cases(pr_distinct()), pr_distinct())
  first <- run()
  for (i in 1:25) expect_identical(run(), first)
})

test_that("simsmok() stops when relapse_cutoff exceeds the columns of pr_relapse", {
  df <- sm_df(3L, 1L, 10L, rn = 0.5)
  expect_error(simsmok(df, matrix(0, 10L, 2L), 3L), "relapse_cutoff")
})


# ---- simsmok_cessation(): health check cessation scenario ----

# One simulant followed for 10 years who is a current smoker in the baseline
# trajectory throughout (smok_quit_yrs 0, smok_dur 20, 21, ..., 29). hc_years
# are the years with hc_eff == 1 (the health check effect is not carried forward).
cess_run <- function(hc_years, pr, rn = rep(0.5, 10L), sex = 1L, qimd = 1L, cutoff = 3L) {
  n  <- 10L
  hc <- integer(n)
  hc[hc_years] <- 1L
  simsmok_cessation(
    smok_status    = rep(4L, n),
    smok_quit_yrs  = rep(0L, n),
    smok_dur       = 20L + 0:(n - 1L),
    sex            = rep(sex, n),
    qimd           = rep(qimd, n),
    new_pid        = c(TRUE, rep(FALSE, n - 1L)),
    hc_eff         = hc,
    relapse_rn     = rn,
    pr_relapse     = pr,
    relapse_cutoff = cutoff
  )
}

test_that("simsmok_cessation(): a quitter who does not relapse stays an ex smoker, also beyond the relapse cut-off", {
  x <- cess_run(2L, matrix(0, 10L, 15L))
  expect_identical(x$smok_status,   c(4L, rep(3L, 9L)))
  expect_identical(x$smok_quit_yrs, 0:9)
  expect_identical(x$smok_dur,      rep(20L, 10L))
})

test_that("simsmok_cessation(): a relapse returns the simulant to the baseline until the next health check", {
  x <- cess_run(c(2L, 7L), matrix(1, 10L, 15L)) # certain relapse
  expect_identical(x$smok_status,   c(4L, 3L, 4L, 4L, 4L, 4L, 3L, 4L, 4L, 4L))
  expect_identical(x$smok_quit_yrs, c(0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L))
  expect_identical(x$smok_dur,      c(20L, 20L, 21L, 23L, 24L, 25L, 25L, 26L, 28L, 29L))
})

# Health check in year 2, so k years since quitting at the relapse draw of year
# 2 + k, which has to read column k (and not the column of the baseline's years
# since quitting, which is 0). Earlier draws never relapse (rn 0.9 > every cell).
for (k in 1:3) {
  test_that(sprintf("simsmok_cessation(): the relapse draw %d year(s) after quitting reads column %d", k, k), {
    pr <- pr_distinct()
    r  <- pr_row(2L, 3L)
    rn <- rep(0.9, 10L)

    rn[2L + k] <- pr[r, k] - eps # relapses unless a cell of a smaller value is read
    x <- cess_run(2L, pr, rn, sex = 2L, qimd = 3L)
    expect_identical(x$smok_status,   c(4L, rep(3L, k), rep(4L, 9L - k)))
    expect_identical(x$smok_quit_yrs, c(0L, 1:k, rep(0L, 9L - k)))
    expect_identical(x$smok_dur,      c(rep(20L, 1L + k), 21L, (22L + k):29L))

    rn[2L + k] <- pr[r, k] + eps # does not relapse unless a cell of a larger value is read
    x <- cess_run(2L, pr, rn, sex = 2L, qimd = 3L)
    expect_identical(x$smok_status,   c(4L, rep(3L, 9L)))
    expect_identical(x$smok_quit_yrs, 0:9)
    expect_identical(x$smok_dur,      rep(20L, 10L))
  })
}

test_that("simsmok_cessation(): a new health check on the quit path does not reset the years since quitting", {
  # health checks within, at the end of, and beyond the relapse cut-off
  x <- cess_run(c(2L, 3L, 5L, 6L, 8L), matrix(0, 10L, 15L))
  expect_identical(x$smok_status,   c(4L, rep(3L, 9L)))
  expect_identical(x$smok_quit_yrs, 0:9)
  expect_identical(x$smok_dur,      rep(20L, 10L))
})

test_that("simsmok_cessation(): the quit path does not carry over to the next simulant", {
  pr <- matrix(0, 10L, 15L)
  pr[10L, ] <- 1 # only women of qimd 5 certainly relapse
  # pid 1 (men, qimd 1): baseline smoker, quits at the health check of year 2 and never relapses
  # pid 2 (women, qimd 5): baseline ex smoker, no health check. Following the path of pid 1 it would relapse
  x <- simsmok_cessation(
    smok_status    = c(rep(4L, 5L), rep(3L, 5L)),
    smok_quit_yrs  = c(rep(0L, 5L), 2L, 3L, 4L, 5L, 6L),
    smok_dur       = c(20L, 21L, 22L, 23L, 24L, rep(15L, 5L)),
    sex            = c(rep(1L, 5L), rep(2L, 5L)),
    qimd           = c(rep(1L, 5L), rep(5L, 5L)),
    new_pid        = c(TRUE, rep(FALSE, 4L), TRUE, rep(FALSE, 4L)),
    hc_eff         = c(0L, 1L, rep(0L, 8L)),
    relapse_rn     = rep(0.5, 10L),
    pr_relapse     = pr,
    relapse_cutoff = 3L
  )
  expect_identical(x$smok_status,   c(4L, 3L, 3L, 3L, 3L, rep(3L, 5L)))
  expect_identical(x$smok_quit_yrs, c(0L, 1L, 2L, 3L, 4L, 2L, 3L, 4L, 5L, 6L))
  expect_identical(x$smok_dur,      c(rep(20L, 5L), rep(15L, 5L)))
})

test_that("simsmok_cessation(): a health check in the first year of a simulant quits with the duration reduced by one", {
  x <- simsmok_cessation(
    smok_status    = rep(4L, 3L),
    smok_quit_yrs  = rep(0L, 3L),
    smok_dur       = c(30L, 31L, 32L),
    sex            = rep(1L, 3L),
    qimd           = rep(1L, 3L),
    new_pid        = c(TRUE, FALSE, FALSE),
    hc_eff         = c(1L, 0L, 0L),
    relapse_rn     = rep(0.5, 3L),
    pr_relapse     = matrix(0, 10L, 15L),
    relapse_cutoff = 3L
  )
  expect_identical(x$smok_status,   c(3L, 3L, 3L))
  expect_identical(x$smok_quit_yrs, c(1L, 2L, 3L))
  expect_identical(x$smok_dur,      c(29L, 29L, 29L))
})

test_that("simsmok_cessation(): the baseline trajectory is returned unchanged when no smoker has a health check effect", {
  status <- c(1L, 1L, 4L, 4L, 3L, 3L, 4L, 4L, 4L, 3L)
  qy     <- c(0L, 0L, 0L, 0L, 1L, 2L, 0L, 0L, 0L, 1L)
  dur    <- c(0L, 0L, 1L, 2L, 2L, 2L, 3L, 4L, 5L, 5L)
  # hc_eff == 1 in years where the baseline is not a current smoker has no effect
  hc     <- c(1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L)
  x <- simsmok_cessation(status, qy, dur, rep(1L, 10L), rep(1L, 10L),
                         c(TRUE, rep(FALSE, 9L)), hc, rep(0, 10L),
                         matrix(1, 10L, 15L), 3L)
  expect_identical(x, list(smok_status = status, smok_quit_yrs = qy, smok_dur = dur))
})

test_that("simsmok_cessation(): repeated calls give identical results (no reads outside pr_relapse)", {
  # relapse in year 3 (column 1), second health check in year 7, relapse in year 8
  rn  <- c(0.9, 0.9, 0.015, 0.9, 0.9, 0.9, 0.9, 0.017, 0.9, 0.9)
  run <- function() cess_run(c(2L, 7L), pr_distinct(), rn, sex = 2L, qimd = 3L)
  first <- run()
  expect_identical(first$smok_status,   c(4L, 3L, 4L, 4L, 4L, 4L, 3L, 4L, 4L, 4L))
  expect_identical(first$smok_quit_yrs, c(0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L))
  expect_identical(first$smok_dur,      c(20L, 20L, 21L, 23L, 24L, 25L, 25L, 26L, 28L, 29L))
  for (i in 1:25) expect_identical(run(), first)
})

test_that("simsmok_cessation() stops when relapse_cutoff exceeds the columns of pr_relapse", {
  expect_error(cess_run(2L, matrix(0, 10L, 2L)), "relapse_cutoff")
})
