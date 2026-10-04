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

# my_qGPO() replaces gamlss.dist::qGPO() and has to reproduce it as it was in
# gamlss.dist 6.1-1, the last release with a correct generalised Poisson. In
# 6.1-11 (CRAN, September 2026) dGPO() returns Poisson densities for
# sigma > 1e-06, so that qGPO() returns Poisson quantiles.
# The reference quantiles in fixtures/qGPO_reference_6.1-1.csv were computed
# with gamlss.dist 6.1-1 by fixtures/make_qGPO_reference.R, they are not
# computed here because the installed gamlss.dist is broken. The points are a
# grid over the ranges of the T2DM duration table (lifecourse_models/
# dm_dur_table.fst: mu 6.39-13.33, sigma 0.182-0.318) and beyond (mu 0.5-30,
# sigma from the Poisson limit to 2, p from 0 to 1 including values near 1)
# plus random draws.

ref <- data.table::fread(test_path("fixtures", "qGPO_reference_6.1-1.csv"))
ref[, q := as.numeric(q)] # Inf is stored for p + 1e-09 >= 1

# The pmf written out, not taken from gamlss.dist
dgpo <- function(x, mu, sigma)
  exp(x * log(mu / (1 + sigma * mu)) + (x - 1) * log(1 + sigma * x) -
        mu * (1 + sigma * x) / (1 + sigma * mu) - lgamma(x + 1))


test_that("my_qGPO reproduces qGPO of gamlss.dist 6.1-1 exactly", {
  expect_identical(my_qGPO(ref$p, ref$mu, ref$sigma), ref$q)
  # the results do not depend on the number of threads
  expect_identical(my_qGPO(ref$p, ref$mu, ref$sigma, n_cpu = 3L), ref$q)
})

test_that("the reference covers the edge cases", {
  expect_gt(nrow(ref), 4000L)
  expect_true(all(c(0, 1e-12, 0.5, 0.999, 1) %in% ref$p))
  expect_true(any(ref$p > 0.999999 & ref$p < 1 & is.finite(ref$q)))        # near 1
  expect_true(any(ref$p < 1 & is.infinite(ref$q)))                         # p + 1e-09 >= 1
  expect_true(any(ref$q == 10000))                                         # max.value is reached
  expect_true(all(c(1e-8, 9.99e-5, 1e-4, 1.01e-4) %in% ref$sigma))         # around the Poisson limit
  expect_true(all(ref$mu >= 0.5 & ref$mu <= 30 & ref$sigma <= 2))
  expect_true(all(c(6.39, 13.33, 0.182, 0.318) %in% c(ref$mu, ref$sigma))) # dm_dur_table ranges
  expect_gt(sum(ref$mu >= 6.39 & ref$mu <= 13.33 & ref$sigma >= 0.182 & ref$sigma <= 0.318), 700L)
})

test_that("my_qGPO is not the Poisson quantile function of gamlss.dist 6.1-11", {
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)
  q <- my_qGPO(p, mu = 8, sigma = 0.25)
  expect_identical(q, c(0, 1, 2, 5, 11, 18, 40, 65)) # as computed by gamlss.dist 6.1-1
  expect_false(isTRUE(all.equal(q, qpois(p, lambda = 8)))) # qpois(p, 8) is 2 5 6 8 10 12 15 18
  expect_true(all(q[p >= 0.75] > qpois(p[p >= 0.75], lambda = 8))) # overdispersed
  # the same for the upper tail in the reference
  up <- ref[sigma >= 0.182 & mu >= 5 & p >= 0.99 & p < 1 - 1e-8 & is.finite(q)]
  expect_gt(nrow(up), 100L)
  expect_true(all(up$q > qpois(up$p, up$mu)))
})

test_that("my_qGPO is the Poisson quantile function for sigma < 1e-04", {
  p <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999)
  expect_identical(my_qGPO(p, mu = 5, sigma = 1e-8), qpois(p, lambda = 5))
  expect_identical(my_qGPO(p, mu = 20, sigma = 9.99e-5), qpois(p, lambda = 20))
})

test_that("my_qGPO inverts the cdf of the generalised Poisson", {
  set.seed(1L)
  n <- 400L
  mu <- runif(n, 0.5, 15); sigma <- runif(n, 0.01, 0.5); p <- runif(n, 0, 0.9999)
  q <- my_qGPO(p, mu, sigma)
  ok <- vapply(seq_len(n), function(i) {
    cdf <- cumsum(dgpo(0:(q[i] + 1), mu[i], sigma[i]))
    # the smallest x with P(X <= x) >= p (the sums are computed in a different way)
    cdf[q[i] + 1] >= p[i] - 1e-12 && (q[i] == 0 || cdf[q[i]] < p[i] + 1e-12)
  }, NA)
  expect_true(all(ok))
  # the parameterisation: mean mu and variance mu * (1 + sigma * mu)^2
  x <- 0:5000
  for (par in list(c(8, 0.25), c(13.33, 0.318), c(2, 0.5))) {
    pmf <- dgpo(x, par[1L], par[2L])
    expect_equal(sum(pmf), 1)
    expect_equal(sum(x * pmf), par[1L])
    expect_equal(sum(x^2 * pmf) - par[1L]^2, par[1L] * (1 + par[2L] * par[1L])^2)
  }
})

test_that("my_qGPO(cdf(x)) is x at the exact values of the cdf", {
  # F(x) is summed in R in the way pGPO() of gamlss.dist 6.1-1 does (cumsum() and
  # sum() accumulate in long double). my_qGPO has to round its cdf in exactly the
  # same way: p = F(x) is the smallest probability with quantile x, and a cdf that
  # is rounded differently, e.g. by fused multiply-adds, gives x + 1.
  set.seed(5L)
  n <- 300L; K <- 100L
  mu <- runif(n, 0.5, 30); sigma <- exp(runif(n, log(1e-3), log(2)))
  wrong <- vapply(seq_len(n), function(i) {
    cdf <- cumsum(dgpo(0:K, mu[i], sigma[i]))
    ok  <- cdf < 1 - 2e-9 & c(TRUE, diff(cdf) > 0) # finite and distinct quantiles
    sum(my_qGPO(cdf[ok], mu[i], sigma[i]) != (0:K)[ok])
  }, 0L)
  expect_identical(sum(wrong), 0L)
  # sigma < 1e-04: the Poisson cdf
  for (par in list(c(5, 1e-8), c(20, 9.99e-5), c(0.5, 1e-6))) {
    cdf <- ppois(0:K, par[1L])
    ok  <- cdf < 1 - 2e-9 & c(TRUE, diff(cdf) > 0)
    expect_identical(my_qGPO(cdf[ok], par[1L], par[2L]), as.numeric((0:K)[ok]))
  }
})

test_that("my_qGPO reaches mu in mean when the probabilities are uniform", {
  # the use in the synthpop: durations of T2DM from uniform ranks
  set.seed(2L)
  n <- 2e5L
  for (par in list(c(6.39, 0.182), c(13.33, 0.318))) {
    q <- my_qGPO(runif(n), par[1L], par[2L])
    se <- sqrt(par[1L] * (1 + par[2L] * par[1L])^2 / n)
    expect_lt(abs(mean(q) - par[1L]), 5 * se)
  }
})

test_that("my_qGPO validates the input like gamlss.dist::qGPO", {
  expect_error(my_qGPO(0.5, mu = 0, sigma = 0.25), "mu must be greater than 0")
  expect_error(my_qGPO(0.5, mu = -1, sigma = 0.25), "mu must be greater than 0")
  expect_error(my_qGPO(0.5, mu = 8, sigma = 0), "sigma must be greater than 0")
  expect_error(my_qGPO(0.5, mu = 8, sigma = -0.25), "sigma must be greater than 0")
  expect_error(my_qGPO(c(0.2, 0.5), mu = c(8, 0), sigma = 0.25), "mu must be greater than 0")
  expect_error(my_qGPO(c(0.2, 0.5), mu = 8, sigma = c(0.25, -1)), "sigma must be greater than 0")
  expect_error(my_qGPO(0.5, mu = NA_real_, sigma = 0.25), "mu must be greater than 0")
  expect_error(my_qGPO(0.5, mu = 8, sigma = NA_real_), "sigma must be greater than 0")
  expect_error(my_qGPO(0.5, mu = Inf, sigma = 0.25), "mu must be greater than 0")
  expect_error(my_qGPO(-0.1, mu = 8, sigma = 0.25), "p must be between 0 and 1")
  expect_error(my_qGPO(1.1, mu = 8, sigma = 0.25), "p must be between 0 and 1")
  expect_error(my_qGPO(NA_real_, mu = 8, sigma = 0.25), "p must be between 0 and 1")
  expect_error(my_qGPO(c(0.1, 0.2, 0.3), mu = c(8, 9), sigma = 0.25), "same length")
  expect_error(my_qGPO(0.5, mu = 8, sigma = 0.25, max_value = -1L), "max_value")
})

test_that("the arguments of my_qGPO", {
  p  <- c(0.05, 0.3, 0.6, 0.95)
  p0 <- p
  q  <- my_qGPO(p, 8, 0.25)
  expect_identical(my_qGPO(p, rep(8, 4), rep(0.25, 4)), q)           # length 1 is recycled
  expect_identical(my_qGPO(p, 8, 0.25, lower_tail = FALSE), my_qGPO(1 - p, 8, 0.25))
  expect_identical(my_qGPO(log(p), 8, 0.25, log_p = TRUE), my_qGPO(exp(log(p)), 8, 0.25))
  expect_identical(p, p0)                                            # the input is not modified
  expect_identical(my_qGPO(c(0.2, 0.99), 8, 0.25, max_value = 20L), c(my_qGPO(0.2, 8, 0.25), 20))
  expect_identical(my_qGPO(c(0.2, 0.9), 8, 0.25, max_value = 0L), c(0, 0))
  expect_identical(my_qGPO(c(1, 1 - 1e-12), 8, 0.25), c(Inf, Inf)) # as qGPO
  expect_identical(my_qGPO(numeric(), 8, 0.25), numeric())
  expect_false(is.unsorted(my_qGPO(seq(0, 0.9999, length.out = 500), 8, 0.25)))
  expect_identical(my_qGPO(c(0.2, 0.9), 8L, 0.25), my_qGPO(c(0.2, 0.9), 8, 0.25)) # integer mu
  expect_type(q, "double")
})
