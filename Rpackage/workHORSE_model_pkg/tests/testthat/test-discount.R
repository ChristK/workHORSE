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

# Expected values: HM Treasury Green Book supplementary guidance on
# discounting (Feb 2026), Annex A, Tables A.1 (3.5%) and A.2 (1.5%, health).

test_that("discount factors reproduce Green Book tables A.1 and A.2", {
  expect_equal(round(discount_factor(3.5, 2021:2024, 2021L), 4),
               c(1, 0.9662, 0.9335, 0.9019))
  expect_equal(round(discount_factor(1.5, 2021:2024, 2021L), 4),
               c(1, 0.9852, 0.9707, 0.9563))
  expect_equal(round(discount_factor(3.5, 2021L + c(10L, 20L, 26L), 2021L), 4),
               c(0.7089, 0.5026, 0.4088))
  expect_equal(round(discount_factor(1.5, 2021L + c(10L, 20L, 26L), 2021L), 4),
               c(0.8617, 0.7425, 0.6790))
})

test_that("discount factors reproduce Annex A tables A.1 and A.2, years 0-30", {
  a1 <- c(1, 0.9662, 0.9335, 0.9019, 0.8714, 0.842, 0.8135, 0.786, 0.7594,
          0.7337, 0.7089, 0.6849, 0.6618, 0.6394, 0.6178, 0.5969, 0.5767,
          0.5572, 0.5384, 0.5202, 0.5026, 0.4856, 0.4692, 0.4533, 0.438,
          0.4231, 0.4088, 0.395, 0.3817, 0.3687, 0.3563)
  a2 <- c(1, 0.9852, 0.9707, 0.9563, 0.9422, 0.9283, 0.9145, 0.901, 0.8877,
          0.8746, 0.8617, 0.8489, 0.8364, 0.824, 0.8118, 0.7999, 0.788,
          0.7764, 0.7649, 0.7536, 0.7425, 0.7315, 0.7207, 0.71, 0.6995,
          0.6892, 0.679, 0.669, 0.6591, 0.6494, 0.6398)
  expect_equal(round(discount_factor(3.5, 2021L + 0:30, 2021L), 4), a1)
  expect_equal(round(discount_factor(1.5, 2021L + 0:30, 2021L), 4), a2)
})

test_that("discount factors are exact and depend only on year - base_year", {
  expect_identical(discount_factor(3.5, 2025L, 2025L), 1)
  expect_equal(discount_factor(3.5, 2022L, 2021L), 1 / 1.035)
  expect_equal(discount_factor(1.5, 2047L, 2021L), 1.015^-26)
  expect_equal(discount_factor(3.5, 2030:2033, 2030L),
               discount_factor(3.5, 2021:2024, 2021L))
  expect_identical(discount_factor(0, 2021:2047, 2021L), rep(1, 27))
  expect_equal(discount(c(100, 100), 3.5, c(2021L, 2023L), 2021L),
               c(100, 100 / 1.035^2))
})

test_that("cumulative present value of a flat stream of 1 a year", {
  expect_equal(round(sum(discount_factor(3.5, 2021:2041, 2021L)), 3), 15.212)
  expect_equal(round(sum(discount_factor(1.5, 2021:2041, 2021L)), 3), 18.169)
  expect_equal(round(sum(discount_factor(3.5, 2021:2047, 2021L)), 3), 17.890)
  expect_equal(round(sum(discount_factor(1.5, 2021:2047, 2021L)), 3), 22.399)
})

test_that("invalid input is rejected", {
  expect_error(discount_factor(3.5, 2020L, 2021L), "before base_year")
  expect_error(discount_factor(-1, 2021L, 2021L), "non-negative")
  expect_error(discount_factor(c(1.5, 3.5), 2021L, 2021L), "single")
  expect_error(discount_factor(NA_real_, 2021L, 2021L), "single")
  expect_error(discount_factor(3.5, 2021L, NA), "base_year")
  expect_error(discount_factor(3.5, 2021L, Inf), "base_year")
})

test_that("deflate() is the inverse of inflate()", {
  x <- c(10, 250, 1e6)
  expect_equal(deflate(inflate(x, 3.5, 2030L, 2021L), 3.5, 2030L, 2021L), x)
  expect_equal(inflate(deflate(x, 2, 2015L, 2021L), 2, 2015L, 2021L), x)
  expect_equal(deflate(1, 3.5, 2031L, 2021L), discount_factor(3.5, 2031L, 2021L))
})

test_that("discount_dt() discounts _cost columns at the cost rate and _utility columns at the QALY rate", {
  dt <- data.table::data.table(
    year = 2021:2024, mc = 1L, net_utility = 1, eq5d = 1,
    policy_cost = 100, healthcare_cost = -50, pops = 10L, nmb_cml = 7
  )
  out <- discount_dt(dt, 3.5, 1.5, 2021L)
  expect_identical(out, dt) # modified in place, returned invisibly
  expect_equal(dt$policy_cost, 100 / 1.035^(0:3))
  expect_equal(dt$healthcare_cost, -50 / 1.035^(0:3))
  expect_equal(dt$net_utility, 1 / 1.015^(0:3))
  expect_identical(dt$eq5d, rep(1, 4))
  expect_identical(dt$pops, rep(10L, 4))
  expect_identical(dt$nmb_cml, rep(7, 4))
})

test_that("discount_dt() with zero rates is a no-op and explicit columns work", {
  dt <- data.table::data.table(year = 2021:2023, net_utility = c(1, 2, 3),
                               x_cost = c(4, 5, 6), eq5d = c(7, 8, 9))
  ref <- data.table::copy(dt)
  discount_dt(dt, 0, 0, 2021L)
  expect_identical(dt, ref)
  discount_dt(dt, 0, 1.5, 2021L, qaly_cols = c("net_utility", "eq5d"))
  expect_equal(dt$eq5d, c(7, 8, 9) / 1.015^(0:2))
  expect_error(discount_dt(dt, 3.5, 1.5, 2022L), "before base_year")
  expect_error(discount_dt(dt, 3.5, 1.5, 2021L, cost_cols = "missing_cost"))
})

test_that("discount_dt() leaves dt unchanged when an argument is invalid", {
  dt <- data.table::data.table(year = 2021:2023, x_cost = c(4, 5, 6),
                               net_utility = c(1, 2, 3))
  ref <- data.table::copy(dt)
  expect_error(discount_dt(dt, 3.5, -1, 2021L), "non-negative")
  expect_identical(dt, ref) # the cost column was not discounted first
  expect_error(discount_dt(dt, 3.5, 1.5, 2022L), "before base_year")
  expect_identical(dt, ref)
})

test_that("discount_dt() takes the years from year_col", {
  dt <- data.table::data.table(yr = 2030:2032, x_cost = 100)
  discount_dt(dt, 3.5, 1.5, 2030L, year_col = "yr")
  expect_equal(dt$x_cost, 100 / 1.035^(0:2))
})
