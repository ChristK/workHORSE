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

# Guards against gamlss.dist d/p/q/r calls binding arguments positionally to
# lower.tail/log.p/log/max.value, and checks the COPD duration statement.

# All calls (as language objects) inside an expression, skipping empty args
calls_in <- function(e, acc = list()) {
  if (is.call(e)) {
    args <- as.list(e)
    for (i in seq_along(args)) {
      if (identical(args[[i]], quote(expr = ))) next
      acc <- calls_in(args[[i]], acc)
    }
    if (is.name(e[[1L]])) acc[[length(acc) + 1L]] <- e
  }
  acc
}

# Every function in the namespace plus the public/private/active methods of R6 generators
pkg_functions <- function(pkg = "workHORSEmisc") {
  ns <- asNamespace(pkg)
  out <- list()
  for (n in ls(ns, all.names = TRUE)) {
    o <- get(n, envir = ns)
    if (is.function(o)) out[[n]] <- o
    if (inherits(o, "R6ClassGenerator"))
      for (m in c("public_methods", "private_methods", "active"))
        for (k in names(o[[m]]))
          if (is.function(o[[m]][[k]])) out[[paste0(n, "$", k)]] <- o[[m]][[k]]
  }
  out
}

test_that("gamlss.dist d/p/q/r calls bind only distribution parameters", {
  gd   <- asNamespace("gamlss.dist")
  seen <- character()
  for (f in names(fns <- pkg_functions())) {
    b <- body(fns[[f]])
    if (is.null(b)) next
    for (cl in calls_in(b)) {
      nm <- as.character(cl[[1L]])
      if (!(grepl("^[dpqr][A-Z]", nm) && exists(nm, envir = gd, inherits = FALSE))) next
      bound <- names(as.list(match.call(get(nm, envir = gd), cl)))[-1L]
      named <- names(as.list(cl))[-1L]
      positional <- setdiff(bound, named[nzchar(named)]) # bound by position
      seen  <- c(seen, nm)
      expect_false(any(positional %in% c("lower.tail", "log.p", "log", "max.value")),
                   label = paste0(f, ": ", paste(deparse(cl, width.cutoff = 500L), collapse = " ")))
    }
  }
  expect_true(all(c("qNBI", "qPIG", "qGEOM") %in% seen)) # the sweep found the calls
  # qGPO of gamlss.dist 6.1-11 returns Poisson quantiles; the package uses my_qGPO()
  expect_false("qGPO" %in% seen)
  # every gamlss.dist function called is imported, not found through the search path
  expect_identical(setdiff(unique(seen), ls(parent.env(asNamespace("workHORSEmisc")))), character())
})

test_that("COPD duration statement in init_prevalence() runs and is geometric", {
  stmt <- Filter(function(s) "qGEOM" %in% all.names(s),
                 as.list(body(workHORSEmisc::init_prevalence))[-1L])
  expect_length(stmt, 1L)
  dt <- data.table::CJ(age = 30:89, sex = c("men", "women"), rep = 1:200)
  dt[, mu := 9 + age / 4]                                 # copd table mu is 9.0-31.4
  dt[, copd_prvl := rep(0:1, length.out = .N)]
  env <- new.env(parent = asNamespace("workHORSEmisc"))   # resolve symbols as the package does
  env$dt <- dt
  dqrng::dqset.seed(42L)
  expect_no_condition(eval(stmt[[1L]], env))
  expect_type(dt$copd_prvl, "integer")
  expect_false(anyNA(dt$copd_prvl))
  expect_true(all(dt[rep(c(TRUE, FALSE), length.out = .N), copd_prvl] == 0L))
  prv <- dt[rep(c(FALSE, TRUE), length.out = .N)]
  expect_gte(min(prv$copd_prvl), 2L)
  # 2 + Geom(mean = mu - 2): mean mu, variance (mu - 2)(mu - 1)
  z <- (mean(prv$copd_prvl) - mean(prv$mu)) /
    sqrt(sum((prv$mu - 2) * (prv$mu - 1))) * nrow(prv)
  expect_lt(abs(z), 4)
})
