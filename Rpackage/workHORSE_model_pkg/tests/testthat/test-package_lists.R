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

# The foreach workers attach every package listed in .packages, and a package
# that cannot be loaded (e.g. qs, archived on CRAN, on R >= 4.6) fails the whole
# parallel run with a FutureLaunchError. Guards the .packages vectors in the
# package and the packages listed in the repo's dependencies.yaml.

# Calls to `name` (as language objects) anywhere inside an expression, also when
# written as pkg::name
find_calls <- function(e, name, acc = list()) {
  if (!is.call(e)) return(acc)
  fn <- e[[1L]]
  if (is.call(fn) && length(fn) == 3L && identical(fn[[1L]], quote(`::`))) fn <- fn[[3L]]
  if (is.name(fn) && identical(as.character(fn), name)) acc[[length(acc) + 1L]] <- e
  args <- as.list(e)
  for (i in seq_along(args)) {
    if (identical(args[[i]], quote(expr = ))) next # empty arg
    acc <- find_calls(args[[i]], name, acc)
  }
  acc
}

# Every function in the namespace plus the public/private/active methods of R6 generators
package_functions <- function(pkg = "workHORSEmisc") {
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

# The .packages vector of every foreach() call in the package, with its location
foreach_packages <- function() {
  out <- list()
  fns <- package_functions()
  for (f in names(fns)) {
    b <- body(fns[[f]])
    if (is.null(b)) next
    for (cl in find_calls(b, "foreach")) {
      pk <- as.list(cl)[[".packages"]]
      if (is.null(pk)) next
      pk <- tryCatch(eval(pk, baseenv()), # a c("a", "b") literal
        error = function(e) stop(f, ": .packages is not a literal vector: ", conditionMessage(e)))
      out[[length(out) + 1L]] <- list(where = f, packages = pk)
    }
  }
  out
}

# Nearest directory above the tests that holds the repo's dependencies.yaml
# (next to global.R); NULL when the package is tested outside the repo
find_repo_root <- function(from = getwd()) {
  d <- normalizePath(from, mustWork = FALSE)
  repeat {
    if (all(file.exists(file.path(d, c("dependencies.yaml", "global.R"))))) return(d)
    up <- dirname(d)
    if (identical(up, d)) return(NULL)
    d <- up
  }
}

test_that("foreach .packages in the package do not list qs and can all be loaded", {
  found <- foreach_packages()
  # the sweep found the calls: one in run_simulation, %do% and %dopar% in write_synthpop
  expect_true(all(c("run_simulation", "SynthPop$write_synthpop") %in%
      vapply(found, `[[`, "", "where")))
  expect_gte(length(found), 3L)
  for (x in found) {
    expect_type(x$packages, "character")
    expect_false("qs" %in% x$packages, label = paste0(x$where, ": .packages lists qs"))
    ok <- vapply(x$packages, requireNamespace, NA, quietly = TRUE)
    expect_true(all(ok),
      label = paste0(x$where, ": cannot load ", paste(x$packages[!ok], collapse = ", ")))
  }
})

test_that("dependencies.yaml does not list qs and every package in it can be loaded", {
  root <- find_repo_root()
  skip_if_not(!is.null(root), "dependencies.yaml not found above the tests")
  deps <- unlist(yaml::read_yaml(file.path(root, "dependencies.yaml")), use.names = FALSE)
  expect_type(deps, "character")
  expect_gt(length(deps), 0L)
  expect_false("qs" %in% deps)
  ok <- vapply(deps, requireNamespace, NA, quietly = TRUE)
  expect_true(all(ok), label = paste0("cannot load ", paste(deps[!ok], collapse = ", ")))
})
