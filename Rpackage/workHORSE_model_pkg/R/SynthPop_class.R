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
`[.SynthPop` = function(x, ...)
  x$pop[...]


#' R6 Class representing a synthetic population
#'
#' A synthpop has a pop data.table that has the synthetic population.
#' @export
SynthPop <-
  R6Class("SynthPop",

          public = list(
            design = NA,
            pop = NA,


            initialize = function(l) {
              stopifnot("Design" %in% class(l))
              self$design <- l

              invisible(self)
            },

            print = function(...) {
              print(self$pop)

              invisible(self)
            }


          )


          )

