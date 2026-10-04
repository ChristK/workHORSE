/* workHORSE is an implementation of the IMPACTncd framework, developed by Chris
 Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
 Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
 funded by NIHR  HTA Project: 16/165/01 - workHORSE: Health Outcomes
 Research Simulation Environment.  The views expressed are those of the
 authors and not necessarily those of the NHS, the NIHR or the Department of
 Health.

 Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos

 workHORSE is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 3 of the License, or (at your option) any later
 version. This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 details. You should have received a copy of the GNU General Public License
 along with this program; if not, see <http://www.gnu.org/licenses/> or write
 to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 Boston, MA 02110-1301 USA. */

#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>
#include <vector>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// The cdf below is accumulated with the very same IEEE double operations that
// R performs, so that the results are identical to the ones of gamlss.dist
// 6.1-1. Compilers must therefore not fuse the multiply-adds, as they do under
// flags such as -march=native.
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

// qGPO

//' Quantile function of the generalised Poisson distribution
//'
//' Replaces \code{gamlss.dist::qGPO}, which is broken in gamlss.dist 6.1-11
//' (\code{dGPO} returns the Poisson density whenever \code{sigma > 1e-06}, so
//' that \code{qGPO} returns Poisson quantiles). \code{my_qGPO} reproduces the
//' results of \code{qGPO} from gamlss.dist 6.1-1, the last release with a
//' correct generalised Poisson, but is much faster. The quantile is the
//' smallest integer \code{x} for which the cdf is at least \code{p}. The cdf is
//' the sum of the pmf up to \code{x}, except for \code{sigma < 1e-04} where it
//' is the Poisson cdf.
//'
//' @param p vector of probabilities.
//' @param mu,sigma vectors of the parameters of the distribution. Both must be
//'   positive and finite. \code{mu} is the mean and the variance is
//'   \code{mu * (1 + sigma * mu)^2}. Vectors of length one are recycled.
//' @param lower_tail,log_p as in \code{\link[stats]{qpois}}.
//' @param max_value the largest quantile that can be returned.
//' @param n_cpu number of OpenMP threads, if the package is built with OpenMP.
//'   It does not affect the results.
//' @return A numeric vector of non-negative whole numbers, \code{Inf} for
//'   \code{p} within 1e-09 of 1, as \code{gamlss.dist::qGPO} does.
//' @export
// [[Rcpp::export]]
NumericVector my_qGPO(const NumericVector& p,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& max_value = 10000,
                      const int& n_cpu = 1)
{
  const R_xlen_t n   = p.length();
  const R_xlen_t nmu = mu.length();
  const R_xlen_t nsg = sigma.length();
  if ((nmu != n && nmu != 1) || (nsg != n && nsg != 1))
    stop("Distribution parameters must be of same length (or of length 1)");
  if (max_value < 0) stop("max_value must be >= 0");

  const double* pp = p.begin();
  const double* pm = mu.begin();
  const double* ps = sigma.begin();
  NumericVector out(n); // holds the probabilities first, then the quantiles
  double* po = out.begin();

  // The tests are negated so that NA and NaN are rejected too
  for (R_xlen_t i = 0; i < nmu; i++)
    if (!(pm[i] > 0.0) || !R_FINITE(pm[i])) stop("mu must be greater than 0 and finite");
  for (R_xlen_t i = 0; i < nsg; i++)
    if (!(ps[i] > 0.0) || !R_FINITE(ps[i])) stop("sigma must be greater than 0 and finite");
  for (R_xlen_t i = 0; i < n; i++)
  {
    double prob = pp[i];
    if (log_p) prob = exp(prob);
    if (!lower_tail) prob = 1.0 - prob;
    if (!(prob >= 0.0 && prob <= 1.0001)) stop("p must be between 0 and 1");
    po[i] = prob;
  }

#ifdef _OPENMP
  const int n_threads = (n_cpu < 1) ? 1 : n_cpu;
#endif

  // No R API calls (and no exceptions) are allowed in the parallel region
#pragma omp parallel num_threads(n_threads)
  {
    std::vector<double> lgam; // lgamma(j + 1), shared by the draws of this thread
#pragma omp for schedule(dynamic, 256)
    for (R_xlen_t i = 0; i < n; i++)
    {
      const double prob = po[i];
      if (prob + 1e-09 >= 1.0) // as in gamlss.dist
      {
        po[i] = R_PosInf;
        continue;
      }

      const double mui = pm[(nmu == 1) ? 0 : i];
      const double sgi = ps[(nsg == 1) ? 0 : i];
      double q = max_value; // the result when the search reaches max_value
      if (sgi < 1e-04)
      {
        // Poisson limit: pGPO() returns the Poisson cdf
        for (int j = 0; ; j++)
        {
          if (prob <= R::ppois(j, mui, 1, 0))
          {
            q = j;
            break;
          }
          if (j == max_value) break;
        }
      }
      else
      {
        // dGPO() evaluated in the same order as in R. R sums the pmf from
        // x = 0 in long double precision and rounds the sum to double
        const double den = 1.0 + sgi * mui;
        const double lth = log(mui / den);
        long double  cum = 0.0L;
        for (int j = 0; ; j++)
        {
          if ((int) lgam.size() == j) lgam.push_back(R::lgammafn(j + 1.0));
          const double x = j;
          const double s = 1.0 + sgi * x;
          const double logL = x * lth + (x - 1.0) * log(s) + (-mui * s) / den - lgam[j];
          cum += exp(logL);
          const double cdf = (double) cum;
          if (cdf != cdf) // NaN, because of absurd values of mu and sigma
          {
            q = R_NaN;
            break;
          }
          if (prob <= cdf)
          {
            q = j;
            break;
          }
          if (j == max_value) break;
        }
      }
      po[i] = q;
    }
  }

  for (R_xlen_t i = 0; i < n; i++)
    if (po[i] != po[i])
    {
      warning("NaNs or NAs were produced");
      break;
    }
  return out;
}

/*** R
# Checked against gamlss.dist 6.1-1 (6.1-11 is broken), run in an R session
# with the library of 6.1-1 first in .libPaths():
# N <- 1e3
# p <- runif(N)
# mu <- runif(N, 0.5, 30)
# sigma <- runif(N, 0.001, 2)
# identical(my_qGPO(p, mu, sigma), gamlss.dist::qGPO(p, mu, sigma))
# library(microbenchmark)
# microbenchmark(my_qGPO(p, mu, sigma), gamlss.dist::qGPO(p, mu, sigma), times = 3)
*/
