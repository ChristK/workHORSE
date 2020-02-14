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
#include <Rmath.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
IntegerVector carry_forward(const IntegerVector& x, const LogicalVector& pid, const int& y) {
  const int n = x.size();
  IntegerVector out = clone(x);
  for (int i = 0; i < n; i++)
  {
    if (!pid[i] && out[i-1] == y) out[i] = y;
  }
  return out;
}

//' @export
// [[Rcpp::export]]
IntegerVector carry_backward(const IntegerVector& x, const LogicalVector& pid) {
  const int n = x.size();
  IntegerVector out = clone(x);
  for (int i = n - 1; i > 0; i--) // Go backwards bat strt from one row before the last
  {
    if (!pid[i] && out[i] > 0) out[i - 1] = out[i] - 1;
    if (i < n && pid[i] && !pid[i + 1] && out[i + 1] > 0 && out[i] == 0) out[i] = out[i + 1] - 1;
    if (out[i] < 0) out[i] = 0;
  }
  return out;
}

//' @export
// [[Rcpp::export]]
LogicalVector mk_new_simulant_markers(const IntegerVector& pid)
{
  // pid should be sorted and same length as x
  const int n = pid.size();
  LogicalVector new_simulant_markers(n);
  new_simulant_markers[0] = true;
  int previous_pid = pid[0];
  // Loop with no conditional branches in the body (therefore branch predictor should get it right almost every time) and minimal memory access by retaining previous_pid.
  for (int i = 1; i < n; i++)
  {
    new_simulant_markers[i] = pid[i] != previous_pid;
    previous_pid = pid[i];
  }
  return new_simulant_markers;
}

//' @export
// [[Rcpp::export]]
LogicalVector identify_longdeads(const IntegerVector& x, const LogicalVector& pid) {
  const int n = x.size();
  LogicalVector out(n);
  for (int i = 0; i < n; i++)
  {
    if (!pid[i] && x[i-1] != 0) out[i] = true; else out[i] = false;
  }
  return out;
}

//' @export
// [[Rcpp::export]]
IntegerVector identify_invitees(const IntegerVector& elig,
                                const NumericVector& prb,
                                const IntegerVector& freq,
                                const LogicalVector& pid
) {
  const int n = elig.size();
  int counter = 1;
  IntegerVector out(n, 0);
  for (int i = 0; i < n; i++)
  {
    if (elig[i] == 1) out[i] = R::rbinom(1.0, prb[i]);
    if (out[i] == 1)
    {
      counter = 1;
      while (counter <=i && counter < freq[i] && !pid[i-counter] && out[i] == 1)
      {
        if (out[i-counter] == 1) out[i] = 0; // had HC recently
        counter++;
      }
    }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
IntegerVector hc_effect(const IntegerVector& x,
                             const double& prb_of_continuation,
                             const LogicalVector& pid)
  {
  const int n = x.size();
  IntegerVector out = clone(x);
  for (int i = 0; i < n; i++)
  {
    if (!pid[i] && out[i-1] == 1) out[i] = R::rbinom(1, prb_of_continuation);
  }
  return out;
}


//quantile implementation for default R method (type 7)
//' @export
// [[Rcpp::export]]
NumericVector fquantile(NumericVector x,
  NumericVector probs,
  bool na_rm = true) {
  if (all(is_na(x)))
  {
    NumericVector out(probs.size(), NA_REAL);
    return(out);
  }
  if (na_rm) x = na_omit(x);
  const int n = x.size();
  NumericVector out(probs.size());
  IntegerVector ii(probs.size());
  NumericVector h(probs.size());
  NumericVector index = 1 + (n - 1) * probs;
  NumericVector lo = floor(index); //floor
  //ceiling
  NumericVector hi = ceiling(index);
  for(int i = 0; i < probs.size(); i++)
  {//catch corner case when index element is int and ceiling = floor
    h[i] = index[i] - lo[i];
  }
  std::sort(x.begin(), x.end());
  out = x[as<IntegerVector>(lo) - 1];
  if (all(is_na(ii)))
  {
    return(out);
  } else
  {
    x = x[as<IntegerVector>(hi) - 1];
    for(int i = 0; i < probs.size(); i++)
    {
      out[i] = (1 - h[i]) * out[i] + h[i] * x[i];
    }
    return(out);
  }
}

//' @export
// [[Rcpp::export]]
List fquantile_byid(NumericVector x,
  NumericVector q,
  StringVector id,
  bool rounding = false,
  bool na_rm = true) {
  // Need to be sorted by id
  const int n = x.size();
  const int m = unique(id).size();
  NumericMatrix z(m, q.size());
  StringVector id_nam(m);
  int start = 0;
  int counter = 0;
  int end = 0;
  int counter_row = 0;

  for (int i = 1; i < n; i++) { // start from 2nd element
    if (id[i] == id[i-1]) counter++;
    else
    {
      start = i - counter - 1;
      end = i - 1;
      counter = 0;
      if (rounding) z.row(counter_row) = round(fquantile(x[seq(start, end)], q, na_rm), 0);
      else z.row(counter_row) = fquantile(x[seq(start, end)], q, na_rm);
      id_nam[counter_row] = id[end];
      counter_row++;
    }
  }
  // take care the last group
  if (rounding) z.row(counter_row) = round(fquantile(x[seq(n-counter-1, n)], q, na_rm), 0);
  else z.row(counter_row) = fquantile(x[seq(n-counter-1, n)], q, na_rm);
  id_nam[counter_row] = id[n-counter-1];

  // return(z);
  const int tt = 1 + q.size();
  List outputList(tt);
  outputList[0] = id_nam;
  for (int i = 1; i < tt; i++) {
    outputList[i] = z(_, i - 1);
  }

  return outputList;
}

//' @export
// [[Rcpp::export]]
NumericVector fbound(const NumericVector &x, NumericVector &a, NumericVector &b) {
  const int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; i++) {
    if (NumericVector::is_na(x[i])) out[i] = NA_REAL;
    else
    {
      if (a[i] > b[i]) {double c = a[i]; a[i] = b[i]; b[i] = c;}; // ensure a < b

      if (x[i] < a[i]) out[i] = a[i];
      else if (x[i] > b[i]) out[i] = b[i];
      else out[i] = x[i];
    }
  }
  return out;
}
