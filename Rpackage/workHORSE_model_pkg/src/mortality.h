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

using namespace Rcpp;
using namespace std;

/**
 * A convenient mortality that does nothing - it never kills.
 */
class NullMortality: public Mortality
{
  public:
    NullMortality(unsigned int count)
      : Mortality(count)
    {
    }

    virtual void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
      _out_mortality[simulant_year_index] = 0;
    }

    virtual void collect_outputs(DataFrame &output_frame, const string disease_name)
    {
      // No relevant outputs
    }
};

class MortalityType1: public Mortality
{
  private:
    const NumericVector _prob;
    const NumericVector _rn;
    const IntegerVector _prevalence;
    const int _code;
  public:
    MortalityType1(const NumericVector &prob, const NumericVector &rn, const IntegerVector &prevalence, const int cod)
      : Mortality(prevalence.size())
      , _prob(prob)
      , _rn(rn)
      , _prevalence(prevalence)
      , _code(cod)
    {
    }

    void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
      if (is_new_simulant)
      {
        if (_prevalence[simulant_year_index] > 0 && _prob[simulant_year_index] > _rn[simulant_year_index])
          _out_mortality[simulant_year_index] = _code;
      }
      else
      {
        if (_out_mortality[simulant_year_index - 1] == _code)
          _out_mortality[simulant_year_index] = _code;
        else if (_prevalence[simulant_year_index] > 0 && _prob[simulant_year_index] > _rn[simulant_year_index])
          _out_mortality[simulant_year_index] = _code;
      }
    }
};

class MortalityType2: public Mortality
{
  private:
    const NumericVector _prob;
    const NumericVector _rn;
    IntegerVector &_prevalence;
    const int _code;
    const int _cured_after_years;
    Disease *_disease;

    void cure(unsigned int simulant_year_index)
    {
      _disease->incidence().cured(simulant_year_index, true);
      _disease->diagnosis().cured(simulant_year_index, true);
    }
  public:
    MortalityType2(const NumericVector &prob, const NumericVector &rn, IntegerVector &prevalence, const int cod, const int cured_after_years)
      : Mortality(prevalence.size())
      , _prob(prob)
      , _rn(rn)
      , _prevalence(prevalence)
      , _code(cod)
      , _cured_after_years(cured_after_years)
    {
    }

    void disease(Disease &value)
    {
      _disease = &value;
    }

    void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
      if (is_new_simulant)
      {
        if (_prevalence[simulant_year_index] > 0 && _prevalence[simulant_year_index] <= _cured_after_years && _prob[simulant_year_index] > _rn[simulant_year_index])
          _out_mortality[simulant_year_index] = _code;
        else
        { // Known: Isn't now dead
          if (_prevalence[simulant_year_index] > _cured_after_years)
            cure(simulant_year_index);
        }
      }
      else
      {
        if (_out_mortality[simulant_year_index - 1] == _code)
          _out_mortality[simulant_year_index] = _code;
        else
        { // Known: Wasn't already dead
          if (_prevalence[simulant_year_index] > 0 && _prevalence[simulant_year_index] <= _cured_after_years && _prob[simulant_year_index] > _rn[simulant_year_index])
            _out_mortality[simulant_year_index] = _code;
          else
          { // Known: Isn't now dead
            if (_prevalence[simulant_year_index] > _cured_after_years)
              cure(simulant_year_index);
          }
        }
      }
    }
};

class MortalityType3: public Mortality
{
private:
  const NumericVector _prob;
  const NumericVector _rn;
  IntegerVector &_prevalence;
  const int _code;

  const Disease *_influencer;
  int _influencer_lag;
  NumericVector _influencer_multiplier;
public:
  MortalityType3(const NumericVector &prob, const NumericVector &rn, IntegerVector &prevalence, const int cod)
    : Mortality(prevalence.size())
  , _prob(prob)
  , _rn(rn)
  , _prevalence(prevalence)
  , _code(cod)
  {
  }

  void set_influencer(const Disease &influencer, const NumericVector &multiplier, const int lag)
  {
    _influencer = &influencer;
    _influencer_multiplier = multiplier;
    _influencer_lag = lag;
  }

  void process(unsigned int simulant_year_index, bool is_new_simulant)
  {
    if (is_new_simulant)
    {
      double used_multiplier = _influencer->incidence().out_prevalence()[simulant_year_index] < _influencer_lag
      ? 1.0
      : _influencer_multiplier[simulant_year_index];
      if (_prevalence[simulant_year_index] > 0 && (_prob[simulant_year_index] * used_multiplier) > _rn[simulant_year_index])
        _out_mortality[simulant_year_index] = _code;
    }
    else
    {
      if (_out_mortality[simulant_year_index - 1] == _code)
        _out_mortality[simulant_year_index] = _code;
      else
      {
        double used_multiplier = _influencer->incidence().out_prevalence()[simulant_year_index] < _influencer_lag
        ? 1.0
        : _influencer_multiplier[simulant_year_index];
        if (_prevalence[simulant_year_index] > 0 && (_prob[simulant_year_index] * used_multiplier) > _rn[simulant_year_index])
          _out_mortality[simulant_year_index] = _code;
      }
    }
  }
};
