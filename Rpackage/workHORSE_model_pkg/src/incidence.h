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
 * A convenient incidence that does nothing.
 * As C++ initialises vectors to 0, this means the results will always be 0.
 */
class NullIncidence: public Incidence
{
  public:
    NullIncidence(unsigned int count)
      : Incidence(count, false)
    {
    }


    void preprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
    }

    void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
    }

    virtual void collect_outputs(DataFrame &output_frame, const string disease_name)
    {
      // No relevant outputs
    }

    virtual bool needs_preprocess() { return false; }
    virtual bool needs_process() { return false; }
};

/**
 * Everyone gets it, every year!
 */
class UniversalIncidence: public Incidence
{
  public:
    UniversalIncidence(unsigned int count)
      : Incidence(count, true)
    {
    }

    void preprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
      _out_prevalence[simulant_year_index] = 1;
    }

    void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
    }

    virtual void collect_outputs(DataFrame &output_frame, const string disease_name)
    {
      // No relevant outputs
    }

    virtual bool needs_process() { return false; }
};

class IncidenceType1: public Incidence
{
  private:
    const NumericVector _prob;
    const IntegerVector _prevalence;
  public:
    IncidenceType1(const NumericVector &prob, const IntegerVector &prevalence)
      : Incidence(clone(prevalence), true)
      , _prob(prob)
      , _prevalence(prevalence)
    {
    }

    void preprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
      if (!is_new_simulant)
        _out_prevalence[simulant_year_index] = _out_prevalence[simulant_year_index - 1];
    }

    void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
      // Carry forward has already happened
      if (_prob[simulant_year_index] == 1)
        _out_prevalence[simulant_year_index]++;
    }
};

class IncidenceType2: public Incidence
{
  private:
    const NumericVector _prob;
    const NumericVector _rn;
    const IntegerVector _prevalence;
    bool _cured;
  public:
    IncidenceType2(const NumericVector &prob, const NumericVector &rn, const IntegerVector &prevalence, const bool can_recur)
      : Incidence(clone(prevalence), can_recur)
      , _prob(prob)
      , _rn(rn)
      , _prevalence(prevalence)
    {
    }

    void cured(unsigned int simulant_year_index, bool value)
    {
      _cured = value;
      _out_prevalence[simulant_year_index] = 0;
    }

    void preprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
      if (is_new_simulant)
      {
        _cured = false;
        // Previous generation of code added 1 to duration of type 2 conditions when they were part of the first processed year, so we do.
        if (_out_prevalence[simulant_year_index] > 0)
          _out_prevalence[simulant_year_index]++;
      }
      else
      {
        if (_out_prevalence[simulant_year_index - 1] > 0 && !_cured)
          _out_prevalence[simulant_year_index] = _out_prevalence[simulant_year_index - 1] + 1;
      }
    }

    void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
      if (_cured && !_can_recur)
        return;

      // Carry forward has already happened
      if (_out_prevalence[simulant_year_index] == 0 && _prob[simulant_year_index] > _rn[simulant_year_index])
        _out_prevalence[simulant_year_index]++;
    }
};

class IncidenceType3: public Incidence
{
  private:
    const NumericVector _prob;
    const NumericVector _rn;
    const IntegerVector _prevalence;

    const Disease *_influencer;
    int _influencer_lag;
    NumericVector _influencer_multiplier;
  public:
    IncidenceType3(const NumericVector &prob, const NumericVector &rn, const IntegerVector &prevalence)
      : Incidence(clone(prevalence), true)
      , _prob(prob)
      , _rn(rn)
      , _prevalence(prevalence)
    {
    }

    void cured(unsigned int simulant_year_index, bool value)
    {
      _out_prevalence[simulant_year_index] = 0;
    }


    void set_influencer(const Disease &influencer, const NumericVector &multiplier, const int lag)
    {
      _influencer = &influencer;
      _influencer_multiplier = multiplier;
      _influencer_lag = lag;
    }

    void preprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
      if ((!is_new_simulant) && _out_prevalence[simulant_year_index - 1] > 0)
        _out_prevalence[simulant_year_index] = _out_prevalence[simulant_year_index - 1] + 1;
    }

    void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
      // Carry forward has already happened
      if (0 == _out_prevalence[simulant_year_index])
      {
        double used_multiplier = _influencer->incidence().out_prevalence()[simulant_year_index] < _influencer_lag
          ? 1.0
          : _influencer_multiplier[simulant_year_index];
        if ((_prob[simulant_year_index] * used_multiplier) > _rn[simulant_year_index])
          _out_prevalence[simulant_year_index] = 1;
      }
    }
};

/**
 * A disease that is present if any of a list of other diseases is present.
 *  For example, a simulant has CVD if they have either stroke or CHD.
 */
class IncidenceType4: public Incidence
{
  private:
    vector<const Disease*> _prerequisites;
  public:
    IncidenceType4(unsigned int simulant_year_count)
      : Incidence(simulant_year_count, true)
    {
    }

    void prerequisites(const vector<const Disease*> &prerequisites)
    {
      _prerequisites = prerequisites;
    }

    void preprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
      if (is_new_simulant)
        _out_prevalence[simulant_year_index] = 0;
      else if (0 == _out_prevalence[simulant_year_index - 1])
        {
        unsigned int previous_simulant_year_index = simulant_year_index - 1;
        int maxSoFar = _out_prevalence[previous_simulant_year_index];
        for (unsigned int prerequisite_index = 0; prerequisite_index < _prerequisites.size(); prerequisite_index++)
          if (_prerequisites[prerequisite_index]->incidence().out_prevalence()[previous_simulant_year_index] > maxSoFar)
            maxSoFar = _prerequisites[prerequisite_index]->incidence().out_prevalence()[previous_simulant_year_index];
        _out_prevalence[simulant_year_index] = 0 == maxSoFar ? 0 : maxSoFar + 1;
      }
      else
        _out_prevalence[simulant_year_index] = _out_prevalence[simulant_year_index - 1] + 1;
    }

    // No processing needed
    void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
    }

    void postprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
      int maxSoFar = _out_prevalence[simulant_year_index];
      for (unsigned int prerequisite_index = 0; prerequisite_index < _prerequisites.size(); prerequisite_index++)
        if (_prerequisites[prerequisite_index]->incidence().out_prevalence()[simulant_year_index] > maxSoFar)
          maxSoFar = _prerequisites[prerequisite_index]->incidence().out_prevalence()[simulant_year_index];
      _out_prevalence[simulant_year_index] = maxSoFar;
    }

    virtual bool needs_preprocess() { return true; }
    virtual bool needs_process() { return false; }
    virtual bool needs_postprocess() { return true; }
};


class IncidenceType5: public Incidence
{
private:
  const NumericVector _prob;
  const NumericVector _rn;
  const IntegerVector _prevalence;

  const Disease *_influencer;
  int _influencer_lag;
  NumericVector _influencer_multiplier;
public:
  IncidenceType5(const NumericVector &prob, const NumericVector &rn, const IntegerVector &prevalence)
    : Incidence(clone(prevalence), true)
  , _prob(prob)
  , _rn(rn)
  , _prevalence(prevalence)
  {
  }

  void set_influencer(const Disease &influencer, const NumericVector &multiplier, const int lag)
  {
    _influencer = &influencer;
    _influencer_multiplier = multiplier;
    _influencer_lag = lag;
  }

  void preprocess(unsigned int simulant_year_index, bool is_new_simulant)
  {
    if ((!is_new_simulant) && _out_prevalence[simulant_year_index - 1] > 0)
      _out_prevalence[simulant_year_index] = _out_prevalence[simulant_year_index - 1] + 1;
  }

  void process(unsigned int simulant_year_index, bool is_new_simulant)
  {
    // Carry forward has already happened
    if (0 == _out_prevalence[simulant_year_index])
    {
      double used_multiplier = _influencer->incidence().out_prevalence()[simulant_year_index] == _influencer_lag
      ? _influencer_multiplier[simulant_year_index]
      : 0.0;
      if ((_prob[simulant_year_index] * used_multiplier) > _rn[simulant_year_index])
        _out_prevalence[simulant_year_index] = 1;
    }
  }
};
