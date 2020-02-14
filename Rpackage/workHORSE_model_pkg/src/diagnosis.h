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
 * A convenient diagnosis that does nothing.
 */
class NullDiagnosis: public Diagnosis
{
  public:
    NullDiagnosis(unsigned int count)
      : Diagnosis(count)
    {
    }

    virtual void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
      _out_diagnosis[simulant_year_index] = 0;
    }

    virtual void collect_outputs(DataFrame &output_frame, const string disease_name)
    {
      // No relevant outputs
    }
};

class DiagnosisType1: public Diagnosis
{
  private:
    const NumericVector _probability;
    const NumericVector _rn;
    const IntegerVector _prevalence;
  public:
    DiagnosisType1(const NumericVector &probability, const NumericVector &rn, const IntegerVector &prevalence, const IntegerVector &diagnosed)
      : Diagnosis(clone(diagnosed))
      , _probability(probability)
      , _rn(rn)
      , _prevalence(prevalence)
    {
    }

    void cured(unsigned int simulant_year_index, bool value)
    {
      _out_diagnosis[simulant_year_index] = 0;
    }

    void process(unsigned int simulant_year_index, bool is_new_simulant)
    {
      if (_prevalence[simulant_year_index] > 0)
      {
        // multiplier to progressively increase prb of diagnosis with year of diagnosis
        double mltp = _prevalence[simulant_year_index] / 5.0;

        if (is_new_simulant)
        {
          if (_out_diagnosis[simulant_year_index] == 0 && _probability[simulant_year_index] * mltp > _rn[simulant_year_index])
            _out_diagnosis[simulant_year_index] = 1;
        }
        else
        {
          if (_out_diagnosis[simulant_year_index - 1] > 0)
            _out_diagnosis[simulant_year_index] = _out_diagnosis[simulant_year_index - 1] + 1;
          else if (_out_diagnosis[simulant_year_index - 1] == 0 && _probability[simulant_year_index] * mltp > _rn[simulant_year_index])
            _out_diagnosis[simulant_year_index] = 1;
        }
      }
    }
};
