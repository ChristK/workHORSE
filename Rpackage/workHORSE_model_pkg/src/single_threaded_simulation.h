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

class SingleTheadedSimulation: public Simulation
{
  private:
    IntegerVector _all_cause_mrtl;
    vector<Processable*> preprocess_these;
    vector<Processable*> process_these;
    vector<Processable*> postprocess_these;

    inline int calculate_cod_for_row(
      const vector<IntegerVector> &cods,
      const NumericVector &rn,
      const int simulant_year_index) const
    {
      const int conditionCount = cods.size();
      // Aggregate the possible causes.
      int possible_causes[conditionCount];
      int next_possible_cause = 0;
      for (int cod_index = 0; cod_index < conditionCount; cod_index++)
        if (cods[cod_index][simulant_year_index] != 0)
          possible_causes[next_possible_cause++] = cods[cod_index][simulant_year_index];
      if (next_possible_cause == 0)
      {
        // No death - the very common case!
        return 0;
      }
      else if (next_possible_cause == 1)
      {
        // Unique cause of death (the commonest case when dying by far) - just put the sole cause into the output and we're done.
        return possible_causes[0];
      }
      else
      {
        // More than one cause of death; select randomly (and equally) between them.
        // TODO: This assumes rn is in the range [0,1) i.e. does not include 1.0!
        return possible_causes[(int)(rn[simulant_year_index] * next_possible_cause)];
      }
    }

    IntegerVector mortality_resolve(
      const vector<IntegerVector> &cods,
      const NumericVector &rn) const
    {
      // pid will be sorted by year. TODO force and check in R side
      const unsigned int n = pids.size();
      IntegerVector out(n);

      int previous_pid = -1;
      for (unsigned int simulant_year_index = 0; simulant_year_index < n; simulant_year_index++)
        if (years[simulant_year_index] >= init_year)  // start logic after init year
        {
          if (previous_pid != pids[simulant_year_index]) // if new simulant
          {
            previous_pid = pids[simulant_year_index];
            out[simulant_year_index] = calculate_cod_for_row(cods, rn, simulant_year_index);
          }
          else // same simulant - which guarantees simulant_year_index > 0 and hence [simulant_year_index - 1] is legal
          {
            if (out[simulant_year_index - 1] != 0)
              out[simulant_year_index] = out[simulant_year_index - 1];
            else
              out[simulant_year_index] = calculate_cod_for_row(cods, rn, simulant_year_index);
          }
        }
      return out;
    }

    void prime_processing_lists(Processable &p)
    {
        if (p.needs_preprocess())
          preprocess_these.push_back(&p);
        if (p.needs_process())
          process_these.push_back(&p);
        if (p.needs_postprocess())
          postprocess_these.push_back(&p);
    }
  public:
    void run()
    {
      // Assemble a list of all and only the calls that need to be made, in the order they need to be made.
      for (auto const& diseasePair : *diseases)
      {
        // Pass these in the order they need to be run: all incidences, then all diagnoses, then all mortalities.
        prime_processing_lists(diseasePair.second->incidence());
        prime_processing_lists(diseasePair.second->diagnosis());
        prime_processing_lists(diseasePair.second->mortality());
      }

      // Run all the phases in order
      const int n = simulant_year_count();
      int previous_pid = -1;
      for (int simulant_year_index = 0; simulant_year_index < n; simulant_year_index++)
      {
        if (years[simulant_year_index] >= init_year)  // start logic after init year
        {
          bool is_new_simulant = previous_pid != pids[simulant_year_index];
          previous_pid = pids[simulant_year_index];
          for (Processable *p : preprocess_these)
            p->preprocess(simulant_year_index, is_new_simulant);
          for (Processable *p : process_these)
            p->process(simulant_year_index, is_new_simulant);
          for (Processable *p : postprocess_these)
            p->postprocess(simulant_year_index, is_new_simulant);
        }
      }

      // Postprocess mortality: where there might be a choice, which one got 'em?
      vector<IntegerVector> mortalities;
      for (auto const& diseasePair : *diseases)
        mortalities.push_back(diseasePair.second->mortality().out_mortality());
      _all_cause_mrtl = mortality_resolve(mortalities, rn_all_cause_mrtl);
    }

    DataFrame output_frame() const
    {
      DataFrame output_frame;
      append_strata(output_frame);

      for (auto const& diseasePair : *diseases)
      {
        // Pass these in the order they need to be run: all incidences, then all diagnoses, then all mortalities.
        diseasePair.second->incidence().collect_outputs(output_frame, diseasePair.first);
        diseasePair.second->diagnosis().collect_outputs(output_frame, diseasePair.first);
        diseasePair.second->mortality().collect_outputs(output_frame, diseasePair.first);
      }
      output_frame.push_back(_all_cause_mrtl, "all_cause_mrtl");
      return output_frame;
    }
};
