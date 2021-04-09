/* IMPACTncd: A decision support tool for primary prevention of NCDs
 Copyright (c) 2015-2021 Chris Kypridemos

 IMPACTncd is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, see <http://www.gnu.org/licenses/>
 or write to the Free Software Foundation, Inc., 51 Franklin Street,
 Fifth Floor, Boston, MA 02110-1301  USA.*/

#include <memory>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

/*
 * This file contains the high-level structure of a simulation:
 * - the notion of a Processable, with different phases of processing;
 * - the different high-level types of Processable;
 * - the notion of a Disease, with incidence, diagnosis, and mortality;
 * - the notion of an Influencer, with the influencing Disease, a lag (in years), and multipliers;
 * - the notion of a simulation of multiple diseases.
 */

/**
 * An abstract superclass for anything that can calculate some state of a simulant in a given year.
 */
class Processable
{
  public:
    virtual ~Processable(){}
    //* simulant_year_index is known to be an appropriate year to process (i.e. >= init_year).  Do whatever the Processable needs to do before anything runs its main phase.
    virtual void preprocess(unsigned int simulant_year_index, bool is_new_simulant) = 0;

    //* simulant_year_index is known to be an appropriate year to process (i.e. >= init_year).  Do whatever the Processable needs to do in its main phase.
    virtual void process(unsigned int simulant_year_index, bool is_new_simulant) = 0;

    //* simulant_year_index is known to be an appropriate year to process (i.e. >= init_year).  Do whatever the Processable needs to do after everything has run its main phase.
    virtual void postprocess(unsigned int simulant_year_index, bool is_new_simulant) = 0;

    //* Use knowledge of the C++ to R naming conventions to add a vector to output_frame for each output variable in the given Processable.
    virtual void collect_outputs(DataFrame &output_frame, const string disease_name) = 0;

    virtual bool needs_preprocess() = 0;
    virtual bool needs_process() = 0;
    virtual bool needs_postprocess() = 0;
};

class Incidence: public Processable
{
  protected:
    IntegerVector _out_prevalence;
    const bool _can_recur;
  public:
    IntegerVector &out_prevalence()
    {
      return _out_prevalence;
    }

    Incidence(int count, bool can_recur)
      : _out_prevalence(count)
      , _can_recur(can_recur)
    {
    }

    /* Note: Do NOT use IntegerVector &out here; clone() doesn't like its output being passing as a reference parameter, and that's used in some calls from subclasses. PJC 2019-04-11 */
    Incidence(IntegerVector out, bool can_recur)
      : _out_prevalence(out)
      , _can_recur(can_recur)
    {
    }


    //* Set to true to mark a disease as cured; it may recur depending on the value of can_recur.
    virtual void cured(unsigned int simulant_year_index, bool value)
    {
      // By default, this information is not used in incidences
    }

    //* Most Incidences need no postprocessing; for convenience, default it.
    virtual void postprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
    }

    virtual bool needs_preprocess() { return true; }
    virtual bool needs_process() { return true; }
    virtual bool needs_postprocess() { return false; }

    virtual void collect_outputs(DataFrame &output_frame, const string disease_name)
    {
      string prevalence_name = disease_name + "_prvl";
      output_frame.push_back(_out_prevalence, prevalence_name);
    }
};

class Diagnosis: public Processable
{
  protected:
    IntegerVector _out_diagnosis;
  public:
    IntegerVector &out_diagnosis()
    {
      return _out_diagnosis;
    }

    Diagnosis(unsigned int count)
      : _out_diagnosis(count)
    {
    }

    Diagnosis(IntegerVector out_diagnosis)
      : _out_diagnosis(out_diagnosis)
    {
    }

    //* Set to true to mark a disease as cured; it may recur depending on the value of can_recur.
    virtual void cured(unsigned int simulant_year_index, bool value)
    {
      // By default, this information is not used in incidences
    }

    virtual void preprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
      // Do nothing
    }

    virtual void postprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
      // Do nothing
    }

    virtual bool needs_preprocess() { return false; }
    virtual bool needs_process() { return true; }
    virtual bool needs_postprocess() { return false; }

    virtual void collect_outputs(DataFrame &output_frame, const string disease_name)
    {
      string diagnosis_name = disease_name + "_dgn";
      output_frame.push_back(_out_diagnosis, diagnosis_name);
    }
};

class Disease;

class Mortality: public Processable
{
  protected:
    IntegerVector _out_mortality;
  public:
    const IntegerVector &out_mortality() const
    {
      return _out_mortality;
    }

    Mortality(unsigned int count)
      : _out_mortality(count)
    {
    }

    Mortality(IntegerVector out_mortality)
      : _out_mortality(out_mortality)
    {
    }

    virtual void disease(Disease &value)
    {
      // Most mortalities don't use this information
    }

    //* Most Mortalities need no preprocessing; for convenience, default it.
    virtual void preprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
    }

    //* Most Mortalities need no postprocessing; for convenience, default it.
    virtual void postprocess(unsigned int simulant_year_index, bool is_new_simulant)
    {
    }

    virtual bool needs_preprocess() { return false; }
    virtual bool needs_process() { return true; }
    virtual bool needs_postprocess() { return false; }

    virtual void collect_outputs(DataFrame &output_frame, const string disease_name)
    {
      string mortality_name = disease_name + "_mrtl";
      output_frame.push_back(_out_mortality, mortality_name);
    }
};

class Disease
{
  private:
    unique_ptr<Incidence> _incidence;
    unique_ptr<Diagnosis> _diagnosis;
    unique_ptr<Mortality> _mortality;
    const string _name;
  public:
    Disease(string name, unique_ptr<Incidence> &incidence, unique_ptr<Diagnosis> &diagnosis, unique_ptr<Mortality> &mortality)
      : _incidence(move(incidence))
      , _diagnosis(move(diagnosis))
      , _mortality(move(mortality))
      , _name(name)
    {
    }

    inline const string name() const
    {
      return _name;
    }

    inline Incidence& incidence() const
    {
      return *_incidence;
    }

    inline Diagnosis& diagnosis() const
    {
      return *_diagnosis;
    }

    inline Mortality& mortality() const
    {
      return *_mortality;
    }
};

class Influencer
{
  private:
    const Disease *_disease;
    int _lag;
    NumericVector _multiplier;
  public:

    Influencer(const Disease &disease, const NumericVector &multiplier, const int lag)
    {
      _disease = &disease;
      _multiplier = multiplier;
      _lag = lag;
    }

    inline const Disease &disease() const
    {
      return *_disease;
    }

    inline const int &lag() const
    {
      return _lag;
    }

    inline const NumericVector &multiplier() const
    {
      return _multiplier;
    }

    const double multiplier_for(const int simulant_year_index) const
    {
      return _disease->incidence().out_prevalence()[simulant_year_index] < _lag
        ? 1.0
        : _multiplier[simulant_year_index];
    }
};

class Simulation
{
  private:
    const DataFrame* _frame;
    string _scenario_name;
    vector<string> _strata_names;
  protected:
    void append_strata(DataFrame &output_frame) const
    {
      for (auto const& stratum_name : _strata_names)
        output_frame.push_back((*_frame)[stratum_name], stratum_name);
    }
  public:
    // pid will be sorted by year. TODO force and check in R side
    IntegerVector pids;
    IntegerVector years;
    int init_year;
    NumericVector rn_all_cause_mrtl;
    unique_ptr<unordered_map<string, unique_ptr<Disease>>> diseases;

    const DataFrame& frame() const { return *_frame; }
    void frame(const DataFrame* value)
    {
      _frame = value;
    }

    Simulation()
      : init_year(0)
    {
    }

    virtual ~Simulation(){}
    void scenario_name(const string &scenario_name)
    {
      _scenario_name = scenario_name;
    }

    const string scenario_name() const
    {
      return _scenario_name;
    }

    int simulant_year_count() const
    {
      return pids.length();
    }

    void add_stratum(const string stratum_name)
    {
      _strata_names.push_back(stratum_name);
    }

    /**
     * Run the simulation to completion. The simulation's concrete subclass is free to take any approach it wishes to fulfil this request:
     * run single-threaded, multi-threaded, on a GPU if it can find one, distribute across a cluster, find a quantum computer, use clairvoyance...
     *
     * Postcondition: run() threw an exception *or* the simulation ran to successful completion and all temporary resources have been cleaned up.
     */
    virtual void run() = 0;

    /**
     * Gather output vectors from all parts of the simulation and return these in a DataFrame.
     *
     * Precondition: Simulation has run to successful completion.
     */
    virtual DataFrame output_frame() const = 0;
};

ostream& operator<< (ostream& os, const Simulation& ss) {
  os << "Simulation("
    << "pids: " << ss.pids.length() << ", "
    << "years: " << ss.years.length() << ", "
    << "init_year: " << ss.init_year << ", "
    << "diseases: " << ss.diseases->size()
    << ')';
  return os;
}
