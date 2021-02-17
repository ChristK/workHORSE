/* workHORSE is an implementation of the IMPACTncd framework, developed by Chris
 Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
 Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
 funded by NIHR  HTA Project: 16/165/01 - workHORSE: Health Outcomes
 Research Simulation Environment.  The views expressed are those of the
 authors and not necessarily those of the NHS, the NIHR or the Department of
 Health.

 Copyright (C) 2018-2021 University of Liverpool, Chris Kypridemos

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


#include <memory>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

/**
 * A simple recursive-descent parser with functions that broadly match the nested list structure.
 */
class SimulationParser
{
private:
  Simulation& _simulation;
  const DataFrame& _frame;

  void friendly_error(const string message)
  {
    stop(message);
  }

  /**
   * Assume the content is a list of the form list(column_name="someName") and return someName
   * if it exists, throw an error if not.
   */
  const char* parse_column_name_or_error(const List &list)
  {
    const char* column_name = "column_name";

    // Check the name "column_name" exists in the list; if not, swear.
    if (!list.containsElementNamed(column_name))
      friendly_error("Cannot find element 'column_name' in a list that's supposed to contain a column name");
    return list[column_name];
  }

  /**
   * Look for innerListName as a named element in outerList. If it doesn't exist, swear. If it does exist, assume the content is a list of the form list(column_name="someName")
   * and return the relevant vector if it exists, throw an error if not.
   */
  const IntegerVector parse_integer_column_or_error(const List &outerList, const char *innerListName, const char *context_name)
  {
    // Check the name exists in the outer list; if not, swear.
    if (!outerList.containsElementNamed(innerListName))
    {
      ostringstream ss;
      ss << "Cannot find element '" << innerListName << "' in list";
      if (0 != context_name)
        ss << " for " << context_name;
      friendly_error(ss.str());
    }

    const char* column_name = parse_column_name_or_error(outerList[innerListName]);

    // Check the column exists in the frame; if not, swear.
    if (!_frame.containsElementNamed(column_name))
    {
      ostringstream ss;
      ss << "Cannot find column '" << column_name << "' in frame";
      friendly_error(ss.str());
    }
    return _frame[column_name];
  }

  /**
   * Look for innerListName as a named element in outerList. If it doesn't exist, swear. If it does exist, assume the content is a list of the form list(column_name="someName")
   * and return the relevant vector if it exists, throw an error if not.
   */
  const NumericVector parse_numeric_column_or_error(const List &outerList, const char *innerListName, const char *context_name)
  {
    // Check the name exists in the outer list; if not, swear.
    if (!outerList.containsElementNamed(innerListName))
    {
      ostringstream ss;
      ss << "Cannot find element '" << innerListName << "' in list";
      if (0 != context_name)
        ss << " for " << context_name;
      friendly_error(ss.str());
    }

    const char* column_name = parse_column_name_or_error(outerList[innerListName]);

    // Check the column exists in the frame; if not, swear.
    if (!_frame.containsElementNamed(column_name))
    {
      ostringstream ss;
      ss << "Cannot find column '" << column_name << "' in frame";
      friendly_error(ss.str());
    }
    return _frame[column_name];
  }

  /**
   * Look for elementName as a named element in list. If it doesn't exist, swear. If it does exist, assume the content is a boolean
   * and return that, throw an error if not.
   */
  const bool parse_boolean_or_error(const List &list, const char *elementName)
  {
    // Check the name exists in the outer list; if not, swear.
    if (!list.containsElementNamed(elementName))
    {
      ostringstream ss;
      ss << "Cannot find element '" << elementName << "' in list";
      friendly_error(ss.str());
    }

    return list[elementName];
  }

  /**
   * Look for elementName as a named element in list. If it doesn't exist, swear. If it does exist, assume the content is an integer
   * and return that, throw an error if not.
   */
  const int parse_int_or_error(const List &list, const char *elementName)
  {
    // Check the name exists in the outer list; if not, swear.
    if (!list.containsElementNamed(elementName))
    {
      ostringstream ss;
      ss << "Cannot find element '" << elementName << "' in list";
      friendly_error(ss.str());
    }

    return list[elementName];
  }

  /**
   * Look for elementName as a named element in list. If it doesn't exist, swear. If it does exist, assume the content is a list of strings
   * and return that as a vector, throw an error if not.
   */
  const vector<string> parse_string_vector_or_error(const List &list, const char *elementName)
  {
    // Check the name exists in the outer list; if not, swear.
    if (!list.containsElementNamed(elementName))
    {
      ostringstream ss;
      ss << "Cannot find element '" << elementName << "' in list";
      friendly_error(ss.str());
    }

    return list[elementName];
  }

  /**
   * Look for elementName as a named element in list. If it doesn't exist, swear. If it does exist, assume the content is a string
   * and return that, throw an error if not.
   */
  const char*const parse_string_or_error(const List &list, const char *elementName)
  {
    // Check the name exists in the outer list; if not, swear.
    if (!list.containsElementNamed(elementName))
    {
      ostringstream ss;
      ss << "Cannot find element '" << elementName << "' in list";
      friendly_error(ss.str());
    }

    return list[elementName];
  }

  Influencer parse_influencer_or_error(const string &influenced_by_disease_name, const List &influenced_by_parameter_list, const char*const disease_name)
  {
    const Disease &influenced_by_disease = *(_simulation.diseases->at(influenced_by_disease_name));
    const NumericVector multiplier = parse_numeric_column_or_error(influenced_by_parameter_list, "multiplier", disease_name);
    const int lag = parse_int_or_error(influenced_by_parameter_list, "lag");
    return Influencer(influenced_by_disease, multiplier, lag);
  }

  const vector<Influencer> parse_influencers_or_error(const List &influenced_by_list, const char*const disease_name)
  {
    const StringVector &influenced_by_disease_names = influenced_by_list.names();
    vector<Influencer> influencers;
    for (int i = 0; i < influenced_by_disease_names.size(); i++)
    {
      const string influenced_by_disease_name = as<string>(influenced_by_disease_names[i]);
      const List &influenced_by_parameter_list = influenced_by_list[i];
      influencers.push_back(parse_influencer_or_error(influenced_by_disease_name, influenced_by_parameter_list, disease_name));
    }
    return influencers;
  }

  unique_ptr<Incidence> parse_incidence_null_or_error(const List &incidence_list, const char*const disease_name)
  {
    return unique_ptr<Incidence>(new NullIncidence(_simulation.simulant_year_count()));
  }

  unique_ptr<Incidence> parse_incidence_universal_or_error(const List &incidence_list, const char*const disease_name)
  {
    return unique_ptr<Incidence>(new UniversalIncidence(_simulation.simulant_year_count()));
  }

  unique_ptr<Incidence> parse_incidence_type_1_or_error(const List &incidence_list, const char*const disease_name)
  {
    const NumericVector incidence = parse_numeric_column_or_error(incidence_list, "incidence", disease_name);
    const IntegerVector prevalence = parse_integer_column_or_error(incidence_list, "prevalence", disease_name);
    return unique_ptr<Incidence>(new IncidenceType1(incidence, prevalence));
  }

  unique_ptr<Incidence> parse_incidence_type_2_or_error(const List &incidence_list, const char*const disease_name)
  {
    const NumericVector incidence = parse_numeric_column_or_error(incidence_list, "incidence", disease_name);
    const NumericVector rn = parse_numeric_column_or_error(incidence_list, "rn", disease_name);
    const IntegerVector prevalence = parse_integer_column_or_error(incidence_list, "prevalence", disease_name);
    const bool can_recur = parse_boolean_or_error(incidence_list, "can_recur");
    return unique_ptr<Incidence>(new IncidenceType2(incidence, rn, prevalence, can_recur));
  }

  unique_ptr<Incidence> parse_incidence_type_3_or_error(const List &incidence_list, const char*const disease_name)
  {
    const NumericVector incidence = parse_numeric_column_or_error(incidence_list, "incidence", disease_name);
    const NumericVector rn = parse_numeric_column_or_error(incidence_list, "rn", disease_name);
    const IntegerVector prevalence = parse_integer_column_or_error(incidence_list, "prevalence", disease_name);
    return unique_ptr<Incidence>(new IncidenceType3(incidence, rn, prevalence));
  }

  unique_ptr<Incidence> parse_incidence_type_4_or_error(const List &incidence_list, const char*const disease_name)
  {
    return unique_ptr<Incidence>(new IncidenceType4(_simulation.simulant_year_count()));
  }

  unique_ptr<Incidence> parse_incidence_type_5_or_error(const List &incidence_list, const char*const disease_name)
  {
    const NumericVector incidence = parse_numeric_column_or_error(incidence_list, "incidence", disease_name);
    const NumericVector rn = parse_numeric_column_or_error(incidence_list, "rn", disease_name);
    const IntegerVector prevalence = parse_integer_column_or_error(incidence_list, "prevalence", disease_name);
    return unique_ptr<Incidence>(new IncidenceType5(incidence, rn, prevalence));
  }


  void parse_incidence_type_3_phase_2(const List &incidence_list, const char*const disease_name)
  {
    // Check the name exists in the outer list; if not, swear.
    if (!incidence_list.containsElementNamed("influenced_by"))
    {
      ostringstream ss;
      ss << "Cannot find element 'influenced_by' in incidence list for " << disease_name;
      friendly_error(ss.str());
    }
    const List &influenced_by_list = incidence_list["influenced_by"];
    vector<Influencer> influencers = parse_influencers_or_error(influenced_by_list, disease_name);
    Disease &me_disease = *(_simulation.diseases->at(disease_name));
    ((IncidenceType3&)me_disease.incidence()).set_influencers(influencers);
  }

  void parse_incidence_type_4_phase_2(const List &incidence_list, const char*const disease_name)
  {
    const vector<string> prerequisite_names = parse_string_vector_or_error(incidence_list, "aggregates");
    vector<const Disease*> prerequisites;
    for (string prerequisite_name : prerequisite_names)
    {
      Disease *prerequisite = _simulation.diseases->at(prerequisite_name).get();
      prerequisites.push_back(prerequisite);
    }
    const Disease &me_disease = *(_simulation.diseases->at(disease_name));
    ((IncidenceType4&)me_disease.incidence()).prerequisites(prerequisites);
  }

  void parse_incidence_type_5_phase_2(const List &incidence_list, const char*const disease_name)
  {
    // Check the name exists in the outer list; if not, swear.
    if (!incidence_list.containsElementNamed("influenced_by"))
    {
      ostringstream ss;
      ss << "Cannot find element 'influenced_by' in incidence list for " << disease_name;
      friendly_error(ss.str());
    }
    const List &influenced_by_list = incidence_list["influenced_by"];
    vector<Influencer> influencers = parse_influencers_or_error(influenced_by_list, disease_name);
    Disease &me_disease = *(_simulation.diseases->at(disease_name));
    ((IncidenceType5&)me_disease.incidence()).set_influencers(influencers);
  }


  unique_ptr<Diagnosis> parse_diagnosis_null_or_error(const List &incidence_list, const IntegerVector &prevalence, const char*const disease_name)
  {
    return unique_ptr<Diagnosis>(new NullDiagnosis(_simulation.simulant_year_count()));
  }

  unique_ptr<Diagnosis> parse_diagnosis_type_1_or_error(const List &incidence_list, const IntegerVector &prevalence, const char*const disease_name)
  {
    const NumericVector probability = parse_numeric_column_or_error(incidence_list, "probability", disease_name);
    const NumericVector rn = parse_numeric_column_or_error(incidence_list, "rn", disease_name);
    const IntegerVector diagnosed = parse_integer_column_or_error(incidence_list, "diagnosed", disease_name);
    return unique_ptr<Diagnosis>(new DiagnosisType1(probability, rn, prevalence, diagnosed));
  }

  unique_ptr<Mortality> parse_mortality_null_or_error(const List &incidence_list, IntegerVector &prevalence, IntegerVector &diagnosed, const char*const disease_name)
  {
    return unique_ptr<Mortality>(new NullMortality(_simulation.simulant_year_count()));
  }

  unique_ptr<Mortality> parse_mortality_type_1_or_error(const List &incidence_list, IntegerVector &prevalence, IntegerVector &diagnosed, const char*const disease_name)
  {
    const NumericVector probability = parse_numeric_column_or_error(incidence_list, "probability", disease_name);
    const NumericVector rn = parse_numeric_column_or_error(incidence_list, "rn", disease_name);
    const int cod = parse_int_or_error(incidence_list, "code");
    return unique_ptr<Mortality>(new MortalityType1(probability, rn, prevalence, cod));
  }

  unique_ptr<Mortality> parse_mortality_type_2_or_error(const List &incidence_list, IntegerVector &prevalence, IntegerVector &diagnosed, const char*const disease_name)
  {
    const NumericVector probability = parse_numeric_column_or_error(incidence_list, "probability", disease_name);
    const NumericVector rn = parse_numeric_column_or_error(incidence_list, "rn", disease_name);
    const int cod = parse_int_or_error(incidence_list, "code");
    const int cured_after_years = parse_int_or_error(incidence_list, "cure");
    return unique_ptr<Mortality>(new MortalityType2(probability, rn, prevalence, cod, cured_after_years));
  }

  unique_ptr<Mortality> parse_mortality_type_3_or_error(const List &incidence_list, IntegerVector &prevalence, IntegerVector &diagnosed, const char*const disease_name)
  {
    const NumericVector probability = parse_numeric_column_or_error(incidence_list, "probability", disease_name);
    const NumericVector rn = parse_numeric_column_or_error(incidence_list, "rn", disease_name);
    const int cod = parse_int_or_error(incidence_list, "code");
    return unique_ptr<Mortality>(new MortalityType3(probability, rn, prevalence, cod));
  }

  void parse_mortality_type_3_phase_2(const List &mortality_list, const char*const disease_name)
  {
    // Check the name exists in the outer list; if not, swear.
    if (!mortality_list.containsElementNamed("influenced_by"))
    {
      ostringstream ss;
      ss << "Cannot find element 'influenced_by' in mortality list for " << disease_name;
      friendly_error(ss.str());
    }
    const List &influenced_by_list = mortality_list["influenced_by"];
    vector<Influencer> influencers = parse_influencers_or_error(influenced_by_list, disease_name);
    Disease &me_disease = *(_simulation.diseases->at(disease_name));
    ((MortalityType3&)me_disease.mortality()).set_influencers(influencers);
  }

  typedef unique_ptr<Incidence>(SimulationParser::*Incidence_function)(const List &, const char*const);
  typedef map<const string, Incidence_function> Incidence_map;
  const Incidence_map incidence_functions
  {
    {"Null", &SimulationParser::parse_incidence_null_or_error},
    {"Universal", &SimulationParser::parse_incidence_universal_or_error},
    {"Type1", &SimulationParser::parse_incidence_type_1_or_error},
    {"Type2", &SimulationParser::parse_incidence_type_2_or_error},
    {"Type3", &SimulationParser::parse_incidence_type_3_or_error},
    {"Type4", &SimulationParser::parse_incidence_type_4_or_error},
	{"Type5", &SimulationParser::parse_incidence_type_5_or_error}
  };
  typedef void(SimulationParser::*Incidence_function_phase_2)(const List &, const char*const);
  typedef map<const string, Incidence_function_phase_2> Incidence_map_phase_2;
  const Incidence_map_phase_2 incidence_phase_2_functions
  {
    {"Type3", &SimulationParser::parse_incidence_type_3_phase_2},
    {"Type4", &SimulationParser::parse_incidence_type_4_phase_2},
	{"Type5", &SimulationParser::parse_incidence_type_5_phase_2}
  };

  typedef unique_ptr<Diagnosis>(SimulationParser::*Diagnosis_function)(const List &, const IntegerVector &, const char*const);
  typedef map<const string, Diagnosis_function> Diagnosis_map;
  const Diagnosis_map diagnosis_functions
  {
    {"Null", &SimulationParser::parse_diagnosis_null_or_error},
    {"Type1", &SimulationParser::parse_diagnosis_type_1_or_error}
  };

  typedef unique_ptr<Mortality>(SimulationParser::*Mortality_function)(const List &, IntegerVector &, IntegerVector &, const char*const);
  typedef map<const string, Mortality_function> Mortality_map;
  const Mortality_map mortality_functions
  {
    {"Null", &SimulationParser::parse_mortality_null_or_error},
    {"Type1", &SimulationParser::parse_mortality_type_1_or_error},
    {"Type2", &SimulationParser::parse_mortality_type_2_or_error},
    {"Type3", &SimulationParser::parse_mortality_type_3_or_error}
  };
  typedef void(SimulationParser::*Mortality_function_phase_2)(const List &, const char*const);
  typedef map<const string, Mortality_function_phase_2> Mortality_map_phase_2;
  const Mortality_map_phase_2 mortality_phase_2_functions
  {
    {"Type3", &SimulationParser::parse_mortality_type_3_phase_2}
  };

  /**
   * Assume incidence_list is a list of parameters for a disease incidence, and construct the corresponding incidence.
   * Return a parsed Incidence on the heap containing the relevant details.
   */
  unique_ptr<Incidence> parse_incidence_or_error(const List &incidence_list, const char*const disease_name)
  {
    const char*const type = parse_string_or_error(incidence_list, "type");
    Incidence_map::const_iterator it = incidence_functions.find(type);
    if (it == incidence_functions.end())
    {
      ostringstream ss;
      ss << "Unknown incidence type '" << type << "' in disease '" << disease_name << "'";
      friendly_error(ss.str());
    }
    Incidence_function to_call = (*it).second;
    return ((*this).*to_call)(incidence_list, disease_name);
  }

  /**
   * Assume incidence_list is a list of parameters for a disease incidence, and fix up the corresponding incidence if required.
   */
  void parse_incidence_phase_2_or_error(const List &incidence_list, const char*const disease_name)
  {
    const char*const type = parse_string_or_error(incidence_list, "type");
    Incidence_map_phase_2::const_iterator it = incidence_phase_2_functions.find(type);
    if (it == incidence_phase_2_functions.end())
      return;
    Incidence_function_phase_2 to_call = (*it).second;
    ((*this).*to_call)(incidence_list, disease_name);
  }

  /**
   * Assume diagnosis_list is a list of parameters for a disease diagnosis, and construct the corresponding diagnosis.
   * Return a parsed Diagnosis on the heap containing the relevant details.
   */
  unique_ptr<Diagnosis> parse_diagnosis_or_error(const List &diagnosis_list, const IntegerVector &prevalence, const char *disease_name)
  {
    const char*const type = parse_string_or_error(diagnosis_list, "type");
    Diagnosis_map::const_iterator it = diagnosis_functions.find(type);
    if (it == diagnosis_functions.end())
    {
      ostringstream ss;
      ss << "Unknown diagnosis type '" << type << "' in disease '" << disease_name << "'";
      friendly_error(ss.str());
    }
    Diagnosis_function to_call = (*it).second;
    return ((*this).*to_call)(diagnosis_list, prevalence, disease_name);
  }

  /**
   * Assume mortality_list is a list of parameters for a disease mortality, and construct the corresponding mortality.
   * Return a parsed Mortality on the heap containing the relevant details.
   */
  unique_ptr<Mortality> parse_mortality_or_error(const List &mortality_list, IntegerVector &prevalence, IntegerVector &diagnosis, const char *disease_name)
  {
    const char*const type = parse_string_or_error(mortality_list, "type");
    Mortality_map::const_iterator it = mortality_functions.find(type);
    if (it == mortality_functions.end())
    {
      ostringstream ss;
      ss << "Unknown mortality type '" << type << "' in disease '" << disease_name << "'";
      friendly_error(ss.str());
    }
    Mortality_function to_call = (*it).second;
    return ((*this).*to_call)(mortality_list, prevalence, diagnosis, disease_name);
  }

  /**
   * Assume mortality_list is a list of parameters for a disease mortality, and fix up the corresponding mortality if required.
   */
  void parse_mortality_phase_2_or_error(const List &mortality_list, const char*const disease_name)
  {
    const char*const type = parse_string_or_error(mortality_list, "type");
    Mortality_map_phase_2::const_iterator it = mortality_phase_2_functions.find(type);
    if (it == mortality_phase_2_functions.end())
      return;
    Mortality_function_phase_2 to_call = (*it).second;
    ((*this).*to_call)(mortality_list, disease_name);
  }

  /**
   * Assume disease_list is a disease, and check for incidence (required) and mortality (optional, default NullMortality).
   * Return a parsed Disease on the heap containing the relevant details.
   */
  unique_ptr<Disease> parse_disease_or_error(const List &disease_list, const char *disease_name)
  {
    // Check incidence exists in the disease_list; if not, swear.
    if (!disease_list.containsElementNamed("incidence"))
    {
      ostringstream ss;
      ss << "Cannot find required element 'incidence' in disease '" << disease_name << "'";
      friendly_error(ss.str());
    }
    unique_ptr<Incidence> incidence(parse_incidence_or_error(disease_list["incidence"], disease_name));

    // By contrast, diagnosis and mortality elements are optional; if missing, use a Null diagnosis or mortality as appropriate.
    unique_ptr<Diagnosis> diagnosis(disease_list.containsElementNamed("diagnosis")
                                      ? parse_diagnosis_or_error(disease_list["diagnosis"], incidence->out_prevalence(), disease_name)
                                        : unique_ptr<Diagnosis>(new NullDiagnosis(_simulation.simulant_year_count())));

    unique_ptr<Mortality> mortality(disease_list.containsElementNamed("mortality")
                                      ? parse_mortality_or_error(disease_list["mortality"], incidence->out_prevalence(), diagnosis->out_diagnosis(), disease_name)
                                        : unique_ptr<Mortality>(new NullMortality(_simulation.simulant_year_count())));

    unique_ptr<Disease> disease(new Disease(disease_name, incidence, diagnosis, mortality));
    disease->mortality().disease(*disease);
    return disease;
  }

  /**
   * Assume disease_list is a disease which is already in the diseases map. Fix up anything that needs fixing up.
   */
  void parse_disease_phase_2_or_error(const List &disease_list, const char *disease_name)
  {
    // At the moment, only Incidences and Mortalities need fixing up, so only parse those.
    if (disease_list.containsElementNamed("incidence"))
      parse_incidence_phase_2_or_error(disease_list["incidence"], disease_name);
    if (disease_list.containsElementNamed("mortality"))
      parse_mortality_phase_2_or_error(disease_list["mortality"], disease_name);
  }

  /**
   * Look for innerListName as a named element in outerList. If it doesn't exist, swear. If it does exist, assume the content is a list of the form list(disease_name=list(...))
   * and return a vector of the contained diseases.
   */
  void parse_diseases_or_error(const List &outerList, const char *innerListName)
  {
    // Check the name exists in the outer list; if not, swear.
    if (!outerList.containsElementNamed(innerListName))
    {
      ostringstream ss;
      ss << "Cannot find element '" << innerListName << "' in list";
      friendly_error(ss.str());
    }

    const List &disease_list = outerList[innerListName];
    const CharacterVector &disease_names = disease_list.names();
    _simulation.diseases = unique_ptr<unordered_map<string, unique_ptr<Disease>>>(new unordered_map<string, unique_ptr<Disease>>(disease_list.size()));
    // Phase 1: Find the diseases and populate what we can.
    for (int i = 0; i < disease_list.size(); i++)
      _simulation.diseases->emplace(as<string>(disease_names[i]), move(parse_disease_or_error(disease_list[i], disease_names[i])));
    // Phase 2: Populate anything that requires all diseases to be known.
    for (int i = 0; i < disease_list.size(); i++)
      parse_disease_phase_2_or_error(disease_list[i], disease_names[i]);
  }

  void parse_strata(const List &outerList, const char *innerListName)
  {
    const vector<string> strata_names = parse_string_vector_or_error(outerList, innerListName);
    for (auto const& stratum_name : strata_names)
      _simulation.add_stratum(stratum_name);
  }
public:
  SimulationParser(const DataFrame& frame, Simulation &simulation)
    : _simulation(simulation)
  , _frame(frame)
  {
    _simulation.frame(&_frame);
  }

  void parse_simulation(const List &simulationStructureList)
  {
    _simulation.pids = parse_integer_column_or_error(simulationStructureList, "pids", "top level");
    _simulation.years = parse_integer_column_or_error(simulationStructureList, "years", "top level");
    _simulation.init_year = parse_int_or_error(simulationStructureList, "init_year");
    _simulation.scenario_name(parse_string_or_error(simulationStructureList, "scenario"));
    _simulation.rn_all_cause_mrtl = parse_numeric_column_or_error(simulationStructureList, "rn_all_cause_mrtl", "top level");
    parse_diseases_or_error(simulationStructureList, "diseases");
    parse_strata(simulationStructureList, "strata_for_outputs");
  }
};
