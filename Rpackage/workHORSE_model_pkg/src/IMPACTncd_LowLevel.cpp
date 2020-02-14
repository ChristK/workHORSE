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



#include <memory>
#include <unordered_map>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;

#include "simulation_structure.h"
#include "incidence.h"
#include "diagnosis.h"
#include "mortality.h"
#include "simulation_parser.h"
#include "single_threaded_simulation.h"

// #include "rng.h"

//' @export
// [[Rcpp::export]]
List run_impactncd_simulation(const List &simulationStructureList, const DataFrame &frame, const List &input_list)
{
  unique_ptr<Simulation> simulation(new SingleTheadedSimulation);
  SimulationParser parser(frame, *simulation);
  parser.parse_simulation(simulationStructureList);
  simulation->run();
  List output_list(input_list);
  output_list.push_back(simulation->output_frame(), simulation->scenario_name());
  return output_list;
}
