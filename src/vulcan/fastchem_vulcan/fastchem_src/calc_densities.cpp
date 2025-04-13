/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2019 Daniel Kitzmann, Joachim Stock
*
* FastChem is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* FastChem is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You find a copy of the GNU General Public License in the main
* FastChem directory under <license.md>. If not, see
* <http://www.gnu.org/licenses/>.
*/


#include "fastchem.h"


#include <algorithm>
#include <vector>
#include <cmath>


#include "input_output_struct.h"


namespace fastchem {


template <class double_type>
unsigned int FastChem<double_type>::calcDensities(FastChemInput& input, FastChemOutput& output)
{
  if (!is_initialized)
    return FASTCHEM_INITIALIZATION_FAILED;
  

  if (is_busy == true)
  {
    std::cout << "FastChem can only be called once per instance!\n";
    return FASTCHEM_IS_BUSY;
  }

  is_busy = true;


  const double temperature = input.temperature;
  const double pressure = input.pressure;


  if ( pressure <= 0.0 || temperature <= 0.0 )
  {
     
     std::cout << "Check input values! Temperature: " << temperature << "  pressure: " << pressure << "\n";
     is_busy = false;
     
     return FASTCHEM_WRONG_INPUT_VALUES;

  }
  
  
  unsigned int status = calcDensity(temperature, pressure, input.use_previous_solution, 
                                    output.number_densities,
                                    output.total_element_density, 
                                    output.mean_molecular_weight,
                                    output.element_conserved,
                                    output.nb_chemistry_iterations);
  
    
  is_busy = false;

  return status;
}




//Solve the chemistry for a single temperature and a single pressure
//Note: this is a private function, that can not be accessed from outside of FastChem
//This function will be called by any public calcDensity function
template <class double_type>
unsigned int FastChem<double_type>::calcDensity(const double temperature, const double pressure, const bool use_previous_solution,
                                                std::vector<double>& number_densities, double& total_element_density, 
                                                double& mean_molecular_weight,
                                                std::vector<unsigned int>& element_conserved,
                                                unsigned int& nb_chemistry_iterations)
{

  for (auto & i : molecules)  i.calcMassActionConstant(temperature);

  //this value will be fixed.
  double_type gas_density = pressure/(CONST_K * temperature);

  
  if (use_previous_solution == true)
  {

   //if we use the previous solution, convert the stored mixing ratios to number densities
   for (auto & i : species)  i->number_density *= gas_density;
  
  }
  else
  {

    //for a fresh start set all species to the minimum value
    for (auto & i : species) i->number_density = options.element_density_minlimit;
    
    //set the initial electron density to 1 (for stability reasons)
    if (e_ != FASTCHEM_UNKNOWN_SPECIES)
      elements[e_].number_density = 1.0;

  }
  


  //call the main FastChem solver  
  bool fastchem_converged = solveFastchem(temperature, gas_density, nb_chemistry_iterations);


  
  if (!fastchem_converged && options.verbose_level >= 1) std::cout << "Convergence problem in FastChem. :(\n";


  //return output
  number_densities.assign(nb_species, 0.0);

  for (size_t i=0; i<nb_species; ++i)
    number_densities[i] = species[i]->number_density; 


  mean_molecular_weight = meanMolecularWeight(gas_density);


  total_element_density = totalElementDensity();


  for (auto & i : elements) i.checkElementConservation(molecules, total_element_density, options.accuracy);
  
  element_conserved.assign(nb_elements, 0);

  for (size_t i=0; i<nb_elements; ++i)
    element_conserved[i] = elements[i].element_conserved;


  unsigned int return_state = FASTCHEM_SUCCESS;

  if (!fastchem_converged) return_state = FASTCHEM_NO_CONVERGENCE;


  //store the mixing ratios in case we want to use them in the future
  for (auto & i : species) i->number_density /= gas_density;

  
  return return_state;
}



template class FastChem<double>;
template class FastChem<long double>;


}


