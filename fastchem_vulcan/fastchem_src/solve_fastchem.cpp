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

#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>


namespace fastchem {


//This is the main FastChem iteration
template <class double_type>
bool FastChem<double_type>::solveFastchem(const double temperature_gas, const double gas_density, unsigned int& nb_iterations)
{

  for (auto & i : elements) i.number_density_maj = 0.0;

  //starting values for contribution of minor species
  for (auto & i : elements) i.calcMinorSpeciesDensities(molecules);


  std::vector<double_type> number_density_old(nb_species, 0.0);

  for (size_t i=0; i<nb_species; ++i)
    number_density_old[i] = species[i]->number_density;


  bool converged = false;
  unsigned int iter_step = 0;
  bool use_backup_solver = false;

  unsigned int max_iter = options.nb_max_fastchem_iter;

  for (iter_step=0; iter_step<max_iter; ++iter_step)
  {
    double_type n_maj = 0.0;
   
    //calculate the element densities in their order
    for (auto it = element_calculation_order.begin(); it<element_calculation_order.end(); it++)
      calculateElementDensities(elements[*it], gas_density, use_backup_solver, n_maj);

   
    for (auto & i : elements) i.calcMinorSpeciesDensities(molecules);
  
    
    if (e_ != FASTCHEM_UNKNOWN_SPECIES)  //only calculate electrons if they are present in the element list
      calculateElectronDensities(elements[e_], number_density_old[e_], gas_density);
 

    //check if n_j_min are small enough, if not use backup solver
    for (auto & i : elements)
      if ( (i.number_density_min + i.number_density_maj > i.epsilon * gas_density) && use_backup_solver == false)
      {
        use_backup_solver = true;
       
        if (options.verbose_level >= 4)
          std::cout << "Too large n_min and n_maj for species " << i.symbol << ". Switching to backup.  Iteration step: " << iter_step << "\n";

        break;
      }


    //convergence check
    if (iter_step > 0)
    {
      converged = true;

      for (size_t i=0; i<nb_species; ++i)
        if (std::fabs((species[i]->number_density - number_density_old[i])) > options.accuracy*number_density_old[i]
             && species[i]->number_density/gas_density > 1.e-155)
        { 
          converged = false;
          break;
        }

    }


    if (converged)
      break;


    //in case the standard FastChem iteration doesn't converge, switch to the backup solver
    if (iter_step == max_iter-1 && !converged && use_backup_solver == false)
    {
      if (options.verbose_level >= 4)
        std::cout << "Standard FastChem iteration failed. Switching to backup. " << "\n";

      use_backup_solver = true;
      max_iter += options.nb_max_fastchem_iter;
    }


    for (size_t i=0; i<nb_species; ++i)
      number_density_old[i] = species[i]->number_density;
  }
  
  nb_iterations = iter_step;
 
  return converged;
}



template class FastChem<double>;
template class FastChem<long double>;

}
