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

#include <string>
#include <vector>
#include <cmath>


namespace fastchem {


//Constructor for the FastChem class
//Requires: parameter file path and the initial verbose level that is used while reading the input files
template <class double_type>
FastChem<double_type>::FastChem(const std::string& model_parameter_file, const unsigned int verbose_level_init) : solver(&options)
{
  options.verbose_level = verbose_level_init;

  bool parameter_file_loaded = false;

  if (model_parameter_file != "")
    parameter_file_loaded = options.readParameterFile(model_parameter_file);


  if (!parameter_file_loaded)
  {
    std::cout << "Error reading parameters\n";
    is_initialized = false;
  }


  if (parameter_file_loaded) init();
}



//Copy constructor
//Could be made more pretty, but this one does the job as well...
template <class double_type>
FastChem<double_type>::FastChem(const FastChem &obj) : solver(&options)
{
  nb_chemical_element_data = obj.nb_chemical_element_data;
  nb_species = obj.nb_species;
  nb_molecules = obj.nb_molecules;
  nb_elements = obj.nb_elements;

  e_ = obj.e_;

  is_initialized = obj.is_initialized;
  is_busy = false;

  chemical_element_data = obj.chemical_element_data;
  elements = obj.elements;
  molecules = obj.molecules;

  element_calculation_order = obj.element_calculation_order;

  for (size_t i=0; i<nb_elements; ++i)
    species.push_back(&elements[i]);

  for (size_t i=0; i<nb_molecules; ++i)
    species.push_back(&molecules[i]);


  //Options object
  options.nb_max_fastchem_iter = obj.options.nb_max_fastchem_iter;
  options.nb_max_bisection_iter = obj.options.nb_max_bisection_iter;
  options.nb_max_neldermead_iter = obj.options.nb_max_neldermead_iter;
  options.nb_max_newton_iter = obj.options.nb_max_newton_iter;

  options.element_density_minlimit = obj.options.element_density_minlimit;
  options.molecule_density_minlimit = obj.options.molecule_density_minlimit;

  options.accuracy = obj.options.accuracy;
  options.newton_err = obj.options.newton_err;

  options.verbose_level = obj.options.verbose_level;
  options.use_scaling_factor = obj.options.use_scaling_factor;
  


  options.chemical_element_file = obj.options.chemical_element_file;
  options.species_data_file = obj.options.species_data_file;
  options.element_abundances_file = obj.options.element_abundances_file;


  //Solver object
  solver.order_anion = obj.solver.order_anion;
  solver.order_cation = obj.solver.order_cation;
}







template class FastChem<double>;
template class FastChem<long double>;


}
