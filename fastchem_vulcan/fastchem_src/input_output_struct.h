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


#ifndef _input_output_struct_h
#define _input_output_struct_h


#include <vector>


namespace fastchem {


struct FastChemInput
{

  double temperature = 0.0; 
  double pressure = 0.0;

  bool use_previous_solution = false;

};



struct FastChemOutput
{

  std::vector<double> number_densities;
  double total_element_density;
  double mean_molecular_weight = 0.0;

  //diagnostic output
  std::vector<unsigned int> element_conserved;
  unsigned int fastchem_flags = 0;
  unsigned int nb_chemistry_iterations = 0;

};



}


#endif