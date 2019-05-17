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


#include "options.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include <cmath>



namespace fastchem {


//Read the FastChem parameter file
template <class double_type>
bool FastChemOptions<double_type>::readParameterFile(const std::string& model_parameter_file)
{
  bool initialization_status = false;

  std::fstream file;

  file.open(model_parameter_file.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to open FastChem parameter file " << model_parameter_file << "\n";

    return initialization_status;
  }


  std::string file_name, line;
  double parameter_value;


  std::getline(file, line);
  file >> file_name;

  if (file_name == "")
    return initialization_status;
  else
    element_abundances_file = file_name;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> file_name;

  if (file_name == "")
    return initialization_status;
  else
    chemical_element_file = file_name;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> file_name;

  if (file_name == "")
    return initialization_status;
  else
    species_data_file = file_name;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    accuracy = parameter_value;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    newton_err = parameter_value;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    nb_max_fastchem_iter = parameter_value;


  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);
  file >> parameter_value;

  if (!(parameter_value > 0))
    return initialization_status;
  else
    nb_max_newton_iter = parameter_value;


  initialization_status = true;


  if (verbose_level >= 3)
  {
    std::cout << "Parameter values read: \n";
    std::cout << "element_abundances " << element_abundances_file << "\n";
    std::cout << "element_file " << chemical_element_file << "\n";
    std::cout << "species_data " << species_data_file << "\n";
    std::cout << "accuracy " << accuracy << "\n";
    std::cout << "newton_err " << newton_err << "\n";
    std::cout << "iter_max " << nb_max_fastchem_iter << "\n";
    std::cout << "iter_max_newton " << nb_max_newton_iter << "\n";
  }


  return initialization_status;
}




template struct FastChemOptions<double>;
template struct FastChemOptions<long double>;


}
