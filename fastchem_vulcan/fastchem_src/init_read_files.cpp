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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include <cmath>



namespace fastchem {


//Read the chemical elements file
template <class double_type>
bool FastChem<double_type>::readElementList()
{
  std::fstream file(options.chemical_element_file.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to open chemical element file " << options.chemical_element_file << "\n";

    return false;
  }


  std::string line;

  while (std::getline(file, line))
  {
    ChemicalElementData<double_type> element;
    std::istringstream input(line);

    input >> element.symbol >> element.name >> element.atomic_weight;

    chemical_element_data.push_back(element);

    chemical_element_data.back().abundance = 0;
  }


  file.close();


  nb_chemical_element_data = chemical_element_data.size();


  return true;
}



template <class double_type>
bool FastChem<double_type>::readSpeciesData()
{
  std::fstream file(options.species_data_file.c_str(), std::ios::in);


  if (file.fail())
  {
    std::cout << "Unable to open species data file " << options.species_data_file << "\n";

    return false;
  }


  std::string line;

  //header
  std::getline(file, line);
  std::getline(file, line);
  std::getline(file, line);



  while (std::getline(file, line))
  {
    std::istringstream input(line);

    std::string symbol, name, element_string, stoichometric_coeff_string;

    input >> symbol >> name;


    if (name == ":")
      name = "";
    else
      input >> element_string;


    std::vector<std::string> elements;
    std::vector<int> stoichometric_coeff;
    int charge = 0;


    while ( (input >> element_string >> stoichometric_coeff_string) && (element_string != "#") )
    {
      elements.push_back(element_string);
      stoichometric_coeff.push_back(std::stoi( stoichometric_coeff_string ));

      if (elements.back() == "e-") charge = -stoichometric_coeff.back();
    }


    std::string line;
    std::getline(file, line);

    std::istringstream ma_input(line);

    std::vector<double_type> mass_action_coeff;


    double ma_coefficient;

    while (ma_input >> ma_coefficient)
      mass_action_coeff.push_back(ma_coefficient);

    addMolecule(name, symbol, elements, stoichometric_coeff, mass_action_coeff, charge);

    //blank separation line
    std::getline(file, line);
  }




  file.close();



  nb_molecules = molecules.size();

  nb_species = nb_molecules + nb_elements;


  return true;
}



//Read the elemental abundances file
template <class double_type>
bool FastChem<double_type>::readElementAbundances()
{

  std::fstream file(options.element_abundances_file.c_str(), std::ios::in);


  if (file.fail())
  {
    std::cout << "Unable to open element abundances file " << options.element_abundances_file << "\n";

    return false;
  }



  std::string line;

  //header
  std::getline(file, line);


  while (std::getline(file, line))
  {
    std::istringstream input(line);

    std::string symbol;
    double abundance;

    input >> symbol >> abundance;


    abundance = std::pow(10., abundance - 12.);


    setElementAbundance(symbol, abundance);

    addAtom(symbol);
  }


  file.close();


  nb_elements = elements.size();


  return true;
}


template class FastChem<double>;
template class FastChem<long double>;


}
