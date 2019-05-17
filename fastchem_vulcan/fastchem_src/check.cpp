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


#include "species_struct.h"


#include <vector>
#include <cmath>


namespace fastchem {


//Check for the number density of elements
template <class double_type>
void Element<double_type>::checkN(const double_type& min_limit, const double_type& gas_density)
{

  if (this->number_density < min_limit) this->number_density = min_limit;

  if (this->number_density > gas_density) this->number_density = gas_density;

}



//Check for the number density of molecules
template <class double_type>
void Molecule<double_type>::checkN(const double_type& min_limit, const double_type& gas_density)
{

  if (this->number_density < min_limit) this->number_density = min_limit;

  if (this->number_density > gas_density) this->number_density = gas_density;

}



//Check for charge conservation
template <class double_type>
bool Element<double_type>::checkChargeConservation(const std::vector< Molecule<double_type> >& molecules, const double_type& accuracy)
{

  //Am I the electron?
  if (this->symbol != "e-") return false;


  //If no ions present, charge conservation is automatically satisfied
  if (molecule_list.size() == 0)
  {
    element_conserved = 1;

    return true;
  }

   
  bool charge_conserved = false;

  //sum up all positive and negative charges in the network
  double_type positive_charge = 0;
  double_type negative_charge = this->number_density;


  for (size_t i=0; i<molecule_list.size(); ++i)
  {
    const unsigned int molecule_index = molecule_list[i];

    if (molecules[molecule_index].stoichometric_vector[index] < 0)
      positive_charge -= molecules[molecule_index].stoichometric_vector[index] * molecules[molecule_index].number_density;

    if (molecules[molecule_index].stoichometric_vector[index] > 0)
      negative_charge += molecules[molecule_index].stoichometric_vector[index] * molecules[molecule_index].number_density;
  }


  if (std::fabs(positive_charge - negative_charge)/std::sqrt(positive_charge*negative_charge) < accuracy)
    charge_conserved = true;
  else
    charge_conserved = false;

  element_conserved = charge_conserved;


  return charge_conserved;
}



//Check for element conservation
template <class double_type>
bool Element<double_type>::checkElementConservation(const std::vector< Molecule<double_type> >& molecules, const double_type total_density, const double_type& accuracy)
{

  //electrons are subject to charge conservation
  if (this->symbol == "e-")
    return checkChargeConservation(molecules, accuracy);


  //sum up the elements contained in each molecule and compare the result to its elemental abundance
  double_type sum = this->number_density;


  for (size_t i=0; i<molecule_list.size(); ++i)
    sum += molecules[molecule_list[i]].stoichometric_vector[index] * molecules[molecule_list[i]].number_density;

  sum /= total_density*epsilon;


  if (std::fabs(sum - 1.0L) < accuracy || molecule_list.size() == 0)
    element_conserved = 1;
  else
    element_conserved = 0;


  return element_conserved;
}



template struct Element<double>;
template struct Element<long double>;
template struct Molecule<double>;
template struct Molecule<long double>;

}
