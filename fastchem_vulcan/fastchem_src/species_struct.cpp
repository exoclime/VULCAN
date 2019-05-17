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

#include <cmath>
#include <iostream>


namespace fastchem {




//Determination of the correction factor n_j_min, that contains the number density of an element j contained in molecules
//with elements that are less abundant than j
//See Eq. (2.23)
template <class double_type>
void Element<double_type>::calcMinorSpeciesDensities(const std::vector< Molecule<double_type> > &molecules)
{
  
  number_density_min = 0.0;
  
  for (size_t jj=0; jj<minor_molecules.size(); ++jj)
  { 
    const unsigned int j = minor_molecules[jj];
    
    number_density_min += (molecules[j].stoichometric_vector[index] + epsilon * molecules[j].sigma) * molecules[j].number_density;
  }
 
}



template <class double_type>
void Element<double_type>::calcEpsilon(const std::vector< Element<double_type> > &elements)
{ 

  double_type element_sum = 0.0;

  for (size_t i=0; i<elements.size(); ++i)
    element_sum += elements[i].abundance;

  epsilon = this->abundance/element_sum;

}


template struct Element<double>;
template struct Element<long double>;
template struct Molecule<double>;
template struct Molecule<long double>;


}
