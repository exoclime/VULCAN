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
#include <string>
#include <vector>



namespace fastchem {


//Set the element abundances for all elements
template <class double_type>
void FastChem<double_type>::setElementAbundance(std::vector<double> abundances)
{

  if (abundances.size() != nb_elements)
  {
    std::cout << "Setting element abundances with an incorrect vector size\n";

    return;
  }


  for (size_t i=0; i<nb_elements; ++i)
  {
    if (i == e_) continue; //the abundance of the electron remains at zero

    chemical_element_data[elements[i].element_data_index].abundance = abundances[i];
    elements[i].abundance = abundances[i];

  }


  reInitialiseFastChem();
}




template class FastChem<double>;
template class FastChem<long double>;

}
