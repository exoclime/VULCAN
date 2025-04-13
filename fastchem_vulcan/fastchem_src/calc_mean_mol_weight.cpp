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

#include <vector>
#include <cmath>



namespace fastchem {


//Calculate the mean molecular weight of the converged network
//See Eq. (2.10)
template <class double_type>
double FastChem<double_type>::meanMolecularWeight(const double gas_density)
{
   double mean_molecular_weight = 0.0;

   
   for (auto & i : species) mean_molecular_weight += i->molecular_weight * i->number_density;

     
   mean_molecular_weight /= gas_density;


   return mean_molecular_weight;
}



template class FastChem<double>;
template class FastChem<long double>;


}
