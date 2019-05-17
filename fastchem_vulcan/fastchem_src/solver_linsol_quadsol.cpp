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


#include "solver.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


namespace fastchem {


//Solver for an element that is not part of other species
//See Eq. (2.32)
template <class double_type>
void FastChemSolver<double_type>::intertSol(Element<double_type>& species, std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules, 
                                            const double_type gas_density)
{

  species.number_density = species.epsilon * gas_density - species.number_density_min - species.number_density_maj;

}


//Analytic solution for linear equation, see Sect. 2.4.2 and Eq. (2.32)
template <class double_type>
void FastChemSolver<double_type>::linSol(Element<double_type>& species, std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules, 
                                         const double_type gas_density)
{
  double_type scaling_factor = 0.0;

  //in case we use the scaling factor, see Appendix A for details
  if (options->use_scaling_factor)
    scaling_factor = solverScalingFactor(species, gas_density);


  if (scaling_factor > 700.0 && options->verbose_level >= 3)
    std::cout << "FastChem: WARNING: Underflow in LinSol for element " << species.symbol << "\n";


  //calculation of coefficient A_j1, see Eq. (2.28)
  const double_type A1 = A1Coeff(species, elements, molecules) + std::exp(-scaling_factor);


  //calculation of coefficient A_j0, see Eq. (2.27)
  const double_type A0 = std::exp(-scaling_factor) * (species.number_density_maj + species.number_density_min - gas_density * species.epsilon);


  //calculation of n_j, Eq. (2.32)
  species.number_density = -A0/A1;

}


//Analytic solution for quadratic equation, see Sect. 2.4.2 and Eq. (2.32)
template <class double_type>
void FastChemSolver<double_type>::quadSol(Element<double_type>& species, std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules, 
                                          const double_type gas_density)
{
  double_type scaling_factor = 0.0;


  //in case we use the scaling factor, see Appendix A for details
  if (options->use_scaling_factor)
    scaling_factor = solverScalingFactor(species, gas_density);


  if (scaling_factor > 700.0 && options->verbose_level >= 3)
    std::cout << "FastChem: WARNING: Underflow in QuadSol for element " << species.symbol << "\n";


  //calculation of coefficient A_j2, see Eq. (2.29)
  const double_type A2 = A2Coeff(species, elements, molecules);

  if (A2<1.e-4900L)
  {
    if (options->verbose_level >= 3) std::cout << "FastChem: Underflow in QuadSol for species " <<  species.symbol << " : switching to LinSol.\n";

    linSol(species, elements, molecules, gas_density);

    return;
  }

  
  //calculation of coefficient A_j1, see Eq. (2.28)
  const double_type A1 = A1Coeff(species, elements, molecules) + std::exp(-scaling_factor);

  //calculation of coefficient A_j0, see Eq. (2.27)
  const double_type A0 = std::exp(-scaling_factor) * (species.number_density_maj + species.number_density_min - gas_density * species.epsilon);
  

  //calculation of n_j, Eq. (2.32)
  const double_type Qj = -0.5 * (A1 + std::sqrt(A1*A1 - 4.*A2*A0));
    
  species.number_density = A0/Qj;

}



template class FastChemSolver<double>;
template class FastChemSolver<long double>;

}



