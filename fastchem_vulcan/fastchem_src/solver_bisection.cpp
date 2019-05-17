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

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>


namespace fastchem {

//Bisection method in one dimension
template <class double_type>
bool FastChemSolver<double_type>::bisectionSolve(Element<double_type>& species, std::vector<double_type>& Aj, const double gas_density)
{
  const unsigned int order = Aj.size() - 1;

  
  auto bisection_function = [&] (const double_type &x)
    {
      double_type f_j = Aj[order]; //Horner scheme

      for (int k = order-1; k >= 1; --k)
        f_j = Aj[k] + x * f_j;

      f_j = Aj[0] + x * f_j;


      return -f_j;
    };


  //initial density interval
  std::vector<double_type> x(2, 0.0);

  x[1] = gas_density;
  x[0] = options->element_density_minlimit;

  unsigned int nb_iterations = options->nb_max_bisection_iter;
  bool converged = false;


  for (unsigned int iter_step = 0; iter_step < nb_iterations; ++iter_step)
  { 
    const double_type x_n = (x[1] - x[0]) * 0.5 + x[0];

    const double_type f_n = bisection_function(x_n);

    if (f_n < 0)
      x[1] = x_n;
    else
      x[0] = x_n;

    //Convergence test. We need to be a little more accurate than the required accuracy.
    //Otherwise FastChem doesn't converge to the desired accuracy.
    if ( std::fabs(x[0] - x[1])/x[1] < options->accuracy * 1e-3 )
    {
      converged = true;
      break;
    }

  }


  //species.number_density = std::exp(x[0]);
  species.number_density = x[0];


  if (!converged && options->verbose_level >= 3)
    std::cout << "Bisection iteration limit reached, result may not be optimal." << "\t" << x[0] << "\t" << x[1]
              << "\t" << std::fabs(std::exp(x[0]) - std::exp(x[1]))/std::exp(x[1]) << "\t" << options->accuracy * 1e-3  << "\n";


  return converged;
}


template class FastChemSolver<double>;
template class FastChemSolver<long double>;

}
