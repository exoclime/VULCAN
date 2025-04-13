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



//Newton's method so solve for element densities
//See Sect. 2.4.2, Paper 1
//The standard version (use_alternative = false) uses the n_min approach to account for minor species
//For use_alternative = true, it will use all species in the law of mass action
template <class double_type>
void FastChemSolver<double_type>::newtonSol(Element<double_type>& species, std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules, 
                                            const double_type gas_density, const bool use_alternative)
{
  double_type scaling_factor = 0.0;
  
  //in case we use the scaling factor, see Appendix A for details
  if (options->use_scaling_factor)
    scaling_factor = solverScalingFactor(species, gas_density); // + std::numeric_limits<double_type>::min_exponent/6.0;


  unsigned int order = 0;

  //calculation of the polynomial coefficients
  std::vector<double_type> Aj; 

  //The standard case using the n_min
  if (!use_alternative)
  {
    order = species.solver_order;

    Aj.assign(order+1, 0.0);

    Aj[0] = std::exp(-scaling_factor) * (species.number_density_maj + species.number_density_min - gas_density * species.epsilon);
    Aj[1] = A1Coeff(species, elements, molecules) + std::exp(-scaling_factor);

    for (size_t k=2; k<order+1; ++k)
      Aj[k] = AmCoeff(species, elements, molecules, k);
  }
  else  //the general case for the alternative version
  {
    for (unsigned int i=0; i<species.molecule_list.size(); ++i)
      if (molecules[species.molecule_list[i]].stoichometric_vector[species.index] > int(order) )
        order = molecules[species.molecule_list[i]].stoichometric_vector[species.index];

    Aj.assign(order+1, 0.0);

    double_type n_exc = 0.0;

    for (size_t i=0; i<molecules.size(); ++i)
      if (molecules[i].stoichometric_vector[species.index] == 0)
        n_exc += molecules[i].sigma * molecules[i].number_density;
  
    n_exc *= species.epsilon;

    Aj[0] = std::exp(-scaling_factor) * (n_exc - gas_density * species.epsilon);
    Aj[1] = AmCoeffAlt(species, elements, molecules, 1) + std::exp(-scaling_factor);

    for (size_t k=2; k<order+1; ++k)
      Aj[k] = AmCoeffAlt(species, elements, molecules, k);
  }

  

  //Newton's method
  bool converged = false;

  //double_type x = gas_density; //initial guess, ensures monotonous convergence.


  double_type x = gas_density;

  if (species.number_density == 0)
    x = gas_density;
  else
    x = species.number_density;


  //one Newton step as lambda function
  auto newton_step = [&] (const double_type &x)
    {
      double_type P_j = Aj[order];        //Horner scheme
      double_type P_j_prime = order*Aj[order];

      for (int k = order-1; k >= 1; --k)
      {
        P_j = Aj[k] + x * P_j;
        P_j_prime = k * Aj[k] + x * P_j_prime;
      }

      P_j = Aj[0] + x * P_j;


      return x - P_j/P_j_prime; //Newton step
    };



  //Newton iteration
  unsigned int mu = 0; 
  for (mu=0; mu<options->nb_max_newton_iter; ++mu)
  {
    double_type x_new = newton_step(x);


    if (std::fabs(x_new - x) < options->newton_err * std::fabs(x_new))  //root found?
    {
      x = x_new;
      converged = true;

      break;
    }


    //prevent x to become negative due to numerical underflow
    if (x_new < 1.e-8*x)
    {
      x_new = 1.e-8*x;
    }


    x = x_new;
  }


  
  //test if root is in (max(0,x*(1-newton_err)),x*(1+newton_err))
  double_type x_lower = std::fmax(0., x * (1. - options->newton_err));
  double_type x_upper = x * (1. + options->newton_err);

  double_type P_j_lower = Aj[order];
  double_type P_j_upper = Aj[order];

  for(int k = order-1; k >=0 ; k--)
  {
    P_j_lower = Aj[k] + x_lower * P_j_lower;
    P_j_upper = Aj[k] + x_upper * P_j_upper;
  }


  if (converged)
    species.number_density = x;

  
  //in case the normal Newton solver does not work, we switch to other solvers
  if (x < 0 || !converged || P_j_lower*P_j_upper > 0.)
  { 
    //if the normal Newton's method does not converge, switch to the alternative version of it
    if (!use_alternative)
    { 
      newtonSol(species, elements, molecules, gas_density, true);

      if (options->verbose_level >= 3)
        std::cout << "FastChem: WARNING: NewtSol failed for species " << species.symbol << " switched to Backup " << x << "\t" << species.number_density << "\n";
    }
    else //if the alternative Newton's method also doesn't work, we use the bisection method
    {
      bisectionSolve(species, Aj, gas_density);


      if (options->verbose_level >= 3)
        std::cout << "FastChem: WARNING: NewtSol Alt failed for species " << species.symbol << " switched to Bisection as backup " << x << "\t" << species.number_density << "\n";
    }
      
    
  }

}  




//Newton's method for the electrons
//Instead of element conservation, solves for charge balance
template <class double_type>
void FastChemSolver<double_type>::newtonSolElectron(Element<double_type>& species, std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules, 
                                                    const double_type gas_density)
{

  //Calculation of the polynomial coefficients
  std::vector<double_type> Aj_cation(order_cation+1, 0.0);
  std::vector<double_type> Aj_anion(order_anion+1, 0.0);

  for (int k=1; k<order_cation+1; ++k)
    Aj_cation[k] = AmCoeffElectron(species, elements, molecules, -k);

  for (int k=1; k<order_anion+1; ++k)
    Aj_anion[k] = AmCoeffElectron(species, elements, molecules, k);


  //Newton's method
  bool converged = false;
  double_type x = order_cation/(1.0 + order_cation) * gas_density; //Initial guess ensures monotonous convergence.

 
  //one Newton step as lambda function
  auto newton_step = [&] (const double_type &x)
    {
      //Horner's method for the anions
      double_type P_anion = Aj_anion[order_anion];
      double_type P_prime_anion = order_anion*Aj_anion[order_anion];

      for (int k = order_anion-1; k >= 1; --k)
      {
        P_anion = Aj_anion[k] + x * P_anion;
        P_prime_anion = k * Aj_anion[k] + x * P_prime_anion;
      }

      //The cations
      double_type P_cation = 0.0;
      double_type P_prime_cation = 0.0;

      for (int k=1; k<order_cation+1; k++)
      {
        P_cation += Aj_cation[k] * std::pow(x, -k);
        P_prime_cation += -k * Aj_cation[k] * std::pow(x, -k-1);
      }

      const double_type P_j = x  + x * P_anion + P_cation;  //this is the charge balance
      const double_type P_j_prime = 1.0 + P_prime_cation + P_prime_anion; //derivative

      return x - P_j/P_j_prime; //Newton step
    };


  //Newton iteration
  for (unsigned int mu=0; mu<options->nb_max_newton_iter; ++mu)
  {
    double_type x_new = newton_step(x);

    if (std::fabs(x_new - x) <= options->newton_err * std::fabs(x_new))  //root found?
    {
      x = x_new;
      converged = true;

      break;
    }


    //prevent x to become negative due to numerical underflow
    if (x_new < 1.e-8*x) x_new = 1.e-8*x;
    
    x = x_new;


    if (std::isnan(x)) break;
  }


  // Test if root is in (max(0,x*(1-newton_err)),x*(1+newton_err))
  const double_type x_lower = std::fmax(0., x * (1. - options->newton_err));
  const double_type x_upper = x * (1. + options->newton_err);


  double_type P_anion_lower = Aj_anion[order_anion];
  double_type P_anion_upper = Aj_anion[order_anion];

  for (int k = order_anion-1; k >= 1; --k)
  {
    P_anion_lower = Aj_anion[k] + x_lower * P_anion_lower;
    P_anion_upper = Aj_anion[k] + x_upper * P_anion_upper;
  }

  double_type P_cation_lower = 0.0;
  double_type P_cation_upper = 0.0;

  for (int k=1; k<order_cation+1; k++)
  {
    P_cation_lower += Aj_cation[k] * std::pow(x_lower, -k);
    P_cation_upper += Aj_cation[k] * std::pow(x_upper, -k);
  }

  const double_type P_j_lower = x_lower  + x_lower * P_anion_lower + P_cation_lower;
  const double_type P_j_upper = x_upper  + x_upper * P_anion_upper + P_cation_upper;



  species.number_density = x; 
  


  //in case something went wrong again, we try to use another backup
  if (x < 0 || !converged || P_j_lower*P_j_upper > 0.)
  {
    const double_type init = std::log(std::fabs(x));
    nelderMeadSolveElectron(species, elements, molecules, init, 0.0);


    if (options->verbose_level >= 3)
      std::cout << "FastChem: WARNING: NewtSol failed for electrons, switching to Nelder-Mead Backup " << x << "\t" << species.number_density << "\n";
  }


}




template class FastChemSolver<double>;
template class FastChemSolver<long double>;

}



