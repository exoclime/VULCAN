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



//Nelder-Mead downhill simplex method in one dimension for the electron density
template <class double_type>
bool FastChemSolver<double_type>::nelderMeadSolveElectron(Element<double_type>& species, std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules, 
                                                          const double_type initial_solution, const double gas_density)
{
  const unsigned int N = 1; //dimension of Nelder-Mead method


  std::vector<double_type> Aj_cation(order_cation+1, 0.0);
  std::vector<double_type> Aj_anion(order_anion+1, 0.0);


  for (int k=1; k<order_cation+1; ++k)
    Aj_cation[k] = AmCoeffElectron(species, elements, molecules, -k);

  for (int k=1; k<order_anion+1; ++k)
    Aj_anion[k] = AmCoeffElectron(species, elements, molecules, k);


  //The function we want to minimise
  auto charge_conservation = [&] (const double_type &log_x)
    {
      const double_type x = std::exp(log_x);
      
      //Anions are calculated with Horner's method
      double_type P_anion = Aj_anion[order_anion];

      for (int k = order_anion-1; k >= 1; --k)
        P_anion = Aj_anion[k] + x * P_anion;

      //Contribution of cations
      //Due to the negative exponents, we can't use Horner's method here
      double_type P_cation = 0.0;

      for (int k=1; k<order_cation+1; k++)
        P_cation += Aj_cation[k] * std::pow(x, -k);


      //charge conservation
      const double_type Pj = - x - (x * P_anion + P_cation);
  

      return Pj;
    };




  //construct initial simplex
  std::vector<double_type> x;

  double_type initial_distance = (1.0 + 0.05) * initial_solution;

  double_type simplex_point = initial_distance;
  x.push_back(simplex_point);

  x.push_back(initial_solution);



  std::vector<double_type> vf(1+1,0);   //values of function evaluated at simplex vertexes

  unsigned int x1 = 0;    //index of best solution
  unsigned int xn = 0;    //index of second worst solution
  unsigned int xnp1 = 0;  //index of worst solution



  const double_type rho = 1.0, chi = 2.0, psi = 0.5, sigma = 0.5; //standard coefficients for Nelder-Mead method


  bool converged = false;


  //downhill simplex method starts
  for (unsigned int iter_step = 0; iter_step < options->nb_max_neldermead_iter; ++iter_step)
  {

    for(size_t i = 0; i < N + 1; ++i)
      vf[i] = std::fabs(charge_conservation(x[i]));
	    

    x1 = 0; xn = 0; xnp1 = 0;

    for (size_t i=0; i<vf.size(); ++i)
    {
      if(vf[i] < vf[x1]) x1 = i;

      if(vf[i] > vf[xnp1]) xnp1 = i;
    }

    xn = xnp1;  //in 1D they are equal


    const double_type xg = x[x1]; //xg: centroid of the N best vertexes; in 1D corresponds to the solution x1


    //check if the function has a root in a delta region around xg
    const double_type delta = xg * options->accuracy*1e-4;

 
    const double_type vf_epsilon_plus = charge_conservation(xg+delta);
    const double_type vf_epsilon_minus = charge_conservation(xg-delta);

    
    //if the function changes sign, we have the solution
    if ((vf_epsilon_minus < 0 && vf_epsilon_plus > 0) || (vf_epsilon_minus > 0 && vf_epsilon_plus < 0)  )
    {
      converged = true;

      break;
    }


    //reflection:
    const double_type xr = (1.0 + rho) * xg - rho*x[xnp1];

    double_type fxr = std::fabs(charge_conservation(xr));


    if (fxr < vf[x1])
    {
      //expansion
      const double_type xe = (1.0 + rho*chi) * xg - rho*chi*x[xnp1];
      const double_type fxe = std::fabs(charge_conservation(xe));

      if (fxe < fxr)
        x[xnp1] = xe;
      else
        x[xnp1] = xr;
    }
    else
    {
      if (fxr < vf[xn])
        x[xnp1] = xr;
      else
      {
        bool perform_shrink = false;

        if (fxr < vf[xnp1])
        {
          //outward contraction
          const double_type xc = (1.0 + psi*rho) * xg - psi*rho*x[xnp1];
          const double_type fxc = std::fabs(charge_conservation(xc));

          if (fxc <= fxr)
            x[xnp1] = xc;
          else
            perform_shrink = true;

        }
        else
        {
          //inward contraction
          const double_type xc = (1.0 - psi) * xg + psi*rho*x[xnp1];
          const double_type fxc = std::fabs(charge_conservation(xc));


          if (fxc < vf[xnp1])
            x[xnp1] = xc;
          else
            perform_shrink = true;
        }

        if (perform_shrink)
        {
          for(size_t i=0; i<x.size(); ++i )
          {

            //total contraction
            if (i != x1)
              x[i] = x[x1] + sigma*(x[i] - x[x1]);

          }

        }

      }

    }


  }//optimisation is finished


  species.number_density = std::exp(x[x1]);

  
  if (!converged && options->verbose_level >= 3)
    std::cout << "Nelder-Mead iteration limit reached, result may not be optimal." << "\t" << x[x1] << "\n";


  return converged;
}






template class FastChemSolver<double>;
template class FastChemSolver<long double>;

}
