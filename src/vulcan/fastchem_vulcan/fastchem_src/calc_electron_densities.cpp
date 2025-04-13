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
#include "solver.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


namespace fastchem {



//Calculation of the electron density
template <class double_type>
void FastChem<double_type>::calculateElectronDensities(Element<double_type>& electron, const double_type& old_number_density, const double_type gas_density)
{

  //Am I the electron? 
  if (electron.symbol != "e-") return;


  //no ions present   
  if (electron.molecule_list.size() == 0)
  {
    
    electron.number_density = 0.0;
    return; 

  } 


  //if we have't determined the maximum order of cations and anions, we do so now
  if (solver.order_anion == -999 && solver.order_cation == -999)
  {
    solver.order_cation = 0;

    for (unsigned int i=0; i<electron.molecule_list.size(); ++i)
      if (molecules[electron.molecule_list[i]].stoichometric_vector[electron.index] < solver.order_cation )
        solver.order_cation = molecules[electron.molecule_list[i]].stoichometric_vector[electron.index];

    solver.order_cation = std::abs(solver.order_cation);


    solver.order_anion = 0;

    for (unsigned int i=0; i<electron.molecule_list.size(); ++i)
      if (molecules[electron.molecule_list[i]].stoichometric_vector[electron.index] > solver.order_anion )
        solver.order_anion = molecules[electron.molecule_list[i]].stoichometric_vector[electron.index];

    solver.order_anion = std::abs(solver.order_anion);
  }



  //for singly-ionised species we use the analytic solution
  if (electron.solver_order == 1)
    calculateSinglyIonElectrons(electron, old_number_density);
  else
    calculateMultIonElectrons(electron, old_number_density, gas_density);
  
}


//Calculation of the electron density for at most singly-ionised species
//Uses the analytical solution from Paper 1, Appendix B
template <class double_type>
void FastChem<double_type>::calculateSinglyIonElectrons(Element<double_type>& electron, const double_type& old_number_density)
{
 
  double_type alpha = 0.0;
  double_type beta = 0.0;

  const unsigned int index = electron.index;

  for (size_t j=0; j<electron.molecule_list.size(); ++j)
  {
    const unsigned int i = electron.molecule_list[j];

    //the anions, Eq. (B3) in Paper 1
    if (molecules[i].stoichometric_vector[index] == 1)
    {
      double_type sum = 0;

      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        unsigned int l = molecules[i].element_indices[k];
  

        if (l != electron.index && molecules[i].stoichometric_vector[l] != 0)
          sum += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density);
      }

      
      
      beta += std::exp(molecules[i].mass_action_constant + sum);
    }
    else if (molecules[i].stoichometric_vector[index] == -1)  //the cations, Eq. (B4) in Paper 1
    {
      double_type sum = 0;

      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        unsigned int l = molecules[i].element_indices[k];
  

        if (l != electron.index && molecules[i].stoichometric_vector[l] != 0)
          sum += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density);
      }

      
      
      alpha += std::exp(molecules[i].mass_action_constant + sum);
    }

  }

  //Eq. (B2) in Paper 1
  double_type electron_density = std::sqrt(alpha/(1.0 + beta));

  elements[e_].number_density = electron_density; 

}




//Calculation of the electron density, based on charge conservation
//This approach is used for multi-ionised species
//First tries to estimate the electron density via Paper 1, Eq. (2.35).
//In case that fails (electron density not sufficiently high enough), it switches to a 1D Newton's method.
//See Sect. 2.4.3 for details.
template <class double_type>
void FastChem<double_type>::calculateMultIonElectrons(Element<double_type>& electron, const double_type& old_number_density, const double_type& gas_density)
{

  electron.number_density = 0.0;


  double_type positive_ion_density = 0;
  double_type negative_ion_density = 0;


  for (size_t i=0; i<electron.molecule_list.size(); ++i)
    if (molecules[electron.molecule_list[i]].stoichometric_vector[e_] > 0)
      negative_ion_density += molecules[electron.molecule_list[i]].stoichometric_vector[e_] * molecules[electron.molecule_list[i]].number_density;
    else
      positive_ion_density -= molecules[electron.molecule_list[i]].stoichometric_vector[e_] * molecules[electron.molecule_list[i]].number_density;


  double_type electron_density = positive_ion_density - negative_ion_density;


  double_type delta = 0.9;
  
  if (electron_density > delta*positive_ion_density)
  {

    electron.number_density = std::sqrt(electron_density * old_number_density);

  }
  else
  {

    //switching to Newton's method
    solver.newtonSolElectron(electron, elements, molecules, gas_density);

  }

}


template class FastChem<double>;
template class FastChem<long double>;

}



