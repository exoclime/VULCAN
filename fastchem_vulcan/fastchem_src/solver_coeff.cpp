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
#include "species_struct.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


namespace fastchem {

template <class double_type>
FastChemSolver<double_type>::FastChemSolver(FastChemOptions<double_type>* options_ptr) 
{

  options = options_ptr; 
  
}


template <class double_type>
double_type FastChemSolver<double_type>::A1Coeff(const Element<double_type>& species, const std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules)
{
  const unsigned int index = species.index;

  //calculation of coefficient A_1, see Eq. (2.28)
  double_type A1 = 0.0;


  for (size_t j=0; j<species.molecule_list.size(); ++j)
  {
    const unsigned int i = species.molecule_list[j];

    
    if (molecules[i].stoichometric_vector[index] == 1 && molecules[i].abundance == species.abundance)
    {
      double_type sum = 0;

      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        unsigned int l = molecules[i].element_indices[k];
  

        if (l != species.index && molecules[i].stoichometric_vector[l] != 0)
          sum += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density);
      }

      const double_type kappa = 1.0 + species.epsilon * molecules[i].sigma;
      
      A1 += std::exp(molecules[i].mass_action_constant + sum) * kappa;
    }
  }

  
  return A1;

}



template <class double_type>
double_type FastChemSolver<double_type>::A2Coeff(const Element<double_type>& species, const std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules)
{
  const unsigned int index = species.index;

  double_type A2 = 0.0;


  for (size_t j=0; j<species.molecule_list.size(); ++j)
  {
    const unsigned int i = species.molecule_list[j];


    if (molecules[i].stoichometric_vector[index] == 2 && molecules[i].abundance == species.abundance)
    {
      double_type sum = 0;

      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        const unsigned int l = molecules[i].element_indices[k];

        if (l != species.index && molecules[i].stoichometric_vector[l] != 0)
          sum += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density);

      }


      const double_type kappa = 2.0 + species.epsilon * molecules[i].sigma; 
      
      A2 += std::exp(molecules[i].mass_action_constant + sum) * kappa;
    }
  }
  
  
  return A2;
}



template <class double_type>
double_type FastChemSolver<double_type>::AmCoeff(const Element<double_type>& species, const std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules, 
                                                 const unsigned int order)
{
  const unsigned int index = species.index;

  double_type Am = 0.0;


  for (size_t j=0; j<species.molecule_list.size(); ++j)
  {
    const unsigned int i = species.molecule_list[j];


    if (molecules[i].stoichometric_vector[index] == int(order) && molecules[i].abundance == species.abundance)
    {
      double_type sum = 0;

      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        const unsigned int l = molecules[i].element_indices[k];

        if (l != species.index && molecules[i].stoichometric_vector[l] != 0)
          sum += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density);

      }


      const double_type kappa = order + species.epsilon * molecules[i].sigma; 
      
      Am += std::exp(molecules[i].mass_action_constant + sum) * kappa;
    }
  }
  
  
  return Am;
}



//Alternative description of the Am coefficients
//Takes the full law of mass action into account, i.e. doesn't stop at minor species as the regular calculation
template <class double_type>
double_type FastChemSolver<double_type>::AmCoeffAlt(const Element<double_type>& species, const std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules, 
                                                    const unsigned int order)
{
  const unsigned int index = species.index;

  double_type Am = 0.0;


  for (size_t j=0; j<species.molecule_list.size(); ++j)
  {
    const unsigned int i = species.molecule_list[j];


    if (molecules[i].stoichometric_vector[index] == int(order))
    {
      double_type sum = 0;

      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        const unsigned int l = molecules[i].element_indices[k];

        if (l != species.index && molecules[i].stoichometric_vector[l] != 0)
          sum += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density);

      }


      const double_type kappa = order + species.epsilon * molecules[i].sigma; 
      
      Am += std::exp(molecules[i].mass_action_constant + sum) * kappa;
    }
  }
  
  
  return Am;
}



//Alternative description of the Am coefficients
//Takes the full law of mass action into account, i.e. doesn't stop at minor species as the regular calculation
template <class double_type>
double_type FastChemSolver<double_type>::AmCoeffElectron(const Element<double_type>& electron, const std::vector< Element<double_type> >& elements, const std::vector< Molecule<double_type> >& molecules, 
                                                         const int order)
{
  const unsigned int index = electron.index;

  double_type Am = 0.0;


  for (size_t j=0; j<electron.molecule_list.size(); ++j)
  {
    const unsigned int i = electron.molecule_list[j];


    if (molecules[i].stoichometric_vector[index] == order)
    {
      double_type sum = 0;

      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        const unsigned int l = molecules[i].element_indices[k];

        if (l != electron.index && molecules[i].stoichometric_vector[l] != 0)
          sum += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density);

      }

      
      Am += std::exp(molecules[i].mass_action_constant + sum) * order;
    }
  }
  
  
  return Am;
}



template class FastChemSolver<double>;
template class FastChemSolver<long double>;

}



