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

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include <cmath>

#include <limits>


namespace fastchem {




template <class double_type>
void FastChem<double_type>::init()
{
  //read input files
  if (    readElementList() == false
       || readElementAbundances() == false
       || readSpeciesData() == false)
  {
    std::cout << "FastChem data file loading failed!\n";

    is_initialized = false;

    return;
  }



  if (std::numeric_limits<double_type>::max_exponent10 > 1000)
  {
    options.element_density_minlimit = 1e-512L;
    options.molecule_density_minlimit = 1e-512L;
  }
  else
  {
    options.element_density_minlimit = 1e-155;
    options.molecule_density_minlimit = 1e-155;
  }



  //setting up the list of species
  for (size_t i=0; i<nb_elements; ++i)
    species.push_back(&elements[i]);

  for (size_t i=0; i<nb_molecules; ++i)
    species.push_back(&molecules[i]);

  nb_species = nb_elements + nb_molecules;

  setMoleculeAbundances();

  determineSolverOrder();

  determineElementCalculationOrder();

  for (auto & i : elements) i.calcEpsilon(elements);


  //some diagnostic output in case you'd like to see it
  if (options.verbose_level >= 4)
  {
    std::cout << "molecule list: \n";
    for (size_t i=0; i<molecules.size(); ++i)
    {
      std::cout << molecules[i].name << "\t" << molecules[i].symbol << std::endl;

      for (size_t j=0; j<molecules[i].mass_action_coeff.size(); ++j)
        std::cout << molecules[i].mass_action_coeff[j] << "\t";

      std::cout << "\n";

      for (size_t j=0; j<molecules[i].stoichometric_vector.size(); ++j)
        std::cout << molecules[i].stoichometric_vector[j] << " ";

      std::cout << "\n";

      std::cout << "elements: ";
      for (size_t j=0; j<molecules[i].element_indices.size(); ++j)
        std::cout << elements[molecules[i].element_indices[j]].symbol << ", index: " << molecules[i].element_indices[j] << "; ";

      std::cout << "\n";

      std::cout << "charge: " << molecules[i].charge << "\n";
    }

    std::cout << "\n element list: \n";
    for (size_t i=0; i<elements.size(); ++i)
    {
      std::cout << elements[i].name << "\t" << elements[i].symbol << "\t" << elements[i].element_data_index << "\t" << elements[i].solver_order << "\t" << elements[i].abundance << "\n";
    }


    std::cout << "calculation order:\n";
      for (size_t i=0; i<element_calculation_order.size(); ++i)

    std::cout << elements[element_calculation_order[i]].symbol << "\t" << elements[element_calculation_order[i]].abundance << "\n";

    std::cout << "\n";
  }


  if (options.verbose_level >= 2)
  {
    std::cout << "considered species:\n";

    for (size_t i=0; i<nb_species; ++i)
      std::cout << species[i]->symbol << "\t" << species[i]->name << "\t" << species[i]->abundance << "\t" << species[i]->molecular_weight << "\n";

    std::cout << "\n";
  }


  e_ = getElementIndex("e-");


  if (options.verbose_level > 0)
    std::cout << "number of species: " << nb_species
              << "  elements: " << nb_elements
              << "  molecules: " << nb_molecules
              << "  chemical elements: " << nb_chemical_element_data << std::endl;


  is_initialized = true;
}



//Reinitializes certain internal data after the element abundances were changed
template <class double_type>
void FastChem<double_type>::reInitialiseFastChem()
{
  //reset the calculation order
  element_calculation_order.resize(0);

  //update the abundances of the molecules
  setMoleculeAbundances();

  //order the elements according to their abundances
  determineElementCalculationOrder();

  //update the solver order for the new abundances
  determineSolverOrder();
  
  //recalculate the epsilons
  for (auto & i : elements) i.calcEpsilon(elements);
}



template <class double_type>
void FastChem<double_type>::setMoleculeAbundances()
{
  
  for (size_t i=0; i<nb_elements; ++i)
    for (size_t j=0; j<nb_elements; ++j)
    {
      if (i==j) continue;

      if (elements[i].abundance == elements[j].abundance)
        elements[j].abundance += std::numeric_limits<double_type>::epsilon()*elements[j].abundance;
    }

  
  //update the definition of the molecular abundances
  for (size_t i=0; i<nb_molecules; ++i)
  {
    molecules[i].abundance = 1e42;

    for (size_t j=0; j<molecules[i].element_indices.size(); ++j)
      if (molecules[i].abundance > elements[molecules[i].element_indices[j]].abundance && elements[molecules[i].element_indices[j]].symbol != "e-")
        molecules[i].abundance = elements[molecules[i].element_indices[j]].abundance;

    //scaled abundances
    molecules[i].abundance_scaled = 1e42;

    for (size_t j=0; j<molecules[i].element_indices.size(); ++j)
    {
      unsigned element_index = molecules[i].element_indices[j];

      if (molecules[i].abundance_scaled > elements[element_index].abundance/molecules[i].stoichometric_vector[element_index]
          && elements[molecules[i].element_indices[j]].symbol != "e-")
        molecules[i].abundance_scaled = elements[element_index].abundance/molecules[i].stoichometric_vector[element_index];
    }
  }

  createMoleculeLists();
}




template <class double_type>
void FastChem<double_type>::createMoleculeLists()
{
  
  //reset the lists
  for (auto & i : elements) i.major_molecules_inc.resize(0);
  for (auto & i : elements) i.major_molecules_exc.resize(0);
  for (auto & i : elements) i.minor_molecules.resize(0);


  
  for (size_t i=0; i<nb_molecules; ++i)
  {

    for (size_t j=0; j<nb_elements; ++j)
    {
       if (elements[j].abundance > molecules[i].abundance)
         elements[j].minor_molecules.push_back(i);
       else
       {

         if (molecules[i].stoichometric_vector[j] == 0)
           elements[j].major_molecules_exc.push_back(i);
         else
           elements[j].major_molecules_inc.push_back(i);
    
       }
    }
  }


  if (options.verbose_level >= 4)
  {
    std::cout << "Molecule lists: \n";

    for (size_t j=0; j<nb_elements; ++j)
    {
      std::cout << "element " << elements[j].symbol << "\n";

      std::cout << "major elements inc:\n";
      for (size_t i=0; i<elements[j].major_molecules_inc.size(); ++i)
        std::cout << molecules[elements[j].major_molecules_inc[i]].symbol << "\n";

      std::cout << "major elements exc:\n";
      for (size_t i=0; i<elements[j].major_molecules_exc.size(); ++i)
        std::cout << molecules[elements[j].major_molecules_exc[i]].symbol << "\n";

      std::cout << "minor elements:\n";
      for (size_t i=0; i<elements[j].minor_molecules.size(); ++i)
        std::cout << molecules[elements[j].minor_molecules[i]].symbol << "\n";
    }

  }
  

}




template class FastChem<double>;
template class FastChem<long double>;


}
