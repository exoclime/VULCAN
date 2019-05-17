
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


#ifndef _species_struct_h
#define _species_struct_h


#include <vector>
#include <string>


namespace fastchem {


//Class for the basic chemical element data
//Those are read in from the elemental abundance file
//and are later cross-correlated with the elements in the network
template <class double_type>
struct ChemicalElementData
{
  std::string symbol;
  std::string name;

  double_type atomic_weight = 0.0;
  double_type abundance = 0.0;   //element abundance
};



//Parent class of all species
template <class double_type>
struct ChemicalSpecies
{
  std::string symbol;
  std::string name;

  double_type molecular_weight = 0.0;
  int charge = 0;

  double_type abundance = 0.0;      //element abundance
  double_type number_density = 0.0;
};



//Class for the molecules
template <class double_type>
struct Molecule : public ChemicalSpecies<double_type>
{
  std::vector<unsigned int> element_indices;
  std::vector<int> stoichometric_vector;

  std::vector<double_type> mass_action_coeff;
  double_type mass_action_constant = 0.0;

  double_type abundance_scaled = 0.0;

  double_type sigma = 0.0;
  double_type sum = 0.0;

  void calcMassActionConstant(const double temperature);
  void checkN(const double_type& min_limit, const double_type& gas_density);
};



//Class for the elements
template <class double_type>
struct Element : public ChemicalSpecies<double_type>
{
  unsigned int element_data_index;
  unsigned int index;

  unsigned int solver_order = 0;

  std::vector<unsigned int> molecule_list;         //contains the list of molecule indices the element is part of
  std::vector<unsigned int> major_molecules_inc;   //contains the list of 'major' molecules that includes the element
  std::vector<unsigned int> major_molecules_exc;   //contains the list of 'major' molecules that doesn't includes the element
  std::vector<unsigned int> minor_molecules;       //contains the list of 'minor' molecules
  unsigned int element_conserved = 0;              //check if element is conserved during calculation, for electrons this is charge conservation

  double_type epsilon = 0.0;
  double_type number_density_maj = 0.0;
  double_type number_density_min = 0.0;

  void calcMinorSpeciesDensities(const std::vector< Molecule<double_type> > &molecules);
  void calcEpsilon(const std::vector< Element<double_type> > &elements);
  
  void checkN(const double_type& min_limit, const double_type& gas_density);
  bool checkElementConservation(const std::vector< Molecule<double_type> >& molecules, const double_type total_density, const double_type& accuracy);
  bool checkChargeConservation(const std::vector< Molecule<double_type> >& molecules, const double_type& accuracy);
};



}


#endif