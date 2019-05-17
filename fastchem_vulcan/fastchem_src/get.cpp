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

#include <string>
#include <vector>



namespace fastchem {


//Query for a basic chemical element index with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the element does not exist
template <class double_type>
unsigned int FastChem<double_type>::getChemicalElementIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<nb_chemical_element_data; ++i)
    if (symbol == chemical_element_data[i].symbol)
    {
      index = i;
      break;
    }


  return index;
}



//Query for a molecule index with chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the molecule does not exist
template <class double_type>
unsigned int FastChem<double_type>::getMoleculeIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<molecules.size(); ++i)
    if (symbol == molecules[i].symbol)
    {
      index = i;
      break;
    }


  return index;
}



//Query for a element's index with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the element does not exist
template <class double_type>
unsigned int FastChem<double_type>::getElementIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<elements.size(); ++i)
    if (symbol == elements[i].symbol)
    {
      index = i;
      break;
    }


  return index;
}



//Query for a species index (both, elements and molecules) with a chemical symbol
//Returns FASTCHEM_UNKNOWN_SPECIES if the species does not exist
template <class double_type>
unsigned int FastChem<double_type>::getSpeciesIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<nb_species; ++i)
    if (symbol == species[i]->symbol)
    {
      index = i;
      break;
    }


  return index;
}



//Query for a species name with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getSpeciesName(const unsigned int species_index)
{
  if (species_index < nb_species)
    return species[species_index]->name;
  else
    return "";
}



//Query for a species symbol with its index
//Returns empty string in case the species does not exist
template <class double_type>
std::string FastChem<double_type>::getSpeciesSymbol(const unsigned int species_index)
{

  if (species_index < nb_species)
    return species[species_index]->symbol;
  else
    return "";

}



//Query for an element name with its index
//Returns empty string in case the element does not exist
template <class double_type>
std::string FastChem<double_type>::getElementName(const unsigned int species_index)
{
  if (species_index < nb_elements)
    return elements[species_index].name;
  else
    return "";
}



//Query for an element symbol with its index
//Returns empty string in case the element does not exist
template <class double_type>
std::string FastChem<double_type>::getElementSymbol(const unsigned int species_index)
{

  if (species_index < nb_elements)
    return elements[species_index].symbol;
  else
    return "";

}



//Query for a the molecular weight of a species with its index
//Returns 0 in case the species does not exist
template <class double_type>
double FastChem<double_type>::getSpeciesMolecularWeight(const unsigned int species_index)
{

  if (species_index < nb_species)
    return species[species_index]->molecular_weight;
  else
    return 0.;

}


//Get the element abundance for a specific element
template <class double_type>
double FastChem<double_type>::getElementAbundance(const unsigned int species_index)
{

  if (species_index < nb_elements)
    return elements[species_index].abundance;
  else
    return 0.;

}



//Get the element abundandes for all elements
template <class double_type>
std::vector<double> FastChem<double_type>::getElementAbundance()
{

  std::vector<double> abundances(nb_elements, 0.0);

  for (size_t i=0; i<nb_elements; ++i)
    abundances[i] = elements[i].abundance;

  return abundances;

}



template class FastChem<double>;
template class FastChem<long double>;

}
