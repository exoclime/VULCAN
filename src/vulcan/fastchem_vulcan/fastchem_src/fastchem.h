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


#ifndef _fastchem_h
#define _fastchem_h

#include <vector>
#include <iostream>
#include <string>

#include "fastchem_constants.h"
#include "species_struct.h"
#include "input_output_struct.h"

#include "options.h"
#include "solver.h"


namespace fastchem {


//FastChem class
template <class double_type>
class FastChem{
  public:
    FastChem(const std::string& model_parameter_file, const unsigned int verbose_level_init);

    FastChem(const FastChem &obj);

    //function calls to calculate number densities
    unsigned int calcDensities(FastChemInput& input, FastChemOutput& output);

    //public query functions
    std::string getSpeciesName(const unsigned int species_index);
    std::string getSpeciesSymbol(const unsigned int species_index);
    unsigned int getSpeciesNumber() {return nb_species;}
    unsigned int getSpeciesIndex(const std::string symbol);

    std::string getElementName(const unsigned int species_index);
    std::string getElementSymbol(const unsigned int species_index);
    unsigned int getElementNumber() {return nb_elements;}

    double getElementAbundance(const unsigned int species_index);
    std::vector<double> getElementAbundance();

    double getSpeciesMolecularWeight(const unsigned int species_index);


    //functions to set internal variables during runtime
    //they will override any read-in values
    void setElementAbundance(std::vector<double> abundances);

    void setVerboseLevel(const unsigned int level) { if (level > 4) options.verbose_level = 4; else options.verbose_level = level;}
    
    void setMaxChemistryIter(const unsigned int nb_steps) {options.nb_max_fastchem_iter = nb_steps;}
    void setMaxNewtonIter(const unsigned int nb_steps) {options.nb_max_newton_iter = nb_steps;}
    void setMaxBisectionIter(const unsigned int nb_steps) {options.nb_max_bisection_iter = nb_steps;}
    void setMaxNelderMeadIter(const unsigned int nb_steps) {options.nb_max_neldermead_iter = nb_steps;}

    void setChemistryAccuracy(const double chem_accuracy) {options.accuracy = chem_accuracy;}
    void setNewtonAccuracy(const double newton_accuracy) {options.newton_err = newton_accuracy;}

    void useScalingFactor(const bool use_switch) {options.use_scaling_factor = use_switch;}


  private:
    FastChemOptions<double_type> options;
    FastChemSolver<double_type> solver;


    unsigned int nb_chemical_element_data = 0;
    unsigned int nb_species = 0;
    unsigned int nb_molecules = 0;
    unsigned int nb_elements = 0;

    unsigned int e_ = FASTCHEM_UNKNOWN_SPECIES; //electron element index


    bool is_initialized = false;
    bool is_busy = false;


    std::vector< ChemicalElementData<double_type> > chemical_element_data;

    std::vector< ChemicalSpecies<double_type>* > species;
    std::vector< Element<double_type> > elements;
    std::vector< Molecule<double_type> > molecules;

    std::vector<unsigned int> element_calculation_order;

    //Initialisation functions
    void init();

    bool readElementList();
    bool readElementAbundances();
    void setElementAbundance(const std::string symbol, const double abundance);
    void setMoleculeAbundances();
  
    bool readSpeciesData();
    void addMolecule(const std::string name, const std::string symbol,
                     const std::vector<std::string> species_elements, const std::vector<int> stoichometric_coeff,
                     const std::vector<double_type> mass_action_coeff, const int charge);
    void addAtom(std::string symbol);

    void reInitialiseFastChem();

    unsigned int determineSolverOrder(const Element<double_type>& species);
    void determineSolverOrder();
    void determineElementCalculationOrder();

    void createMoleculeLists();


    //Internal query functions
    unsigned int getChemicalElementIndex(const std::string symbol);
    unsigned int getElementIndex(const std::string symbol);
    unsigned int getMoleculeIndex(const std::string symbol);


    //Functions for the calculations of the number densities
    unsigned int calcDensity(const double temperature, const double pressure, const bool use_previous_solution,
                             std::vector<double>& number_densities, double& total_element_density, 
                             double& mean_molecular_weight,
                             std::vector<unsigned int>& element_conserved,
                             unsigned int& nb_chemistry_iterations);


    bool solveFastchem(const double temperature_gas, const double gas_density, unsigned int& nb_iterations);

    void calculateElementDensities(Element<double_type>& species, const double_type gas_density,
                                   bool use_backup_solver, double_type& n_major);
    double_type calculateMoleculeDensities(Element<double_type>& species, const double_type gas_density);
    
    void calculateElectronDensities(Element<double_type>& species, const double_type& old_number_density, const double_type gas_density);
    void calculateSinglyIonElectrons(Element<double_type>& electron, const double_type& old_number_density);
    void calculateMultIonElectrons(Element<double_type>& electron, const double_type& old_number_density, const double_type& gas_density);
    
    
    double totalElementDensity();
    double meanMolecularWeight(const double gas_density);
};



}

#endif
