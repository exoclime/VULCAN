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

#include "../fastchem_src/fastchem.h"
#include "../fastchem_src/input_output_struct.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cmath>


//forward declaration
bool read_config_file(std::string file_path, std::string& fastchem_options_file, std::string& atmosphere_file,
                      std::string& chem_output_file, std::string& monitor_output_file,
                      unsigned int& verbose_level, bool& output_mixing_ratios);




int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    std::cout << "config file command line parameter missing!\n";

    return 1;
  }


  std::string config_file_name = argv[1];

  std::string fastchem_options_file;
  std::string atmosphere_file;
  std::string chem_output_file;
  std::string monitor_output_file;
  unsigned int verbose_level = 5;
  bool output_mixing_ratios = false;


  if (!read_config_file(config_file_name, fastchem_options_file,atmosphere_file,chem_output_file, monitor_output_file, verbose_level, output_mixing_ratios))
    return 1;



  fastchem::FastChem<long double> fastchem(fastchem_options_file, verbose_level);


  std::vector<double> temperature;
  std::vector<double> pressure;


  //read the input data
  std::fstream file(atmosphere_file.c_str(), std::ios::in);

  std::string line;

  //read in the comment header
  //std::getline(file, line);

  while(std::getline(file, line))
  {
    std::stringstream line_stream(line);

    double temperature_in;
    double pressure_in;


    if (!(line_stream >> pressure_in >> temperature_in)) continue;
    //if (!(line_stream >> temperature_in >> pressure_in)) continue;


    pressure.push_back(pressure_in*1.e+6);  //pressure in dyn cm-2
    temperature.push_back(temperature_in);
  }

  file.close();



  unsigned int nb_grid_points = pressure.size();



  std::vector<double> total_density(pressure);

  for (unsigned int i=0; i<nb_grid_points; i++)
    total_density[i] /= fastchem::CONST_K * temperature[i];

  
  //std::cout << "\n" << "Read in p-T structure: \n";

  //for (unsigned int i=0; i<nb_grid_points; i++)
  //  std::cout << i << "\t" << pressure[i] << "\t" << temperature[i] << std::endl;

  //std::cout << "\n";


  std::vector< std::vector<double> > densities;
  densities.resize(nb_grid_points);

  std::vector<double> total_element_density(nb_grid_points, 0.0);
  std::vector<double> mean_molecular_weights(nb_grid_points, 0.0);


  std::vector< std::vector<unsigned int> > element_conserved;
  element_conserved.resize(nb_grid_points);

  std::vector<unsigned int> nb_chemistry_iterations (nb_grid_points, 0);
  std::vector<unsigned int> fastchem_flags (nb_grid_points, 0);


  //set terminal output of FastChem
  fastchem.setVerboseLevel(verbose_level);


  //call FastChem point-wise
  for (unsigned int i=0; i<nb_grid_points; i++)
  {
    fastchem::FastChemInput input;
    fastchem::FastChemOutput output;

    input.temperature = temperature[i];
    input.pressure = pressure[i];

    if (i==0) input.use_previous_solution = false;
    else input.use_previous_solution = false;

    fastchem_flags[i] = fastchem.calcDensities(input, output);


    if (fastchem_flags[i] == fastchem::FASTCHEM_INITIALIZATION_FAILED)
    {
       std::cout << "FastChem Initialization failed!\n";
       return 1.0;
    }

    
    densities[i] = output.number_densities;
    total_element_density[i] = output.total_element_density;
    mean_molecular_weights[i] = output.mean_molecular_weight;
    element_conserved[i] = output.element_conserved;
    nb_chemistry_iterations[i] = output.nb_chemistry_iterations;

  }



  file.open(chem_output_file.c_str(), std::ios::out);

  unsigned int nb_species = fastchem.getSpeciesNumber();


  file << std::setw(16) << std::left << "P" << "\t"
       << std::setw(16) << std::left << "T" << "\t"
       << std::setw(16) << std::left << "n_<tot>" << "\t"
       << std::setw(16) << std::left << "n_g" << "\t"
       << std::setw(16) << std::left << "m";
  for (unsigned int i=0; i<nb_species; i++)
    file << "\t" << std::setw(16) << std::left << fastchem.getSpeciesSymbol(i);

  file << "\n";


  for (unsigned int i=0; i<nb_grid_points; i++)
  {
    file << std::setprecision(10) << std::scientific
         << pressure[i]/1.e6 << "\t"
         << temperature[i] << "\t"
         << total_element_density[i] << "\t"
         << total_density[i] << "\t"
         << mean_molecular_weights[i];

    for (unsigned int j=0; j<nb_species; j++)
      if (!output_mixing_ratios) file << "\t" << densities[i][j];
      else file << "\t" << densities[i][j] /total_density[i];

    file << "\n";
  }

  file.close();



  file.open(monitor_output_file.c_str(), std::ios::out);

  unsigned int nb_elements = fastchem.getElementNumber();


  file << std::setw(10) << std::left << "#grid point" << "\t"
       << std::setw(10) << std::left << "c_iterations" << "\t"
       << std::setw(16) << std::left << "c_convergence" << "\t"
       << std::setw(16) << std::left << "P" << "\t"
       << std::setw(16) << std::left << "T(k)" << "\t"
       << std::setw(16) << std::left << "n_<tot> (cm-3)" << "\t"
       << std::setw(16) << std::left << "n_g (cm-3)" << "\t"
       << std::setw(16) << std::left << "m(u)";
  for (unsigned int i=0; i<nb_elements; i++)
    file << "\t" << std::setw(5) << std::left << fastchem.getElementSymbol(i);

  file << "\n";


  std::vector<std::string> output_flags {"fail", "ok"};
  std::vector<std::string> convergence_flags {"yes", "no conv", "init fail"};


  for (unsigned int i=0; i<nb_grid_points; i++)
  {
    std::string c_conv;

    if (fastchem_flags[i] == fastchem::FASTCHEM_NO_CONVERGENCE)
      c_conv = output_flags[0];
    else
      c_conv = output_flags[1];




    file << std::setw(10) << i << "\t"
         << std::setw(10) << nb_chemistry_iterations[i] << "\t"
         << std::setw(16) << c_conv << "\t";


    file << std::setprecision(10) << std::scientific
         << pressure[i]/1.e6 << "\t" << temperature[i] << "\t"
                                     << total_element_density[i] << "\t"
                                     << total_density[i] << "\t"
                                     << mean_molecular_weights[i];

    for (unsigned int j=0; j<nb_elements; j++)
      file << "\t" << std::setw(5) << output_flags[element_conserved[i][j]];

    file << "\n";
  }

  file.close();


  file.open("output/chem_species.dat", std::ios::out);


  for (unsigned int j=0; j<nb_species; j++)
    file << fastchem.getSpeciesSymbol(j) << "_i=" << j+1+5 << ";\n";

  file.close();

  std::cout << "Equilibrium calculation (Fastchem) finished. " << std::endl;

  return 0;
}
