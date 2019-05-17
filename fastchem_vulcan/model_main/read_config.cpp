/*
* This file is part of the FastChem code (https://github.com/exoclime/fastchem).
* Copyright (C) 2018 Daniel Kitzmann, Joachim Stock
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


#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <functional>

#include <cmath>


bool read_config_file(std::string file_path, std::string& fastchem_options_file, std::string& atmosphere_file,
                      std::string& chem_output_file, std::string& monitor_output_file,
                      unsigned int& verbose_level, bool& output_mixing_ratios)
{
  std::fstream file;

  file.open(file_path.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to read config file: " << file_path << "\n";
    return false;
  }


  std::string line;

  //FastChem parameter file
  std::getline(file, line);

  std::getline(file, line);
  fastchem_options_file = line;

  if (line == "")
  {
    std::cout << "Unable to read FastChem parameter file location from: " << file_path.c_str() << "\n";

    return false;
  }


  //p-T file
  std::getline(file, line);
  std::getline(file, line);

  std::getline(file, line);
  atmosphere_file = line;

  if (line == "")
  {
    std::cout << "Unable to read p-T file location from: " << file_path.c_str() << "\n";

    return false;
  }


  //chemistry output file
  std::getline(file, line);
  std::getline(file, line);

  std::getline(file, line);
  chem_output_file = line;

  if (line == "")
  {
    std::cout << "Unable to read chemistry output file location from: " << file_path.c_str() << "\n";

    return false;
  }


  //monitor file
  std::getline(file, line);
  std::getline(file, line);

  std::getline(file, line);
  monitor_output_file = line;

  if (line == "")
  {
    std::cout << "Unable to read monitor file location from: " << file_path.c_str() << "\n";

    return false;
  }


  //verbose option
  std::getline(file, line);
  std::getline(file, line);

  std::getline(file, line);
  if (line == "")
    std::cout << "Unable to read verbose option from: " << file_path.c_str() << " , using value of 4 from here on" << "\n";
  else
    verbose_level = std::stoi(line);


  //output option
  std::getline(file, line);
  std::getline(file, line);

  std::getline(file, line);
  if (line == "MR")
    output_mixing_ratios = true;
  else
    output_mixing_ratios = false;


  if (verbose_level == 4)
  std::cout << "Read from config file " << file_path.c_str() << " : \n"
            << "FastChem parameter file: " << fastchem_options_file << "\n"
            << "Temperature-pressure file: " << atmosphere_file << "\n"
            << "Chemistry output file: " << chem_output_file << "\n"
            << "Monitor file: " << monitor_output_file << "\n"
            << "Verbose level: " << verbose_level << "\n"
            << "Output mixing ratios: " << output_mixing_ratios << "\n"
            << "\n";


  return true;
}
