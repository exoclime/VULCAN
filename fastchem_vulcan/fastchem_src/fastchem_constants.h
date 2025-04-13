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


#ifndef _fastchem_constants_h
#define _fastchem_constants_h


namespace fastchem {


//FastChem constants
const unsigned int FASTCHEM_UNKNOWN_SPECIES = 9999999;

const unsigned int FASTCHEM_SUCCESS = 0;
const unsigned int FASTCHEM_NO_CONVERGENCE = 1;
const unsigned int FASTCHEM_INITIALIZATION_FAILED = 2;
const unsigned int FASTCHEM_IS_BUSY = 3;
const unsigned int FASTCHEM_WRONG_INPUT_VALUES = 4;


//Physical constants
const double CONST_K = 1.3806504e-16;    //Boltzmann's constant in erg K-1
const double CONST_AMU = 1.66055e-24;    //Atomic mass unit in g


}

#endif
