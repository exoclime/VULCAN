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


#include "fastchem_constants.h"
#include "species_struct.h"

#include <cmath>



namespace fastchem {


//Calculates the mass action constant, see Eq. (2.9)
//Change this function if you want to implement your own parametrisation
template <class double_type>
void Molecule<double_type>::calcMassActionConstant(const double temperature)
{
  double_type thermal_energy = 1.0e-6 * CONST_K * temperature;
  double_type log_K;
  double_type log_C, log_H, log_N, log_O, log_Ti, log_V;

  if (temperature <= 1000.0) {
	log_K = mass_action_coeff[0]/2.* std::pow(temperature, -2.)
	- mass_action_coeff[1]*( (1+std::log(temperature))/temperature )
	- mass_action_coeff[2]*(1-std::log(temperature))
	+ 1./2*mass_action_coeff[3]*temperature
	+ 1./6*mass_action_coeff[4]*std::pow(temperature, 2.)
	+ 1./12*mass_action_coeff[5]*std::pow(temperature, 3.)
	+ 1./20*mass_action_coeff[6]*std::pow(temperature, 4.)
	- mass_action_coeff[8]/temperature + mass_action_coeff[9];

   log_C = 6.495031470E+02/2.* std::pow(temperature, -2.)
- (-9.649010860E-01)*( (1+std::log(temperature))/temperature )
- 2.504675479E+00*(1-std::log(temperature))
+ 1./2*(-1.281448025E-05)*temperature
+ 1./6*1.980133654E-08*std::pow(temperature, 2.)
+ 1./12*(-1.606144025E-11)*std::pow(temperature, 3.)
+ 1./20*5.314483411E-15*std::pow(temperature, 4.)
- 8.545763110E+04/temperature + 4.747924288E+00;

   log_H = -2.5E+00*(1-std::log(temperature))
	- 2.547370801E+04/temperature -4.466828530E-01;

   log_N = -2.5E+00*(1-std::log(temperature))
	- 5.610463780E+04/temperature +4.193905036E+00;

   log_O = -7.953611300E+03/2.* std::pow(temperature, -2.)
- 1.607177787E+02*( (1+std::log(temperature))/temperature )
- 1.966226438E+00*(1-std::log(temperature))
+ 1./2*1.013670310E-03*temperature
- 1./6*1.110415423E-06*std::pow(temperature, 2.)
+ 1./12*6.517507500E-10*std::pow(temperature, 3.)
+ 1./20*(-1.584779251E-13)*std::pow(temperature, 4.)
- 2.840362437E+04/temperature + 8.404241820E+00;
	
   log_Ti =   -4.570179400E+04/2.* pow(temperature, -2.)
	-  6.608092020E+02*( (1+log(temperature))/temperature )
	- 4.295257490E-01*(1-log(temperature))
	+ 1./2*3.615029910E-03*temperature
	- 1./6*3.549792810E-06*pow(temperature, 2.)
	+ 1./12*1.759952494E-09*pow(temperature, 3.)
	- 1./20*3.052720871E-13*pow(temperature, 4.)
	- 5.270947930E+04/temperature + 2.026149738E+01;
   
   log_V = -5.535376020E+04/2.* pow(temperature, -2.)
	- 5.593338510E+02*( (1+log(temperature))/temperature )
	- 2.675543482E+00*(1-log(temperature))
	- 1./2*6.243049630E-03*temperature
	+ 1./6*1.565902337E-05*pow(temperature, 2.)
	- 1./12*1.372845314E-08*pow(temperature, 3.)
	+ 1./20*4.168388810E-12*pow(temperature, 4.)
	- 5.820664360E+04/temperature + 9.524567490E+00;
	} else {
		
   log_K =  mass_action_coeff[10]/2.* std::pow(temperature, -2.)
	- mass_action_coeff[11]*( (1+std::log(temperature))/temperature )
	- mass_action_coeff[12]*(1-std::log(temperature))
	+ 1./2*mass_action_coeff[13]*temperature
	+ 1./6*mass_action_coeff[14]*std::pow(temperature, 2.)
	+ 1./12*mass_action_coeff[15]*std::pow(temperature, 3.)
	+ 1./20*mass_action_coeff[16]*std::pow(temperature, 4.)
	- mass_action_coeff[18]/temperature + mass_action_coeff[19];

   log_C = -1.289136472E+05/2.* std::pow(temperature, -2.)
   - 1.719528572E+02*( (1+std::log(temperature))/temperature )
   - 2.646044387E+00*(1-std::log(temperature))
   + 1./2*(-3.353068950E-04)*temperature
   + 1./6*1.742092740E-07*std::pow(temperature, 2.)
   + 1./12*(-2.902817829E-11)*std::pow(temperature, 3.)
   + 1./20*1.642182385E-15*std::pow(temperature, 4.)
   - 8.410597850E+04/temperature + 4.130047418E+00;

   log_H = 6.078774250E+01/2.* std::pow(temperature, -2.)
   - (-1.819354417E-01)*( (1+std::log(temperature))/temperature )
   - 2.500211817E+00*(1-std::log(temperature))
   + 1./2*(-1.226512864E-07)*temperature
   + 1./6*(3.732876330E-11)*std::pow(temperature, 2.)
   + 1./12*(-5.687744560E-15)*std::pow(temperature, 3.)
   + 1./20*3.410210197E-19*std::pow(temperature, 4.)
   - 2.547486398E+04/temperature -4.481917770E-01;

   log_N = 8.876501380E+04/2.* std::pow(temperature, -2.)
  	- (-1.071231500E+02)*( (1+std::log(temperature))/temperature )
  	- (2.362188287E+00)*(1-std::log(temperature))
  	+ 1./2*2.916720081E-04*temperature
  	+ 1./6*(-1.729515100E-07)*std::pow(temperature, 2.)
  	+ 1./12*(4.012657880E-11)*std::pow(temperature, 3.)
  	+ 1./20*(-2.677227571E-15)*std::pow(temperature, 4.)
  	- (5.697351330E+04)/temperature + 4.865231506E+00;

  log_O = 2.619020262E+05/2.* std::pow(temperature, -2.)
  	- (-7.298722030E+02)*( (1+std::log(temperature))/temperature )
  	- 3.317177270E+00*(1-std::log(temperature))
  	+ 1./2*(-4.281334360E-04)*temperature
  	+ 1./6*1.036104594E-07*std::pow(temperature, 2.)
  	+ 1./12*(-9.438304330E-12)*std::pow(temperature, 3.)
  	+ 1./20*2.725038297E-16*std::pow(temperature, 4.)
  	- (3.392428060E+04)/temperature + -6.679585350E-01;
  
  log_Ti =  -1.704786714E+05/2.* pow(temperature, -2.)
	- 1.073852803E+03*( (1+log(temperature))/temperature )
	- 1.181955014E+00*(1-log(temperature))
	+ 1./2*2.245246352E-04*temperature
	+ 1./6*3.091697848E-07*pow(temperature, 2.)
	+ 1./12*-5.740027280E-11*pow(temperature, 3.)
	+ 1./20*2.927371014E-15*pow(temperature, 4.)
	- 4.978069910E+04/temperature + 1.740431368E+01;
  
  log_V = 1.200390300E+06/2.* pow(temperature, -2.)
  	+ 5.027005300E+03*( (1+log(temperature))/temperature )
  	- 1.058830594E+01*(1-log(temperature))
  	- 1./2*5.044326100E-03*temperature
  	+ 1./6*1.488547375E-06*pow(temperature, 2.)
  	- 1./12*1.785922508E-10*pow(temperature, 3.)
  	+ 1./20*8.113013866E-15*pow(temperature, 4.)
  	- 9.170740910E+04/temperature -4.768336320E+01;
}
   
	// in the order of the element_abundances.dat file
	unsigned int index_C = 0;
	unsigned int index_H = 1;
	unsigned int index_N = 3;
	unsigned int index_O = 4;
	unsigned int index_Ti = 8;
	unsigned int index_V = 9;

  log_K -= stoichometric_vector[index_C]*log_C + stoichometric_vector[index_H]*log_H + stoichometric_vector[index_N]*log_N + stoichometric_vector[index_O]*log_O + stoichometric_vector[index_Ti]*log_Ti+ stoichometric_vector[index_V]*log_V;
  mass_action_constant = log_K - sigma * std::log(thermal_energy);
  
}


template struct Molecule<double>;
template struct Molecule<long double>;

}
