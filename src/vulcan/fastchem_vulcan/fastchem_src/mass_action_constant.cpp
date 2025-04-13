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
//#include <iostream>


namespace fastchem {


//Calculates the mass action constant, see Eq. (2.9)
//Change this function if you want to implement your own parametrisation
template <class double_type>
void Molecule<double_type>::calcMassActionConstant(const double temperature)
{
  double_type thermal_energy = 1.0e-6 * CONST_K * temperature;
  double_type log_K;
  double_type log_C, log_H, log_N, log_O, log_Ti, log_V, log_Cl, log_S, log_P, log_Si, log_Km, log_Na, log_Mg, log_F;
  double_type log_Ca, log_Fe;
  double_type log_He, log_e; 
  // Km: potassium since log_K has already been used
  // He: for reactions with He+, HeH+ etc.
  
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
 
   log_Cl = 2.276215854E+04/2.* std::pow(temperature, -2.)
	+ 2.168413293E+02*( (1+std::log(temperature))/temperature )
	- 2.745185115E+00*(1-std::log(temperature))
	+ 1./2*2.451101694E-03*temperature
	- 1./6*5.458011990E-06*std::pow(temperature, 2.)
	+ 1./12*4.417986880E-09*std::pow(temperature, 3.)
	- 1./20*1.288134004E-12*std::pow(temperature, 4.)
	- 1.501357068E+04/temperature + 3.102963457E+00;
	   
   log_Km = 9.665143930E+00/2.* std::pow(temperature, -2.)
	+ 1.458059455E-01*( (1+std::log(temperature))/temperature )
	- 2.500865861E+00*(1-std::log(temperature))
	- 1./2*2.601219276E-06*temperature
	+ 1./6*4.187306580E-09*std::pow(temperature, 2.)
	- 1./12*3.439722110E-12*std::pow(temperature, 3.)
	+ 1./20*1.131569009E-15*std::pow(temperature, 4.)
	- 9.959493490E+03/temperature + 5.035822260E+00;

   log_Na = 0.000000000E+00/2.* std::pow(temperature, -2.)
	- 0.000000000E+00*( (1+std::log(temperature))/temperature )
	- 2.500000000E+00*(1-std::log(temperature))
	+ 1./2*0.000000000E+00*temperature
	+ 1./6*0.000000000E+00*std::pow(temperature, 2.)
	+ 1./12*0.000000000E+00*std::pow(temperature, 3.)
	+ 1./20*0.000000000E+00*std::pow(temperature, 4.)
	- 1.218382949E+04/temperature + 4.244028180E+00;

   log_Mg = 0.000000000E+00/2.* std::pow(temperature, -2.)
	- 0.000000000E+00*( (1+std::log(temperature))/temperature )
	- 2.500000000E+00*(1-std::log(temperature))
	+ 1./2*0.000000000E+00*temperature
	+ 1./6*0.000000000E+00*std::pow(temperature, 2.)
	+ 1./12*0.000000000E+00*std::pow(temperature, 3.)
	+ 1./20*0.000000000E+00*std::pow(temperature, 4.)
	- 1.694658761E+04/temperature + 3.634330140E+00;
	   
   log_S = -3.174841820E+02/2.* std::pow(temperature, -2.)
	+ 1.924704923E+02*( (1+std::log(temperature))/temperature )
	- 4.686825930E+00*(1-std::log(temperature))
	- 1./2*5.841365600E-03*temperature
	+ 1./6*7.538533520E-06*std::pow(temperature, 2.)
	- 1./12*4.863586040E-09*std::pow(temperature, 3.)
	+ 1./20*1.256976992E-12*std::pow(temperature, 4.)
	- 3.323592180E+04/temperature - 5.718523969E+00;
                 
   log_P = 5.040866570E+01/2.* std::pow(temperature, -2.)
	+ 7.639418650E-01*( (1+std::log(temperature))/temperature )
	- 2.504563992E+00*(1-std::log(temperature))
	- 1./2*1.381689958E-05*temperature
	+ 1./6*2.245585515E-08*std::pow(temperature, 2.)
	- 1./12*1.866399889E-11*std::pow(temperature, 3.)
	+ 1./20*6.227063395E-15*std::pow(temperature, 4.)
	- 3.732421910E+04/temperature + 5.359303481E+00;
   
   log_Si = 9.836140810E+01/2.* std::pow(temperature, -2.)
	- 1.546544523E+02*( (1+std::log(temperature))/temperature )
	- 1.876436670E+00*(1-std::log(temperature))
	+ 1./2*1.320637995E-03*temperature
	- 1./6*1.529720059E-06*std::pow(temperature, 2.)
	+ 1./12*8.950562770E-10*std::pow(temperature, 3.)
	- 1./20*1.952873490E-13*std::pow(temperature, 4.)
	- 5.263510310E+04/temperature + 9.698288880E+00;

   log_F = 1.137409088E+03/2.* std::pow(temperature, -2.)
	+ 1.453392797E+02*( (1+std::log(temperature))/temperature )
	- 4.077403610E+00*(1-std::log(temperature))
	- 1./2*4.303360140E-03*temperature
	+ 1./6*5.728897740E-06*std::pow(temperature, 2.)
	- 1./12*3.819312900E-09*std::pow(temperature, 3.)
	+ 1./20*1.018322509E-12*std::pow(temperature, 4.)
	- 9.311110120E+03/temperature -3.558982650E+00;
	
   log_Ca = 0.000000000E+00/2.* std::pow(temperature, -2.)
	- 0.000000000E+00*( (1+std::log(temperature))/temperature )
	- 2.500000000E+00*(1-std::log(temperature))
	+ 1./2*0.000000000E+00*temperature
	+ 1./6*0.000000000E+00*std::pow(temperature, 2.)
	+ 1./12*0.000000000E+00*std::pow(temperature, 3.)
	+ 1./20*0.000000000E+00*std::pow(temperature, 4.)
	- 2.063892786E+04/temperature + 4.384548330E+00;
   
  	log_Fe = 6.790822660E+04/2.* std::pow(temperature, -2.)
  	+ 1.197218407E+03*( (1+std::log(temperature))/temperature )
  	- 9.843393310E+00*(1-std::log(temperature))
  	- 1./2*1.652324828E-02*temperature
  	+ 1./6*1.917939959E-05*std::pow(temperature, 2.)
  	- 1./12*1.149825371E-08*std::pow(temperature, 3.)
  	+ 1./20*2.832773807E-12*std::pow(temperature, 4.)
  	- 5.466995940E+04/temperature -3.383946260E+01;
	
    log_He = 0* std::pow(temperature, -2.)
  	- 0*( (1+std::log(temperature))/temperature )
  	- 2.5*(1-std::log(temperature))
  	+ 1./2*0*temperature
  	+ 1./6*0*std::pow(temperature, 2.)
  	+ 1./12*0*std::pow(temperature, 3.)
  	+ 1./20*0*std::pow(temperature, 4.)
  	+ 7.453750000E+02/temperature + 9.287239740E-01;
			   
	log_e = 0/2.* std::pow(temperature, -2.)
	- 0*( (1+std::log(temperature))/temperature )
	- 2.500000000E+00*(1-std::log(temperature))
	+ 1./2*0*temperature
	+ 1./6*0*std::pow(temperature, 2.)
	+ 1./12*0*std::pow(temperature, 3.)
	+ 1./20*0*std::pow(temperature, 4.)
	+ 7.453750000E+02/temperature -1.172081224E+01;		   
				   
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
   
  log_Cl = -1.697932930E+05/2.* std::pow(temperature, -2.)
	- 6.081726460E+02*( (1+std::log(temperature))/temperature )
	- 2.128664090E+00*(1-std::log(temperature))
	+ 1./2*1.307367034E-04*temperature
	- 1./6*2.644883596E-08*std::pow(temperature, 2.)
	+ 1./12*2.842504775E-12*std::pow(temperature, 3.)
	- 1./20*1.252911731E-16*std::pow(temperature, 4.)
	- 9.934387400E+03/temperature + 8.844772103E+00;

  log_Km =  -3.566422360E+06/2.* std::pow(temperature, -2.)
	- 1.085289825E+04*( (1+std::log(temperature))/temperature )
	+ 1.054134898E+01*(1-std::log(temperature))
	+ 1./2*8.009801350E-03*temperature
	- 1./6*2.696681041E-06*std::pow(temperature, 2.)
	+ 1./12*4.715294150E-10*std::pow(temperature, 3.)
	- 1./20*2.976897350E-14*std::pow(temperature, 4.)
	+ 5.875337010E+04/temperature + 9.738551240E+01;
  
  log_Na =  9.525723380E+05/2.* std::pow(temperature, -2.)
	+ 2.623807254E+03*( (1+std::log(temperature))/temperature )
	- 5.162596620E+00*(1-std::log(temperature))
	- 1./2*1.210218586E-03*temperature
	+ 1./6*2.306301844E-07*std::pow(temperature, 2.)
	- 1./12*1.249597843E-11*std::pow(temperature, 3.)
	+ 1./20*7.226771190E-16*std::pow(temperature, 4.)
	- 2.912963564E+04/temperature -1.519717061E+01;
  
  log_Mg =  -5.364831550E+05/2.* std::pow(temperature, -2.)
	- 1.973709576E+03*( (1+std::log(temperature))/temperature )
	+ 3.633776900E-01*(1-std::log(temperature))
	+ 1./2*2.071795561E-03*temperature
	- 1./6*7.738051720E-07*std::pow(temperature, 2.)
	+ 1./12*1.359277788E-10*std::pow(temperature, 3.)
	- 1./20*7.766898397E-15*std::pow(temperature, 4.)
	- 4.829188110E+03/temperature + 2.339104998E+01;

  log_S =  -4.854244790E+05/2.* std::pow(temperature, -2.)
	- 1.438830408E+03*( (1+std::log(temperature))/temperature )
	- 1.258504116E+00*(1-std::log(temperature))
	+ 1./2*3.797990430E-04*temperature
	+ 1./6*1.630685864E-09*std::pow(temperature, 2.)
	- 1./12*9.547095850E-12*std::pow(temperature, 3.)
	+ 1./20*8.041466646E-16*std::pow(temperature, 4.)
	- 2.334995270E+04/temperature + 1.559554855E+01;

  log_P =  1.261794642E+06/2.* std::pow(temperature, -2.)
	+ 4.559838190E+03*( (1+std::log(temperature))/temperature )
	- 8.918079310E+00*(1-std::log(temperature))
	- 1./2*4.381401460E-03*temperature
	+ 1./6*1.454286224E-06*std::pow(temperature, 2.)
	- 1./12*2.030782763E-10*std::pow(temperature, 3.)
	+ 1./20*1.021022887E-14*std::pow(temperature, 4.)
	- 6.541723960E+04/temperature -3.915974795E+01;

  log_Si =  -6.169298850E+05/2.* std::pow(temperature, -2.)
	- 2.240683927E+03*( (1+std::log(temperature))/temperature )
	+ 4.448619320E-01*(1-std::log(temperature))
	+ 1./2*1.710056321E-03*temperature
	- 1./6*4.107714160E-07*std::pow(temperature, 2.)
	+ 1./12*4.558884780E-11*std::pow(temperature, 3.)
	+ 1./20*-1.889515353E-15*std::pow(temperature, 4.)
	- 3.953558760E+04/temperature + 2.679668061E+01;

  log_F =  1.473506226E+04/2.* std::pow(temperature, -2.)
	- 8.149927360E+01*( (1+std::log(temperature))/temperature )
	- 2.444371819E+00*(1-std::log(temperature))
	+ 1./2*2.120210026E-05*temperature
	- 1./6*4.546918620E-09*std::pow(temperature, 2.)
	+ 1./12*5.109528730E-13*std::pow(temperature, 3.)
	- 1./20*2.333894647E-17*std::pow(temperature, 4.)
	- 8.388374650E+03/temperature + 5.478710640E+00;
  
  log_Ca = 7.547341240E+06/2.* std::pow(temperature, -2.)
  + 2.148642662E+04*( (1+std::log(temperature))/temperature )
  - 2.530849567E+01*(1-std::log(temperature))
  - 1./2*1.103773705E-02*temperature
  + 1./6*2.293249636E-06*std::pow(temperature, 2.)
  - 1./12*1.209075383E-10*std::pow(temperature, 3.)
  - 1./20*4.015333268E-15*std::pow(temperature, 4.)
  - 1.585862323E+05/temperature -1.609512955E+02;          
  
  log_Fe =  -1.954923682E+06/2.* std::pow(temperature, -2.)
	- 6.737161100E+03*( (1+std::log(temperature))/temperature )
	+ 5.486410970E+00*(1-std::log(temperature))
	+ 1./2*4.378803450E-03*temperature
	- 1./6*1.116286672E-06*std::pow(temperature, 2.)
	+ 1./12*1.544348856E-10*std::pow(temperature, 3.)
	- 1./20*8.023578182E-15*std::pow(temperature, 4.)
	- 7.137370060E+03/temperature + 6.504979860E+01;
  
  log_He = 0* std::pow(temperature, -2.)
	- 0*( (1+std::log(temperature))/temperature )
	- 2.500000000E+00*(1-std::log(temperature))
	+ 0.000000000E+00*temperature
	+ 0.000000000E+00*std::pow(temperature, 2.)
	+ 1./12*0*std::pow(temperature, 3.)
	+ 1./20*0*std::pow(temperature, 4.)
	+ 7.453750000E+02/temperature + 9.287239740E-01;  
  
  log_e = 0/2.* std::pow(temperature, -2.)
	- 0*( (1+std::log(temperature))/temperature )
	- 2.500000000E+00*(1-std::log(temperature))
	+ 1./2*0*temperature
	+ 1./6*0*std::pow(temperature, 2.)
	+ 1./12*0*std::pow(temperature, 3.)
	+ 1./20*0*std::pow(temperature, 4.)
	+ 7.453750000E+02/temperature -1.172081224E+01;	 
 
   		 
}
   
	// in the order of the element_abundances.dat file
	unsigned int index_C = 0;
	unsigned int index_H = 1;
	unsigned int index_N = 3;
	unsigned int index_O = 4;
	unsigned int index_Ti = 8;
	unsigned int index_V = 9;
	
	unsigned int index_P = 5;
	unsigned int index_S = 6;
	unsigned int index_Si = 7;
	unsigned int index_Cl = 10;
	unsigned int index_Km = 11; 
	unsigned int index_Na = 12;
	unsigned int index_Mg = 13;
	unsigned int index_F = 14;
	unsigned int index_Ca = 15;
	unsigned int index_Fe = 16;
	
	unsigned int index_He = 2;
	unsigned int index_e = 17;
		
  log_K -= stoichometric_vector[index_C]*log_C + stoichometric_vector[index_H]*log_H + stoichometric_vector[index_N]*log_N + stoichometric_vector[index_O]*log_O + stoichometric_vector[index_Ti]*log_Ti+ stoichometric_vector[index_V]*log_V;
  log_K -= stoichometric_vector[index_P]*log_P + stoichometric_vector[index_S]*log_S + stoichometric_vector[index_Si]*log_Si + stoichometric_vector[index_Cl]*log_Cl + stoichometric_vector[index_Km]*log_Km + stoichometric_vector[index_Na]*log_Na + stoichometric_vector[index_Mg]*log_Mg + stoichometric_vector[index_F]*log_F + stoichometric_vector[index_Ca]*log_Ca + stoichometric_vector[index_Fe]*log_Fe;
  log_K -= stoichometric_vector[index_He]*log_He;
  log_K -= stoichometric_vector[index_e]*log_e;
  mass_action_constant = log_K - sigma * std::log(thermal_energy);
  
}


template struct Molecule<double>;
template struct Molecule<long double>;

}
