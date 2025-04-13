# Generating Rayleigh scattering cross cestions
import numpy as np

# 1./(lmd/1.e7) converts the wavength (nm) to the wavenumber, nu (cm^-1) in Daniel's note 

def cross(lmd, n_ref, nr, K):
    '''
    The cross sections times n_ref ^2 from Daniel's note
    lmd: wavelength (nm in the code; cm in the note); n_ref: a reference number density (not the actual # density); 
    nr(lmd): the refractive index; K(lmd): the King factor
    '''
    cross_ns = 24*np.pi**3/(n_ref**2 *(lmd/1.e7)**4.) * ( (nr**2-1.)/(nr**2+2.) )**2 *K 

    return cross_ns
     

    
lmd_array = np.arange(0.1,800.1, 0.1) 
n_indx, King, n_ref = {}, {}, {}
n_indx['H2'] = lambda lmd: 1.358E-4 * (1. + 7.52E-3* (lmd/1.e3)**(-2) ) + 1.
King['H2'] = lambda lmd: 1.
n_ref['H2'] = 2.65163E19 

n_indx['O2'] = lambda lmd: (20564.8 + 2.480899E13/(4.09E9 - (1./(lmd/1.e7))**2.) )*1.e-8 +1. 
n_ref['O2'] = 2.68678E19 
King['O2'] = lambda lmd: 1.09 + 1.385E-11* (1./(lmd/1.e7))**2. + 1.448E-20*(1./(lmd/1.e7))**4.

n_indx['N2'] = lambda lmd: (5677.465 + 3.1881874E14/(1.44E10 - (1./(lmd/1.e7))**2. )  )*1.e-8 +1.
n_ref['N2'] = 2.546899E19 
King['N2'] = lambda lmd: 1.034 + 3.17E-12*(1./(lmd/1.e7)) 

n_indx['CO2'] = lambda lmd: ( 5799.25/(128908.9**2-(1./(lmd/1.e7))**2) + 120.05/(89223.8**2-(1./(lmd/1.e7))**2) + 5.3334/(75037.5**2-(1./(lmd/1.e7))**2) +4.3244/(67837.7**2-(1./(lmd/1.e7))**2) + 0.1218145E-6/(2418.136**2-(1./(lmd/1.e7))**2)  )*1.1427e3 +1. 
n_ref['CO2'] = 2.546899E19 
King['CO2'] = lambda lmd: 1.1364 + 25.3E-12 * (1./(lmd/1.e7))**2.

n_indx['He'] = lambda lmd: ( 2283. + (1.8102E13)/(1.5342E10 - (1./(lmd/1.e7))**2. ) )*1.E-8 + 1
n_ref['He'] = 2.546899E19 
King['He'] = lambda lmd: 1.


# Choose the species
sp = 'He'

ost = '#lambda (nm)  cross(cm^2)\n'

for lm in lmd_array:
    ost += "{0:<14s}".format(str(lm)) + "{:.3E}".format(cross(lm, n_ref[sp], n_indx[sp](lm) , King[sp](lm))) +'\n'

ost = ost[:-1]
with open(sp + '_rayleigh.txt', 'w') as f: f.write(ost)