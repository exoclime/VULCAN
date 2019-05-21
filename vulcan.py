#!/usr/bin/env python2

# ==============================================================================
# This is the main file of VULCAN: the chemical kinetics code.                      
# Copyright (C) 2016 Shang-Min Tsai (Shami)                                      
#
# To run VULCAN simply run this file with python.                              
# 
# Requirements:
# - Python and the following packages:
#   numpy 
#   scipy
#   Sympy
#   matplotlib 
#   PIL/Pillow (optional)
#
# - Following files in the same directory:
#   build_atm.py
#   chem_funs.py
#   CHO_network.txt (or other text file for the chemical network)
#   op.py
#   prepipe.py
#   store.py
#   vulcan.py
#   vulcan_cfg.py
#   
# - Subdirectories:
#   /atm/ - for the atmospheric input files
#   /output/ - the default directory for the binary output files  
#   /plot/ - for plots
#   /thermo/ - for thermodynamic data
#
# ============================================================================== 
#     VULCAN is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     VULCAN is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You find a copy of the GNU General Public License in the main VULCAN
#     directory under <GPL_license.txt>. If not, see <http://www.gnu.org/licenses/>.
# ==============================================================================
# ============ \V/ ========== Live long and prosper! ========== \V/ ============


# import public modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.legend as lg
import scipy
import scipy.optimize as sop
import time, timeit, os, sys

# import VULCAN modules
import store, build_atm, op
try: import chem_funs
except: 
    raise IOError ('\nThe module "chem_funs" does not exist.\nPlease run prepipe.py first to create the module...')
     
# import the configuration inputs
import vulcan_cfg
from phy_const import kb, Navo


# Setting the current working directory to the script location
abspath = os.path.abspath(sys.argv[0])
dname = os.path.dirname(abspath)
os.chdir(dname)

from chem_funs import ni, nr  # number of species and reactions in the network
np.set_printoptions(threshold='nan')  # print all for debuging


### sholud these be here?
species = chem_funs.spec_list
compo = np.genfromtxt(vulcan_cfg.com_file,names=True,dtype=None)
compo_row = list(compo['species'])

### creat the instances for storing the variables and parameters
data_var = store.Variables()
data_atm = store.AtmData()
data_para = store.Parameters()

# record starting CPU time
data_para.start_time = time.time()
make_atm = build_atm.Atm()

# for plotting and printing
output = op.Output()

# construct pico
data_atm = make_atm.f_pico(data_atm)
# construct Tco and Kzz 
data_atm =  make_atm.load_TPK(data_atm, output)
# for reading rates
rate = op.ReadRate()

# read-in network and calculating forward rates
data_var = rate.read_rate(data_var, data_atm)
# reversing rates
data_var = rate.rev_rate(data_var, data_atm)
# removing rates
data_var = rate.remove_rate(data_var)

ini_abun = build_atm.InitialAbun()
# initialing y and ymix (the number density and the mixing ratio of every species)
data_var = ini_abun.ini_y(data_var, data_atm)

# storing the initial total number of atmos
data_var = ini_abun.ele_sum(data_var)

# calculating mean molecular weight, dz, and dzi
data_atm = make_atm.f_mu_dz(data_var, data_atm)

# ============== Execute VULCAN  ==============
# time-steping in the while loop until conv() returns True or count > count_max 

# setting the solver to the desinated one in vulcan_cfg
solver = eval("op." + vulcan_cfg.ode_solver + "()") 
quasi = op.QuasiSteady(solver, output)

# Running the integration loop
integ = op.Integration(solver, output) 
integ(data_var, data_atm, data_para, make_atm)
quasi(data_var, data_atm, data_para, make_atm)

output.save_out(data_var, data_atm, data_para, dname)
