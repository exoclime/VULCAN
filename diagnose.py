# Diagnostic tool for plotting diffdf and k[i] at each level

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.legend as lg
import vulcan_cfg
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False
import os, sys
import pickle

import chem_funs, op
from chem_funs import nr, re_wM_dict

# Setting the 2nd input argument as the filename of vulcan output   
vul_data = 'output/st99-nz60-cap1e6e0-2nd-Earth.vul'

# setting the numerical solver to the desinated one in vulcan_cfg
solver_str = vulcan_cfg.ode_solver
solver = getattr(op, solver_str)()

# the number of fastest reactions to print out
top_num = 10

with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)
species = data['variable']['species'] 

rate_list = []
max_re_list = []  
for re in range(1,nr+1):
    rate = data['variable']['k'][re].astype(float)
    for sp in re_wM_dict[re][0]: # [0]: reactants; [1]: prodcuts
        if sp == 'M': rate *= data['atm']['n_0']
        else: rate *= data['variable']['ymix'][:,species.index(sp)] 
    rate_list.append(rate)
    max_re_list.append(np.amax(rate))

# 1+ is to shift the index to match starting with R1  
re_sort_indx = 1 + np.argsort(max_re_list)[::-1]
#rate_sort = np.sort(max_re_list)[::-1]
top_re = re_sort_indx[0:top_num]
for re in top_re:
    print (re)
    if re % 2 == 1: print (data['variable']['Rf'][re] + ' max rate: ' + "{:.2e}".format(max_re_list[re-1]))
    else: print('The reverse of ' + str(data['variable']['Rf'][re-1]))
     
#plt.figure()
#plt.plot(data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][:-1]/1.e5, color=colors[color_index], label=sp_lab, lw=1.5)