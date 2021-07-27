# Diagnostic tool for output largest rates

import sys
sys.path.insert(0, '../') # including the upper level of directory for the path of modules

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
from chem_funs import nr, re_wM_dict, re_dict


vul_data = 'output/.vul'

# setting the numerical solver to the desinated one in vulcan_cfg
solver_str = vulcan_cfg.ode_solver
solver = getattr(op, solver_str)()

# the number of fastest reactions to print out
top_num = 100

# the species involved
diag_sp = 'CH4'

# the region to look at
# the bottom layer 
n_bot = 25
# the top layer
n_top = -1

with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)
species = data['variable']['species'] 

rate_list = []
max_re_list = [] 
total_re_list = [] 
for re in range(1,nr+1):
    rate = data['variable']['k'][re][n_bot:n_top].astype(float)
    for sp in re_wM_dict[re][0]: # [0]: reactants; [1]: prodcuts
        if sp == 'M': rate *= data['atm']['n_0'][n_bot:n_top]
        #else: rate *= data['variable']['y'][:,species.index(sp)] 
        elif sp != 'S8_l_s': rate *= data['variable']['y'][:,species.index(sp)][n_bot:n_top] 
        rate_list.append(rate)
    rate_list.append(rate)
    max_re_list.append(np.amax(rate))
    total_re_list.append(rate)

# 1+ is to shift the index to match starting with R1  
re_sort_indx = 1 + np.argsort(max_re_list)[::-1]
#rate_sort = np.sort(max_re_list)[::-1]
top_re = re_sort_indx[0:top_num]
for re in top_re:
    if diag_sp in re_dict[re][0] or diag_sp in re_dict[re][1]:
        print (re)
        if re % 2 == 1: 
            print (data['variable']['Rf'][re] + ' max rate: ' + "{:.2e}".format(max_re_list[re-1]))
            print (total_re_list[re-1])
        else: print('The reverse of ' + str(data['variable']['Rf'][re-1]))
     