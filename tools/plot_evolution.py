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


# Setting the vulcan output file to read in    
vul_data = '../output/rtol005-continue-1e7conden-evap-cap-Earth.vul'
#vul_data2 = '../output/COini-noPhoto-HD189.vul'

# Setting the list of species to plot
plot_spec = ['CH4', 'CO', 'H2O', 'H']
plot_spec = ['CH4','H2O', 'OH', 'O2', 'O3']
# Setting the output plot filename        
plot_name = 'evolution'

# Setting the pressure level (cgs) to plot 
plot_p = 1e5


plot_dir = '../' + vulcan_cfg.plot_dir
colors = ['c','b','g','r','m','y','k','orange','pink','grey','darkred','darkblue','salmon','chocolate','steelblue','plum','hotpink']

tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$','N2O':'N$_2$O','O3':'O$_3$' }

with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)
species = data['variable']['species']

#plot_spec = species[20:30]

# Find the index of pco closest to p_ana
p_indx1 = min( range(len(data['atm']['pco'])), key=lambda i: abs(data['atm']['pco'][i]-plot_p))

try: 
    with open(vul_data2, 'rb') as handle:
        data2 = pickle.load(handle)
        species2 = data2['variable']['species']
        p_indx2 = min( range(len(data2['atm']['pco'])), key=lambda i: abs(data2['atm']['pco'][i]-plot_p))
except NameError: pass

color_index = 0
for sp in plot_spec:
    if sp in tex_labels: sp_lab = tex_labels[sp]
    else: sp_lab = sp
    plt.plot(data['variable']['t_time'], np.array(data['variable']['y_time'])[:,p_indx1,species.index(sp)]/float(data['atm']['n_0'][p_indx1]),  color=colors[color_index], label=sp_lab) 
    
    try: plt.plot(data2['variable']['t_time'], np.array(data2['variable']['y_time'])[:,p_indx2,species2.index(sp)]/float(data2['atm']['n_0'][p_indx2]),  color=colors[color_index], ls='--') 
    except NameError: pass

    color_index += 1

plt.gca().set_xscale('log')       
plt.gca().set_yscale('log') 
#plt.gca().invert_yaxis() 
#plt.xlim((1.E-12, 1.))
plt.ylim((1e-22, 2.))
plt.legend(frameon=0, prop={'size':13}, loc='best')

plt.xlabel("Time(s)", fontsize=12)
#plt.ylabel("Pressure (bar)")
plt.ylabel("Mixing Ratio", fontsize=12)
#plt.title('Earth (CIRA equator in January 1986)', fontsize=14)
plt.savefig(plot_dir + plot_name + '.png')
plt.savefig(plot_dir + plot_name + '.eps')
if vulcan_cfg.use_PIL == True:
    plot = Image.open(plot_dir + plot_name + '.png')
    plot.show()
else: plt.show()
