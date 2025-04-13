'''
This script reads VULCAN output (.vul) files using pickle and plot the species volumn mixing ratios as a function of pressure, with the initial abundances (typically equilibrium) shown in dashed lines.
Plots are saved in the folder assigned in vulcan_cfg.py, with the default plot_dir = 'plot/'.
'''

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

# swtich for plot
if '-h' in sys.argv: use_height = True
else: use_height = False 

# Setting the 2nd input argument as the filename of vulcan output   
vul_data = sys.argv[1]
# Setting the 3rd input argument as the species names to be plotted (separated by ,)
plot_spec = sys.argv[2]
# Setting the 4th input argument as the output eps filename        
plot_name = sys.argv[3]

plot_dir = '../' + vulcan_cfg.plot_dir
# Checking if the plot folder exsists
if not os.path.exists(plot_dir):
    print ('The plotting directory assigned in vulcan_cfg.py does not exist.')
    print( 'Directory ' , plot_dir,  " created.")
    os.mkdir(plot_dir)

# taking user input species and splitting into separate strings and then converting the list to a tuple
plot_spec = tuple(plot_spec.split(','))
nspec = len(plot_spec)

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180),(255, 127, 14),(44, 160, 44),(214, 39, 40),(148, 103, 189),(140, 86, 75), (227, 119, 194),(127, 127, 127),(188, 189, 34),(23, 190, 207),\
(174, 199, 232),(255, 187, 120),(152, 223, 138),(255, 152, 150),(197, 176, 213),(196, 156, 148),(247, 182, 210),(199, 199, 199),(219, 219, 141),(158, 218, 229)] 
# 


# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

# tex labels for plotting
tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$' }



with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)

color_index = 0
vulcan_spec = data['variable']['species']
for color_index,sp in enumerate(plot_spec):
    if color_index == len(tableau20): # when running out of colors
        tableau20.append(tuple(np.random.rand(3)))
    
    if sp in tex_labels: sp_lab = tex_labels[sp]
    else: sp_lab = sp  
    
    #plt.plot(data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][:-1]/1.e5, color=tableau20[color_index], label=sp_lab, lw=1.5)
    if use_height == False:
        plt.plot(data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['pco']/1.e6, color=tableau20[color_index], label=sp_lab, lw=1.5)
    else: 
        plt.plot(data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][1:]/1.e5, color=tableau20[color_index], label=sp_lab, lw=1.5)
    #plt.plot(data['variable']['y_ini'][:,vulcan_spec.index(sp)]/data['atm']['n_0'], data['atm']['pco']/1.e6, color=tableau20[color_index], ls=':', lw=1.5) # plotting the initial (equilibrium) abundances


if use_height == False:
    plt.gca().set_yscale('log') 
    plt.gca().invert_yaxis() 
    plt.ylim((data['atm']['pco'][0]/1e6,data['atm']['pco'][-1]/1e6))
    plt.ylabel("Pressure (bar)")
else:
    plt.ylim((data['atm']['zmco'][0]/1e5,data['atm']['zmco'][-1]/1e5)) 
    plt.xlabel("Mixing Ratio")  
    
#plt.title('T1400')
   
plt.gca().set_xscale('log')       
plt.xlim((1.E-12, 1.e-2))
plt.legend(frameon=0, prop={'size':12}, loc='best')
# handles, labels = plt.gca().get_legend_handles_labels()
# display = range(len(sp_list))
# #Create custom artists
# art0 = plt.Line2D((0,0),(0,0), ls='None')
# Artist1 = plt.Line2D(range(10),range(10), color='black')
# Artist2 = plt.Line2D((0,1),(0,0), color='black', ls='--',lw=1.5)
# plt.legend([Artist1,Artist2],['Equilibrium','Kinetics'], frameon=False, prop={'size':12}, loc='best')

plt.savefig(plot_dir + plot_name + '.png')
plt.savefig(plot_dir + plot_name + '.eps')
if vulcan_cfg.use_PIL == True:
    plot = Image.open(plot_dir + plot_name + '.png')
    plot.show()
else: plt.show()

