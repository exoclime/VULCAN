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


#vul_data = 'output/photo_Moses_HD209_nominalKzz.vul'
#vul_data2 = 'output/res01-800nm-photo_Moses_HD209_nominalKzz.vul'

# vul_data = '../output/newBC-100nmcut-cap_condenR1e-6-Paul-ini-Earth.vul'
# vul_data2 = '../output/test-Earth.vul'
# vul_data3 = '../output/no-Coldtrap-Earth.vul'

vul_data = '../output/EQini-noPhoto-HD189.vul'
vul_data2 = '../output/fastchem-test-HD189.vul'

# Setting the 2rd input argument as the species names to be plotted (separated by ,)
plot_spec = sys.argv[1]
# Setting the 3th input argument as the output eps filename        
plot_name = sys.argv[2]

plot_dir = '../'+ vulcan_cfg.plot_dir

# taking user input species and splitting into separate strings and then converting the list to a tuple
plot_spec = tuple(plot_spec.split(','))
nspec = len(plot_spec)

colors = ['c','b','g','r','m','y','darkgray','orange','pink','grey','darkred','darkblue','salmon','chocolate','steelblue','plum','hotpink']

tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH'}


with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)
with open(vul_data2, 'rb') as handle:
  data2 = pickle.load(handle)
  
# with open(vul_data3, 'rb') as handle:
#   data3 = pickle.load(handle)

color_index = 0
vulcan_spec = data['variable']['species']
vulcan_spec2 = data2['variable']['species']

for sp in plot_spec:
    if color_index == len(colors): # when running out of colors
        colors.append(tuple(np.random.rand(3)))
    if sp in tex_labels: sp_lab = tex_labels[sp]
    else: sp_lab = sp
    
    #plt.plot(data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][:-1]/1.e5, color=colors[color_index], label=sp_lab, alpha=0.9)
    plt.plot(data['variable']['y_ini'][:,vulcan_spec.index(sp)]/data['atm']['n_0'], data['atm']['pco']/1.e6, color=colors[color_index], ls='-',lw=1.2, alpha=0.9)
    if sp in data2['variable']['species']:
        #plt.plot(data2['variable']['ymix'][:,vulcan_spec2.index(sp)], data2['atm']['zco'][:-1]/1.e5, color=colors[color_index], ls='--',lw=1.2, alpha=0.9)
        plt.plot(data2['variable']['y_ini'][:,vulcan_spec2.index(sp)]/data2['atm']['n_0'], data2['atm']['pco']/1.e6, color=colors[color_index], ls=':',lw=1.2, alpha=0.9)
    # if sp in data3['variable']['species']:
    #     plt.plot(data3['variable']['ymix'][:,vulcan_spec2.index(sp)], data3['atm']['zco'][:-1]/1.e5, color=colors[color_index], ls='-.',lw=1.2, alpha=0.9)

    color_index +=1
      
plt.gca().set_xscale('log')       
plt.gca().set_yscale('log') 
plt.gca().invert_yaxis() 
plt.xlim((1.E-16, 2.))
#plt.ylim((1.E3,1.E-8))
plt.legend(frameon=0, prop={'size':12}, loc=3)
plt.xlabel("Mixing Ratio")
#plt.ylabel("Pressure (bar)")
plt.ylabel("Height (km)")
plt.title('Earth')
plt.savefig(plot_dir + plot_name + '.png')
plt.savefig(plot_dir + plot_name + '.eps')
if vulcan_cfg.use_PIL == True:
    plot = Image.open(plot_dir + plot_name + '.png')
    plot.show()
else: plt.show()

