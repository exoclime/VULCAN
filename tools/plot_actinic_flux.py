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
       
plot_name = 'Earth-photosphere_sp_sum'
vul_data = '../output/v321-100nmcut-cap_condenR1e-6-Paul-ini-Earth.vul'
plot_dir = '../' + vulcan_cfg.plot_dir

color_index = 0
colors = ['c','b','g','r','m','y','chocolate','orange','pink','grey','darkred','salmon','steelblue','hotpink','k']

tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$','O3':'O$_3$' }


with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)

# hight to plot
height_list = [10, 20, 30, 40, 50, 60]

# photosphere of each species
tau_sp= {}

photo_sp = ['O2', 'O3', 'H2O'] # , 'H2', 'CH4', 'CO2'   
#, 'N2', 'C2H4', 'C2H2','C2H6', 'CH3', 'CO2', 'HCO','HCN'  ,'CH3CHO','NO' , 'NO2'
#scat_sp = ['H2','He']
scat_sp = ['N2', 'O2']

bins = data['variable']['bins']
dz = data['atm']['dz']
nz = len(dz)
vulcan_spec = data['variable']['species']


for zz in height_list:
    # Find the index of zco
    z_indx = min( range(len(data['atm']['zco'])), key=lambda i: abs(data['atm']['zco'][i]-zz*1e5))    
    plt.plot(bins, data['variable']['aflux'][z_indx], c=colors[color_index], lw=1.2, label=str(int(data['atm']['zco'][z_indx]/1e5))  )
    color_index += 1

#plt.plot(bins, photosph_sum, c='plum', lw=3, alpha=0.6)


#plt.plot(data['variable']['bins'], photosph_abs, c='red')
           
plt.gca().set_yscale('log') 
#plt.gca().invert_yaxis() 
plt.ylim((1.E9,1.E15))
plt.xlim((150,400))
plt.legend(frameon=0, prop={'size':12}, loc='best')
plt.xlabel("wavelength (nm)")
plt.ylabel("actinic flux ()")
#plt.title('Photosphere')
plt.savefig(plot_dir + plot_name + '.png')
plt.savefig(plot_dir + plot_name + '.pdf')
if vulcan_cfg.use_PIL == True:
    plot = Image.open(plot_dir + plot_name + '.png')
    plot.show()
else: plt.show()

