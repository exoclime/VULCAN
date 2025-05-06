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
import matplotlib.ticker as mtick

data_list = ['../output/helios_output/output/NEQ_HD189/0_TOA_flux_eclipse.dat', '../output/helios_output/output/NEQ_HD189/0_TOA_flux_eclipse.dat']
label_list = ['EQ', 'NEQ']

plot_name = sys.argv[1]

plot_dir = '../' + vulcan_cfg.plot_dir

# taking user input species and splitting into separate strings and then converting the list to a tuple
#plot_spec = tuple(plot_spec.split(','))
#nspec = len(plot_spec)

colors = ['r', 'orange', 'g','b','c','m','pink','grey','darkred','darkblue','salmon','chocolate','steelblue','plum','hotpink']

rJ = 7.1492E4 #km
r_sun = 6.957E5 #km
emi = {}

#fig, ax = plt.subplots()

for index, data in enumerate(data_list):
    emi[index] = np.genfromtxt(data,dtype=float,skip_header=3)
    plt.plot(emi[index][:,1], emi[index][:,6], color=colors[index], label=label_list[index], lw=1, alpha=0.5)
        
plt.gca().set_xscale('log') 
plt.xlim((emi[0][:,1][0],20.))
#plt.xlim((0.0001,10.))
#plt.ylim((0,0.002))
#plt.xticks([1,5,10])
#plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%i'))
#plt.gca().xaxis.set_minor_formatter(mtick.NullFormatter())
#plt.minorticks_off()
#plt.ylim((data['atm']['pco'][0][0]/1e6,data['atm']['pco'][0][-1]/1e6))
plt.legend(frameon=0, prop={'size':13}, loc=4)
plt.xlabel(r"Wavelength ($\mu$m)") # "\n" new line
#plt.ylabel("% R$_p$/R$_s$$^2$")
#plt.title(tex_labels[sp])
plt.savefig(plot_dir + plot_name + '.png')
plt.savefig(plot_dir + plot_name + '.pdf')
if vulcan_cfg.use_PIL == True:
    plot = Image.open(plot_dir + plot_name + '.png')
    plot.show()
else: plt.show()

