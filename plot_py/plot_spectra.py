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

#data_list = ['transit_spectra_Kzz/spectrum_EQ_lr.dat','transit_spectra_Kzz/spectrum_vxExp1E5-KzzE8_lr.dat','transit_spectra_Kzz/spectrum_vxExp1E5-KzzE10_lr.dat','transit_spectra_Kzz/spectrum_vxExp1E5-KzzE12_lr.dat']
data_list = ['../output/transit_spectra/vulcan_Moses_Venot/spectrum_EQ_HD189_lr.dat', '../output/transit_spectra/vulcan_Moses_Venot/spectrum_HD189-vulcan_lr.dat','../output/transit_spectra/vulcan_Moses_Venot/spectrum_HD189-Moses_lr.dat',\
'../output/transit_spectra/vulcan_Moses_Venot/spectrum_HD189-Venot_lr.dat']

# Wasp43-b
#data_list = ['transit_spectra/spectrum_4X-wasp43b-vx0-wKzz_lr.dat', 'transit_spectra/spectrum_4X-wasp43b-vx-wKzz_lr.dat']

#label_list = [r'K$_{zz}$ = 0', '         10$^8$', '         10$^{10}$', '         10$^{12}$']
label_list = ['U = 0', '      10 (m/s)', '      100 (m/s)', '      1000 (m/s)', '      5000 (m/s)']
label_list = ['EQ', 'wo H2O', 'wo CH4']
label_list = ['EQ', 'Vulcan', 'Moses', 'Venot']

plot_name = sys.argv[1]

plot_dir = '../' + vulcan_cfg.plot_dir

# taking user input species and splitting into separate strings and then converting the list to a tuple
#plot_spec = tuple(plot_spec.split(','))
#nspec = len(plot_spec)

colors = ['r', 'orange', 'g','b','c','m','pink','grey','darkred','darkblue','salmon','chocolate','steelblue','plum','hotpink']

rJ = 6.9911E4 #km
r_sun = 6.957E5 #km
tran = {}

#fig, ax = plt.subplots()

for index, data in enumerate(data_list):
    tran[index] = np.genfromtxt(data,dtype=float,skip_header=1, names = ['lambda','R'])
    
    plt.plot(tran[index]['lambda'], 100*(tran[index]['R']/r_sun)**2, color=colors[index], label=label_list[index], lw=1, alpha=0.5)
        
plt.gca().set_xscale('log') 
plt.xlim((1,6.01))
#plt.xticks([1,5,10])
plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%i'))
plt.gca().xaxis.set_minor_formatter(mtick.NullFormatter())
#plt.minorticks_off()
#plt.ylim((data['atm']['pco'][0][0]/1e6,data['atm']['pco'][0][-1]/1e6))
plt.legend(frameon=0, prop={'size':13}, loc=4)
plt.xlabel(r"Wavelength ($\mu$m)") # "\n" new line
plt.ylabel("% R$_p$/R$_s$$^2$")
#plt.title(tex_labels[sp])
plt.savefig(plot_dir + plot_name + '.png')
plt.savefig(plot_dir + plot_name + '.pdf')
if vulcan_cfg.use_PIL == True:
    plot = Image.open(plot_dir + plot_name + '.png')
    plot.show()
else: plt.show()

