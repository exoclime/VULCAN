import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.legend as lg
#import vulcan_cfg
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False
import os, sys
import pickle
#from phy_const import au, r_sun
       
plot_name = 'stellar-flux'
# vul_data = 'output/HCO_new_H2O_HD189_Moses_nominalKzz.vul'
# vul_data2 = 'output/test.vul'

gj436 = np.genfromtxt('sflux-GJ436.txt', names=['lambda','flux'], dtype=None, skip_header=1)
gj876 = np.genfromtxt('sflux-GJ876.txt', names=['lambda','flux'], dtype=None, skip_header=1)
gj551 = np.genfromtxt('sflux-GJ551.txt', names=['lambda','flux'], dtype=None, skip_header=1)
gj1214 = np.genfromtxt('sflux-GJ1214.txt', names=['lambda','flux'], dtype=None, skip_header=1)

sun = np.genfromtxt('VPL_solar.txt', names=['lambda','flux'], dtype=None, skip_header=1)
moses_HD189 = np.genfromtxt('sflux-HD189_Moses11.txt', names=['lambda','flux'], dtype=None, skip_header=1)

#Eridani = np.genfromtxt('atm/flux-HD189_Moses11.txt', names=['lambda','flux'], dtype=None, skip_header=1)
#actinc_flux = data['variable']['aflux'][-1]


plt.figure()
plt.plot(sun['lambda'],sun['flux']  , c='orange', label='Sun' , alpha = 0.7)
plt.plot(gj876['lambda'],gj876['flux'] , c='r' ,label='GJ876', alpha = 0.7)
plt.plot(gj436['lambda'],gj436['flux'] , c='g' ,label='GJ436', alpha = 0.7)
plt.plot(gj551['lambda'],gj551['flux'] , c='purple' ,label='GJ551', alpha = 0.7)
plt.plot(gj1214['lambda'],gj1214['flux'] , c='b' ,label='GJ1214', alpha = 0.7)
#plt.plot(moses_HD189['lambda'],moses_HD189['flux'] , c='b' ,label='HD189', alpha = 0.7)

plt.gca().set_yscale('log') 
#plt.gca().invert_yaxis() 
plt.xlim((0,500.))
#plt.ylim((1e6,1.1e18))
plt.legend(frameon=0, prop={'size':12}, loc=1)
plt.xlabel("wavelength (nm)")
plt.ylabel("Flux (ergs s-1 nm-1 cm-2)")
plt.savefig('../../plot/' + plot_name + '.png')
plt.savefig('../../plot/' + plot_name + '.eps')
#if vulcan_cfg.use_PIL == True:
plot = Image.open('../../plot/' + plot_name + '.png')
plot.show()
#else: plt.show()

