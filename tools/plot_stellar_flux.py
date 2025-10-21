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
from phy_const import au, r_sun
       
plot_name = 'Stellar-flux-compare'
#vul_data = 'output/HCO_new_H2O_HD189_Moses_nominalKzz.vul'
#vul_data2 = 'output/test.vul'

# with open(vul_data, 'rb') as handle:
#   data = pickle.load(handle)
# with open(vul_data2, 'rb') as handle:
#   data2 = pickle.load(handle)

#kelt9 = np.genfromtxt('atm/Kelt9_flux_phoneix.txt', names=['lambda','flux'], dtype=None)
sun = np.genfromtxt('../atm/stellar_flux/VPL_solar.txt', names=['lambda','flux'], dtype=None, skip_header=1)
GJ436 = np.genfromtxt('../atm/stellar_flux/sflux-GJ436.txt', names=['lambda','flux'], dtype=None, skip_header=1)
eri = np.genfromtxt('../atm/stellar_flux/h_hd22049_uvsum_1x_51620_etc.txt', names=['lambda','flux'], dtype=None, skip_header=1)
              
Gueymard_2013 = np.genfromtxt('sun_Gueymard.txt', names=['lambda','flux'], dtype=['float','float'], skip_header=0)
Gueymard_2018 = np.genfromtxt('Gueymard_2018.txt', names=['lambda','flux'], dtype=['float','float'], skip_header=12)


#actinc_flux = data['variable']['aflux'][-1]

#moses_HD189 = np.genfromtxt('atm/HD189_Moses11.txt', names=['lambda','flux'], dtype=None, skip_header=1)

plt.figure()
plt.plot(sun['lambda'],sun['flux']  , c='orange', label='Sun VPL' , lw=1.2, alpha = 0.7)

# scaled back to the surface of the Sun
#plt.plot(Gueymard_2013['lambda'], Gueymard_2013['flux']*1e7/1e4 *(au/r_sun)**2, c='b', label='Sun Gueymard' , lw=0.8, alpha = 0.6)
plt.plot(Gueymard_2018['lambda'], Gueymard_2018['flux']*1e7/1e4 *(au/r_sun)**2, c='purple', label='Sun Gueymard18' , lw=0.9, alpha = 0.7)


output = open('Gueymard2018_solar.txt', "w")

ost = '# WL(nm)     Flux(ergs/cm**2/s/nm)  Synthetic/composite Extraterrestrial Spectrum, April 2018\n'
for n,lda in enumerate(Gueymard_2018['lambda']):
    ost += '{:<10.1f}'.format(lda) + '{:>15.3E}'.format(Gueymard_2018['flux'][n] *1e7/1e4 *(au/r_sun)**2 ) + '\n'
ost = ost[:-1]
output.write(ost)
output.close()

#plt.plot(GJ436['lambda'],GJ436['flux'] , c='r' ,label='GJ436', alpha = 0.7, lw=0.7 )

au = 1.4959787E13  # cm
r_sun = 6.957E10 # cm

#plt.plot(eri['lambda']/10.,eri['flux']*10. *(662449.475*au/r_sun*0.735)**2 , c='b' ,label='$\epsilon$ Eri', alpha = 0.7, lw=0.7 )


#plt.plot(kelt9['lambda'][::10],kelt9['flux'][::10]*kelt9['lambda'][::10] *(2.362*r_sun/(620.*9.4607E17))**2 ,label='KELT9' ) #*kelt9['lambda']  *(2.362*r_sun/620./au)**2
#
#plt.plot(data['variable']['bins'], actinc_flux)
#plt.plot(moses_HD189['lambda'][::10],moses_HD189['flux'][::10]  ,label='HD189' )
# plt.plot(data['variable']['bins'], data['variable']['aflux'][-1], c='g', label='old')
# plt.plot(data2['variable']['bins'], data2['variable']['aflux'][-1], c='red')
     
plt.gca().set_yscale('log') 
#plt.gca().invert_yaxis() 
plt.xlim((1,600.))
plt.ylim(ymin = 20)
#plt.ylim((1e6,1.1e18))
plt.legend(frameon=0, prop={'size':14}, loc='best')
plt.xlabel("wavelength (nm)")
plt.ylabel("Flux (ergs s-1 nm-1 cm-2)")
plt.savefig('../plot/' + plot_name + '.png')
plt.savefig('../plot/' + plot_name + '.pdf')
if vulcan_cfg.use_PIL == True:
    plot = Image.open('../plot/' + plot_name + '.png')
    plot.show()
else: plt.show()

