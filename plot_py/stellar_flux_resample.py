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
from scipy import interpolate
       
plot_name = 'stellar-flux-resample'

# *(vulcan_cfg.r_star*r_sun/(au*vulcan_cfg.orbit_radius) )**2 
#sflux_raw = np.genfromtxt('../atm/stellar_flux/Gueymard_solar.txt', dtype=float, skip_header=1, names = ['lambda','flux'])
sflux_raw = np.genfromtxt('../atm/stellar_flux/sflux-epseri.txt', dtype=float, skip_header=1, names = ['lambda','flux'])


# Stellar flux at TOA; not converting to actinic flux yet *1/(hc/ld)
# for values outside the boundary => fill_value = 0
inter_sflux = interpolate.interp1d(sflux_raw['lambda'], sflux_raw['flux'], bounds_error=False, fill_value=0)
inter_quard_sflux = interpolate.interp1d(sflux_raw['lambda'], sflux_raw['flux'], bounds_error=False, fill_value=0, kind='slinear')


dbin = 0.1
#bins = np.arange(2.89,505.1,dbin)
bins = np.arange(115.,283,dbin)
last_bin = bins[-1]
sflux_top = np.zeros(len(bins))

for n, ld in enumerate(bins):
    # define the next bin in the new uniform grid
    if ld != bins[-1]: next_ld = bins[n+1]
    
    if ld in sflux_raw['lambda'] or ld == last_bin or sflux_raw['lambda'][np.searchsorted(sflux_raw['lambda'],ld,side='right')] - sflux_raw['lambda'][np.searchsorted(sflux_raw['lambda'],ld,side='right')-1] > dbin: # when bin coincide with raw_bin or dbin is smaller than the width of raw_bin
        sflux_top[n] = inter_sflux(ld)

    else:
        
        # finding the index for the left & right pts in the raw data
        raw_left_indx = np.searchsorted(sflux_raw['lambda'],ld,side='right')
        raw_right_indx = np.searchsorted(sflux_raw['lambda'],next_ld,side='right') - 1
        flux_left = inter_sflux(ld)
        flux_right = inter_sflux(next_ld)
        
        # trapezoid integral
        bin_flux = (flux_left + sflux_raw['flux'][raw_left_indx])*0.5*(sflux_raw['lambda'][raw_left_indx]-ld)
        # print (ld)
        # print (flux_left)
        # print (sflux_raw['flux'][raw_left_indx])
        # print ( sflux_raw['lambda'][raw_left_indx]-ld )
        # print (bin_flux/((sflux_raw['lambda'][raw_left_indx]-ld)) )
        bin_flux += (flux_right + sflux_raw['flux'][raw_right_indx])*0.5*(next_ld-sflux_raw['lambda'][raw_right_indx])

    
        if raw_right_indx - raw_left_indx > 0:
            #print (ld)
            #print (raw_right_indx - raw_left_indx)
            for raw_i in range(raw_left_indx,raw_right_indx):
                bin_flux += (sflux_raw['flux'][raw_i] + sflux_raw['flux'][raw_i+1])*0.5*(sflux_raw['lambda'][raw_i+1] - sflux_raw['lambda'][raw_i])
              
        bin_flux = bin_flux/dbin
        
        #print (( sflux_raw['lambda'][raw_right_indx]-sflux_raw['lambda'][raw_left_indx] + sflux_raw['lambda'][raw_left_indx]-ld + sflux_raw['lambda'][raw_right_indx]-next_ld ))
        sflux_top[n] = bin_flux
        
# Check for eneergy conservation
# finding the index for the left & right pts in the raw data
raw_left_indx = np.searchsorted(sflux_raw['lambda'],bins[0],side='right')
raw_right_indx = np.searchsorted(sflux_raw['lambda'],bins[-1],side='right')-1
sum_orgin, sum_bin = 0, 0

sum_bin2 = 0
for n in range(raw_left_indx,raw_right_indx):
    sum_orgin += 0.5*(sflux_raw['flux'][n] + sflux_raw['flux'][n+1]) * (sflux_raw['lambda'][n+1]- sflux_raw['lambda'][n])

if sflux_raw['lambda'][raw_left_indx]-bins[0] <0: print ('left index outside of bound')
if bins[-1]-sflux_raw['lambda'][raw_right_indx] <0: print ('right index outside of bound')
sum_orgin += 0.5 *(inter_sflux(bins[0])+sflux_raw['flux'][raw_left_indx])* (sflux_raw['lambda'][raw_left_indx]-bins[0])
sum_orgin += 0.5 *(inter_sflux(bins[-1])+sflux_raw['flux'][raw_right_indx])* (bins[-1]-sflux_raw['lambda'][raw_right_indx])   

sflux_inter, sflux_inter2 = np.zeros(len(bins)), np.zeros(len(bins))
for n, ld in enumerate(bins):
    sflux_inter[n] = inter_sflux(ld) 
    sflux_inter2[n] = inter_quard_sflux(ld)
    
sum_old = dbin * np.sum(sflux_inter)
sum_old -= dbin*0.5*(sflux_inter[0]+sflux_inter[-1])

sum_old2 = dbin * np.sum(sflux_inter2)
sum_old2 -= dbin*0.5*(sflux_inter2[0]+sflux_inter2[-1])

sum_bin = dbin * np.sum( sflux_top )
sum_bin -= dbin*0.5*(sflux_top[0]+sflux_top[-1])    
    
print ("conservation check:")
print (sum_bin/sum_orgin)
print (sum_old/sum_orgin)
#print (sum_old2/sum_orgin)
plt.figure()
plt.plot(sflux_raw['lambda'], sflux_raw['flux'], c='orange', label='original' , lw=1, alpha = 0.7)                                              
plt.plot(bins, sflux_top, label='resample' , lw=1, alpha = 0.7, c='b')

plt.plot(bins, sflux_inter, label='interpolation' , lw=1, alpha = 0.7, c='grey')
#plt.plot(bins, sflux_inter2, label='interpolation2' , lw=1, alpha = 0.7, c='green')

plt.gca().set_yscale('log') 
#plt.gca().invert_yaxis() 
plt.xlim((125,150.))
#plt.ylim(ymin = 20)
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

new_str = '# WL(nm)    Flux(ergs/cm**2/s/nm)\n'
for i,ld in enumerate(bins):
    new_str += '{:<8.3f}'.format(float(ld)) + "{:>12.2E}".format(sflux_top[i]) + '\n'
new_str = new_str[:-1]
with open('epseri_01nm.txt', 'w+') as f: f.write(new_str)  