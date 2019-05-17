import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.legend as lg
#import vulcan_cfg
from PIL import Image
from scipy.signal import savgol_filter
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
from scipy import interpolate

import os, sys
import pickle
       
plot_name = 'cross_sections'
#vul_data = 'output/test.vul'
plot_dir = '../plot/'

scat_cross = {}
for n, sp in enumerate(['H2','He', 'N2','O2','CO2']):
    scat_cross[sp] = np.genfromtxt(sp+'_rayleigh.txt',dtype=float,skip_header=1, names = ['lambda','cross'])

#scat_cross['H2'] = np.genfromtxt('H2_rayleigh.txt',dtype=float,skip_header=1, names = ['lambda','cross'])       
bins = scat_cross['H2']['lambda']

# smooth test
o2_smooth = savgol_filter(scat_cross['O2']['cross'], 1001, 3, mode='wrap')
o2_segments = np.append(scat_cross['O2']['cross'][0:750], scat_cross['O2']['cross'][2300:]) 
o2_bins = np.append(bins[0:750], bins[2300:]) 

n2_segments = np.append(scat_cross['N2']['cross'][0:205], scat_cross['N2']['cross'][1400:]) 
n2_bins = np.append(bins[0:205], bins[1400:]) 

co2_segments = np.append(scat_cross['CO2']['cross'][0:205], scat_cross['CO2']['cross'][2200:]) 
co2_bins = np.append(bins[0:205], bins[2200:]) 


he_segments = np.append(scat_cross['He']['cross'][0:250], scat_cross['He']['cross'][1300:]) 
he_bins = np.append(bins[0:250], bins[1300:]) 

o2_test = InterpolatedUnivariateSpline(o2_bins, o2_segments)
o2_test2 = interpolate.splrep(o2_bins, o2_segments, s=10000)

o2_fit = interpolate.interp1d(o2_bins, np.log10(o2_segments),kind='cubic')
n2_fit = interpolate.interp1d(n2_bins, np.log10(n2_segments),kind='quadratic')
co2_fit = interpolate.interp1d(co2_bins, np.log10(co2_segments),kind='cubic')

he_fit = interpolate.interp1d(he_bins, np.log10(he_segments),kind='cubic')

ost = '#lambda (nm)  cross(cm^2)\n'
for lm in bins:
   ost += "{0:<14s}".format(str(lm)) + "{:.3E}".format(float(10**co2_fit(lm))) +'\n'

ost = ost[:-1]
with open('CO2_smooth.txt', 'w') as f: f.write(ost)




plt.figure()   
plt.plot(scat_cross['H2']['lambda'], scat_cross['H2']['cross'], label='H2')

plt.plot(scat_cross['He']['lambda'], scat_cross['He']['cross'], label='He', c='purple')
plt.plot(bins, 10**he_fit(bins), c='purple', lw=1.5, ls='--')

plt.plot(bins, scat_cross['N2']['cross'], c='green', label='N2')
plt.plot(bins, 10**n2_fit(bins), c='green', lw=1.5, ls='--')
plt.plot(bins, scat_cross['O2']['cross'], c='red', label='O2')
plt.plot(bins, 10**o2_fit(bins), c='red', lw=1.5, ls='--')
plt.plot(bins, scat_cross['CO2']['cross'], c='black', label='CO2')
plt.plot(bins, 10**co2_fit(bins), c='black', lw=1.5, ls='--')
#plt.plot(data['variable']['bins'], photosph_abs, c='red')
           
plt.gca().set_yscale('log') 
#plt.gca().invert_yaxis() 
plt.xlim(xmax=300)
#plt.ylim((1.E3,1.E-8))
plt.legend(frameon=0, prop={'size':14}, loc='best')
plt.xlabel("wavelength (nm)")
plt.ylabel("cross-sections (cm2)")
plt.savefig(plot_dir + plot_name + '.png')
plt.savefig(plot_dir + plot_name + '.eps')
#if vulcan_cfg.use_PIL == True:
plot = Image.open(plot_dir + plot_name + '.png')
plot.show()
#else: plt.show()

