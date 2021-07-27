import numpy as np
import scipy

### PHYS CONSTANTS
# Planck constant times the light speed
hc = 1.98644568E-9 # erg.nm
au = 1.4959787E13  # cm
r_sun = 6.957E10 # cm
###

# VPL SUN: W/cm^2/micron
new_str = '# solar flux at the "surface of Sun" from VPL\n# WL(nm)\t Flux(ergs/cm**2/s/nm)\n'

with open('vpl_sun_original.txt') as f:
    for line in f.readlines():
        if not line.startswith("#") and line.split():
            li = line.split()
            new_str += '{:<12}'.format(float(li[0])*1.e3) + '\t' + "{:.2E}".format(float(li[1])*1.e4 * (au/r_sun)**2.  ) + '\n'   
            
   
with open('VPL_solar.txt', 'w+') as f: f.write(new_str)   
