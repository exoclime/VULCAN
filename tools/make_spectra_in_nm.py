import numpy as np
import scipy

au = 1.4959787E13  # cm
r_sun = 6.957E10 # cm

# Epsilon Eridani is 10.475 light years away and with 0.735 solar radius
# GJ876 is 15.2 light years away and has 0.3761 solar radius

new_str = '# WL(nm)    Flux(ergs/cm**2/s/nm)\n'

with open('../atm/stellar_flux/Gueymard_solar.txt') as f:
    for line in f.readlines():
        if not line.startswith("#") and line.split():
            li = line.split()
            if float(li[0]) < 115.:
                new_str += '{:<8.3f}'.format(float(li[0])) + "{:>12.2E}".format(float(li[1])) + '\n'
            else: break
                       	
#with open('h_epseri_uvsum_spc.txt') as f: # wl was in Angstroms in the file
with open('../atm/stellar_flux/epseri_01nm.txt') as f: # wl was in Angstroms in the file
    for line in f.readlines():
        if not line.startswith("#") and line.split(): 
            li = line.split()
            #wl = float(li[0])*0.1
            #flux = float(li[1])*10. *(10.475*63241*au/r_sun*0.735)**2     # eps Eradian is 10.475 light years away
            wl = float(li[0])
            flux = float(li[1])
            if flux > 0:
                new_str += '{:<8.3f}'.format(float(wl)) + "{:>12.2E}".format(flux) + '\n'

with open('../atm/stellar_flux/Gueymard_solar.txt') as f:
    for line in f.readlines():
        if not line.startswith("#") and line.split():
            li = line.split()
            if float(li[0]) >= 283.:
                new_str += '{:<8.3f}'.format(float(li[0])) + "{:>12.2E}".format(float(li[1]) *0.1) + '\n'
            else: pass

    
#with open('sflux-HD189_Moses11.txt', 'w+') as f: f.write(new_str)   
with open('../atm/stellar_flux/sflux-HD189_Moses11.txt', 'w+') as f: f.write(new_str)  