# import sys
# sys.path.insert(0, '../') # including the upper level of directory for the path of modules

import numpy as np 
import os, sys
import pickle
       
output = open('../output/mix_table/vulcan_test6_photo.txt', "w")
vul = '../output/ISSI_test6.vul'

with open(vul, 'rb') as handle:
  vul = pickle.load(handle)
  

bins = vul['variable']['bins']

ost = '{:<8s}'.format('(bar)')  + '{:>9s}'.format('Wavelength (nm)') + '\n'
ost += '{:<8s}'.format('Pressure')
for lda in bins: ost += '{:>11.0f}'.format(lda)
ost += '\n' 

for n, p in enumerate(vul['atm']['pco']):
    ost += '{:<8.3E}'.format(p/1e6)
    for i in range(len(bins)):
        ost += '{:>11.3e}'.format(vul['variable']['aflux'][n][i])
    ost += '\n'

ost = ost[:-1]
output.write(ost)
output.close()

# Print the photodissociation rate
output = open('../output/mix_table/vulcan_test6_J-rate.txt', "w")
photo_sp = ['H2O', 'CH4', 'CH3', 'CO', 'H2', 'C2H2', 'CO2', 'C2H4', 'C2H6', 'OH', 'HCO', 'H2CO']


ost = '{:<8s}'.format('(bar)') + '  photodissociation rate J (s−1) for each species\n'
ost += '{:<8s}'.format('Pressure')
for sp in photo_sp: ost += '{:>11s}'.format(sp) 
ost +='\n'

for n, p in enumerate(vul['atm']['pco']):
    ost += '{:<8.3E}'.format(p/1e6)
    for sp in photo_sp: ost += '{:>11.3e}'.format(vul['variable']['J_sp'][(sp,0)][n])
    ost += '\n'
ost = ost[:-1]
output.write(ost)
output.close()