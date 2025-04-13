import numpy as np

#species list
spec_list = ['H', 'H2O', 'OH', 'H2', 'O', 'CH', 'C', 'CH2', 'CH3', 'CH4', 'C2', 'C2H2', 'C2H', 'C2H3', 'C2H4', 'C2H5', 'C2H6', 'CO', 'CO2', 'CH2OH', 'H2CO', 'HCO', 'CH3O', 'CH3OH', 'CH3CO', 'O2', 'H2CCO', 'HCCO', 'He']
#spec_list = ['H']

ost = ''
for sp in spec_list:
    with open(sp+'.txt') as f:
        ost += sp + '\n'
        lines = f.read()
        ost += lines + '\n' + '\n'
        #print (lines)
       
with open('issi_nasa9.txt', "w") as fout:      # This removes the file contents
    fout.write(ost)