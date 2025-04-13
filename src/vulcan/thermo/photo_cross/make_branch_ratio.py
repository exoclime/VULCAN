# -*- coding: utf-8 -*-
"""
Making the braching ratio files
"""
import numpy as np
from matplotlib import pyplot as plt

molecule = 'SO2'
# number of branches
num_br = 2

# 907    [ SO2 -> SO + O                      ]              SO2     1
# 909    [ SO2 -> S + O2                      ]              SO2     2

outstr = '# Branching ratios for '+molecule+' -> (1)SO + O  (2)S+ O2  nm  from PhiDRates\n'
outstr += '# lambda '
for i in range(1,num_br+1): outstr += "{:>12}".format(', br_ratio_' + str(i))
outstr += '\n'

# Lambda  Total   SO/O     S/O2     SO2band  SO2+    

phid = np.genfromtxt('../../../../../../Bern/python/photo_cross_data/phidrates_data/'+molecule+'.txt',dtype=float,skip_header=1, names = ['lambda', 'tot', '1', '2'], usecols = [0,1, 2,3])

for n,lmd in enumerate(phid['lambda']):
    if lmd/10 > 0:
        dis_tot = phid['1'][n] + phid['2'][n] #+ phid['3'][n] #+ phid['4'][n]
        #print ("{:<9.3e}".format(lmd/10.) + "," + "{:>9.3f}".format(float(phid['1'][n])/dis_tot) + "," + "{:>9.3f}".format(float(phid['2'][n])/dis_tot)+ "," + "{:>9.3f}".format(float(phid['5'][n])/dis_tot) +'\n')
        
        if dis_tot > 0:
            outstr +=  "{:<9.3e}".format(lmd/10.) + "," + "{:>9.3f}".format(float(phid['1'][n])/dis_tot) + "," + "{:>9.3f}".format(float(phid['2'][n])/dis_tot) +'\n'  #+ "," + "{:>9.3f}".format(float(phid['3'][n])/dis_tot)+ "," + "{:>9.3f}".format(float(phid['4'][n])/dis_tot) + '\n'


# outstr +=  "{:<9.3e}".format(121.) + "," + "{:>9.3f}".format(0.) + "," + "{:>9.3f}".format(0.41) + "," + "{:>9.3f}".format(0.26) + "," + "{:>9.3f}".format(0.16)+ "," + "{:>9.3f}".format(0.17) + '\n'
# outstr +=  "{:<9.3e}".format(121.6) + "," + "{:>9.3f}".format(0.15) + "," + "{:>9.3f}".format(0.31) + "," + "{:>9.3f}".format(0.25) + "," + "{:>9.3f}".format(0.26)+ "," + "{:>9.3f}".format(0.03) + '\n'
# outstr +=  "{:<9.3e}".format(121.6) + "," + "{:>9.3f}".format(0.15) + "," + "{:>9.3f}".format(0.31) + "," + "{:>9.3f}".format(0.25) + "," + "{:>9.3f}".format(0.26)+ "," + "{:>9.3f}".format(0.03) + '\n'

# outstr +=  "{:<9.3e}".format(210.) + "," + "{:>9.3f}".format(0.) + "," + "{:>9.3f}".format(1.) + '\n'
# outstr +=  "{:<9.3e}".format(220.) + "," + "{:>9.3f}".format(0.1) + "," + "{:>9.3f}".format(0.9) + '\n'
# outstr +=  "{:<9.3e}".format(300.) + "," + "{:>9.3f}".format(0.1) + "," + "{:>9.3f}".format(0.9) + '\n'
# outstr +=  "{:<9.3e}".format(306.) + "," + "{:>9.3f}".format(1.-0.875) + "," + "{:>9.3f}".format(0.875) + '\n'
# outstr +=  "{:<9.3e}".format(307.) + "," + "{:>9.3f}".format(1.-0.844) + "," + "{:>9.3f}".format(0.844) + '\n'
# outstr +=  "{:<9.3e}".format(308.) + "," + "{:>9.3f}".format(1.-0.76) + "," + "{:>9.3f}".format(0.76) + '\n'
# outstr +=  "{:<9.3e}".format(309.) + "," + "{:>9.3f}".format(1.-0.616) + "," + "{:>9.3f}".format(0.616) + '\n'
# outstr +=  "{:<9.3e}".format(310.) + "," + "{:>9.3f}".format(1.-0.443) + "," + "{:>9.3f}".format(0.443) + '\n'
# outstr +=  "{:<9.3e}".format(311.) + "," + "{:>9.3f}".format(1.-0.298) + "," + "{:>9.3f}".format(0.298) + '\n'
# outstr +=  "{:<9.3e}".format(312.) + "," + "{:>9.3f}".format(1.-0.208) + "," + "{:>9.3f}".format(0.208) + '\n'
# outstr +=  "{:<9.3e}".format(313.) + "," + "{:>9.3f}".format(1.-0.162) + "," + "{:>9.3f}".format(0.162) + '\n'
# outstr +=  "{:<9.3e}".format(314.) + "," + "{:>9.3f}".format(1.-0.143) + "," + "{:>9.3f}".format(0.143) + '\n'
# outstr +=  "{:<9.3e}".format(315.) + "," + "{:>9.3f}".format(1.-0.136) + "," + "{:>9.3f}".format(0.136) + '\n'
# outstr +=  "{:<9.3e}".format(316.) + "," + "{:>9.3f}".format(1.-0.133) + "," + "{:>9.3f}".format(0.133) + '\n'
# outstr +=  "{:<9.3e}".format(317.) + "," + "{:>9.3f}".format(1.-0.129) + "," + "{:>9.3f}".format(0.129) + '\n'
# outstr +=  "{:<9.3e}".format(318.) + "," + "{:>9.3f}".format(1.-0.123) + "," + "{:>9.3f}".format(0.123) + '\n'
# outstr +=  "{:<9.3e}".format(319.) + "," + "{:>9.3f}".format(1.-0.116) + "," + "{:>9.3f}".format(0.116) + '\n'
# outstr +=  "{:<9.3e}".format(320.) + "," + "{:>9.3f}".format(1.-0.109) + "," + "{:>9.3f}".format(0.109) + '\n'
# outstr +=  "{:<9.3e}".format(321.) + "," + "{:>9.3f}".format(1.-0.101) + "," + "{:>9.3f}".format(0.101) + '\n'
# outstr +=  "{:<9.3e}".format(322.) + "," + "{:>9.3f}".format(1.-0.095) + "," + "{:>9.3f}".format(0.095) + '\n'
# outstr +=  "{:<9.3e}".format(323.) + "," + "{:>9.3f}".format(1.-0.089) + "," + "{:>9.3f}".format(0.089) + '\n'
# outstr +=  "{:<9.3e}".format(324.) + "," + "{:>9.3f}".format(1.-0.085) + "," + "{:>9.3f}".format(0.085) + '\n'
# outstr +=  "{:<9.3e}".format(325.) + "," + "{:>9.3f}".format(1.-0.082) + "," + "{:>9.3f}".format(0.082) + '\n'
# outstr +=  "{:<9.3e}".format(326.) + "," + "{:>9.3f}".format(1.-0.080) + "," + "{:>9.3f}".format(0.080) + '\n'
# outstr +=  "{:<9.3e}".format(121.6) + "," + "{:>9.3f}".format(0.42) + "," + "{:>9.3f}".format(0.48) + "," + "{:>9.3f}".format(0.03) + "," + "{:>9.3f}".format(0.07) + '\n'

# for H2CO -> H + HCO    
# def hco_2(lmd):
#     return 557.95835182 - 7.31994058026*lmd + 0.03553521598*lmd**2 -7.54849718e-5*lmd**3 +5.91001021e-8*lmd**4


# for ld in np.arange(250.,339., 1.):
#     outstr +=  "{:<9.3e}".format(ld) + "," + "{:>9.3f}".format(1.-hco_2(ld))+ ","  + "{:>9.3f}".format(hco_2(ld)) +'\n'
    

# 314    0.143
# 315    0.136
# 316    0.133
# 317    0.129
# 318    0.123
# 319    0.116
# 320    0.109
# 321    0.101
# 322    0.095
# 323    0.089
# 324    0.085
# 325    0.082
# 326    0.080
# 327    0.079
# 328    0.078
# 330    0.080
# 335    0.080
# 340    0.080



with open(molecule+'/'+ molecule+"_branch.csv","w") as of:
    outstr = outstr[:-1]
    of.write(outstr)
print (molecule+' done')