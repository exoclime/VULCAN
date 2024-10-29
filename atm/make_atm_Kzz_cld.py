import numpy as np

line_1 = '#(dyne/cm2) (K)          (cm2/s)    (cm-3)   (um)   (cm2)'
line_2 = 'Pressure        Temp     Kzz   nd   rm   sig2'

fname = 'atm_GJ1214b_Kzz.txt'


data = np.loadtxt(fname,skiprows=2)

p = data[:,0]
T = data[:,1]
Kzz = data[:,2]

nz = len(p)

# Set up aerosol component

nd = np.zeros(nz)
rm = np.zeros(nz)
sig2 = np.zeros(nz)


for i in range(nz):
  if (p[i]/1e6 < 1e-1):
    nd[i] = 1e5
    rm[i] = 1e-2
    sig2[i] = 1.0

fout = 'atm_GJ1214b_Kzz_aer.txt'
f = open(fout, 'w')

f.write(line_1 + '\n')
f.write(line_2 + '\n')
for i in range(nz):
  f.write(str(p[i]) + ' ' + str(T[i]) + ' ' + str(Kzz[i]) + ' ' + str(nd[i]) + ' ' + str(rm[i]) + ' ' + str(sig2[i]) + '\n')



