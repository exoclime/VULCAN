'''
Aerosol opacity module for VULCAN (Elspeth K.H. Lee October 2024)

See:
 
YY

for details

Uses the miepython package as the Mie theory calculator.
Follows similar layout to gCMCRT method to calculate the aerosol opacity.

'''

import numpy as np
import matplotlib.pylab as plt
import miepython.miepython as miepython

import vulcan_cfg

def read_nk(sp,wl_nm):

  # Number of wavelengths
  nwl = len(wl_nm)

  # Convert to um from nm
  wl = np.zeros(nwl)
  wl[:] = wl_nm[:] * 1e-3

  # Read nk file
  nk_data = np.loadtxt('./nk/'+sp+'.txt',skiprows=5)
  wl_file = nk_data[:,0]
  n_file = nk_data[:,1]
  k_file = nk_data[:,2]

  # Allocate space for n and k
  n = np.zeros(nwl)
  k = np.zeros(nwl)

  # Linearly interpolate to VULCAN wavelength grid 
  # assume constant values outside interpolation range
  n = np.interp(wl,wl_file,n_file)
  k = np.interp(wl,wl_file,k_file)

  return n, k

def calc_mie(nz,wl_nm,n,k,nd,rm,sig2):

  # Number of wavelengths
  nwl = len(wl_nm)

  # Convert to cm from nm
  wl_cm = np.zeros(nwl)
  wl_cm[:] = wl_nm[:] * 1e-7

  aer_abs = np.zeros((nz,nwl))
  aer_scat = np.zeros((nz,nwl))
  aer_g = np.zeros((nz,nwl))

  for i in range(nz):

    if ((rm[i] < 1.0e-4) or (nd[i] < 1e-10)):
      aer_abs[i,:] = 0.0
      aer_scat[i,:] = 0.0 
      aer_g[i,:] = 0.0
      continue

    if (vulcan_cfg.idist == 1):
      # Single particle size calculation

      r_cm = rm[i] * 1e-4
      xsec = np.pi * r_cm**2

      x = (2.0*np.pi * r_cm)/wl_cm
      m = np.empty(nwl,dtype=np.complex128)
      m.real = n
      m.imag = -k

      qext, qsca, qback, g = miepython.mie(m,x)


    aer_abs[i,:] = (qext[:] - qsca[:])  * xsec 
    aer_scat[i,:] = qsca[:] * xsec
    aer_g[i,:] = g[:]

  return aer_abs, aer_scat, aer_g

# N=500
# m=1.5
# x = np.linspace(0.1,20,N)  # also in microns

# qext, qsca, qback, g = miepython_jit.mie(m,x)

# plt.figure(figsize=(8,4.5))
# plt.plot(x,g)
# plt.show()