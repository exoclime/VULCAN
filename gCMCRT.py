'''
1D plane-parallel Monte Carlo radiative-transfer for VULCAN (Elspeth K.H. Lee May 2024)

See:
 
https://ui.adsabs.harvard.edu/abs/2024RNAAS...8...96L/abstract 

for details

Uses a 1D Lucy path length method to calculate the mean intensity, which can then be used to calculate dissociation rates
Useful for when high accuracy is required.
We use JIT to accelerate the calculations - this can make the first iteration slow, but later iterations fast.

NOTE: Modifications will be required when non-isotropic scattering is required e.g. aerosol opacity + asymmetry factor etc
Default is Rayleigh scattering
'''

from random import seed, random
import numpy as np
from numba import jit, config, int32, float64, prange
from numba.experimental import jitclass

## Flag to disable JIT
config.DISABLE_JIT = False

## Definitions of packet properties required by numba
pac_def = [
    ('flag', int32),      
    ('id', int32),
    ('cost', float64),
    ('nzp', float64),
    ('zc', int32),
    ('zp', float64),
    ('tau_p', float64),
    ('tau', float64),
    ('iscat',int32),
]

## Class that defines packet properties
@jitclass(pac_def)
class pac:
  def __init__(self, flag, id_ph, cost, nzp, zc, zp, tau_p, tau, iscat):
    self.flag = flag
    self.id = id_ph
    self.cost = cost
    self.nzp = nzp
    self.zc = zc
    self.zp = zp
    self.tau_p = tau_p
    self.tau = tau
    self.iscat = iscat

## Function that integrates a packet through the 1D plane-parallel grid
@jit(nopython=True, cache=True)
def tauint_1D_pp(ph, nlay, z, sig_ext, l, Jdot):

  # Initial tau is at zero
  ph.tau = 0.0

  # Do loop until tau is equal or more than sampled tau
  while (ph.tau < ph.tau_p):

    # Calculate dsz, the distance to the next vertical level
    if (ph.nzp > 0.0):
      # Packet traveling upward, find distance to upper level
      dsz = (z[ph.zc+1]-ph.zp)/ph.nzp
      zoffset = 1
    elif (ph.nzp < 0.0):
      # Packet traveling downward, find distance to lower level
      dsz = (z[ph.zc]-ph.zp)/ph.nzp
      zoffset = -1
    else:
      # Packet traveling directly in z plane
      # Return, packet does not move in z direction
      break
    
    # Calculate optical depth to level edge
    taucell = dsz * sig_ext[ph.zc,l]

    # Check if packet ends path in this layer
    if ((ph.tau + taucell) >= ph.tau_p):

      # Packet stops in this cell - move distance then exit loop
      d1 = (ph.tau_p-ph.tau)/sig_ext[ph.zc,l]
      ph.zp += d1 * ph.nzp

      # Update estimator
      Jdot[ph.zc,l] += d1 # Mean intensity estimator

      # tau of packet is now the sampled tau
      ph.tau = ph.tau_p

    else:

      #Packet continues to level edge - update position, cell index and tau
      ph.zp += (dsz + 1.0e-12) * ph.nzp

      # Update estimator
      Jdot[ph.zc,l] += dsz # Mean intensity estimator

      # Apply integer offset to cell number
      ph.zc += zoffset

      # Add tau of cell to tau counter of packet
      ph.tau += taucell

      # Check is packet has exited the domain
      if ((ph.zc > nlay-1) or (ph.zp >= z[-1])):
        # Packet has exited the domain at the top of atmosphere
        ph.flag = 1
        break
      elif((ph.zc < 0) or (ph.zp <= z[0])):
        # Packet has hit the lower part od the domain (surface)
        ph.flag = -2
        break

  return

## Function for incident stellar radiation properties
@jit(nopython=True, cache=True)
def inc_stellar(ph, nlay, z, mu_z):

  ph.cost = -mu_z   # Initial direction is zenith angle
  ph.nzp = ph.cost  # give nzp cost
  ph.zc = nlay-1    # Initial layer number is top of atmosphere
  ph.zp = z[-1] - 1.0e-12 # Initial position is topic atmosphere minus a little bit

  return

## Function to scatter a packet into a new direction
@jit(nopython=True, cache=True)
def scatter(ph, g):

  # Find the type of scattering the packet undergoes 
  # (1 = isotropic, 2 = Rayleigh, 3 = Henyey-Greenstein 
  # 4 = Two-Term Henyey-Greenstein, 5 = Draine 2003)
  match ph.iscat:

    case 1:
      # Isotropic scattering
      ph.cost = 2.0 * random() - 1.0
      ph.nzp = ph.cost
      return
    
    case 2:
      # Rayleigh scattering via direct spherical coordinate sampling
      # Assumes non-polarised incident packet
      q = 4.0*random() - 2.0
      u = np.cbrt(-q + np.sqrt(1.0 + q**2))
      bmu = u - 1.0/u

    case 3:
      # Sample from single HG function
      if (abs(g) > 1e-4):
        # Avoid division by zero - is isotropic as g -> 0 anyway.
        g2 = g**2

        bmu = ((1.0 + g2) - \
          ((1.0 - g2) / (1.0 - g + 2.0 * g * random()))**2) \
          / (2.0*g)
      else:
        # g ~ 0, Isotropic scattering
        ph.cost = 2.0 * random() - 1.0
        ph.nzp = ph.cost
        return
    case 4:
      # Sample from two-term HG function following Cahoy et al. (2010)
      # Check if near isotropic and sample isotropic
      if (abs(g) < 1e-4):
        # g ~ 0, Isotropic scattering
        ph.cost = 2.0 * random() - 1.0
        ph.nzp = ph.cost
        return
      
      gb = -g/2.0
      alph = 1.0 - gb**2

      zeta1 = random()
      if (zeta1 < alph):
        # sample forward direction
        g2 = g**2
        bmu = ((1.0 + g2) - \
          ((1.0 - g2) / (1.0 - g + 2.0 * g * random()))**2) \
          / (2.0*g)
      else:
        # sample backward direction
        gb2 = gb**2
        bmu = ((1.0 + gb2) - \
          ((1.0 - gb2) / (1.0 - gb + 2.0 * gb * random()))**2) \
          / (2.0*gb)
    case 5:
      # Sample from Draine (2003) phase function (= Cornette and Shanks (1992) when alpha = 1)
      # Use the analytical approach of Jendersie & d’Eon (2023)

      # First calculate big G from small g
      alpha = 1.0
      a1 = 1.0/2.0  + 5.0/(6.0*alpha) - (25.0/81.0) * g**2
      b1 = (125.0/729.0)*g**3 + 5.0/(9.0*alpha) * g
      G1 = np.cbrt(np.sqrt(a1**3 + b1**2) + b1) - np.cbrt(np.sqrt(a1**3 + b1**2) - b1) + (5.0/9.0) * g
      G2 = G1**2
      G4 = G1**4

      # Sample Draine (2003) function analytically (Jendersie & d’Eon 2023)
      zeta1 = random()
      t0 = alpha - alpha*G2
      t1 = alpha*G4 - alpha
      t2 = -3.0*(4.0*(G4 - G2) + t1*(1.0 + G2))
      t3 = G1*(2.0*zeta1 - 1.0)
      t4 = 3.0*G2*(1.0+t3) + alpha*(2.0+G2*(1.0+(1.0+2.0*G2)*t3))
      t5 = t0*(t1*t2+t4**2)+t1**3
      t6 = t0*4.0*(G4 - G2)
      t7 = np.cbrt(t5 + np.sqrt(t5**2 - t6**3))
      t8 = 2.0*(t1+t6/t7+t7)/t0
      t9 = np.sqrt(6.0*(1.0+G2) + t8)
      bmu = G1/2.0 + (1.0/(2.0*G1) - 1.0/(8.0*G1)*(np.sqrt(6.0*(1.0+G2) -  t8 + 8.0*t4/(t0*t9)) - t9)**2)

    case _ :
      # Invalid iscat value
      print('Invalid iscat chosen: ' + str(ph.iscat) + ' doing isotropic')
      ph.cost = 2.0 * random() - 1.0
      ph.nzp = ph.cost
      return
    
  # Now apply rotation if non-isotropic scattering
  # Change direction of packet given by sampled direction
 
  # Calculate change in direction in grid reference frame
  if (bmu >= 1.0):
    # Packet directly forward scatters - no change in direction
    ph.cost = ph.cost
  elif (bmu <= -1.0):
    # Packet directly backward scatters - negative sign for cost
    ph.cost = -ph.cost
  else:
    # Packet scatters according to sampled cosine and current direction
    # Save current cosine direction of packet and calculate sine
    costp = ph.cost
    sintp = 1.0 - costp**2
    if (sintp <= 0.0):
      sintp = 0.0
    else:
      sintp = np.sqrt(sintp)

    # Randomly decide if scattered in +/- quadrant direction and find new cosine direction
    sinbt = np.sqrt(1.0 - bmu**2)
    ri1 = 2.0 * np.pi * random()
    if (ri1 > np.pi):
      cosi3 = np.cos(2.0*np.pi - ri1)
      # Calculate new cosine direction
      ph.cost = costp * bmu + sintp * sinbt * cosi3
    else: #(ri1 <= pi)
      cosi1 = np.cos(ri1)
      ph.cost = costp * bmu + sintp * sinbt * cosi1

    #Give nzp the cost value
    ph.nzp = ph.cost

  return

## Function to scatter from a surface (Assumed Lambertian)
@jit(nopython=True, cache=True)
def scatter_surf(ph, z):
    
    ph.cost = np.sqrt(random()) # Sample cosine from Lambertian cosine law
    ph.nzp = ph.cost            # Give nzp cost value
    ph.zc = 0                   # Initial layer is lowest (0)
    ph.zp = z[0] + 1.0e-12      # Initial position is lowest altitude plus a little bit

    return

## Main gCMCRT function for wavelength and packet loop
@jit(nopython=True, cache=True, parallel=True)
def gCMCRT_main(nit, Nph, nlay, nwl, all_sp, ph_sp, ray_sp, nd_sp, \
  cross_sp, cross_sca_sp, nd_aer, cross_aer, cross_sca_aer, g_aer, Iinc, mu_z, dze):

  # Number of levels
  nlev = nlay + 1

  # Initialise Jdot - mean intensity estimator
  Jdot = np.zeros((nlay, nwl))

  # Initialise packet energy array
  e0dt = np.zeros(nwl)

  # Initialise surface albedo
  surf_alb = np.zeros(nwl)

  # Initialise asymmetry factor - need to alter this if non-isotropic
  g = np.zeros((nlay, nwl))

  # Initialise single scattering albedo
  alb = np.zeros((nlay, nwl))

  # Initialise gas phase scattering probability
  P_gas = np.zeros((nlay, nwl))

  # Initialise extinction and scattering opacity arrays
  sig_ext_tot = np.zeros((nlay, nwl))
  sig_abs_tot = np.zeros((nlay, nwl))
  sig_sca_tot = np.zeros((nlay, nwl))

  sig_abs_gas = np.zeros((nlay, nwl))
  sig_sca_gas = np.zeros((nlay, nwl))

  sig_abs_aer = np.zeros((nlay, nwl))
  sig_sca_aer = np.zeros((nlay, nwl))

  # Reconstruct altitude grid from dze
  z = np.zeros(nlev)
  for i in range(nlay):
    z[i+1] = z[i] + dze[i]

  # Set number of cross section species and rayleigh
  n_cross = len(ph_sp)
  n_ray = len(ray_sp)

  # Regular python loop
  #for l in range(nwl):
  # Parallel loop
  for l in prange(nwl):

    # Find the total absorption opacity from photocross sections
    for i in range(n_cross):
      sig_abs_gas[:,l] += nd_sp[:,all_sp.index(ph_sp[i])] * cross_sp[i,l] # Absorption opacity [cm-1]

    # Find the total rayleigh opacity from Rayleigh species cross sections
    for i in range(n_ray):
      sig_sca_gas[:,l] += nd_sp[:,all_sp.index(ray_sp[i])] * cross_sca_sp[i,l]   # Scattering opacity [cm-1]

    # Find absorption opacity of aerosol
    sig_abs_aer[:,l] = nd_aer[:] * cross_aer[:,l]

    # Find scattering opacity of aerosol
    sig_sca_aer[:,l] = nd_aer[:] * cross_sca_aer[:,l]

    # Total absorption and scattering is gas + aerosol
    sig_abs_tot[:,l] = sig_abs_gas[:,l] + sig_abs_aer[:,l]
    sig_sca_tot[:,l] = sig_sca_gas[:,l] + sig_sca_aer[:,l]

    # Total extinction = absorption + scattering
    sig_ext_tot[:,l] = sig_abs_tot[:,l] + sig_sca_tot[:,l]

    # Find total scattering albedo - MCRT only works well up to around 0.98 ssa, so limit to that
    alb[:,l] = np.minimum(sig_sca_tot[:,l]/sig_ext_tot[:,l],0.95)

    # Asymmetry factor for this wavelength 
    g[:,l] = g_aer[:,l]

    # Probability of gas phase scattering
    P_gas[:,l] = sig_sca_gas[:,l]/sig_sca_tot[:,l]

    # Energy carried by each packet for this wavelength
    e0dt[l] = (mu_z * Iinc[l])/float(Nph)

    # Initialise random seed for this wavelength 
    iseed = int(l + l**2 + l/2)
    seed(iseed)

    # Packet number loop
    for n in range(Nph):

      # Initialise packet variables (janky python way)
      flag = 0
      id_ph = l*Nph + n

      # Initialise photon packet with initial values
      ph = pac(flag, id_ph, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 2)

      # Place photon at top of atmosphere with negative zenith angle
      inc_stellar(ph, nlay, z, mu_z)

      # Scattering loop
      while(ph.flag == 0):

        # Sample a random optical depth of the photon
        ph.tau_p = -np.log(random())

        # Move the photon through the grid given by the sampled tau
        tauint_1D_pp(ph, nlay, z, sig_ext_tot, l, Jdot)

        # Check status of the photon after moving
        match ph.flag:
          case 1:
            # Photon exited top of atmosphere, exit loop
            break
          case -2: 
            # Photon hit surface surface, test surface albedo
            if (random() < surf_alb[l]):
              # Packet scattered off surface (Lambertian assumed)
              scatter_surf(ph, z)
            else:
              # Packet was absorbed by the surface
              break
          case 0:
            # Photon still in grid, test atmospheric albedo
            if (random() < alb[ph.zc,l]):

              # Check if gas or aerosol scattering and allocate iscat
              if (random() < P_gas[ph.zc,l]):
                ph.iscat = 2
              else:
                ph.iscat = 5

              # Packet get scattered into new direction
              scatter(ph, g[ph.zc,l])
            else:
              # Packet was absorbed, exit loop
              break
          case _:
            # Photon has an invalid flag, raise error
            print(str(ph.id) + ' has invalid flag: ' + str(ph.flag))

    ## Scale estimators to dze
    Jdot[:,l] = (e0dt[l]*Jdot[:,l])/dze[:]

  return Jdot
