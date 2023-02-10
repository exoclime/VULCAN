'''Generate a petitRadTrans spectrum from a VULCAN run.

SYNTAX:
python prt_spectrum.py  ../output_path/results.vul  molecules_in_spectrum  plotname

EXAMPLE:
python prt_spectrum.py   ../wasp39b_helios-sIC_SNCHO_100x/wasp39b_helios.vul   H2O,OH,CO,CO2,CH4,CH3,HCN,C2H2,C2H4,CN,CH,SO2,SH,H2S,Na,K   wasp39_IC_SNCHO_100x_all


NOTES:

  Not all molecules listed in the above example are included in pRT's
  default opacity tarball. You may need to download them manually.

  pRT sometimes uses rather specific names for its opacity folders
  (e.g., 'CO_all_iso_HITEMP'). The script below tries to be slightly
  clever about this, but you may need to edit the dictionary
  'prt_longnames' below if there are other, nonstandard names.



HISTORY:
--------
2023-02-08 10:44 IJMC: Created by adapting 'plot_vulcan.py' script.

'''

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

#IJMC additions:
import os
import shutil
import pandas as pd
from petitRADTRANS import Radtrans
from petitRADTRANS.retrieval import util as prt_util

# swtich for plot
if '-h' in sys.argv: use_height = True
else: use_height = False 

# Setting the 2nd input argument as the filename of vulcan output   
vul_data = sys.argv[1]
# Setting the 3rd input argument as the species names to be plotted (separated by ,)
plot_spec = sys.argv[2]
# Setting the 4th input argument as the output eps filename        
if len(sys.argv)<4:
    saveplot=False
else:
    saveplot=True
    plot_name = sys.argv[3]
    plot_dir = '../' + vulcan_cfg.plot_dir
    # Checking if the plot folder exsists
    if not os.path.exists(plot_dir):
        print ('The plotting directory assigned in vulcan_cfg.py does not exist.')
        print( 'Directory ' , plot_dir,  " created.")
        os.mkdir(plot_dir)

# taking user input species and splitting into separate strings and then converting the list to a tuple
plot_spec = tuple(plot_spec.split(','))
nspec = len(plot_spec)

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180),(255, 127, 14),(44, 160, 44),(214, 39, 40),(148, 103, 189),(140, 86, 75), (227, 119, 194),(127, 127, 127),(188, 189, 34),(23, 190, 207),\
(174, 199, 232),(255, 187, 120),(152, 223, 138),(255, 152, 150),(197, 176, 213),(196, 156, 148),(247, 182, 210),(199, 199, 199),(219, 219, 141),(158, 218, 229)] 
# 


# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

# tex labels for plotting
tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$' }



with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)


  
########################################
# Shenanigans to import *this* CFG file (saved as cfg_xxxxx.txt):
########################################
_load_directory = os.path.split(vul_data)[0]
_filenames = os.listdir(_load_directory)
_cfg_file = _filenames[np.nonzero([fn[0:4]=='cfg_' for fn in _filenames])[0][0]]


shutil.copyfile(_load_directory+'/'+_cfg_file, 'cfg_file.py')
import cfg_file as cfg
os.remove('cfg_file.py')

########################################
# Now prepare petitRadTrans!
########################################
wmin, wmax = 1, 25
P0 = 0.01
atmo = 'H2'
xscale = 'log'
oneatatime = True

rsun_cgs = 69550800000.0
c_cgs = 29979245800.

if atmo=='H2':
    rayleigh_species = ['H2', 'He']
    continuum_opacities = ['H2-H2', 'H2-He']
else:
    print('Only H2-dominated atmospheres have been implemented for this script. But please consider submitted an update to GitHub!')

# Convert plot_spec molecule names to pRT-style opacity names:

# CO_all_iso_Chubb      CO_all_iso_HITEMP    
# H2O_Exomol H2O_HITEMP
# K_allard   K_burrows   K_lor_cut 
# Na_allard  Na_burrows  Na_lor_cut 
# TiO_all_Exomol          TiO_all_Plez 
# VO_Plez      
pRT_longnames = dict(CO='CO_all_iso_HITEMP', H2O='H2O_HITEMP', K='K_allard', Na='Na_allard', TiO='TiO_all_Exomol', VO='VO_Plez')
pRT_species = ['Al','Al+','AlH','AlO','C2H2','C2H4','CH4','CO','CO2','Ca','Ca+','CaH','CrH','Fe','Fe+','FeH','H2O','H2S','HCN','K','Li','Mg','Mg+','MgH','MgO','NH3','Na','NaH','O','O+','O2','O3','OH','PH3','SH','Si','Si+','SiO','SiO2','Ti','Ti+','TiO','V','V+','VO']
pRT_species += ['SO2', 'SO3', 'CH3', 'C2', 'CH', 'CN']



line_species_prt =  []
line_species_vul = []
for species in plot_spec:
    if species in pRT_longnames.keys():
        line_species_prt.append(pRT_longnames[species])
        line_species_vul.append(species)
    elif species in pRT_species:
        line_species_prt.append(species)
        line_species_vul.append(species)
    else:
        print('<WARNING>')
        print('Input species %s is probably not among petitRadTrans opacity files.' % species )
        print('This may not work out.')
        print('<\WARNING>')
        line_species_vul.append(species)
        line_species_prt.append(species)

        
atmosphere = Radtrans(line_species = line_species_prt,
      rayleigh_species = rayleigh_species,
      continuum_opacities = continuum_opacities,
      wlen_bords_micron = [wmin, wmax])

pressures = data['atm']['pco']/1e6
psort = np.argsort(pressures)
pressures = pressures[psort]
temperature = data['atm']['Tco'][psort]
atmosphere.setup_opa_structure(pressures)

vulcan_species = data['variable']['species']
vulcan_vmr = data['variable']['ymix'] # Volume Mixing Ratios

# Convert VULCAN's volume mixing ratios into pRT's MASS-mixing-ratios:
compositions = pd.read_csv('../thermo/all_compose.txt', delim_whitespace=True)
c_inds = np.zeros(len(vulcan_species), dtype=int)
for ii,species in enumerate(vulcan_species):
    vec = species==compositions.species.values
    if vec.any():
        c_inds[ii] = vec.nonzero()[0][0]
    else:
        c_inds[ii] = -1

mmw = compositions.mass.values[c_inds]
mmw_profile = (vulcan_vmr * mmw).sum(1) / vulcan_vmr.sum(1)
vulcan_mmr = vulcan_vmr * mmw / mmw_profile.reshape(mmw_profile.size, 1)

mass_fractions = {}
#import pdb
#pdb.set_trace()
for sp,sv in zip(line_species_prt,line_species_vul):
    if sv in vulcan_species:
        ind = vulcan_species.index(sv)
        mass_fractions[sp] = vulcan_mmr[:,ind][psort] # sort by increasing pressure
    else:
        print('<WARNING>')
        print('Species %s was input but has no abundances in the VULCAN network.' % sv )
        print('Setting to zero abundance in pRT model... what a waste!')
        print('</WARNING>')
        mass_fractions[sp] = np.zeros(pressures.shape) # sort by increasing pressure

if atmo=='H2':
    h2ind = vulcan_species.index('H2')
    heind = vulcan_species.index('He')
    mass_fractions['H2'] = vulcan_mmr[:,h2ind][psort] # sort by increasing pressure    
    mass_fractions['He'] = vulcan_mmr[:,heind][psort] # sort by increasing pressure    
    

    
atmosphere.calc_transm(temperature, mass_fractions, cfg.gs, mmw_profile, R_pl=cfg.Rp, P0_bar=P0)

wlen_micron = c_cgs/atmosphere.freq/1e-4
rprs2 = (atmosphere.transm_rad / (cfg.r_star * rsun_cgs))**2
if oneatatime:
    rprss = np.zeros((len(line_species_prt), rprs2.size), dtype=float)
    for ii,spec in enumerate(line_species_prt):
        mf = dict(H2=mass_fractions['H2'], He=mass_fractions['He'])
        for jj,spec in enumerate(line_species_prt):
            if ii==jj:
                mf[spec] = mass_fractions[spec]
            else:
                mf[spec] = np.zeros(mass_fractions[spec].size)
        atmosphere.calc_transm(temperature, mf, cfg.gs, mmw_profile, R_pl=cfg.Rp, P0_bar=P0)
        rprss[ii] = (atmosphere.transm_rad / (cfg.r_star * rsun_cgs))**2
    
########################################
########################################


output = pd.DataFrame(dict(wave=wlen_micron))
output['rprs2'] = rprs2
for jj,spec in enumerate(line_species_prt):
    output['rprs2_only_'+spec] = rprss[jj]

output.to_csv(vul_data.replace('.vul', '_transmission.csv'), index=False)


if saveplot:
    plt.figure()
    plt.plot(wlen_micron, rprs2*1e6)
    plt.gca().set_xscale(xscale)
    plt.ylabel('Transit Depth [ppm]', fontsize=14)
    plt.xlabel('Wavelength [micron]', fontsize=14)
    plt.minorticks_on()
    plt.title(cfg.output_dir + cfg.out_name, fontsize=12)
    plt.tight_layout()
    plt.xlim(wlen_micron.min(), wlen_micron.max())
    
    plt.savefig(plot_dir + plot_name + '_transmission.png')
    plt.savefig(plot_dir + plot_name + '_transmission.eps')

    
    plt.figure()
    plt.plot(wlen_micron, rprs2*1e6, '-k', label='All species')
    for jj,spec in enumerate(line_species_prt):
        plt.plot(wlen_micron, rprss[jj]*1e6, label='only %s' % spec  )
    plt.gca().set_xscale(xscale)
    plt.ylabel('Transit Depth [ppm]', fontsize=14)
    plt.xlabel('Wavelength [micron]', fontsize=14)
    plt.legend(loc=4)
    plt.minorticks_on()
    plt.title(cfg.output_dir + cfg.out_name, fontsize=12)
    plt.tight_layout()
    plt.xlim(wlen_micron.min(), wlen_micron.max())

    plt.savefig(plot_dir + plot_name + '_transmission_1by1.png')
    plt.savefig(plot_dir + plot_name + '_transmission_1by1.eps')
    
    

########################################
########################################
