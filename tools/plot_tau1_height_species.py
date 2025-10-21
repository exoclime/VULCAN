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
from scipy import interpolate
 
bin_down = True
bin_array = np.arange(10,500,1.)
       
plot_name = 'photosphere-Earth'
vul_data = '../output/Earth-rtol01.vul'
plot_dir = '../' + vulcan_cfg.plot_dir

plot_all_sp = False
plot_sp = [ 'O2', 'O3','N2', 'H2O', 'CO2', 'CH4'] # costomized what species to plot

color_index = 0
# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180),(255, 127, 14),(44, 160, 44),(214, 39, 40),(148, 103, 189),(140, 86, 75), (227, 119, 194),(127, 127, 127),(188, 189, 34),(23, 190, 207),\
(174, 199, 232),(255, 187, 120),(152, 223, 138),(255, 152, 150),(197, 176, 213),(196, 156, 148),(247, 182, 210),(199, 199, 199),(219, 219, 141),(158, 218, 229)] 

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

# tex labels for plotting
tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$','N2O':'N$_2$O','O3':'O$_3$' }

with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)

# photosphere
tau1 = 1.

# photosphere of each species
tau_sp, tau_scat = {}, {}

photo_sp = data['variable']['cross'].keys() # all species contribute to absorption
scat_sp = data['variable']['cross_scat'].keys() # all species contribute to scattering
# storing the main absorbers
main_sp, spectrum_sp = set(), [] # main_sp is a set


bins = data['variable']['bins']
dz = data['atm']['dz']
nz = len(dz)
vulcan_spec = data['variable']['species']

tau_sum = np.zeros( (nz+1, len(bins)) )



for sp in photo_sp: tau_sp[sp] = np.zeros( (nz+1, len(bins)) )
tau_scat = np.zeros( (nz+1, len(bins)) )
            
for j in range(nz-1,-1,-1):
    for sp in photo_sp: # scat_sp are not necessary photo_sp, e.g. He
        #absorption of speecies sp at level j
        if sp in data['variable']['cross_T'].keys():
            # cross_T is 2D
            cross_sp = data['variable']['cross_T'][sp][j]
        else: cross_sp = data['variable']['cross'][sp]
        
        tau_sp[sp][j] += data['variable']['y'][j,vulcan_spec.index(sp)] * dz[j] * cross_sp
        #absorption of all speecies sp at level j
        tau_sum[j] += tau_sp[sp][j]
        
        #adding the layer above for species j
        tau_sp[sp][j] += tau_sp[sp][j+1]
    
    for sp in scat_sp:
        tau_scat[j] += data['variable']['y'][j,vulcan_spec.index(sp)] * dz[j] * data['variable']['cross_scat'][sp]
        tau_sum[j] += data['variable']['y'][j,vulcan_spec.index(sp)] * dz[j] * data['variable']['cross_scat'][sp]
    
    # adding the layer above only at the end of species loop    
    tau_scat[j] += tau_scat[j+1] 
    tau_sum[j] += tau_sum[j+1] 

       
photosph = []
photosph_abs = []
photosph_sp = {}
photosph_scat = []

photosph_sum = []

for n in range(len(bins)):
    photosph.append( data['atm']['zco'][np.argmin(np.abs(data['variable']['tau'][:,n]-tau1)) ]/1.e5 )
    #photosph_abs.append( data['atm']['pco'][np.argmin(np.abs(data['variable']['tau_abs'][:,n]-tau1)) ]/1.e6 )
    photosph_scat.append( data['atm']['zco'][np.argmin(np.abs(tau_scat[:,n]-tau1)) ]/1.e5 )
    
    photosph_sum.append( data['atm']['zco'][np.argmin(np.abs(tau_sum[:,n]-tau1)) ]/1.e5 )
    
    
    
    
# finding out the main absorber
print ('zero tau:')
for n in range(len(bins)): #
    tau_bot = 0
    tau_one = 0 # where tau = 1
    max_sp = 'NA'
    for species in photo_sp:
        if tau_sp[species][0,n] > tau_bot: 
            tau_bot = tau_sp[species][0,n]
            max_sp = species
        
    main_sp.add(max_sp)
    spectrum_sp.append(max_sp)
    if tau_bot == 0: print (bins[n])

print ('Dominant species at the bottom:')
print (spectrum_sp[::10])
    
    
    
    
    
if 'NA' in main_sp:
    print ('Warning! No absorbers at all in some bins!')
    main_sp.remove('NA')

#main_sp.add('N2')
#main_sp.add('NH3')

# flag to plot all main_sp
if plot_all_sp == True: plot_sp = main_sp
    
for sp in plot_sp:
    photosph_sp[sp] = []
    for n in range(len(bins)):
        photosph_sp[sp].append( data['atm']['zco'][np.argmin(  np.abs(tau_sp[sp][:,n] -tau1)) ]/1.e5 )

# Plotting
plt.figure()
if bin_down == False: plt.plot(bins, photosph, c='k', lw=1.75, alpha=0.75, label=r'$\tau$=1')
else:
    tau_inter = interpolate.interp1d(bins, photosph)
    tau_tot = np.empty(len(bin_array))
    for _,ld in enumerate(bin_array): tau_tot[_] = float(tau_inter(ld))
    plt.plot(bin_array, tau_tot, c='k', lw=1.5, alpha=0.75, label=r'$\tau$=1')


for sp in plot_sp:
    if sp in tex_labels: tex_lab = tex_labels[sp]
    else: tex_lab = sp
    
    tau_sp = np.array(photosph_sp[sp])
    # above 1000 bar
    if bin_down == False: plt.plot(bins[tau_sp<1000.], tau_sp[tau_sp<1000.], color=tableau20[color_index], label = tex_lab, alpha=0.75 )
    else:
        tau_inter = interpolate.interp1d(bins, tau_sp)
        tau_sp = np.empty(len(bin_array))
        for _,ld in enumerate(bin_array): tau_sp[_] = float(tau_inter(ld))
        
        plt.plot(bin_array[tau_sp<1000.], tau_sp[tau_sp<1000.], color=tableau20[color_index], label = tex_lab, alpha=0.75 )
    # if sp == 'N2':
    #     plt.plot(bins, photosph_sp[sp], color=tableau20[color_index], label = tex_lab, alpha=0.6, lw=5 )
    
    color_index += 1

plt.plot(bins, photosph_scat, color=tableau20[color_index], label = 'Rayleigh', alpha=0.6)

#plt.plot(bins, photosph_sum, c='plum', lw=3, alpha=0.6)
#plt.plot(data['variable']['bins'], photosph_abs, c='red')
           
#plt.gca().set_yscale('log') 
#plt.gca().invert_yaxis() 
plt.ylim((0,120))
plt.xlim((20,360))
plt.legend(frameon=0, prop={'size':12}, loc='best')
plt.xlabel("wavelength (nm)")
plt.ylabel("Height (km)")
#plt.title('Photosphere')
plt.savefig(plot_dir + plot_name + '.png')
plt.savefig(plot_dir + plot_name + '.pdf')
#if vulcan_cfg.use_PIL == True:
plot = Image.open(plot_dir + plot_name + '.png')
plot.show()
#else: plt.show()

