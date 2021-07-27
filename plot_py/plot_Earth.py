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


# Setting the 2nd input argument as the filename of vulcan output   
vul_data = sys.argv[1]
# Setting the 3rd input argument as the species names to be plotted (separated by ,)
plot_spec = sys.argv[2]
# Setting the 4th input argument as the output eps filename        
plot_name = sys.argv[3]

plot_dir = '../' + vulcan_cfg.plot_dir

# taking user input species and splitting into separate strings and then converting the list to a tuple
plot_spec = tuple(plot_spec.split(','))
nspec = len(plot_spec)

colors = ['c','b','g','r','m','y','k','orange','pink','grey','darkred','darkblue','salmon','chocolate','steelblue','plum','hotpink']

tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$','N2O':'N$_2$O','O3':'O$_3$' }

# Data from MASSIE & HUNTEN (1981)
z_MH = np.arange(0,66,5)
z_MH[0] = 1
n2o_MH = 1e-9* np.array([316., 316., 316., 304, 251., 177., 112., 51., 16, 5.4, 2.1, 1.1, 0.62, 0.35])  
ch4_MH = 1e-6* np.array([1.54,1.54,1.54, 1.48,1.28,1.09,0.93,0.76,0.60,0.47,0.36,0.27,0.22,0.18])
o3_MH =  1e-6* np.array([0.3,0.32,0.05,0.33,2.2,5.7,6.5,6.2,4.7,2.9,1.6,0.86,0.56,0.41])

# By eye... Data from Sen et al (1998)
z_sen = np.arange(20,40)
no_sen = 1e-7/30.*np.array([0.1, 0.15, 0.2, 0.3,0.5,0.7,0.85,1., 1.5,2,3,5,7,8,9,9.1,9.3,9.5,9.8,9.9])

# Data from U.S. Standard Atmosphere 1976
z_us = np.arange(0,16.1,2)
z_us[0] = 1.
z_us = np.append(z_us,np.arange(18,33,2)) 
z_us = np.append(z_us,np.arange(35,81,5))

h2o_us = 1e-6/18*100*np.array([3700,2843,1268,554,216,43.2,11.3,3.3,3.3,3.3,4.5,7.2,11.6,18.6,18.2,17.6,16.8,15.4,12.2,11.1,7.6,4.9,3.8,2.3,1.4,1.1,0.6])

with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)

with open('../output/Earth-S-rtol05-st085.vul', 'rb') as handle:
  data2 = pickle.load(handle)

color_index = 0
vulcan_spec = data['variable']['species']
vulcan_spec2 = data2['variable']['species']

for sp in plot_spec:
    if color_index == len(colors): # when running out of colors
        colors.append(tuple(np.random.rand(3)))
    
    if sp in tex_labels: sp_lab = tex_labels[sp]
    else: sp_lab = sp  
    
    plt.plot(data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][:-1]/1.e5, color=colors[color_index], label=sp_lab)
    #plt.plot(data['variable']['y_ini'][:,vulcan_spec.index(sp)]/data['atm']['n_0'], data['atm']['pco']/1.e6, color=colors[color_index], ls=':', lw=1.5)
    
    plt.plot(data2['variable']['ymix'][:,vulcan_spec2.index(sp)], data2['atm']['zco'][:-1]/1.e5, color=colors[color_index], ls='--')
    
    if sp == 'H2O': plt.plot(data['atm']['sat_p']['H2O']/data['atm']['pco'], data['atm']['zco'][:-1]/1.e5, color=colors[color_index], ls='-.', label='H2O saturation')
        
        
        
    if sp == 'O3':
        plt.scatter(o3_MH, z_MH, marker='o', color=colors[color_index], facecolor= 'None', alpha=0.7)
        plt.errorbar(o3_MH, z_MH, xerr=np.vstack((o3_MH*0.9,o3_MH*9)), color=colors[color_index], linestyle='None', alpha=0.7)
    elif sp == 'N2O':
        plt.scatter(n2o_MH, z_MH, marker='o', color=colors[color_index], facecolor= 'None', alpha=0.7)
        plt.errorbar(n2o_MH, z_MH, xerr=np.vstack((n2o_MH*0.9,n2o_MH*9)), color=colors[color_index], linestyle='None', alpha=0.7)
    elif sp == 'CH4':
        plt.scatter(ch4_MH, z_MH, marker='o', color=colors[color_index], facecolor= 'None', alpha=0.7)
        plt.errorbar(ch4_MH, z_MH, xerr=np.vstack((ch4_MH*0.9,ch4_MH*9)), color=colors[color_index], linestyle='None', alpha=0.7)
    
    elif sp == 'NO':
        plt.scatter(no_sen, z_sen, marker='o', color=colors[color_index], facecolor= 'None', alpha=0.7)
        plt.errorbar(no_sen, z_sen, xerr=np.vstack((no_sen*0.9,no_sen*9)), color=colors[color_index], linestyle='None', alpha=0.7)
    
    elif sp == 'H2O':
        plt.scatter(h2o_us, z_us, marker='o', color=colors[color_index], facecolor= 'None', alpha=0.7)
        plt.errorbar(h2o_us, z_us, xerr=np.vstack((h2o_us*0.9,h2o_us*9)), color=colors[color_index], linestyle='None', alpha=0.7)
        
    
    color_index +=1


plt.gca().set_xscale('log')       
#plt.gca().set_yscale('log') 
#plt.gca().invert_yaxis() 
plt.xlim((1.E-14, 1.e-2))
plt.ylim((0, 80.))
#plt.ylim((1.E3,data['atm']['pco'][-1]/1e6))
plt.legend(frameon=0, prop={'size':13}, loc='best')
# handles, labels = plt.gca().get_legend_handles_labels()
# display = range(len(sp_list))
# #Create custom artists
# art0 = plt.Line2D((0,0),(0,0), ls='None')
# Artist1 = plt.Line2D(range(10),range(10), color='black')
# Artist2 = plt.Line2D((0,1),(0,0), color='black', ls='--',lw=1.5)
# plt.legend([Artist1,Artist2],['Equilibrium','Kinetics'], frameon=False, prop={'size':12}, loc='best')

plt.xlabel("Mixing Ratio", fontsize=16)
#plt.ylabel("Pressure (bar)")
plt.ylabel("Height (km)", fontsize=16)
plt.title('Earth (CIRA equator in January 1986)', fontsize=14)
plt.savefig(plot_dir + plot_name + '.png')
plt.savefig(plot_dir + plot_name + '.eps')
if vulcan_cfg.use_PIL == True:
    plot = Image.open(plot_dir + plot_name + '.png')
    plot.show()
else: plt.show()

# Plot concentration
plt.figure()
plt.plot(data['variable']['y'][:,vulcan_spec.index('H2O')], data['atm']['zco'][:-1]/1.e5, color='b', label='H2O')    
plt.plot(data['variable']['y'][:,vulcan_spec.index('OH')], data['atm']['zco'][:-1]/1.e5, color='g', label='OH')
plt.plot(data['variable']['y'][:,vulcan_spec.index('HO2')], data['atm']['zco'][:-1]/1.e5, color='r', label='HO2')

plt.gca().set_xscale('log')  
plt.legend(frameon=0, prop={'size':12}, loc='best')


plt.xlabel("Num density")
#plt.ylabel("Pressure (bar)")
plt.ylabel("Height (km)")
#plt.title('HD189733b')
plt.savefig(plot_dir + plot_name + '-concentration.png')
plt.savefig(plot_dir + plot_name + '-concentration.eps')
# if vulcan_cfg.use_PIL == True:
#     plot = Image.open(plot_dir + plot_name + '-concentration.png')
#     plot.show()
# else: plt.show()

plt.figure()
#plt.plot(data['atm']['sat_p']['H2O']/data['atm']['pco'], data['atm']['zco'][:-1]/1.e5, color='b')
    
plt.plot(data['variable']['y'][:,vulcan_spec.index('HO2')], data['atm']['zco'][:-1]/1.e5, color='b')

plt.xlabel("Saturation mixing ratio", fontsize=12)
plt.ylim((0, 80.))
#plt.xlim((1.E-14, 1.e-2))
plt.ylabel("Height (km)", fontsize=12)

plt.savefig(plot_dir + plot_name + '-saturation.png')
plt.savefig(plot_dir + plot_name + '-saturation.eps')
if vulcan_cfg.use_PIL == True:
    plot = Image.open(plot_dir + plot_name + '-concentration.png')
    plot.show()
else: plt.show()