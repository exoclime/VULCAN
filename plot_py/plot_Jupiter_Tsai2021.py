'''
This script reads VULCAN output (.vul) files using pickle and plot the species volumn mixing ratios as a function of pressure, with the initial abundances (typically equilibrium) shown in dashed lines.
Plots are saved in the folder assigned in vulcan_cfg.py, with the default plot_dir = 'plot/'.
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
from phy_const import kb, Navo

# swtich for plot
if '-h' in sys.argv: use_height = True
else: use_height = False 

# Setting the 3rd input argument as the species names to be plotted (separated by ,)
plot_spec = sys.argv[1]
# Setting the 4th input argument as the output eps filename        
plot_name = sys.argv[2]

#vul_data = '../output/Jupiter_rtol025.vul'

vul_data = '../output/Jupiter_rtol005.vul'




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
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O', 'H2O_l_s':'H2O(s)', 'NH3_l_s':'NH3(s)'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$', 'H2O_l_s':'H$_2$O ice','NH3_l_s':'NH$_3$ ice' ,'C6H6':'C$_6$H$_6$','C3H3':'C$_3$H$_3$','C3H2':'C$_3$H$_2$','C4H5':'C$_4$H$_5$'}



with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)

color_index = 0
vulcan_spec = data['variable']['species']

m_h2o = 18./Navo
m_nh3 = 17./Navo 
n_mol_h2o = 4./3*np.pi*data['atm']['r_p']['H2O_l_s']**3 *data['atm']['rho_p']['H2O_l_s'] /m_h2o
n_mol_nh3 = 4./3*np.pi*data['atm']['r_p']['NH3_l_s']**3 *data['atm']['rho_p']['NH3_l_s'] /m_nh3

p_sp = {} # to store matplot obj
fig, ax = plt.subplots()


# in the form of N pts... (xerr1_low,xerr2_low,xerr3_low,...), (xerr1_high, x2err_high,...), (yerr1_low,yerr2_low), (yerr1_high, y2err_high)
# note that y is inverted!
# Jupiter_obs = {'C2H2': [ (1.4E-6,0.25/1e3), (3E-8,30./1e3),(4E-6,1e-5),  (1E-6,1.5e-8,3E-6), (1E-6,2e-8,5E-6), (0,0.025,0),(0,0.2,0)],  'C2H6': [ (9E-6,1/1e3),(2.6E-6,10/1e3),(2E-6,1E-6), (1E-5,3E-6), (0,0),(0,0) ],\
# 'C2H4': [(4E-10,7/1e3),(2E-6,7E-3/1e3), (3E-10,1.E-6),(5E-10,5E-7),(0.005,0),(0.09,0)], 'CH4':[(2E-4,1.3e-7),(3E-5,4e-6),(1E-4,2E-5),(2E-4,5E-5),(0.7E-7,0),(2E-7,0)], 'CO':[(3e-9,50/1e3),(3e-9,50/1e3),(2E-9,2E-9),(1e-9,1e-9),(0,0),(0,0)]          }


#Observations for Jupiter: name, [pressures], [mixing ratios], [pressure-error-low], [pressure-error-high], [mix-error-low], [mix-error-high]
# obs = {'CH4-drossart':[1.03e-05, 9.24e-06, 8.27e-06, 7.44e-06, 6.68e-06, 5.98e-06, 5.37e-06, 4.82e-06, 4.32e-06,\
#  3.87e-06, 3.48e-06, 3.12e-06, 2.8e-06, 2.51e-06, 2.26e-06, 2.03e-06, 1.82e-06, 1.63e-06, 1.46e-06, 1.31e-06,\
#   1.18e-06, 1.06e-06, 9.48e-07, 8.49e-07, 7.63e-07, 6.85e-07, 6.13e-07, 5.5e-07, 4.94e-07, 4.44e-07, 3.98e-07,\
#    3.57e-07, 3.2e-07, 2.87e-07, 2.58e-07, 2.32e-07, 2.08e-07, 1.86e-07, 1.67e-07, 1.5e-07, 1.34e-07, 1.2e-07, 1.08e-07, 9.7e-08],\
#    [0.00108781, 0.0010624699999999998, 0.0010371299999999998, 0.00101179, 0.00098645, 0.00096111, 0.0009339600000000001, 0.00090862,\
#     0.00087966, 0.0008525099999999999, 0.00082536, 0.00079821, 0.00076925, 0.0007421, 0.00071314, 0.00068599, 0.00065703, 0.00062988,\
#      0.00060092, 0.00057196, 0.00054481, 0.00051766, 0.00049051, 0.00046336, 0.00043801999999999996, 0.00041268, 0.00038733999999999996,\
#       0.000362, 0.00033847, 0.00031675, 0.00029322, 0.0002715, 0.00025159, 0.00023168, 0.00021177000000000001, 0.00019548, 0.000177742,0.00016109, 0.00014570500000000002, 0.000131587, 0.000118012, 0.000105342, 9.375799999999999e-05, 8.3441e-05]],\
#        ['CH4-Festou',[5e-6],[2.5e-5],[3e-6],[7e-6],[1.2e-5],[5e-5]],['CH4-Yelle',[2e-7],[1.5e-4],[1e-7],[4e-7],[1e-4],[2e-4]],\
#        ['C2H2-Fouchet',[4e-3],[3.62e-8],[2e-3],[8e-3],[3.28e-8],[4.22e-8]],['C2H2-Moses',[2.5e-4,2.0e-3],[1.4e-6,1.5e-7],\
#        [1.25e-4,1e-3],[5e-4,4e-3],[6.0e-7,1.1e-7],[2.2e-6,1.9e-7]],['C2H2-Kim',[1e-4],[1e-6],[1e-5],[1e-3],[1e-7],[1e-5]],\
#        ['C2H4-Romani',[5e-6,2.2e-6],[5.5e-7,1.1e-6],[2.5e-6,1.1e-6],[1e-5,4.4e-6],[2.75e-7,5.5e-7],[1.1e-6,2.2e-6]],\
#        ['C2H4-Bezard',[1e-3],[6e-10],[5e-4],[2e-3],[3.9e-10],[1.02e-9]],\
#        ['C2H6-Fouchet',[1e-3,1e-2],[8.62e-6,2.24e-6],[5e-4,5e-3],[2e-3,2e-2],[6.896e-6,1.724e-6],[1.034e-5,2.672e-6]],\
#        ['C2H6-Moses',[3.5e-3,7e-3],[4e-6,2.7e-6],[1.75e-3,3.5e-3],[7e-3,1.4e-2],[3e-6,2e-6],[5e-6,3.4e-6]],\
#        ['C2H6-Yelle',[5e-3],[4.65e-6],[4e-4],[1e-2],[2.8e-6],[6.5e-6]],['C2H6-Kim',[1e-5],[5e-6],[1e-6],[1e-4],[2.5e-6],[1e-5]],\
#        ['C4H2-Fouchet',[5.25e-4],[1.26e-10],[2.5e-4],[1e-3]],['C4H2-Moses',[1e-3],[1.8e-10],[5e-4],[2e-3]]]}


# Gladstone+1996
obsYung = [['C2H2',[1.3e-2,1.0e-2,1.0e-2,1.5e-3,1.3e-2,1.3e-2,1.3e-2],[1.5e-8,1.0e-7,3.0e-8,1.0e-7,2.6e-8,2.3e-8,9.0e-8],[1e-3,5e-3,5e-3,1e-4,1e-3,6e-3,6e-3],[6e-2,1.5e-2,1.5e-2,4e-3,1e-1,1e-1,1e-1],[7e-9,9e-8,2e-8,7e-8,1.1e-8,1e-8,7e-8],[2.3e-8,1.1e-7,4e-8,1.3e-7,4.1e-8,3.6e-8,1.1e-7]],['C2H4',[5e-3,1.3e-2],[5e-10,7e-9],[7e-4,6e-3],[3e-2,1e-1],[2e-10,4e-9],[8e-10,1e-8]],['C2H6',[2.1e-2,5e-6,4e-3,4e-3,1.5e-3,5e-3,1.8e-2,1.7e-2,1.7e-2,5e-3],[1.9e-6,2.5e-6,6e-6,2.5e-6,5.5e-6,3.75e-6,1.1e-6,2.6e-6,5e-6,4e-6],[3e-3,2.5e-6,1e-3,1e-3,1e-4,7e-4,2e-3,6e-3,6e-3,7e-4],[6e-2,1e-5,2e-2,2e-2,4e-3,3e-2,1e-1,1e-1,1e-1,3e-2],[8e-7,1e-6,1e-6,1e-6,4e-6,1.5e-6,5e-7,1.1e-6,4e-6,1.6e-6],[3e-6,5e-6,1.2e-5,4e-6,7e-6,6e-6,1.7e-6,4.1e-6,6e-6,6.4e-6]]]


CH4_drossart = {'Pbar': [1.03e-05, 9.24e-06, 8.27e-06, 7.44e-06, 6.68e-06, 5.98e-06, 5.37e-06, 4.82e-06, 4.32e-06,3.87e-06, 3.48e-06, 3.12e-06, 2.8e-06, 2.51e-06, 2.26e-06, 2.03e-06, 1.82e-06, 1.63e-06, 1.46e-06, 1.31e-06,\
  1.18e-06, 1.06e-06, 9.48e-07, 8.49e-07, 7.63e-07, 6.85e-07, 6.13e-07, 5.5e-07, 4.94e-07, 4.44e-07, 3.98e-07,3.57e-07, 3.2e-07, 2.87e-07, 2.58e-07, 2.32e-07, 2.08e-07, 1.86e-07, 1.67e-07, 1.5e-07, 1.34e-07, 1.2e-07, 1.08e-07, 9.7e-08],\
"mix": [0.00108781, 0.0010624699999999998, 0.0010371299999999998, 0.00101179, 0.00098645, 0.00096111, 0.0009339600000000001, 0.00090862,0.00087966, 0.0008525099999999999, 0.00082536, 0.00079821, 0.00076925, 0.0007421, 0.00071314, 0.00068599, 0.00065703, 0.00062988,\
     0.00060092, 0.00057196, 0.00054481, 0.00051766, 0.00049051, 0.00046336, 0.00043801999999999996, 0.00041268, 0.00038733999999999996,0.000362, 0.00033847, 0.00031675, 0.00029322, 0.0002715, 0.00025159, 0.00023168, 0.00021177000000000001, 0.00019548, 0.000177742,\
     0.00016109, 0.00014570500000000002, 0.000131587, 0.000118012, 0.000105342, 9.375799999999999e-05, 8.3441e-05]} 

# Be'zard 2002
CO_bevard = [1e-9, 6.]

# Moses+2005 given by Paul
C2H2_Moses = {'Pbar':[2.5e-4,2.0e-3], 'mix':[1.4e-6,1.5e-7], 'PerrMin':np.array([1.25e-4,1e-3]), 'PerrMax':np.array([5e-4,4e-3]), 'MixMin':np.array([6.0e-7,1.1e-7]), 'MixMax':np.array([2.2e-6,1.9e-7]) }
# Kim+2010 given by Paul
C2H2_Kim = {'Pbar':[1e-4], 'mix':[1e-6], 'PerrMin':np.array([1e-5]),'PerrMax':np.array([1e-3]),'MixMin':np.array([1e-7]),'MixMax':np.array([1e-5])}
# Gladstone+ 1996 given by Paul
C2H2_Gladstone = {'Pbar':[1.3e-2,1.0e-2,1.0e-2,1.5e-3,1.3e-2,1.3e-2,1.3e-2],'mix':[1.5e-8,1.0e-7,3.0e-8,1.0e-7,2.6e-8,2.3e-8,9.0e-8], 'PerrMin':np.array([1e-3,5e-3,5e-3,1e-4,1e-3,6e-3,6e-3]),'PerrMax':np.array([6e-2,1.5e-2,1.5e-2,4e-3,1e-1,1e-1,1e-1]), 'MixMin':np.array([7e-9,9e-8,2e-8,7e-8,1.1e-8,1e-8,7e-8]), 'MixMax':np.array([2.3e-8,1.1e-7,4e-8,1.3e-7,4.1e-8,3.6e-8,1.1e-7])}

# Romani+2008 and Be'zard 2001a given by Paul
C2H4_Romani = {'Pbar':[5e-6,2.2e-6],'mix':[5.5e-7,1.1e-6], 'PerrMin':np.array([2.5e-6,1.1e-6]),'PerrMax':np.array([1e-5,4.4e-6]),'MixMin':np.array([2.75e-7,5.5e-7]),'MixMax':np.array([1.1e-6,2.2e-6])}
C2H4_Bezard = {'Pbar':[1e-3],'mix':[6e-10], 'PerrMin':np.array([5e-4]),'PerrMax':np.array([2e-3]), 'MixMin':np.array([3.9e-10]), 'MixMax':np.array([1.02e-9])}

# Yelle   given by Paul
C2H6_Gladstone = {'Pbar':[2.1e-2,5e-6,4e-3,4e-3,1.5e-3,5e-3,1.8e-2,1.7e-2,1.7e-2,5e-3],'mix':[1.9e-6,2.5e-6,6e-6,2.5e-6,5.5e-6,3.75e-6,1.1e-6,2.6e-6,5e-6,4e-6], 'PerrMin':np.array([3e-3,2.5e-6,1e-3,1e-3,1e-4,7e-4,2e-3,6e-3,6e-3,7e-4]),'PerrMax':np.array([6e-2,1e-5,2e-2,2e-2,4e-3,3e-2,1e-1,1e-1,1e-1,3e-2]), 'MixMin':np.array([8e-7,1e-6,1e-6,1e-6,4e-6,1.5e-6,5e-7,1.1e-6,4e-6,1.6e-6]),'MixMax':np.array([3e-6,5e-6,1.2e-5,4e-6,7e-6,6e-6,1.7e-6,4.1e-6,6e-6,6.4e-6])}
C2H6_Fouchet = {'Pbar':[1e-3,1e-2],'mix':[8.62e-6,2.24e-6],'PerrMin':np.array([5e-4,5e-3]),'PerrMax':np.array([2e-3,2e-2]),'MixMin':np.array([6.896e-6,1.724e-6]),'MixMax':np.array([1.034e-5,2.672e-6])}  



for color_index,sp in enumerate(plot_spec):
    if color_index == len(tableau20): # when running out of colors
        tableau20.append(tuple(np.random.rand(3)))
    
    if sp in tex_labels: sp_lab = tex_labels[sp]
    else: sp_lab = sp  
    
    if use_height == False:
        
        if sp == 'H2O_l_s': 
            #particle_mass = 4/3*np.pi*data['atm']['r_p']['H2O_l_s']**3 *data['atm']['rho_p']['H2O_l_s']
            p_sp[sp], = plt.plot(data['variable']['y'][:,vulcan_spec.index(sp)]/Navo*18., data['atm']['pco']/1.e6, color=tableau20[color_index], label=sp_lab, alpha=0.9, ls='-.')

        elif sp == 'NH3_l_s': 
            #particle_mass = 4/3*np.pi*data['atm']['r_p']['NH3_l_s']**3 *data['atm']['rho_p']['NH3_l_s']
            p_sp[sp], = plt.plot(data['variable']['y'][:,vulcan_spec.index(sp)]/Navo*17., data['atm']['pco']/1.e6, color=tableau20[color_index], label=sp_lab, alpha=0.9, ls='-.')
        
        else:
        
            p_sp[sp], = plt.plot(data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['pco']/1.e6, color=tableau20[color_index], lw=1.5) # , label=sp_lab
            #plt.plot(data['variable']['y_ini'][:,vulcan_spec.index(sp)]/data['atm']['n_0'], data['atm']['pco']/1.e6, color=tableau20[color_index], ls=':', lw=1.5)
            
            # Obs. data
            if sp == 'CH4': plt.scatter(CH4_drossart['mix'][::4],CH4_drossart['Pbar'][::4], edgecolors=tableau20[color_index],facecolors= 'None', marker='s', alpha=0.75, label='Drossart et al.(1999)' )
            elif sp == 'CO': 
                plt.scatter(CO_bevard[0], CO_bevard[1], edgecolors=tableau20[color_index],facecolors= 'None', marker='s', alpha=0.75 ,label='Bézard et al.(2002)') 
                #plt.errorbar(CO_bevard[0], CO_bevard[1], xerr=0.2e-9, color=tableau20[color_index], linestyle='None', alpha=0.75) # ms=80
            elif sp == 'C2H2':     
                #plt.scatter(C2H2_Moses['mix'], C2H2_Moses['Pbar'], edgecolors=tableau20[color_index],facecolors= 'None', marker='s', alpha=0.75 ) 
                plt.errorbar(C2H2_Moses['mix'], C2H2_Moses['Pbar'], xerr=[C2H2_Moses['mix']-C2H2_Moses['MixMin'],C2H2_Moses['MixMax']-C2H2_Moses['mix']], yerr=[C2H2_Moses['Pbar']-C2H2_Moses['PerrMin'],C2H2_Moses['PerrMax']-C2H2_Moses['Pbar']], color=tableau20[color_index], fmt='o', alpha=0.75,capsize=2, mfc='none', label='Moses et al.(2005)' ) # mfc='none' for open circles
                plt.errorbar(C2H2_Gladstone['mix'], C2H2_Gladstone['Pbar'], xerr=[C2H2_Gladstone['mix']-C2H2_Gladstone['MixMin'],C2H2_Gladstone['MixMax']-C2H2_Gladstone['mix']], yerr=[C2H2_Gladstone['Pbar']-C2H2_Gladstone['PerrMin'],C2H2_Gladstone['PerrMax']-C2H2_Gladstone['Pbar']], color=tableau20[color_index], fmt='.', alpha=0.75,capsize=2, mfc='none', label='Gladstone et al.(1999)')
                #plt.errorbar(C2H2_Kim['mix'], C2H2_Kim['Pbar'], xerr=[C2H2_Kim['mix']-C2H2_Kim['MixMin'],C2H2_Kim['MixMax']-C2H2_Kim['mix']], yerr=[C2H2_Kim['Pbar']-C2H2_Kim['PerrMin'],C2H2_Kim['PerrMax']-C2H2_Kim['Pbar']], color=tableau20[color_index], fmt='s', alpha=0.75,capsize=2, mfc='none')
            
            elif sp == 'C2H4':    
                plt.errorbar(C2H4_Romani['mix'], C2H4_Romani['Pbar'], xerr=[C2H4_Romani['mix']-C2H4_Romani['MixMin'],C2H4_Romani['MixMax']-C2H4_Romani['mix']], yerr=[C2H4_Romani['Pbar']-C2H4_Romani['PerrMin'],C2H4_Romani['PerrMax']-C2H4_Romani['Pbar']], color=tableau20[color_index], fmt='.', alpha=0.75,capsize=2, mfc='none', label='Romani et al.(2008)')
                plt.errorbar(C2H4_Bezard['mix'], C2H4_Bezard['Pbar'], xerr=[C2H4_Bezard['mix']-C2H4_Bezard['MixMin'],C2H4_Bezard['MixMax']-C2H4_Bezard['mix']], yerr=[C2H4_Bezard['Pbar']-C2H4_Bezard['PerrMin'],C2H4_Bezard['PerrMax']-C2H4_Bezard['Pbar']], color=tableau20[color_index], fmt='o', alpha=0.75,capsize=2, mfc='none', label='Bézard et al.(2002)')
            
            elif sp == 'C2H6':    
                plt.errorbar(C2H6_Gladstone['mix'], C2H6_Gladstone['Pbar'], xerr=[C2H6_Gladstone['mix']-C2H6_Gladstone['MixMin'],C2H6_Gladstone['MixMax']-C2H6_Gladstone['mix']], yerr=[C2H6_Gladstone['Pbar']-C2H6_Gladstone['PerrMin'],C2H6_Gladstone['PerrMax']-C2H6_Gladstone['Pbar']], color=tableau20[color_index], fmt='.', alpha=0.75,capsize=2, mfc='none', label='Gladstone et al.(1999)')
                plt.errorbar(C2H6_Fouchet['mix'], C2H6_Fouchet['Pbar'], xerr=[C2H6_Fouchet['mix']-C2H6_Fouchet['MixMin'],C2H6_Fouchet['MixMax']-C2H6_Fouchet['mix']], yerr=[C2H6_Fouchet['Pbar']-C2H6_Fouchet['PerrMin'],C2H6_Fouchet['PerrMax']-C2H6_Fouchet['Pbar']], color=tableau20[color_index], fmt='o', alpha=0.75,capsize=2, mfc='none', label='Fouchet et al.(2000)')
                

    else: 
        plt.plot(data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['zco'][1:]/1.e5, color=tableau20[color_index], label=sp_lab, lw=1.5)
    # # plotting the initial (equilibrium) abundances
        

                
                

if use_height == False:
    plt.gca().set_yscale('log') 
    plt.gca().invert_yaxis() 
    plt.ylim((data['atm']['pco'][0]/1e6,data['atm']['pco'][-1]/1e6))
    plt.ylabel("Pressure (bar)",fontsize=12)
else:
    plt.ylim((data['atm']['zmco'][0]/1e5,data['atm']['zmco'][-1]/1e5)) 


plt.xlabel("Mixing Ratio / Cloud Density (g/cm$^3$)",fontsize=12)  
plt.xlabel("Mixing Ratio",fontsize=12)     

# plt.plot(data['atm']['sat_p']['H2O']/data['atm']['pco'], data['atm']['pco']/1.e6, color='k', label='H$_2$O\nsaturation', alpha=0.7, ls=':')
# plt.plot(data['atm']['sat_p']['NH3']/data['atm']['pco'], data['atm']['pco']/1.e6, color='k', label='NH$_3$\nsaturation', alpha=0.7, ls='--')
   
plt.gca().set_xscale('log')       
plt.xlim((1.E-16, 1.))
#plt.xlim((5.E-11, 1.e-2))
#
# plt.xlim((1.E-30, 1e-2))

#plt.tick_params(axis='x',which='majpr',bottom=True)
plt.minorticks_on
ax.minorticks_on()

leg1 = ax.legend([p_sp[sp] for sp in plot_spec], [tex_labels[sp] for sp in plot_spec],frameon=0, prop={'size':11}, loc=1)
# Add second legend for the species
leg2 = ax.legend(frameon=0, prop={'size':9.5}, loc=4)
# Manually add the first legend back! WTF
ax.add_artist(leg1)

# plt.legend(frameon=0, prop={'size':12}, loc=1)
# handles, labels = plt.gca().get_legend_handles_labels()
# display = range(len(sp_list))
# #Create custom artists
# art0 = plt.Line2D((0,0),(0,0), ls='None')
# Artist1 = plt.Line2D(range(10),range(10), color='black')
# Artist2 = plt.Line2D((0,1),(0,0), color='black', ls='--',lw=1.5)
# plt.legend([Artist1,Artist2],['Equilibrium','Kinetics'], frameon=False, prop={'size':12}, loc='best')

plt.savefig(plot_dir + plot_name + '.png')
#plt.savefig(plot_dir + plot_name + '.pdf')
if vulcan_cfg.use_PIL == True:
    plot = Image.open(plot_dir + plot_name + '.png')
    plot.show()
else: plt.show()

