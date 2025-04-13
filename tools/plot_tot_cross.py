# -*- coding: utf-8 -*-
"""
@author: 
"""
import sys
sys.path.insert(0, '../') # including the upper level of directory for the path of modules

import vulcan_cfg
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
import csv, ast
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False


plot_dir = '../plot/'
plot_name = 'S-H'
    
#molecule='CO'
sp_list = ['SH','H2S'] 
lmd_list = np.arange(10,400.05,0.5)

 
# H/OH     H2/O1D   O/H/H    OH+H     O+H2     H+OH     H2O+ 
# sCH2/H2  CH3/H    CH2/H/H  CH4+     CH3+H    CH2+H2   CH+H2/H  H+CH3    CH/H2/H  

# C2H2:  C2H2+    C2H+H    C2/H2    C2H/H    C2H2**
#CO2: Lambda  Total   sCO/O1D  CO2+     CO+O     O+CO     C+O2     sCO/O    tCO/O
#C2H4: Lambda  Total   C2H4+    C2H2+H2  C2H3+H   C2H2/2H  C2H2/H2
#C2H6: Lambda  Total   CH3/CH3  C2H5/H   1CH2/CH4 C2H6+    C2H4/H2   
#CH3CHO: Lambda  Total   CH4/CO   CH3/HCO  tCH3CHO    
#CH3OH:  Lambda  Total   CH3/OH   CH3OH+   CH3O+H   H2CO+H2  H2CO/H2 
#N2: Lambda  Total   N/N      N+N      N2+    
#NH3: Lambda  Total   NH2-H    sNH-H2   NH-H-H   NH3+     NH2+H    NH+H2    N+H2/H   H+NH2
#O2: Lambda  Total   O/O      O/O1D    O+O      O1S/O1S  O2+       
# number of branches

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180),(255, 127, 14),(44, 160, 44),(214, 39, 40),(148, 103, 189),(140, 86, 75), (227, 119, 194),(127, 127, 127),(188, 189, 34),(23, 190, 207),\
(174, 199, 232),(255, 187, 120),(152, 223, 138),(255, 152, 150),(197, 176, 213),(196, 156, 148),(247, 182, 210),(199, 199, 199),(219, 219, 141),(158, 218, 229)] 
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

n_branch = 3

sp_wavelen, sp_br_ratio = {}, {}


tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$', 'CH3CHO':'CH$_3$CHO','O3':'O$_3$','NO3':'NO$_3$' }

cross_raw, cross, disso = {}, {}, {}
color_indx = 0
for sp in sp_list:
    if sp in tex_labels: tex_sp = tex_labels[sp]
    else: tex_sp = sp
    cross_raw[sp] = np.genfromtxt('../'+ vulcan_cfg.cross_folder+sp+'_cross.csv',dtype=float,delimiter=',',skip_header=1, names = ['lambda','cross','disso','ion'])
    inter_cross = interpolate.interp1d(cross_raw[sp]['lambda'], cross_raw[sp]['cross'], bounds_error=False, fill_value=0)
    inter_disso = interpolate.interp1d(cross_raw[sp]['lambda'], cross_raw[sp]['disso'], bounds_error=False, fill_value=0)
    cross[sp] = np.zeros(len(lmd_list))
    disso[sp] = np.zeros(len(lmd_list))
    for _,ld in enumerate(lmd_list): 
        cross[sp][_] = inter_cross(ld)
        disso[sp][_] = inter_disso(ld) 
    
    plt.plot( lmd_list, cross[sp], c=tableau20[color_indx], alpha=0.8, label = tex_sp, lw=1.25)
    plt.plot( lmd_list, disso[sp], c=tableau20[color_indx], alpha=0.8, ls = '--', lw=1.5)
    color_indx += 1
    
# name_list = np.array(['lmd','cross']); name_list=np.append(name_list, np.array(list(range(1,n_branch+1)))); name_list = list(name_list)
# phid = np.genfromtxt('phidrates_data/'+molecule+".txt", skip_header=1, names=['lmd','cross','1','2','3']   ) #names=['lmd','cross','2','1','3','_','_','_','_','_','4']
# leiden = np.genfromtxt('Leiden_csv/'+molecule+'_cross.csv',dtype=float,delimiter=',',skip_header=1, names = ['lambda','cross','disso'])
# exomol = np.genfromtxt('ExoMol_highT/'+molecule+'_181-231_T1630K.txt',dtype=float, skip_header=1, names = ['lambda','cross'])

# # Read-in branching ratios from
# with open(network) as f:
#     all_lines = f.readlines()
#     br_read = False
#     for line_indx, line in enumerate(all_lines):
#         if line.startswith("# braching info start"):
#             br_read = True
#         elif line.startswith("# braching info end"):
#             br_read = False
#
#         if br_read == True and not line.startswith("#"):
#             # read in the quantum yields of photolysis reactions
#             sp_list = line.partition(':')
#             species = sp_list[0]
#             lists = sp_list[-1]
#             wavelen_yield = lists.partition(';')
#             # wavelen_yield is tuple of string in wavelength seitch, ;, Q yield e.g. ('[165.]', ';', '[(1.,0),(0,1.)]')
#             sp_wavelen[species] = ast.literal_eval(wavelen_yield[0].strip())
#             sp_br_ratio[species] = ast.literal_eval(wavelen_yield[-1].strip())
    
# # constructing the effecective cross section of each branch
# leiden_cross = {}
# wl_num = len(sp_wavelen[molecule])
# leiden_abs = leiden['cross']
# leiden_dis = leiden['disso']
# bins = leiden['lambda']

# for lmd in bins:
#
#     if wl_num == 0:
#         for nbr in range(1, n_branch+1):
#             leiden_cross[nbr] = leiden_dis *sp_br_ratio[molecule][0][nbr-1]
#     else: # wl_num == 1 or wl_num >= 2
#         for nbr in range(1, n_branch+1):
#             # the first wavelength region
#             leiden_cross[nbr] = leiden_dis *sp_br_ratio[molecule][0][nbr-1]* (bins<sp_wavelen[molecule][0])
#             # the last wavelength region
#             leiden_cross[nbr] += leiden_dis *sp_br_ratio[molecule][-1][nbr-1]* (bins>=sp_wavelen[molecule][-1])
#
#             if wl_num >= 2:
#                 for region_num in range(1,wl_num):
#                     leiden_cross[nbr] += leiden_dis *sp_br_ratio[molecule][region_num][nbr-1] * (np.array(bins>=sp_wavelen[molecule][region_num-1]) & np.array(bins<sp_wavelen[molecule][region_num]))
#
#
# color_list = ['r', 'g', 'b', 'purple', 'c', 'orange']

# plt.plot( phid['lmd']/10., phid['cross'], c='grey', alpha=0.6, label = '$\sigma_{abs}$ phidrates')
# plt.plot( leiden['lambda'], leiden['cross'], c='k', alpha=0.6, ls='--', lw=1.2,label = '$\sigma_{abs}$ Leiden')
# plt.plot( leiden['lambda'], leiden['disso'], c='k', alpha=0.6, ls=':', lw=1.3,label = '$\sigma_{dis}$ Leiden')
#
# plt.plot( exomol['lambda'], exomol['cross'], c='r', alpha=0.6, lw=1.2,label = '$\sigma_{abs}$ExoMol 1700K')


plt.yscale('log')
#plt.title(molecule)
plt.xlim((10,200))
plt.xlim((lmd_list[0],lmd_list[-1]))
#plt.ylim(bottom=1e-24)
plt.ylim((1e-23,1e-15))
plt.legend(frameon=0, prop={'size':12}, loc=3)
plt.savefig(plot_dir + plot_name + '_tot_cross.png')
plt.savefig(plot_dir + plot_name + '_tot_cross.eps')
plot = Image.open(plot_dir + plot_name + '_tot_cross.png')
plot.show()
 

# plt.figure()
# # plotting to verify
# #plt.plot( phid['lmd']/10., phid['cross'], c='grey', alpha=0.6, lw=0.9, label = '$\sigma_{abs}$ PHIDRATES')
# #plt.plot( leiden['lambda'], leiden['cross'], c='k', alpha=0.6, ls='--', lw=1.1, label = '$\sigma_{abs}$ Leiden')
# #plt.plot( leiden['lambda'], leiden['disso'], c='k', alpha=0.6, ls=':', lw=1.2, label = '$\sigma_{diss}$ Leiden')
#
#
# for nbr in range(1, n_branch+1):
#     if nbr == 1: brname = 'OH + H'
#     if nbr == 2: brname = r'H2 + $^1$O'
#     if nbr == 3: brname = 'O + H + H'
#
#     # nan is used to prevent plotting the line connecting to zero
#     plt.plot( phid['lmd']/10., [float('nan') if x==0 else x for x in phid[str(nbr)]], alpha=0.6, c=color_list[nbr-1], label=brname, lw=0.8)
#     plt.plot( bins, [float('nan') if x==0 else x for x in leiden_cross[nbr]], alpha=0.6, ls='--', lw=1.5, c=color_list[nbr-1])
#
# plt.yscale('log')
# plt.title(tex_labels[molecule] + r' + h$\nu$')
#
# plt.xlim((12,600))
# #plt.xlim(right=700.)
# plt.ylim(bottom=1e-24)
# plt.legend(frameon=0, prop={'size':10}, loc=3)
# plt.savefig(plot_dir + plot_name + '_cross.png')
# #plt.savefig(plot_dir +  plot_name + '_cross.eps')
# plot = Image.open(plot_dir + plot_name + '_cross.png')
# plot.show()
    

# H2O
    # if nbr == 1: brname = 'OH + H'
    # if nbr == 2: brname = r'H2 + $^1$O'
    # if nbr == 3: brname = 'O + H + H'
    # if nbr == 3: brname = 'O + H + H'
# CH4
    # if nbr == 1: brname = 'CH3 + H'
    # if nbr == 2: brname = r'$^1$CH2 + H2'
    # if nbr == 3: brname = '$^1$CH2 + H + H'
    # if nbr == 4: brname = 'CH + H2 + H'
# C2H2 
    # if nbr == 1: brname = r'C$_2$H + H'
    # if nbr == 2: brname = r'C$_2$ + H$_2$'
# CO2
    # if nbr == 1: brname = r'CO + O'
    # if nbr == 2: brname = r'$^1$CO + O'



# C2H4
    # if nbr == 1: brname = r'C2H2 + H$_2$'
    # if nbr == 2: brname = r'C2H2 + 2H'
    # if nbr == 3: brname = r'C2H3 + H'
    
# C2H3
    # if nbr == 1: brname = r'C2H2 + H$_2$'
    # if nbr == 2: brname = r'C2H2 + 2H'
    # if nbr == 3: brname = r'C2H3 + H'
# C2H6
    # if nbr == 1: brname = r'$C_2H_4$ + H$_2$'
    # if nbr == 2: brname = r'$C_2H_4$ + 2H'
    # if nbr == 3: brname = r'$C_2H_2$ + 2H$_2$'
    # if nbr == 4: brname = r'$CH_4$ + $^1CH_2$'
    # if nbr == 5: brname = r'$CH_3$ + $CH_3$'
# CH3CHO    
    # if nbr == 1: brname = r'CH$_4$ + CO'
    # if nbr == 2: brname = r'CH$_3$ + HCO'
# CH3OH
    # if nbr == 1: brname = r'CH$_4$ + CO'
    # if nbr == 2: brname = r'CH$_3$ + HCO'
# NH3
    # if nbr == 1: brname = r'NH$_2$ + H'
    # if nbr == 2: brname = r'NH + 2H'
 
# for molecule in mol_list:
#     L_csv_open=np.genfromtxt ('./leiden_csv/'+molecule+"_cross.csv", delimiter=",")
#     L_wave = L_csv_open[:,0]
#     L_sigma = L_csv_open[:,1]
#     plt.plot(L_wave, L_sigma,label='leiden')
#
#     csv_open=np.genfromtxt ('./phidrates_csv/'+molecule+"_cross.csv", delimiter=",")
#     wave = csv_open[:,0]
#     sigma = csv_open[:,1]
#     plt.plot(wave, sigma)
#     plt.yscale('log')
#     plt.legend(loc='left')
#
#     plt.show()