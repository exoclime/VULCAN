# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 10:29:09 2017

@author: 
"""
import numpy as np
from matplotlib import pyplot as plt
import csv, ast
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False

network = 'New-NCHO_earth_photo_network.txt'
plot_dir = '../plot/'

    
#molecule='CO'
#mol_list = ['H2O'] # 'CH4','CH3','C2H2','C2H4','C2H6','CO','CO2','H2','H2O','H2CO','HCO','HCN','N2','NH3','NO','O2', 'OH', 'NO2','HO2','H2O2'
molecule = 'SH'
# H/OH     H2/O1D   O/H/H    OH+H     O+H2     H+OH     H2O+ 
# sCH2/H2  CH3/H    CH2/H/H  CH4+     CH3+H    CH2+H2   CH+H2/H  H+CH3    CH/H2/H  

#O2: Lambda  Total   O/O      O/O1D    O+O      O1S/O1S  O2+       
# number of branches
n_branch = 1

sp_wavelen, sp_br_ratio = {}, {}


tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$', 'CH3CHO':'CH$_3$CHO','O3':'O_3','NO3':'NO_3' }

T_list = [500, 750, 1000, 1500, 2000, 2500, 3000]
cross = {}
for T in T_list:
    cross[T] = np.genfromtxt('../thermo/photo_cross/'+molecule+'/'+molecule+'_cross_'+str(T)+'K.csv',delimiter=',', skip_header=1, names = ['lambda','cross'])

leiden = np.genfromtxt('../thermo/photo_cross/'+molecule+'/'+molecule+'_cross.csv',delimiter=',', skip_header=1, names = ['lambda','cross'])
#for molecule in mol_list:
    # input from phiDrates database

# name_list = np.array(['lmd','cross']); name_list=np.append(name_list, np.array(list(range(1,n_branch+1)))); name_list = list(name_list)
# phid = np.genfromtxt('phidrates_data/'+molecule+".txt", skip_header=1, names=['lmd','cross','1','2','3']   ) #names=['lmd','cross','2','1','3','_','_','_','_','_','4']
# leiden = np.genfromtxt('Leiden_csv/'+molecule+'_cross.csv',dtype=float,delimiter=',',skip_header=1, names = ['lambda','cross','disso'])
#
# exomol = np.genfromtxt('ExoMol_highT/'+molecule+'_181-231_T1630K.txt',dtype=float, skip_header=1, names = ['lambda','cross'])

color_list = ['r', 'g', 'b', 'purple', 'c', 'orange']    
for _,T in enumerate(T_list):
    plt.plot( cross[T]['lambda'], cross[T]['cross'], alpha=0.8, label = str(T), lw=1.5)
    
    
plt.plot( leiden['lambda'], leiden['cross'], c='k', alpha=0.6, label = 'reference')


plt.yscale('log')
plt.title(molecule)
plt.xlim((260,450))
plt.ylim(bottom=1e-29)
plt.legend(frameon=0, prop={'size':12}, loc=3)
plt.savefig(plot_dir + '_cross.png')
#plt.savefig(plot_dir + molecule + '_tot_cross.eps')
plot = Image.open(plot_dir + '_cross.png')
plot.show()
 

# plt.figure()
# # plotting to verify
# #plt.plot( phid['lmd']/10., phid['cross'], c='grey', alpha=0.6, lw=0.9, label = '$\sigma_{abs}$ PHIDRATES')
# #plt.plot( leiden['lambda'], leiden['cross'], c='k', alpha=0.6, ls='--', lw=1.1, label = '$\sigma_{abs}$ Leiden')
# #plt.plot( leiden['lambda'], leiden['disso'], c='k', alpha=0.6, ls=':', lw=1.2, label = '$\sigma_{diss}$ Leiden')
#
#
#
#
# # H2O
#     # if nbr == 1: brname = 'OH + H'
#     # if nbr == 2: brname = r'H2 + $^1$O'
#     # if nbr == 3: brname = 'O + H + H'
#     # if nbr == 3: brname = 'O + H + H'
# # CH4
#     # if nbr == 1: brname = 'CH3 + H'
#     # if nbr == 2: brname = r'$^1$CH2 + H2'
#     # if nbr == 3: brname = '$^1$CH2 + H + H'
#     # if nbr == 4: brname = 'CH + H2 + H'
# # C2H2
#     # if nbr == 1: brname = r'C$_2$H + H'
#     # if nbr == 2: brname = r'C$_2$ + H$_2$'
# # CO2
#     # if nbr == 1: brname = r'CO + O'
#     # if nbr == 2: brname = r'$^1$CO + O'
#
#
#
# # C2H4
#     # if nbr == 1: brname = r'C2H2 + H$_2$'
#     # if nbr == 2: brname = r'C2H2 + 2H'
#     # if nbr == 3: brname = r'C2H3 + H'
#
# # C2H3
#     # if nbr == 1: brname = r'C2H2 + H$_2$'
#     # if nbr == 2: brname = r'C2H2 + 2H'
#     # if nbr == 3: brname = r'C2H3 + H'
# # C2H6
#     # if nbr == 1: brname = r'$C_2H_4$ + H$_2$'
#     # if nbr == 2: brname = r'$C_2H_4$ + 2H'
#     # if nbr == 3: brname = r'$C_2H_2$ + 2H$_2$'
#     # if nbr == 4: brname = r'$CH_4$ + $^1CH_2$'
#     # if nbr == 5: brname = r'$CH_3$ + $CH_3$'
# # CH3CHO
#     # if nbr == 1: brname = r'CH$_4$ + CO'
#     # if nbr == 2: brname = r'CH$_3$ + HCO'
# # CH3OH
#     # if nbr == 1: brname = r'CH$_4$ + CO'
#     # if nbr == 2: brname = r'CH$_3$ + HCO'
# # NH3
#     # if nbr == 1: brname = r'NH$_2$ + H'
#     # if nbr == 2: brname = r'NH + 2H'
#
# # for molecule in mol_list:
# #     L_csv_open=np.genfromtxt ('./leiden_csv/'+molecule+"_cross.csv", delimiter=",")
# #     L_wave = L_csv_open[:,0]
# #     L_sigma = L_csv_open[:,1]
# #     plt.plot(L_wave, L_sigma,label='leiden')
# #
# #     csv_open=np.genfromtxt ('./phidrates_csv/'+molecule+"_cross.csv", delimiter=",")
# #     wave = csv_open[:,0]
# #     sigma = csv_open[:,1]
# #     plt.plot(wave, sigma)
# #     plt.yscale('log')
# #     plt.legend(loc='left')
# #
# #     plt.show()