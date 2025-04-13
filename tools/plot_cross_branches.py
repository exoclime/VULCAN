# -*- coding: utf-8 -*-
"""
Created on 
@author: 
"""
import sys
sys.path.insert(0, '../') # including the upper level of directory for the path of modules
import vulcan_cfg
import numpy as np
from matplotlib import pyplot as plt
import csv, ast
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False
from scipy import interpolate
 
network = '_photo_network.txt'
plot_dir = '../plot/'


sp = 'CO2'
# number of branches
num_br = 2
# determin the reolution of bins in nm
bins = np.arange(30,240,0.2)

tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$', 'CH3CHO':'CH$_3$CHO','O3':'O_3','NO3':'NO_3' }

cross_raw, ratio_raw = {}, {}
cross_J, inter_ratio = {}, {}

cross_raw[sp] = np.genfromtxt('../'+ vulcan_cfg.cross_folder+sp+'/'+sp+'_cross.csv',dtype=float,delimiter=',',skip_header=1, names = ['lambda','cross','disso'])
ratio_raw[sp] = np.genfromtxt('../'+ vulcan_cfg.cross_folder+sp+'/'+sp+'_branch.csv',dtype=float,delimiter=',',skip_header=1, names = True)

cross =  np.zeros(len(bins))
for i in range(1,num_br+1): # fill_value extends the first and last elements for branching ratios
    cross_J[i] = np.zeros(len(bins))
    
    inter_cross = interpolate.interp1d(cross_raw[sp]['lambda'], cross_raw[sp]['cross'], bounds_error=False, fill_value=0) 
    inter_cross_J = interpolate.interp1d(cross_raw[sp]['lambda'], cross_raw[sp]['disso'], bounds_error=False, fill_value=0)
    br_key = 'br_ratio_' + str(i)
    inter_ratio[i] = interpolate.interp1d(ratio_raw[sp]['lambda'], ratio_raw[sp][br_key], bounds_error=False, fill_value=(ratio_raw[sp][br_key][0],ratio_raw[sp][br_key][-1]))
    # interpolate
    
    for n, ld in enumerate(bins):
        cross[n] = inter_cross(ld)
        cross_J[i][n] = inter_cross_J(ld) * inter_ratio[i](ld)                     
                
plt.plot( bins[cross>=0], cross[cross>=0], c='k', alpha=0.7, label = tex_labels[sp])
for i in range(1,num_br+1):
    if i == 1:
        plt.plot( bins[cross_J[i]>=0], cross_J[i][cross_J[i]>=0],  alpha=0.7, label=r'CO + O')
    elif i == 2:
        plt.plot( bins[cross_J[i]>=0], cross_J[i][cross_J[i]>=0],  alpha=0.7, label=r'$^1$CO + O')
    # elif i == 3:
    #     plt.plot( bins, cross_J[i],  alpha=0.8, label=r'$^1$CH$_2$ + H + H')
    # elif i == 4:
    #     plt.plot( bins, cross_J[i],  alpha=0.8, label=r'CH + H$_2$ + H')

plt.yscale('log')
plt.title(sp)
plt.xlabel('wavelength (nm)')
plt.ylabel(r'cross sections (cm$^2$)')
plt.xlim((bins[0],bins[-1]))
plt.ylim(bottom=1e-24)
plt.legend(frameon=0, prop={'size':12}, loc='best')
plt.savefig(plot_dir + sp + '_cross_branches.png')
plt.savefig(plot_dir + sp + '_cross_branches.pdf')
plot = Image.open(plot_dir + sp + '_cross_branches.png')
plot.show()
 

plt.figure('T')
T_list = [195,230,300,420,500,585,700,800,1160]
cross_T = {}
# high-T 423 573 1630
for _,tt in enumerate(T_list):
    cross_T[tt] = np.genfromtxt('../'+ vulcan_cfg.cross_folder+sp+'/'+sp+'_cross_'+str(tt)+'K.csv',dtype=float,delimiter=',',skip_header=1, names = ['lambda','cross'])

    inter_cross = interpolate.interp1d(cross_T[tt]['lambda'], cross_T[tt]['cross'], bounds_error=False, fill_value=0)

    cross_T_inter = np.zeros(len(bins))
    for nn,ld in enumerate(bins):
        cross_T_inter[nn] = inter_cross(ld)

    plt.plot( bins[cross_T_inter>0], cross_T_inter[cross_T_inter>0], color = plt.cm.RdYlGn(0.9 - _/(len(T_list)+1)), alpha=0.8, label = str(tt) + ' K')

plt.yscale('log')
plt.title(sp)
plt.xlim((bins[0],bins[-1]))
plt.ylim(bottom=1e-22)
plt.legend(frameon=0, prop={'size':12}, loc=3)
plt.savefig(plot_dir + sp + '_cross_T.png')
plt.savefig(plot_dir + sp + '_cross_T.eps')
plot = Image.open(plot_dir + sp + '_cross_T.png')
plot.show()

# plt.figure()
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
# plt.xlim((0,200))
# #plt.xlim(right=700.)
# plt.ylim(bottom=1e-23)
# plt.legend(frameon=0, prop={'size':10}, loc=3)
# plt.savefig(plot_dir + molecule + '_cross.png')
# plt.savefig(plot_dir +  molecule + '_cross.eps')
# plot = Image.open(plot_dir + molecule + '_cross.png')
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