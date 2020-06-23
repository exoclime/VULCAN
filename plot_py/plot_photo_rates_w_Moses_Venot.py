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
from phy_const import kb

atm = 'hd189_day_photo'
#atm = 'hd209_day_photo'

vul_data = 'output/db02-HD189.vul'
vul_data2 = 'output/no-JCH3-test-v118-HD189.vul'


venot_HD189 = np.genfromtxt('output/venot/venot_HD189_steady.dat',names=True)

venot_data_ch4 = 'output/venot/taux_loss_CH4.dat'
venot_data_nh3 = 'output/venot/taux_loss_NH3.dat'
venot_ch4_rates = np.genfromtxt(venot_data_ch4, skip_header=6, names = ['hight','P', 'total', 'R11','R12','R13','R14']) # P in mbar
venot_nh3_rates = np.genfromtxt(venot_data_nh3, skip_header=3, names = ['hight','P', 'total','R17']) # P in mbar


# Setting the input argument as the species names to be plotted (separated by ,)
#plot_spec = sys.argv[1]
# Setting the input argument as the output eps filename        
plot_name = sys.argv[1]


plot_dir = vulcan_cfg.plot_dir

colors = ['c','b','g','r','m','y','darkblue','orange','pink','silver','darkred','salmon','chocolate','steelblue','plum','hotpink','k']

tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$' }


with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)
  
with open(vul_data2, 'rb') as handle:
  data2 = pickle.load(handle)


 
# Read Moses 2011 for Tiso and const-Kzz 
JM_data={}
if atm == 'hd189':
    JM_f='atm/JM189_K9_nopho.inp'
    with open(JM_f, "r") as f:
        for line in f:
            if "%" in line:
                name = line[2:-1].strip()
                JM_data.update({name: "" })
            else:
                JM_data[name] = JM_data[name] + " " + line[:-1].strip() 
    for _ in JM_data:
        JM_data[_] = np.fromstring(JM_data[_], dtype=float, sep=" ")
        
elif atm == 'hd209_K9':
    JM_f='atm/JM209_K9_nopho.inp'
    with open(JM_f, "r") as f:
        for line in f:
            if "%" in line:
                name = line[2:-1].strip()
                JM_data.update({name: "" })
            else:
                JM_data[name] = JM_data[name] + " " + line[:-1].strip() 
    for _ in JM_data:
        JM_data[_] = np.fromstring(JM_data[_], dtype=float, sep=" ")
        
elif atm == 'hd189_day_photo':
    JM, JM_labels = {}, {}
    JM[0] = np.genfromtxt('atm/JM/reorder0_HD189_wtday.txt', names=True)
    for i in  [1,2,3,4,9,10,11,12,13]:
        JM[i] = np.genfromtxt('atm/JM/photo_rates_HD189_wtday/photo_'+str(i)+'.txt',names=True)
        
        JM_labels[i] = np.genfromtxt('atm/JM/photo_rates_HD189_wtday/photo_'+str(i)+'.txt', dtype=str)[0]
 
elif atm == 'hd209_day_photo':
    JM, JM_labels = {}, {}
    JM[0] = np.genfromtxt('atm/JM/TPK_HD209_dayave.txt', names=True)
    for i in  [1,2,3,4,9,10,11]:
        JM[i] = np.genfromtxt('atm/JM/photo_rates/photo_'+str(i)+'.txt',names=True)
        
        JM_labels[i] = np.genfromtxt('atm/JM/photo_rates/photo_'+str(i)+'.txt', dtype=str)[0]
 

 
color_index = 0
#vulcan_spec = data['variable']['species']

# Index w Moses
# H2O-1: 86
# H2O-2: 88
# H2O-3: 87
# CH4-1: 5
# CH4-2: 6
# CH4-3: 7
# CH4-4: 9
# CH3-1: 3
# CH3-2: 8
# CO: 89
# H2: 1
# NH3-1: 105
# NH3 -> NH + 2H: 106 

moses_index = {('H2O',1):'K86',('H2O',2):'K88',('H2O',3):'K87', ('CH4',1):'K5', ('CH4',2):'K6', ('CH4',3):'K7', ('CH4',4):'K9', ('CH3',1):'K3', ('CH3',2):'K8'\
,('NH3',1):'K105',('NH3',2):'K106', ('H2',1):'K1',('CO',1):'K89',('CO2',1):'K90',('CO2',2):'K91',('N2',1):'K103',('HCN',1):'K109',('H2CO',1):'K94',('H2CO',2):'K93',('HCO',1):'K92'\
,('C2H2',1):'K10',('C2H2',2):'K11',('C2H4',1):'K13',('C2H4',2):'K14',('C2H4',3):'K15',('C2H6',1):'K17',('C2H6',2):'K18',('C2H6',3):'K19',('C2H6',4):'K20',('C2H6',5):'K21',('OH',1):'K85',\
('NO',1):'K128', ('NO2',1):'K129'    }

venot_index = {('CH4',1): 'R11', ('CH4',2): 'R12', ('CH4',3): 'R13', ('CH4',4): 'R14', ('NH3',1): 'R17' }

vulcan_sp = 'CH4'
branch = 1
re = moses_index[(vulcan_sp,branch)]

vulcan_TP = interpolate.interp1d(data['atm']['pco'], data['atm']['Tco'], assume_sorted = False, bounds_error=False,fill_value=(data['atm']['Tco'][np.argmin(data['atm']['pco'])], data['atm']['Tco'][np.argmax(data['atm']['pco'])] )  )
venot_inter_Tco = vulcan_TP(venot_HD189['P']*1e3)
venot_inter_n_0 = venot_HD189['P']*1e3/(venot_inter_Tco*kb)

plt.plot(data['variable']['J_sp'][(vulcan_sp,0)], data['atm']['pco']/1.e6, color='k', label='total (this work)')



plt.plot(data['variable']['J_sp'][(vulcan_sp,branch)], data['atm']['pco']/1.e6, color='r', label='CH$_3$ + H (this work)', alpha=0.5)

#plt.plot(data2['variable']['J_sp'][(vulcan_sp,branch)], data2['atm']['pco']/1.e6, color=colors[color_index], ls=':', lw=1.5, alpha=0.8)

 
 
    #plt.plot(EQ_data['variable']['ymix'][:,vulcan_spec.index(sp)], EQ_data['atm']['pco']/1.e6, color=colors[color_index], ls='--', lw=1.5)
    

            
if atm == 'hd189_day_photo':
        for i in  [1,2,3,4,9,10,11,12,13]:
            re = moses_index[(vulcan_sp,1)]
            if re in JM_labels[i]:
                plt.plot(JM[i][re], JM[0]['PRESSURE']/1.e3, color='r', ls='--', label='CH$_3$ + H (M11)', alpha=0.5)
        
                JM_tot = JM[i][re]
        #plt.plot(venot_nh3_rates[venot_index[(vulcan_sp,branch)]]/(venot_inter_n_0*venot_HD189[vulcan_sp]), venot_nh3_rates['P']/1.e3, color='k', ls='-.', label='CH$_3$ + H (V12)')   
        plt.plot(venot_ch4_rates[venot_index[(vulcan_sp,branch)]]/(venot_inter_n_0*venot_HD189[vulcan_sp]), venot_ch4_rates['P']/1.e3, color='r', ls='-.', label='CH$_3$ + H (V12)', alpha=0.5)   
        
        plt.plot(data['variable']['J_sp'][(vulcan_sp,2)], data['atm']['pco']/1.e6, color='g', label='$^1$CH$_2$ + H2 (this work)', alpha=0.5)
        
        for i in  [1,2,3,4,9,10,11,12,13]:
            re = moses_index[(vulcan_sp,2)]
            if re in JM_labels[i]:
                plt.plot(JM[i][re], JM[0]['PRESSURE']/1.e3, color='g', ls='--', label='$^1$CH$_2$ + H2 (M11)')
                JM_tot += JM[i][re]
                
        plt.plot(venot_ch4_rates[venot_index[(vulcan_sp,2)]]/(venot_inter_n_0*venot_HD189[vulcan_sp]), venot_ch4_rates['P']/1.e3, color='g', ls='-.', label='$^1$CH$_2$ + H2 (V12)')   
        
        ### J-3
        plt.plot(data['variable']['J_sp'][(vulcan_sp,3)], data['atm']['pco']/1.e6, color='b', label='$^1$CH$_2$ + H + H (this work)')
        for i in  [1,2,3,4,9,10,11,12,13]:
            re = moses_index[(vulcan_sp,3)]
            if re in JM_labels[i]:
                plt.plot(JM[i][re], JM[0]['PRESSURE']/1.e3, color='b', ls='--', label='$^1$CH$_2$ + H + H (M11)')
                JM_tot += JM[i][re]
                
        plt.plot(venot_ch4_rates[venot_index[(vulcan_sp,3)]]/(venot_inter_n_0*venot_HD189[vulcan_sp]), venot_ch4_rates['P']/1.e3, color='b', ls='-.', label='$^1$CH$_2$ + H + H (V12)')   
        
        ### J-4
        plt.plot(data['variable']['J_sp'][(vulcan_sp,4)], data['atm']['pco']/1.e6, color='m', label='CH + H2 + H (this work)')
        for i in  [1,2,3,4,9,10,11,12,13]:
            re = moses_index[(vulcan_sp,4)]
            if re in JM_labels[i]:
                plt.plot(JM[i][re], JM[0]['PRESSURE']/1.e3, color='m', ls='--', label='CH + H2 + H (M11)')
                JM_tot += JM[i][re]
                
        plt.plot(venot_ch4_rates[venot_index[(vulcan_sp,4)]]/(venot_inter_n_0*venot_HD189[vulcan_sp]), venot_ch4_rates['P']/1.e3, color='m', ls='-.', label='CH + H2 + H (V12)')   
        
        
        plt.plot(JM_tot, JM[0]['PRESSURE']/1.e3, color='k', ls='--', label='CH + H2 + H (M11)')
        plt.plot((venot_ch4_rates[venot_index[(vulcan_sp,1)]]+venot_ch4_rates[venot_index[(vulcan_sp,2)]]+venot_ch4_rates[venot_index[(vulcan_sp,3)]]+venot_ch4_rates[venot_index[(vulcan_sp,4)]])/(venot_inter_n_0*venot_HD189[vulcan_sp]), venot_ch4_rates['P']/1.e3, color='k', ls='-.', label='CH + H2 + H (V12)')   
        
        
elif atm == 'hd209_day_photo':
        for i in  [1,2,3,4,9,10]:
            if re in JM_labels[i]:
                plt.plot(JM[i][re], JM[0]['PRESSURE']/1.e3, color=colors[color_index], ls='--')
           
         
plt.gca().set_xscale('log')       
plt.gca().set_yscale('log') 
plt.gca().invert_yaxis() 
plt.xlim((1.e-6,0.05))
#plt.title('J-'+vulcan_sp+'--'+str(branch))
#plt.title('NH$_3$ photodissociation')
plt.ylim((data['atm']['pco'][0]/1.e6,data['atm']['pco'][-1]/1.e6))
plt.legend(frameon=0, prop={'size':10}, loc=4)
plt.xlabel('CH$_4$ photodissociation rate (cm$^{-1}$)')
plt.ylabel("Pressure (bar)")
plt.savefig(plot_dir + plot_name + '.png')
plt.savefig(plot_dir + plot_name + '.eps')
if vulcan_cfg.use_PIL == True:
    plot = Image.open(plot_dir + plot_name + '.png')
    plot.show()
else: plt.show()

