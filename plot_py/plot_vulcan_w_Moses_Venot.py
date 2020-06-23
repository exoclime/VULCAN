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
#vul_data = sys.argv[1]
# Setting the 3rd input argument as the species names to be plotted (separated by ,)
plot_spec = sys.argv[1]
# Setting the 4th input argument as the output eps filename        
plot_name = sys.argv[2]

vul_data = '../output/HD189.vul'
vul_data2 = '../output/HD189-HNO.vul'

venot_HD189 = np.genfromtxt('../output/venot/venot_HD189_steady.dat',names=True)


plot_dir = vulcan_cfg.plot_dir

# taking user input species and splitting into separate strings and then converting the list to a tuple
plot_spec = tuple(plot_spec.split(','))
nspec = len(plot_spec)

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180),(255, 127, 14),(44, 160, 44),(214, 39, 40),(148, 103, 189),(140, 86, 75), (227, 119, 194),(127, 127, 127),(188, 189, 34),(23, 190, 207),\
(174, 199, 232),(255, 187, 120),(152, 223, 138),(255, 152, 150),(197, 176, 213),(196, 156, 148),(247, 182, 210),(199, 199, 199),(219, 219, 141),(158, 218, 229)] 

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)
tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
,'CH2OH':'CH$_2$OH','N2':'N$_2$','NH3':'NH$_3$', 'NO2':'NO$_2$','HCN':'HCN','NO':'NO', 'NO2':'NO$_2$','CH3CHO': 'CH$_3$CHO'}


with open(vul_data, 'rb') as handle:
  data = pickle.load(handle)
  
with open(vul_data2, 'rb') as handle:
  data2 = pickle.load(handle)
#
# with open('output/hd189_EQ_moses.vul', 'rb') as handle:
#   EQ_data = pickle.load(handle)  


atm = 'hd189_day_photo'
#atm = 'hd209_day_photo'
 
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
    JM[0] = np.genfromtxt('../atm/JM/reorder0_HD189_wtday.txt', names=True)
    for i in  range(1,10):
        JM[i] = np.genfromtxt('../atm/JM/reorder'+str(i)+'_HD189_wtday.txt',names=True)
        JM_labels[i] = np.genfromtxt('../atm/JM/reorder'+str(i)+'_HD189_wtday.txt', dtype=str)[0]
 
elif atm == 'hd209_day_photo':
    JM, JM_labels = {}, {}
    JM[0] = np.genfromtxt('../atm/JM/TPK_HD209_dayave.txt', names=True)
    for i in  range(1,7):
        JM[i] = np.genfromtxt('../atm/JM/reorder'+str(i)+'_HD209_dayave.txt',names=True)
        JM_labels[i] = np.genfromtxt('../atm/JM/reorder'+str(i)+'_HD209_dayave.txt', dtype=str)[0]
        
# plt.figure('TP')
# plt.plot(data['atm']['Tco'], data['atm']['pco']/1.e6,  c='k')
# #plt.plot(venot_HD189['T'], venot_HD189['P']/1.e3, color=tableau20[color_index], ls='-.', c='g', lw=1.5)
# plt.plot(JM[0]['TEMPERATURE'], JM[0]['PRESSURE']/1.e3, ls='--', lw=1.5, c='b')
# plt.gca().set_yscale('log')
# plt.gca().invert_yaxis()
# plt.savefig(plot_dir + 'TP' + '.png')
# plt.savefig(plot_dir + 'TP' + '.eps')
# if vulcan_cfg.use_PIL == True:
#     plot = Image.open(plot_dir + 'TP' + '.png')
#     plot.show()
# else: plt.show()

plt.figure()
color_index = 0
vulcan_spec = data['variable']['species']

#plot_spec = vulcan_spec[50:60]

for color_index,sp in enumerate(plot_spec):
    if sp in tex_labels: sp_lab = tex_labels[sp]
    else: sp_lab = sp
    if color_index == len(tableau20): # when running out of colors
        tableau20.append(tuple(np.random.rand(3)))
    plt.plot(data['variable']['ymix'][:,vulcan_spec.index(sp)], data['atm']['pco']/1.e6, color=tableau20[color_index], label=sp_lab, alpha=0.7, lw=1.2)
    
    #plt.plot(data['variable']['y_ini'][:,vulcan_spec.index(sp)] / data['atm']['n_0'], data['atm']['pco']/1.e6, color=tableau20[color_index], ls=':', alpha=0.6, lw=1.2)
    if sp in venot_HD189.dtype.names:
        plt.plot(venot_HD189[sp], venot_HD189['P']/1.e3, color=tableau20[color_index], ls='-.', alpha=0.7, lw=1.5)
        pass
    
    plt.plot(data2['variable']['ymix'][:,vulcan_spec.index(sp)], data2['atm']['pco']/1.e6, color=tableau20[color_index], ls=':', lw=1.8, alpha=0.8)

    #plt.plot(EQ_data['variable']['ymix'][:,vulcan_spec.index(sp)], EQ_data['atm']['pco']/1.e6, color=tableau20[color_index], ls='--', lw=1.5)
    
    if sp == 'He': sp = 'HE'
    if sp == 'CH2': sp = '3CH2'
    if sp == 'CH2_1': sp = '1CH2'
    if sp == 'O_1': sp = '(1)O'

    if not (atm == 'hd189_day_photo' or atm == 'hd209_day_photo'):

        if sp in JM_data.keys():

            plt.plot(JM_data[sp]/JM_data['DENSITY'], JM_data['PRESSURE']/1.e3, color=tableau20[color_index], ls='--', lw=1.5)
        else:
            print (sp + ' not in JM.')

    elif atm == 'hd189_day_photo':
            for i in  range(1,10):
                if sp in JM_labels[i]:
                    plt.plot(JM[i][sp]/JM[0]['DENSITY'], JM[i]['PRESSURE']/1.e3, color=tableau20[color_index], ls='--', lw=1.5, alpha=0.7)

    elif atm == 'hd209_day_photo':
            for i in  range(1,7):
                if sp in JM_labels[i]:
                    plt.plot(JM[i][sp]/JM[0]['DENSITY'], JM[i]['PRESSURE']/1.e3, color=tableau20[color_index], ls='--', lw=1.5)
           
         
plt.gca().set_xscale('log')       
plt.gca().set_yscale('log') 
plt.gca().invert_yaxis() 
plt.xlim((1e-25,0.999))
#plt.xlim((1.e-14,0.99))
plt.ylim((data['atm']['pco'][0]/1.e6,data['atm']['pco'][-1]/1e6))
plt.legend(frameon=0, prop={'size':14}, loc=3)
plt.title('HD189733b\ncompared with Moses et al.(2011) and Venot et al.(2012)', fontsize=10)
#plt.title('HD189733b')
#plt.title('HD189733b resolution test')
#plt.title('HD189')
plt.xlabel("Mixing Ratio", fontsize=14)
plt.ylabel("Pressure (bar)", fontsize=14)
plt.savefig('../' + plot_dir + plot_name + '.png')
plt.savefig('../' +plot_dir + plot_name + '.eps')
if vulcan_cfg.use_PIL == True:
    plot = Image.open('../' +plot_dir + plot_name + '.png')
    plot.show()
else: plt.show()





# debug
# re_list = []
# for re in range(1,627):
#     err = data['variable']['k'][re]/ data2['variable']['k'][re]
#     if max(err) >1.01 or min(err) < 0.99:
#         print (re)
#         print (err)
#         re_list.append(re)
    