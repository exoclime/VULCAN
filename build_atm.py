import numpy as np
from numpy import polynomial
import scipy
from scipy import interpolate
import scipy.optimize as sop
import subprocess
import pickle
from shutil import copyfile 

import vulcan_cfg
from phy_const import kb, Navo, r_sun, au
from vulcan_cfg import nz
import chem_funs
from chem_funs import ni, nr  # number of species and reactions in the network
species = chem_funs.spec_list

### read in the basic chemistry data
with open(vulcan_cfg.com_file, 'r') as f:
    columns = f.readline() # reading in the first line
    num_ele = len(columns.split())-2 # number of elements (-2 for removing "species" and "mass") 
type_list = ['int' for i in range(num_ele)]
type_list.insert(0,'U20'); type_list.append('float')
compo = np.genfromtxt(vulcan_cfg.com_file,names=True,dtype=type_list)
# dtype=None in python 2.X but Sx -> Ux in python3
compo_row = list(compo['species'])
### read in the basic chemistry data


class InitialAbun(object):
    """
    Calculating the appropriate initial mixing ratios with the assigned elemental abundance
    """
    
    def __init__(self):
        self.ini_m = [0.9,0.1,0.,0.,0] # initial guess
        #self.EQ_ini_file = vulcan_cfg.EQ_ini_file
        
        self.atom_list = vulcan_cfg.atom_list

    def abun_lowT(self, x):
        """
        calculating the initial mixing ratios of the following 5 molecules (with CH4) 
        satisfying the assigned elemental abundance
        x1:H2 x2:H2O x3:CH4 x4:He x5:NH3
        """
        O_H, C_H, He_H, N_H = vulcan_cfg.O_H, vulcan_cfg.C_H, vulcan_cfg.He_H, vulcan_cfg.N_H
        x1,x2,x3,x4,x5 = x
        f1 = x1+x2+x3+x4-1.
        f2 = x2 - (2*x1+2*x2+4*x3+3*x5)*O_H
        f3 = x3 - (2*x1+2*x2+4*x3+3*x5)*C_H
        f4 = x4 - (2*x1+2*x2+4*x3+3*x5)*He_H
        f5 = x5 - (2*x1+2*x2+4*x3+3*x5)*N_H
        return f1,f2,f3,f4,f5
        
    def abun_highT(self, x):
        """
        calculating the initial mixing ratios of the following 4 molecules (with CO)
        satisfying the assigned elemental abundance
        x1:H2 x2:H2O x3:CO x4:He x5:N2
        """
        O_H, C_H, He_H, N_H = vulcan_cfg.O_H, vulcan_cfg.C_H, vulcan_cfg.He_H, vulcan_cfg.N_H
        x1,x2,x3,x4,x5 = x
        f1 = x1+x2+x3+x4-1.
        f2 = x2+x3 - (2*x1+2*x2)*O_H
        f3 = x3 - (2*x1+2*x2)*C_H
        f4 = x4 - (2*x1+2*x2)*He_H
        f5 = x5*2 - (2*x1+2*x2)*N_H
        return f1,f2,f3,f4,f5
    
    
    def ini_mol(self):
        if vulcan_cfg.ini_mix == 'const_lowT':
            return np.array(sop.fsolve(self.abun_lowT, self.ini_m))
        # somehow const_highT is not stable at high P...
        # elif vulcan_cfg.ini_mix == 'const_highT':
            # return np.array(sop.fsolve(self.abun_highT, self.ini_m))
            
    def ini_fc(self, data_var, data_atm):
        # reading-in the default elemental abundances from Lodders 2009
        # depending on including ion or not (whether there is e- in the fastchem elemental abundance dat)
        tmp_str = ""
        solar_ele = 'fastchem_vulcan/input/solar_element_abundances.dat'
        if vulcan_cfg.use_ion == True:
            copyfile('fastchem_vulcan/input/parameters_ion.dat', 'fastchem_vulcan/input/parameters.dat')
        else:
            copyfile('fastchem_vulcan/input/parameters_wo_ion.dat', 'fastchem_vulcan/input/parameters.dat')
            
        with open(solar_ele ,'r') as f:
            new_str = ""
            ele_list = list(vulcan_cfg.atom_list)
            ele_list.remove('H')
            
            fc_list = ['C', 'N', 'O', 'S', 'P', 'Si', 'Ti','V','Cl','K','Na','Mg','F','Ca','Fe']
            
            if vulcan_cfg.use_solar == True: 
                new_str = f.read() # read in as a string
                print ("Initializing with the default solar abundance.")
                
            else: # using costomized elemental abundances
                print ("Initializing with the customized elemental abundance:")
                print ("{:4}".format('H') + str('1.'))
                for line in f.readlines():   
                        li = line.split()
                        sp = li[0].strip()
                        
                        if sp in ele_list:
                            # read-in vulcan_cfg.sp_H
                            sp_abun = getattr(vulcan_cfg, sp+'_H')
                            fc_abun = 12. + np.log10(sp_abun)
                            line = sp + '\t' + "{0:.4f}".format(fc_abun) + '\n'
                            print ("{:4}".format(sp) + "{0:.4E}".format(sp_abun))

                        elif sp in fc_list: # other elements included in fastchem but not in VULCAN 
                            sol_ratio = li[1].strip()
                            # print (sp + ":  " + str(sol_ratio))
                            if hasattr(vulcan_cfg, 'fastchem_met_scale'): # vulcan_cfg.fastchem_met_scale exists 
                                met_scale = vulcan_cfg.fastchem_met_scale
                            else:
                                met_scale = 1.
                                print ("fastchem_met_scale not specified in vulcan_cfg. Using solar metallicity for other elements not included in vulcan.")
                            
                            new_ratio = float(sol_ratio) + np.log10(met_scale)
                            line = sp + '\t' + "{0:.4f}".format(new_ratio) + '\n'
                            
                        new_str += line
                
            # make the new elemental abundance file
            with open('fastchem_vulcan/input/element_abundances_vulcan.dat', 'w') as f: f.write(new_str)
            
        # write a T-P text file for fast_chem
        with open('fastchem_vulcan/input/vulcan_TP/vulcan_TP.dat' ,'w') as f:
            ost = '#p (bar)    T (K)\n'   
            for n, p in enumerate(data_atm.pco): # p in bar in fast_chem
                ost +=  '{:.3e}'.format(p/1.e6) + '\t' + '{:.1f}'.format(data_atm.Tco[n])  + '\n'
            ost = ost[:-1]
            f.write(ost)
        
        try: subprocess.check_call(["./fastchem input/config.input"], shell=True, cwd='fastchem_vulcan/') # check_call instead of call can catch the error 
        except: print ('\n FastChem cannot run properly. Try compile it by running make under /fastchem_vulcan\n'); raise
           
    def ini_y(self, data_var, data_atm): 
        # initial mixing ratios of the molecules
        
        ini_mol = self.ini_mol()  
        ini = np.zeros(ni)
        y_ini = data_var.y
        gas_tot = data_atm.M
        charge_list = [] # list of charged species excluding echarge_list
        
        if vulcan_cfg.ini_mix == 'EQ':
        
            self.ini_fc(data_var, data_atm)
            fc = np.genfromtxt('fastchem_vulcan/output/vulcan_EQ.dat', names=True, dtype=None, skip_header=0)
            for sp in species:
                if sp in fc.dtype.names:
                    y_ini[:,species.index(sp)] = fc[sp]*gas_tot # this also changes data_var.y because the address of y array has passed to y_ini
                
                else: print (sp + ' not included in fastchem.')
                
                if vulcan_cfg.use_ion == True:
                    if compo[compo_row.index(sp)]['e'] != 0: charge_list.append(sp)
            
            # remove the fc output
            subprocess.call(["rm vulcan_EQ.dat"], shell=True, cwd='fastchem_vulcan/output/')
                             
        elif vulcan_cfg.ini_mix == 'vulcan_ini':
            print ("Initializing with compositions from the prvious run " + vulcan_cfg.vul_ini)
            with open(vulcan_cfg.vul_ini, 'rb') as handle:
              vul_data = pickle.load(handle) 
            
            #y_ini = np.copy(vul_data['variable']['y'])
            #data_var.y = np.copy(y_ini)
            for sp in species:
                if sp in vul_data['variable']['species']:
                    y_ini[:,species.index(sp)] = vul_data['variable']['y'][:,vul_data['variable']['species'].index(sp)]
                else: print (sp + " not included in the prvious run.") 
            
            #if vulcan_cfg.use_ion == True: charge_list = vul_data['variable']['charge_list']
            if vulcan_cfg.use_ion == True:
                for sp in species: 
                    if compo[compo_row.index(sp)]['e'] != 0: charge_list.append(sp)
        
        elif vulcan_cfg.ini_mix  == 'table':
            table = np.genfromtxt(vulcan_cfg.vul_ini, names=True, dtype=None, skip_header=1) 
            if not len(data_atm.pco) == len(table['Pressure']): 
                print ("Warning! The initial profile has different layers than the current setting...")
                raise
            for sp in species:  
                data_var.y[:,species.index(sp)] = data_atm.n_0 * table[sp]

    
        elif vulcan_cfg.ini_mix == 'const_mix':
            print ("Initializing with constant (well-mixed): " + str(vulcan_cfg.const_mix))
            for sp in vulcan_cfg.const_mix.keys():
                y_ini[:,species.index(sp)] = gas_tot* vulcan_cfg.const_mix[sp] # this also changes data_var.y
            if vulcan_cfg.use_ion == True:
                for sp in species: 
                    if compo[compo_row.index(sp)]['e'] != 0: charge_list.append(sp)
                
        else:
            for i in range(nz):
                
                if vulcan_cfg.ini_mix == 'const_lowT':
                    y_ini[i,:] = ini
                    y_ini[i,species.index('H2')] = ini_mol[0]*gas_tot[i]; y_ini[i,species.index('H2O')] = ini_mol[1]*gas_tot[i]; y_ini[i,species.index('CH4')] = ini_mol[2]*gas_tot[i]
                    y_ini[i,species.index('NH3')] = ini_mol[4]*gas_tot[i]
                    # assign rest of the particles to He
                    y_ini[i,species.index('He')] = gas_tot[i] - np.sum(y_ini[i,:])

                else:
                    raise IOError ('\nInitial mixing ratios unknown. Check the setting in vulcan_cfg.py.')
        
        if vulcan_cfg.use_condense == True:
            for sp in vulcan_cfg.condense_sp:
                data_atm.sat_mix[sp] = data_atm.sat_p[sp]/data_atm.pco
                
                # fixed 2022
                data_atm.sat_mix[sp] = np.minimum(1., data_atm.sat_mix[sp])
                
                if sp == 'H2O': data_atm.sat_mix[sp] *= vulcan_cfg.humidity
                
                if vulcan_cfg.use_ini_cold_trap == True:
                    
                    if  vulcan_cfg.ini_mix != 'table' and vulcan_cfg.ini_mix != 'vul_ini':
                        # the level where condensation starts    
                        conden_bot = np.argmax( data_atm.n_0*data_atm.sat_mix[sp] <= data_var.y[:,species.index(sp)] )
                        # conden_status: ture if the partial p >= the saturation p
                        sat_rho = data_atm.n_0 * data_atm.sat_mix[sp]
                        conden_status = data_var.y[:,species.index(sp)] >= sat_rho

                        # take the min between the mixing ratio and the saturation mixing ratio
                        data_var.y[:,species.index(sp)] = np.minimum(data_atm.n_0 * data_atm.sat_mix[sp], data_var.y[:,species.index(sp)])

                        if list(data_var.y[conden_status,species.index(sp)]): # if it condenses
                            min_sat = np.amin(data_atm.sat_mix[sp][conden_status]) # the mininum value of the saturation p within the saturation region
                            conden_min_lev = np.where(data_atm.sat_mix[sp] == min_sat)[0][0]
                            
                            data_atm.conden_min_lev = conden_min_lev
                            
                            print ( sp + " condensed from nz = " + str(conden_bot) + " to the minimum level nz = "+ str(conden_min_lev) + " (cold trap)") 
                            #data_var.y[conden_min_lev:,species.index(sp)] = (y_ini[conden_min_lev,species.index(sp)]/data_atm.n_0[conden_min_lev]) *data_atm.n_0[conden_min_lev:]
                            data_var.y[conden_min_lev:,species.index(sp)] = data_atm.sat_mix[sp][conden_min_lev] * data_atm.n_0[conden_min_lev:]  
                           
        # re-normalisation 
        # TEST
        # Excluding the non-gaseous species
        if vulcan_cfg.use_condense == True:
            exc_conden = [_ for _ in range(ni) if species[_] not in vulcan_cfg.non_gas_sp]
            ysum = np.sum(y_ini[:,exc_conden], axis=1).reshape((-1,1))
        else:
            ysum = np.sum(y_ini, axis=1).reshape((-1,1))
        
        data_var.y_ini = np.copy(y_ini)
        data_var.ymix = y_ini/ysum
        
        if vulcan_cfg.use_ion == True: 
            # if the charge_list is empty (no species with nonzero charges include)
            if not charge_list: 
                print ( "vulcan_cfg.use_ion = True but the network with ions is not supplied.\n" )
                raise IOError("vulcan_cfg.use_ion = True but the network with ions is not supplied.\n")
            else:
                if 'e' in charge_list: charge_list.remove('e') 
                data_var.charge_list = charge_list
   
        return data_var
        


    def ele_sum(self, data_var):
        
        for atom in self.atom_list:
            data_var.atom_ini[atom] = np.sum([compo[compo_row.index(species[i])][atom] * data_var.y[:,i] for i in range(ni)])
            data_var.atom_loss[atom] = 0.
            data_var.atom_conden[atom] = 0.
            
        return data_var


class Atm(object):
    
    def __init__(self):
        self.gs = vulcan_cfg.gs # gravity
        self.P_b = vulcan_cfg.P_b
        self.P_t = vulcan_cfg.P_t
        self.type = vulcan_cfg.atm_type
        self.use_Kzz = vulcan_cfg.use_Kzz
        self.Kzz_prof = vulcan_cfg.Kzz_prof
        self.const_Kzz = vulcan_cfg.const_Kzz
        self.use_vz = vulcan_cfg.use_vz
        self.vz_prof = vulcan_cfg.vz_prof
        self.const_vz = vulcan_cfg.const_vz 
        self.use_settling = vulcan_cfg.use_settling 
        self.non_gas_sp = vulcan_cfg.non_gas_sp      
        
    def f_pico(self, data_atm):
        '''calculating the pressure at interface'''
        
        pco = data_atm.pco
        
        # construct pico
        pco_up1 = np.roll(pco,1)
        pi = (pco * pco_up1)**0.5
        pi[0] = pco[0]**1.5 * pco[1]**(-0.5)
        pi = np.append(pi,pco[-1]**1.5 * pco[-2]**(-0.5))
        
        # store pico
        data_atm.pico = pi
        #data_atm.pco = pco
             
        return data_atm
    
    
    def load_TPK(self, data_atm):
        
        PTK_fun = {}
        
        # IF switches for TP types
        if self.type == 'isothermal': 
            data_atm.Tco = np.repeat(vulcan_cfg.Tiso,nz)
            #data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
            #data_atm.vz = np.repeat(self.const_vz,nz-1)
            
        elif self.type == 'analytical': 
            
            # plotting T-P on the fly                               
            para_atm = vulcan_cfg.para_anaTP 
            
            # return the P-T function
            PTK_fun['pT'] = lambda pressure: self.TP_H14(pressure, *para_atm)        
            data_atm.Tco = PTK_fun['pT'](data_atm.pco)
            #data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
            #data_atm.vz = np.repeat(self.const_vz,nz-1)
        
        # for atm_type = 'file' and also Kzz_prof = 'file
        elif self.type == 'file':
            
            if self.Kzz_prof == 'file':
                atm_table = np.genfromtxt(vulcan_cfg.atm_file, names=True, dtype=None, skip_header=1)
                p_file, T_file, Kzz_file = atm_table['Pressure'], atm_table['Temp'], atm_table['Kzz']
            
            else:     
                atm_table = np.genfromtxt(vulcan_cfg.atm_file, names=True, dtype=None, skip_header=1)
                p_file, T_file = atm_table['Pressure'], atm_table['Temp']

            if max(p_file) < data_atm.pco[0] or min(p_file) > data_atm.pco[-1]:
                print ('Warning: P_b and P_t assgined in vulcan.cfg are out of range of the input.\nConstant extension is used.')
            
            PTK_fun['pT'] = interpolate.interp1d(p_file, T_file, assume_sorted = False, bounds_error=False,\
             fill_value=(T_file[np.argmin(p_file)], T_file[np.argmax(p_file)] )  )
            # store Tco in data_atm
            try:
                data_atm.Tco = PTK_fun['pT'](data_atm.pco)
                
            # for SciPy earlier than v0.18.0
            except ValueError:
                PTK_fun['pT'] = interpolate.interp1d(p_file, T_file, assume_sorted = False, bounds_error=False, fill_value=T_file[np.argmin(p_file)] )  
                data_atm.Tco = PTK_fun['pT'](data_atm.pco)
            
            
            if self.use_Kzz == True and self.Kzz_prof == 'file': 
                
                PTK_fun['pK'] = interpolate.interp1d(p_file, Kzz_file, assume_sorted = False, bounds_error=False, fill_value=(Kzz_file[np.argmin(p_file)], Kzz_file[np.argmax(p_file)]) )
                # store Kzz in data_atm
                try:
                    data_atm.Kzz = PTK_fun['pK'](data_atm.pico[1:-1])
                # for SciPy earlier than v0.18.0
                except ValueError:
                    PTK_fun['pK'] = interpolate.interp1d(p_file, Kzz_file, assume_sorted = False, bounds_error=False, fill_value=Kzz_file[np.argmin(p_file)] )
                    data_atm.Kzz = PTK_fun['pK'](data_atm.pico[1:-1])
            
            elif self.Kzz_prof == 'const': data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
        
        elif self.type == 'vulcan_ini':
            print ("Initializing PT from the prvious run " + vulcan_cfg.vul_ini)
            with open(vulcan_cfg.vul_ini, 'rb') as handle:
              vul_data = pickle.load(handle) 
            
            data_atm.Tco = vul_data['atm']['Tco']
            
        elif self.type == 'table':
            print ("Initializing PT from the prvious run " + vulcan_cfg.vul_ini)
            table = np.genfromtxt(vulcan_cfg.vul_ini, names=True, dtype=None, skip_header=1) 
            if not len(data_atm.pco) == len(table['Pressure']): 
                print ("Warning! The initial profile has different layers than the current setting...")
                raise    
            data_atm.pco = table['Pressure']     
            data_atm.Tco = table['Temp']
              
        else: raise IOError ('\n"atm_type" cannot be recongized.\nPlease trassign it in vulcan_cfg.')
        
        # IF switches for Kzz types
        if self.Kzz_prof == 'const': data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
        elif self.Kzz_prof == 'JM16': # Kzz profiles assumed in Moses et al.2016
            data_atm.Kzz = 1e5 * (300./(data_atm.pico[1:-1]*1e-3))**0.5
            data_atm.Kzz = np.maximum(vulcan_cfg.K_deep, data_atm.Kzz)
        elif self.Kzz_prof == 'Pfunc': # Kzz profiles assumed in Tsai 2020
            data_atm.Kzz = vulcan_cfg.K_max * (vulcan_cfg.K_p_lev*1e6 /(data_atm.pico[1:-1]))**0.4
            data_atm.Kzz = np.maximum(vulcan_cfg.K_max, data_atm.Kzz)
        
        elif self.Kzz_prof == 'file': pass # already defined within atm_type = 'file     
        else: raise IOError ('\n"Kzz_prof" (the type of Kzz profile) cannot be recongized.\nPlease assign it as "file" or "const" or "JM16" in vulcan_cfg.')
        
        # IF switches for Vz types
        if self.vz_prof == 'const': data_atm.vz = np.repeat(self.const_vz,nz-1)
        elif self.vz_prof == 'file': 
            inter_vz = interpolate.interp1d( atm_table['Pressure'], atm_table['vz'], assume_sorted = False, bounds_error=False, fill_value=0 )
            data_atm.vz =  inter_vz(data_atm.pico[1:-1])
        else: raise IOError ('\n"vz_prof" cannot be recongized.\nPlease assign it as "file" or "const" in vulcan_cfg.')
        
        
        
                        
        if self.use_Kzz == False:
            # store Kzz in data_atm
            data_atm.Kzz = np.zeros(nz-1)
        if self.use_vz == False: 
            data_atm.vz = np.zeros(nz-1)   
        
        # TEST 
        # moved to calculating g
        
        # if self.use_settling == True:
        # # TESTing settling velocity
        # # based on L. D. Cloutman: A Database of Selected Transport Coefficients for Combustion Studies (Table 1.)
        #     if vulcan_cfg.atm_base == 'N2':
        #         na = 1.52; a = 1.186e-5; b = 86.54
        #     elif vulcan_cfg.atm_base == 'H2':
        #         na = 1.67; a = 1.936e-6; b = 2.187
        #     elif vulcan_cfg.atm_base == 'CO2':
        #         print ("NO CO2 viscosity yet! (using N2 instead)")
        #         na = 1.52; a = 1.186e-5; b = 86.54
        #     elif vulcan_cfg.atm_base == 'H2O':
        #         na = 1.5; a = 1.6e-5; b = 0
        #     elif vulcan_cfg.atm_base == 'O2':
        #         na = 1.46; a = 2.294e-5; b = 164.4
        #
        #     dmu = a * data_atm.Tco**na /(b + data_atm.Tco) # g cm-1 s-1 dynamic viscosity
        #
        #     for sp in vulcan_cfg.non_gas_sp:
        #         try:
        #             rho_p = data_atm.rho_p[sp]
        #             r_p = data_atm.r_p[sp]
        #
        #        # if sp == 'H2O_l_s':
        #        #      rho_p = data_atm.rho_p_h2o
        #        #      r_p = data_atm.r_p_h2o
        #        #  elif sp == 'H2SO4_l':
        #        #      rho_p = 1.8302
        #        #      r_p = data_atm.r_p_h2so4
        #        #  elif sp == 'H2SO4_l':
        #        #      rho_p = 1.8302
        #        #      r_p = data_atm.r_p_h2so4
        #
        #         except: print (sp + " has not been prescribed size and density!");raise
        #
        #         # Calculating the setteling (terminal) velocity
        #         data_atm.vs[:,species.index(sp)] = -1. *(2./9*rho_p * r_p**2 * data_atm.g / dmu[1:])
        
        # calculating and storing M(the third body)
        data_atm.M = data_atm.pco/(kb*data_atm.Tco)
        data_atm.n_0 = data_atm.M.copy()
        
        # moved to f_mu_dz()
        # # plot T-P profile
        # if vulcan_cfg.plot_TP == True: output.plot_TP(data_atm)
        #
        # # print warning when T exceeds the valid range of Gibbs free energy (NASA polynomials)
        # if np.any(np.logical_or(data_atm.Tco < 200, data_atm.Tco > 6000)): print ('Temperatures exceed the valid range of Gibbs free energy.\n')
        
        return data_atm
        
        
    # T(P) profile in Heng et al. 2014 (126)
    def TP_H14(self, pco, *args_analytical):
        
        # convert args_analytical tuple to a list so we can modify it
        T_int, T_irr, ka_0, ka_s, beta_s, beta_l = list(args_analytical) 
        
        g = vulcan_cfg.gs
        P_b = vulcan_cfg.P_b 
     
        # albedo(beta_s) also affects T_irr
        albedo = (1.0-beta_s)/(1.0+beta_s)
        T_irr *= (1-albedo)**0.25    
        eps_L = 3./8; eps_L3=1./3; ka_CIA=0
        m = pco/g; m_0 = P_b/g; ka_l = ka_0 + ka_CIA*m/m_0
        term1 = T_int**4/4*(1/eps_L + m/(eps_L3*beta_l**2)*(ka_0 + ka_CIA*m/(2*m_0) ) )
        term2 = (1/(2*eps_L) + scipy.special.expn(2,ka_s*m/beta_s)*(ka_s/(ka_l*beta_s)- (ka_CIA)*m*beta_s/(eps_L3*ka_s*m_0*beta_l**2) ) )
        term3 = ka_0*beta_s/(eps_L3*ka_s*beta_l**2)*(1./3 - scipy.special.expn(4,ka_s*m/beta_s))
        term4 = 0. #related to CIA
        T = (term1 + T_irr**4/8*(term2 + term3 + term4) )**0.25
    
        return T
    
        
    def mol_mass(self, sp):
        ''' calculating the molar mass of each species'''
        return compo['mass'][compo_row.index(sp)]

    def mean_mass(self, var, atm, ni):
        mu = np.zeros(nz)
        for i in range(ni):
            mu += self.mol_mass(species[i]) * var.ymix[:,i]
        atm.mu = mu
        return atm        
        
    def f_mu_dz(self, data_var, data_atm, output): # Initilising mean molecular weight and dz 
            
        dz, zco = data_atm.dz, data_atm.zco # pressure defined at interfaces
        Tco, pico = data_atm.Tco.copy(), data_atm.pico.copy()
        Hp = data_atm.Hp
        
        if vulcan_cfg.rocky == False and self.P_b >= 1e6: # if the lower BC greater than 1bar for gas giants
            # Find the index of pico closest to 1bar
            pref_indx = min( range(nz+1), key=lambda i: abs(np.log10(pico[i])-6.))
        else: pref_indx = 0
        # print ("g_s starts from " + str(pico[pref_indx]/1e6) + " bar")
        
        # updating and storing mu
        data_atm = self.mean_mass(data_var, data_atm, ni)
        mu = data_atm.mu
        gs = self.gs
        gz = data_atm.g
        # updating and storing mu
        data_atm = self.mean_mass(data_var, data_atm, ni)
        
        for i in range(pref_indx,nz):
            if i == pref_indx:
                gz[i] = gs
                Hp[i] = kb*Tco[i]/(mu[i]/Navo*gs)    
            else:
                gz[i] = gs * (vulcan_cfg.Rp/(vulcan_cfg.Rp+ zco[i]))**2
                Hp[i] = kb*Tco[i]/(mu[i]/Navo*gz[i])
            dz[i] = Hp[i] * np.log(pico[i]/pico[i+1]) # pico[i+1] has a lower P than pico[i] (higer height)
            zco[i+1] = zco[i] + dz[i] # zco is set zero at 1bar for gas giants

        # for pref_indx != zero 
        if not pref_indx == 0:
            for i in range(pref_indx-1,-1,-1):
                gz[i] = gs * (vulcan_cfg.Rp/(vulcan_cfg.Rp + zco[i+1]))**2
                Hp[i] = kb*Tco[i]/(mu[i]/Navo*gz[i])
                dz[i] = Hp[i] * np.log(pico[i]/pico[i+1]) 
                zco[i] = zco[i+1] - dz[i] # from i+1 propogating down to i
            
        zmco = 0.5*(zco + np.roll(zco,-1))
        zmco = zmco[:-1]
        dzi = 0.5*(dz + np.roll(dz,1))
        dzi = dzi[1:]
        # for the j grid, dzi[j] from the grid above and dz[j-1] from the grid below
        
        # for the molecular diffsuion
        if vulcan_cfg.use_moldiff == True:
            Ti = 0.5*(Tco + np.roll(Tco,-1))
            data_atm.Ti = Ti[:-1]
            Hpi = 0.5*(Hp + np.roll(Hp,-1))
            data_atm.Hpi = Hpi[:-1]
            
        # updating and storing dz and dzi
        #data_atm.dz = dz
        data_atm.dzi = dzi
        #data_atm.zco = zco
        data_atm.zmco = zmco
        data_atm.g, data_atm.gs = gz, gs 
        data_atm.pref_indx = pref_indx
                    
        if self.use_settling == True:
        # TESTing settling velocity
        # based on L. D. Cloutman: A Database of Selected Transport Coefficients for Combustion Studies (Table 1.)
            if vulcan_cfg.atm_base == 'N2':
                na = 1.52; a = 1.186e-5; b = 86.54
            elif vulcan_cfg.atm_base == 'H2':
                na = 1.67; a = 1.936e-6; b = 2.187
            elif vulcan_cfg.atm_base == 'CO2':
                print ("NO CO2 viscosity yet! (using N2 instead)")
                na = 1.52; a = 1.186e-5; b = 86.54
            elif vulcan_cfg.atm_base == 'H2O':
                na = 1.5; a = 1.6e-5; b = 0
            elif vulcan_cfg.atm_base == 'O2':
                na = 1.46; a = 2.294e-5; b = 164.4
                
            dmu = a * data_atm.Tco**na /(b + data_atm.Tco) # g cm-1 s-1 dynamic viscosity
            
            for sp in vulcan_cfg.non_gas_sp:
                try:
                    rho_p = data_atm.rho_p[sp]
                    r_p = data_atm.r_p[sp]
                 
                except: print (sp + " has not been prescribed size and density!");raise
                 
                # Calculating the setteling (terminal) velocity
                gi = 0.5*(data_atm.g + np.roll(data_atm.g,-1))
                gi = gi[:-1]        
                data_atm.vs[:,species.index(sp)] = -1. *(2./9*rho_p * r_p**2 * gi / dmu[1:])
        
        # plot T-P profile
        if vulcan_cfg.plot_TP == True: output.plot_TP(data_atm)
        # print warning when T exceeds the valid range of Gibbs free energy (NASA polynomials)
        if np.any(np.logical_or(data_atm.Tco < 200, data_atm.Tco > 6000)): print ('Temperatures exceed the valid range of Gibbs free energy.\n')
            
        return data_atm
        
    def read_sflux(self, var, atm):
        '''reading in stellar stpactal flux at the stellar surface and converting it to the flux on the planet to the uniform grid using trapezoidal integral'''
        atm.sflux_raw = np.genfromtxt(vulcan_cfg.sflux_file, dtype=float, skip_header=1, names = ['lambda','flux'])
        
        # for values outside the boundary => fill_value = 0
        bins = var.bins
        
        dbin1 = vulcan_cfg.dbin1
        dbin2 = vulcan_cfg.dbin2
        
        inter_sflux = interpolate.interp1d(atm.sflux_raw['lambda'], atm.sflux_raw['flux']* (vulcan_cfg.r_star*r_sun/(au*vulcan_cfg.orbit_radius) )**2, bounds_error=False, fill_value=0)
        for n, ld in enumerate(var.bins):
            var.sflux_top[n] = inter_sflux(ld) 
            if ld == vulcan_cfg.dbin_12trans: var.sflux_din12_indx = n
            # not converting to actinic flux yet *1/(hc/ld)
        
        # Check for energy conservation
        # finding the index for the left & right pts that match var.bins in the raw data
        raw_flux = atm.sflux_raw['flux']* (vulcan_cfg.r_star*r_sun/(au*vulcan_cfg.orbit_radius) )**2
        raw_left_indx = np.searchsorted(atm.sflux_raw['lambda'],bins[0],side='right')        
        raw_right_indx = np.searchsorted(atm.sflux_raw['lambda'],bins[-1],side='right')-1

        sum_orgin = 0        
        # for checking the trapezoidal error in energy conservation
        for n in range(raw_left_indx,raw_right_indx):
            sum_orgin += 0.5*(raw_flux[n] + raw_flux[n+1]) * (atm.sflux_raw['lambda'][n+1]- atm.sflux_raw['lambda'][n])
        sum_orgin += 0.5 *(inter_sflux(bins[0])+raw_flux[raw_left_indx])* (atm.sflux_raw['lambda'][raw_left_indx]-bins[0])
        sum_orgin += 0.5 *(inter_sflux(bins[-1])+raw_flux[raw_right_indx])* (bins[-1]-atm.sflux_raw['lambda'][raw_right_indx])
        
        # dbin_12trans outside the bins
        if not 'sflux_din12_indx' in vars(var).keys(): var.sflux_din12_indx = -1
        
        sum_bin = dbin1 * np.sum(var.sflux_top[:var.sflux_din12_indx])
        sum_bin -= dbin1 *0.5*(var.sflux_top[0]+var.sflux_top[var.sflux_din12_indx-1])
        sum_bin += dbin2 * np.sum(var.sflux_top[var.sflux_din12_indx:])
        sum_bin -= dbin2 *0.5*(var.sflux_top[var.sflux_din12_indx]+var.sflux_top[-1])
         
        print ("The stellar flux is interpolated onto uniform grid of " +str(vulcan_cfg.dbin1) + " (<" +str(vulcan_cfg.dbin_12trans)+" nm) and "+str(vulcan_cfg.dbin2)\
        + " (>="+str(vulcan_cfg.dbin_12trans)+" nm)" + " and conserving " + "{:.2f}".format(100* sum_bin/sum_orgin)+" %" + " energy." )
        #print (str(100* sum_old/sum_orgin)+" %" )

        
    
    def mol_diff(self, atm):
        '''
        choosing the formulea of molecular diffusion for each species
        then constucting Dzz(z) 
        '''
        Tco = atm.Tco
        n_0 = atm.n_0 
        
        # using the value defined on the interface
        Tco_i = np.delete((Tco + np.roll(Tco,1))*0.5, 0)
        n0_i = np.delete((n_0 + np.roll(n_0,1))*0.5, 0)
        
        if vulcan_cfg.use_moldiff == False:
            for i in range(len(species)):
                # this is required even without molecular weight
                atm.ms[i] = compo[compo_row.index(species[i])][-1]
            return
        
        if vulcan_cfg.atm_base == 'H2':
            Dzz_gen = lambda T, n_tot, mi: 2.2965E17*T**0.765/n_tot *( 16.04/mi*(mi+2.016)/18.059 )**0.5 # from Moses 2000a
            
            # scaling with (15.27) in [Aeronomy part B by Banks & Kockarts(1973)]
            # *( m_ref/mi*(mi+ m_base)/m_tot )**0.5 (m_ref is the molecular mass of the ref-base e.g. CH4 in CH4-H2 in Moses 2000a)
            
            # # thermal diffusion factor (>0 means (heavier) components diffuse toward colder rigions)
            if 'H' in species: atm.alpha[species.index('H')] = -0.1 # simplified from Moses 2000a
            if 'He' in species:  atm.alpha[species.index('He')] = 0.145
            for sp in species:
                if self.mol_mass(sp) > 4.: atm.alpha[species.index(sp)] = 0.25
            
        elif vulcan_cfg.atm_base == 'N2': # use CH4-N2 in Aeronomy [Banks ] as a reference to scale by the molecular mass
            Dzz_gen = lambda T, n_tot, mi: 7.34E16*T**0.75/n_tot *( 16.04/mi*(mi+28.014)/44.054 )**0.5
            
            # # thermal diffusion factor (>0 means (heavier) components diffuse toward colder rigions)
            if 'H' in species: atm.alpha[species.index('H')] = -0.25
            if 'H2' in species:  atm.alpha[species.index('H2')] = -0.25
            if 'He' in species:  atm.alpha[species.index('He')] = -0.25
            if 'Ar' in species:  atm.alpha[species.index('Ar')] = 0.17
            
        elif vulcan_cfg.atm_base == 'O2': # use CH4-O2 in Aeronomy [Banks ] as a reference to scale by the molecular mass
            Dzz_gen = lambda T, n_tot, mi: 7.51E16*T**0.759/n_tot *( 16.04/mi*(mi+32)/48.04 )**0.5
            
            # # thermal diffusion factor (>0 means (heavier) components diffuse toward colder rigions)
            if 'H' in species: atm.alpha[species.index('H')] = -0.25
            if 'H2' in species:  atm.alpha[species.index('H2')] = -0.25
            if 'He' in species:  atm.alpha[species.index('He')] = -0.25
            if 'Ar' in species:  atm.alpha[species.index('Ar')] = 0.17
            
        elif vulcan_cfg.atm_base == 'CO2': # use H2-CO2 in Hu seager as a reference to scale by the molecular mass
            Dzz_gen = lambda T, n_tot, mi: 2.15E17*T**0.750/n_tot *( 2.016/mi*(mi+44.001)/46.017 )**0.5
            
            # # thermal diffusion factor (>0 means (heavier) components diffuse toward colder rigions)
            if 'H' in species: atm.alpha[species.index('H')] = -0.25
            if 'H2' in species:  atm.alpha[species.index('H2')] = -0.25
            if 'He' in species:  atm.alpha[species.index('He')] = -0.25
            if 'Ar' in species:  atm.alpha[species.index('Ar')] = 0.17
            
        else: raise IOError ('\n Unknow atm_base!')
        
        for i in range(len(species)):
            # input should be float or in the form of nz-long 1D array
            atm.Dzz[:,i] = Dzz_gen(Tco_i, n0_i, self.mol_mass(species[i]))
            
            # constructing the molecular weight for every species
            # this is required even without molecular weight
            atm.ms[i] = compo[compo_row.index(species[i])][-1]
        
        # setting the molecuar diffusion of the non-gaseous species to zero
        for sp in [_ for _ in vulcan_cfg.non_gas_sp if _ in species]: atm.Dzz[:,species.index(sp)] = 0
                
    
    def BC_flux(self, atm):
        '''
        Reading-in the boundary conditions of constant flux (cm^-2 s^-1) at top/bottom
        '''
        # read in the const top BC
        if vulcan_cfg.use_topflux == True: 
            print ("Using the prescribed constant top flux.")
            with open (vulcan_cfg.top_BC_flux_file) as f:
                for line in f.readlines():
                    if not line.startswith("#") and line.strip():
                        li = line.split()                   
                        atm.top_flux[species.index(li[0])] = li[1]
        
        # read in the const bottom BC
        if vulcan_cfg.use_botflux == True: 
            print ("Using the prescribed constant bottom flux.")
            with open (vulcan_cfg.bot_BC_flux_file) as f:
                for line in f.readlines():
                    if not line.startswith("#") and line.strip():
                        li = line.split()                   
                        atm.bot_flux[species.index(li[0])] = li[1]
                        atm.bot_vdep[species.index(li[0])] = li[2]
                        
        # using fixed-mixing-ratio BC          
        if vulcan_cfg.use_fix_sp_bot == True: 
            print ("Using the prescribed fixed bottom mixing ratios.")
            with open (vulcan_cfg.bot_BC_flux_file) as f:
                for line in f.readlines():
                    if not line.startswith("#") and line.strip():
                        li = line.split()                   
                        atm.bot_fix_sp[species.index(li[0])] = li[3]
    
    # TEST condensation
    def sp_sat(self, atm):
        '''
        For all the species in vulcan_cfg.condense_sp, pre-calculating the  
        saturation varpor pressure (in dyne/cm2) and storing in atm.sat_p. 
        '''
        # the list that the data has been coded 
        sat_sp_list = ['H2O','NH3','H2SO4','S2','S8' ,'C','H2S' ]
        
        for sp in vulcan_cfg.condense_sp:
            if sp not in sat_sp_list: raise IOError ( "No saturation vapor data for " +sp + ". Check the sp_sat function in build_atm.py" )
            
            T = np.copy(atm.Tco)
            
            if sp == "H2O":
                # T is in C
                T -= 273.
                # from Seinfeld & Pandis 2006, P in mbar in the book
                a_water = (6.107799961, 4.436518521E-1, 1.428945805E-2, 2.650648471E-4, 3.031240396E-6, 2.034080948E-8, 6.136820929E-11)
                a_ice = (6.109177956, 5.034698970E-1, 1.886013408E-2, 4.176223716E-4, 5.824720280E-6, 4.838803174E-8, 1.838826904E-10)

                # saturate_p_1 = (T<0)*( a_ice[0] + a_ice[1]*T + a_ice[2]*T**2 + a_ice[3]*T**3 + a_ice[4]*T**4 + a_ice[5]*T**5 + a_ice[6]*T**6 ) +\
                #  (T>0)*(a_water[0] + a_water[1]*T + a_water[2]*T**2 + a_water[3]*T**3 + a_water[4]*T**4 + a_water[5]*T**5 + a_water[6]*T**6)
                
                # ice from Ackerman&Marley (2001)
                c0 = 6111.5; c1 = 23.036; c2 = -333.7; c3 = 279.82
                # water from Ackerman&Marley (2001)
                w0 = 6112.1; w1 = 18.729; w2 = -227.3; w3 = 257.87
    
                saturate_p = (T<0)*( c0 * np.exp( (c1*T + T**2/c2)/(T + c3) ) ) + (T>0)*(w0 * np.exp( (w1*T + T**2/w2)/(T + w3) ) )
                
                atm.sat_p[sp] = saturate_p
                #atm.sat_p[sp] = c0 * np.exp( (c1*T + T**2/c2)/(T + c3) )
                #atm.sat_p[sp] = moses_ice
            
            elif sp == "NH3":
                # from Weast (1971) in bar
                c0 = 10.53; c1 = -2161.0;  c2 = -86596.0
                saturate_p = np.exp(c0 + c1/T + c2/T**2)
                atm.sat_p[sp] = saturate_p * 1.e6
            
            elif sp == "H2SO4":
                # change to Kulmala later
                p_ayers = np.e**(-10156./T + 16.259) # in atm
                atm.sat_p[sp] = p_ayers * 1.01325*1e6 # arm to cgs
                
            elif sp == "S2":
                atm.sat_p[sp] = np.zeros(nz)
                # from Zahnle 2017 (refitted from Lyons 2008)
                atm.sat_p[sp][T<413] = np.exp(27. - 18500./T[T<413]) *1e6 # in bar => cgs
                atm.sat_p[sp][T>=413] = np.exp(16.1 - 14000./T[T>=413]) *1e6
            
            elif sp == "S4":
                atm.sat_p[sp] = np.zeros(nz)
                # from Lyons 2008
                atm.sat_p[sp] = 10**(6.0028 -6047.5/T) *1.01325e6 # atm to vgs
            
            elif sp == "S8":
                atm.sat_p[sp] = np.zeros(nz)
                # from Zahnle 2017 (refitted from Lyons 2008)
                atm.sat_p[sp][T<413] = np.exp(20. - 11800./T[T<413]) *1e6 # in bar => cgs
                atm.sat_p[sp][T>=413] = np.exp(9.6 - 7510./T[T>=413]) *1e6
            
            elif sp == "C":
                atm.sat_p[sp] = np.zeros(nz)
                # from NIST (dyna cm^-2)
                a = 3.27860E+01
                b = -8.65139E+04
                c = 4.80395E-01 
                atm.sat_p[sp] = np.exp(a+b/(atm.Tco +c) )
                
            elif sp == "H2S":
                # from Giauque and Blue(1936) in cmHg (adoped in Atreya's book)
                h2s_ice_log10 = -1329./T + 9.28588 - 0.0051263*T # 164.9 <= T <= 187.6
                h2s_l_log10 = -1145./T + 7.94746 - 0.00322*T     # 187.6 < T <= 213.2
                
                saturate_p = 10**( (T <= 187.6)*h2s_ice_log10 + (T > 187.6)*h2s_l_log10 )
                atm.sat_p[sp] = saturate_p * 0.001333 * 1.e6

            
if __name__ == "__main__":
    print("This module stores classes for constructing atmospheric structure \
    and initializing its chemical composition from the desinated elemental abudance.")