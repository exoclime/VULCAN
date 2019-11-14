import numpy as np
from numpy import polynomial
import scipy
from scipy import interpolate
import scipy.optimize as sop
import subprocess
import pickle

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
        #with open('fastchem_vulcan/chemistry/elements/element_abundances_lodders.dat' ,'r') as f:
        with open('fastchem_vulcan/input/element_abundances_lodders.dat' ,'r') as f:
            new_str = ""
            ele_list = list(vulcan_cfg.atom_list)
            ele_list.remove('H')
            
            if vulcan_cfg.use_solar == True: 
                new_str = f.read() # read in as a string
                print ("Initializing with the default solar abundance.")
                
            else: # using costomized elemental abundances
                print ("Initializing with the customized elemental abundance:")
                print ("{:4}".format('H') + str('1.'))
                for line in f.readlines():   
                        li = line.split()
                        if li[0] in ele_list:
                            sp = li[0].strip()
                            # read-in vulcan_cfg.sp_H
                            sp_abun = getattr(vulcan_cfg, sp+'_H')
                            fc_abun = 12. + np.log10(sp_abun)
                            line = sp + '\t' + "{0:.4f}".format(fc_abun) + '\n'
                            print ("{:4}".format(sp) + "{0:.4E}".format(sp_abun))
                        new_str += line
            
            # make the new elemental abundance file
            #with open('fastchem_vulcan/chemistry/elements/element_abundances_vulcan.dat', 'w') as f: f.write(new_str)
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
      
        if vulcan_cfg.ini_mix == 'EQ':
        
            self.ini_fc(data_var, data_atm)
            fc = np.genfromtxt('fastchem_vulcan/output/vulcan_EQ.dat', names=True, dtype=None, skip_header=0)
            neutral_sp = [sp for sp in species if sp not in vulcan_cfg.excit_sp]

            for sp in neutral_sp:
                if sp in fc.dtype.names:
                    y_ini[:,species.index(sp)] = fc[sp]*gas_tot
                else: print (sp + ' not included in fastchem.')
            
            # remove the fc output
            subprocess.call(["rm vulcan_EQ.dat"], shell=True, cwd='fastchem_vulcan/output/')
        
        elif vulcan_cfg.ini_mix == 'fc_precal':
            
            pre_fc = 'fastchem_vulcan/output/vulcan_EQ_pre.dat'
            print ('\n Using the precalculated fastchem output: '+ pre_fc)
            
            fc = np.genfromtxt(pre_fc, names=True, dtype=None, skip_header=0)   
            for sp in species:
                y_ini[:,species.index(sp)] = fc[sp]*gas_tot
        
        elif vulcan_cfg.ini_mix == 'vulcan_ini':
            with open(vulcan_cfg.vul_ini, 'rb') as handle:
              vul_data = pickle.load(handle) 
            
            y_ini = np.copy(vul_data['variable']['y'])
        
        elif vulcan_cfg.ini_mix == 'const_mix':
            for sp in vulcan_cfg.const_mix.keys():
                y_ini[:,species.index(sp)] = gas_tot* vulcan_cfg.const_mix[sp]

            
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
            for sp in vulcan_cfg.condesne_sp:
                data_atm.sat_mix[sp] = data_atm.sat_p[sp]/data_atm.pco
                
                #if not vulcan_cfg.ini_mix == 'vulcan_ini':
                    # TEST Cold trap
                    # np.argmax stops at the first TRUE value and returns its index
                conden_lev = np.argmax( data_atm.n_0*data_atm.sat_mix[sp] <= data_var.y[:,species.index(sp)] )
                if not conden_lev == 0: 
                    print ( sp + " condensed from nz = " + str(conden_lev) + " (cold trap)")
                    data_var.y[conden_lev:,species.index(sp)] = 0
                    data_var.y[:,species.index(sp)] = np.minimum(data_atm.n_0 * data_atm.sat_mix[sp], data_var.y[:,species.index(sp)])
                else:
                    print ( sp + " should be condensed from the surface!")
        
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
        
        return data_var
        


    def ele_sum(self, data_var):
        
        for atom in self.atom_list:
            data_var.atom_ini[atom] = np.sum([compo[compo_row.index(species[i])][atom] * data_var.y[:,i] for i in range(ni)])
            data_var.atom_loss[atom] = 0.
            data_var.atom_conden[atom] = 0.
            
        return data_var


class Atm(object):
    
    def __init__(self):
        self.g = vulcan_cfg.g # gravity
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
        
        if self.type == 'isothermal': 
            data_atm.Tco = np.repeat(vulcan_cfg.Tiso,nz)
            data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
            data_atm.vz = np.repeat(self.const_vz,nz-1)
            
        elif self.type == 'analytical': 
            
            # plotting T-P on the fly                               
            para_atm = vulcan_cfg.para_anaTP 
            
            # return the P-T function
            PTK_fun['pT'] = lambda pressure: self.TP_H14(pressure, *para_atm)        
            data_atm.Tco = PTK_fun['pT'](data_atm.pco)
            data_atm.Kzz = np.repeat(self.const_Kzz,nz-1)
            data_atm.vz = np.repeat(self.const_vz,nz-1)
            
        elif self.type == 'file':
            
            if self.Kzz_prof == 'const':     
                atm_table = np.genfromtxt(vulcan_cfg.atm_file, names=True, dtype=None, skip_header=1)
                p_file, T_file = atm_table['Pressure'], atm_table['Temp']
            
            elif self.Kzz_prof == 'file':
                atm_table = np.genfromtxt(vulcan_cfg.atm_file, names=True, dtype=None, skip_header=1)
                p_file, T_file, Kzz_file = atm_table['Pressure'], atm_table['Temp'], atm_table['Kzz']
            else: raise IOError ('\n"Kzz_prof" (the type of Kzz profile) cannot be recongized.\nPlease assign it as "file" or "const" in vulcan_cfg.')

            if self.vz_prof == 'const': data_atm.vz = np.repeat(self.const_vz,nz-1)
            elif self.vz_prof == 'file': vz_file =  atm_table['vz']


            if max(p_file) < data_atm.pco[0] or min(p_file) > data_atm.pco[-1]:
                print ('Warning! P_b and P_t assgined in vulcan.cfg are out of range of the file.\nConstant extension will be used.')
            
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
        
        else: raise IOError ('\n"atm_type" cannot be recongized.\nPlease trassign it in vulcan_cfg.')
                        
        if self.use_Kzz == False:
            # store Kzz in data_atm
            data_atm.Kzz = np.zeros(nz-1)
        if self.use_vz == False: 
            data_atm.vz = np.zeros(nz-1)   
        
        # if self.use_settling == True:
        # TESTing settling velocity
        # def vs(self, var, atm):
            # dmu = 1.729e-4 # the dynamics viscosity of air
            # for sp in self.non_gas_sp:
            #     if sp == 'H2O_l_s':
            #         cs = 1. # the slip correction factor
            #         rho_p = 1.
            #         r_p = 2.e-4
            #         vs = -2./9*r_p**2* rho_p*self.g *cs/dmu # negative for downward
            #         data_atm.vs[:,species.index(sp)] = vs
            
#             # Using Gao 2018 (50)
#             g = self.g
#             r_p = 1.e-4
#             R_uni = kb*Navo # the universal gas const
#
#             Ti = 0.5*(data_atm.Tco + np.roll(data_atm.Tco,-1))
#             Ti = Ti[:-1]
#             mui = 0.5*(data_atm.mu + np.roll(data_atm.mu,-1))
#             mui = mui[:-1]
#             pi = data_atm.pico[1:-1]
#
#             data_atm.mui = mui
#             data_atm.Ti = Ti
#             data_atm.pi = pi
#
#             # mass density of the air
#             rho_a = pi/Ti/(R_uni/mui)
#             for sp in self.non_gas_sp:
#                 if sp == 'H2O_l_s':
#                     rho_p = 1.
#                     vs = 0.5*rho_p*g*r_p / rho_a *(np.pi*mui /(2*R_uni*Ti))**0.5
#                     data_atm.vs[:,species.index(sp)] = vs
        
        
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
        
        g = vulcan_cfg.g
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
        ''' calculating the molar mass of each species?'''
        return compo['mass'][compo_row.index(sp)]

    def mean_mass(self, var, atm, ni):
        mu = np.zeros(nz)
        for i in range(ni):
            mu += self.mol_mass(species[i]) * var.ymix[:,i]
        atm.mu = mu
        return atm        
        
    def f_mu_dz(self, data_var, data_atm, output): # Initilising mean molecular weight and dz 
            
        dz, zco = np.empty(nz), np.zeros(nz+1) # pressure defined at interfaces
        Tco, pico = data_atm.Tco.copy(), data_atm.pico.copy()
        g = self.g
        
        # updating and storing mu
        data_atm = self.mean_mass(data_var, data_atm, ni)
        # updating and storing Hp
        data_atm.Hp = kb*Tco/(data_atm.mu/Navo*g) 

        for i in range(nz):
            dz[i] = data_atm.Hp[i] * np.log(pico[i]/pico[i+1])
            zco[i+1] = zco[i] + dz[i]

        zmco = 0.5*(zco + np.roll(zco,-1))
        zmco = zmco[:-1]
        dzi = 0.5*(dz + np.roll(dz,1))
        dzi = dzi[1:]
        # for the j grid, dzi[j] from the grid above and dz[j-1] from the grid below
        
        # updating and storing dz and dzi
        data_atm.dz = dz
        data_atm.dzi = dzi
        data_atm.zmco = zmco
        
        if self.use_settling == True:
            # Using Gao 2018 (50)
            r_p = 1.e-4
            R_uni = kb*Navo # the universal gas const
        
            Ti = 0.5*(Tco + np.roll(Tco,-1))
            Ti = Ti[:-1]
            mui = 0.5*(data_atm.mu + np.roll(data_atm.mu,-1))
            mui = mui[:-1]
            pi = pico[1:-1]
        
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
        dbin = vulcan_cfg.dbin
        inter_sflux = interpolate.interp1d(atm.sflux_raw['lambda'], atm.sflux_raw['flux']* (vulcan_cfg.r_star*r_sun/(au*vulcan_cfg.orbit_radius) )**2, bounds_error=False, fill_value=0)
        for n, ld in enumerate(var.bins):
            var.sflux_top[n] = inter_sflux(ld) 
            # not converting to actinic flux yet *1/(hc/ld)
            
        # Stellar flux at TOA; not converting to actinic flux yet *1/(hc/ld)
        # for values outside the boundary => fill_value = 0
        # inter_sflux = interpolate.interp1d(atm.sflux_raw['lambda'], atm.sflux_raw['flux'], bounds_error=False, fill_value=0)
#         dbin = vulcan_cfg.dbin
#         bins = var.bins
#         last_bin = bins[-1]
#
#         for n, ld in enumerate(var.bins):
#             # define the next bin in the new uniform grid
#             if ld != bins[-1]: next_ld = bins[n+1]
#
#             if ld in atm.sflux_raw['lambda'] or ld == last_bin or atm.sflux_raw['lambda'][np.searchsorted(atm.sflux_raw['lambda'],ld,side='right')] - atm.sflux_raw['lambda'][np.searchsorted(atm.sflux_raw['lambda'],ld,side='right')-1] > dbin: # when bin coincide with raw_bin or dbin is smaller than the width of raw_bin
#                 var.sflux_top[n] = inter_sflux(ld)
#
#             else:
#                 # finding the index for the left & right pts in the raw data
#                 raw_left_indx = np.searchsorted(atm.sflux_raw['lambda'],ld,side='right')
#                 raw_right_indx = np.searchsorted(atm.sflux_raw['lambda'],next_ld,side='right') - 1
#                 flux_left = inter_sflux(ld)
#                 flux_right = inter_sflux(next_ld)
#
#                 # trapezoid integral
#                 bin_flux = (flux_left + atm.sflux_raw['flux'][raw_left_indx])*0.5*(atm.sflux_raw['lambda'][raw_left_indx]-ld)
#                 bin_flux += (flux_right + atm.sflux_raw['flux'][raw_right_indx])*0.5*(next_ld-atm.sflux_raw['lambda'][raw_right_indx])
#
#                 if raw_right_indx - raw_left_indx > 0:
#                     for raw_i in range(raw_left_indx,raw_right_indx):
#                         bin_flux += (atm.sflux_raw['flux'][raw_i] + atm.sflux_raw['flux'][raw_i+1])*0.5*(atm.sflux_raw['lambda'][raw_i+1] - atm.sflux_raw['lambda'][raw_i])
#
#                 bin_flux = bin_flux/dbin
#                 var.sflux_top[n] = bin_flux
#
#         var.sflux_top *= (vulcan_cfg.r_star*r_sun/(au*vulcan_cfg.orbit_radius) )**2
#
#         # the old direct interpolation
#         sflux_inter = np.zeros(len(bins))
#         for n, ld in enumerate(bins):
#             sflux_inter[n] = inter_sflux(ld)
#         sflux_inter *= (vulcan_cfg.r_star*r_sun/(au*vulcan_cfg.orbit_radius) )**2
        
        # Check for energy conservation
        # finding the index for the left & right pts that match var.bins in the raw data
        raw_flux = atm.sflux_raw['flux']* (vulcan_cfg.r_star*r_sun/(au*vulcan_cfg.orbit_radius) )**2
        raw_left_indx = np.searchsorted(atm.sflux_raw['lambda'],bins[0],side='right')
        raw_right_indx = np.searchsorted(atm.sflux_raw['lambda'],bins[-1],side='right')-1
        #sum_orgin, sum_bin = 0, 0
        sum_orgin = 0
        for n in range(raw_left_indx,raw_right_indx):
            sum_orgin += 0.5*(raw_flux[n] + raw_flux[n+1]) * (atm.sflux_raw['lambda'][n+1]- atm.sflux_raw['lambda'][n])
        sum_orgin += 0.5 *(inter_sflux(bins[0])+raw_flux[raw_left_indx])* (atm.sflux_raw['lambda'][raw_left_indx]-bins[0])
        sum_orgin += 0.5 *(inter_sflux(bins[-1])+raw_flux[raw_right_indx])* (bins[-1]-atm.sflux_raw['lambda'][raw_right_indx])

        
        sum_bin = dbin * np.sum(var.sflux_top)
        sum_bin -= dbin*0.5*(var.sflux_top[0]+var.sflux_top[-1])
        
        
        #sum_old = dbin * np.sum(sflux_inter)
        #sum_old -= dbin*0.5*(sflux_inter[0]+sflux_inter[-1])
        
        print ("The stellar flux is interpolated onto uniform grid of " +str(vulcan_cfg.dbin)+ " nm and conserving " + "{:.2f}".format(100* sum_bin/sum_orgin)+" %" + " energy." )
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
            Dzz_gen = lambda T, n_tot, mi: 2.2965E17*T**0.765/n_tot *( 16.04/mi*(mi+2.016)/18.059 )**0.5
            
            # scaling with (15.27) in [Aeronomy part B by Banks & Kockarts(1973)]
            # *( m_ref/mi*(mi+ m_base)/m_tot )**0.5 (m_ref is the molecular mass of the ref-base e.g. CH4 in CH4-H2 in Moses 2000a)
            
        elif vulcan_cfg.atm_base == 'N2': # use CH4-N2 in Aeronomy [Banks ] as a reference to scale by the molecular mass
            Dzz_gen = lambda T, n_tot, mi: 7.34E16*T**0.75/n_tot *( 16.04/mi*(mi+28.014)/44.054 )**0.5
        elif vulcan_cfg.atm_base == 'O2': # use CH4-O2 in Aeronomy [Banks ] as a reference to scale by the molecular mass
            Dzz_gen = lambda T, n_tot, mi: 7.51E16*T**0.759/n_tot *( 16.04/mi*(mi+32)/48.04 )**0.5
        elif vulcan_cfg.atm_base == 'CO2': # use H2-CO2 in Hu seager as a reference to scale by the molecular mass
            Dzz_gen = lambda T, n_tot, mi: 2.15E17*T**0.750/n_tot *( 2.016/mi*(mi+44.001)/46.017 )**0.5
        else: raise IOError ('\n Unknow atm_base!')
        
        for i in range(len(species)):
            # input should be float or in the form of nz-long 1D array
            atm.Dzz[:,i] = Dzz_gen(Tco_i, n0_i, self.mol_mass(species[i]))
            
            # constructing the molecular weight for every species
            # this is required even without molecular weight
            atm.ms[i] = compo[compo_row.index(species[i])][-1]
        
        # setting the molecuar diffusion of the non-gaseous species to zero
        for sp in [_ for _ in vulcan_cfg.non_gas_sp if _ in species]: atm.Dzz[:,species.index(sp)] = 0
        
        # no exception needed!?    
        # exceptions
        # if vulcan_cfg.atm_base == 'H2':
#             atm.Dzz[:,species.index('H2')] = np.zeros(nz-1)
#         else: raise IOError ('\n Unknow atm_base!')
        
        
        # # thermal diffusion for H and H2
        # atm.alpha[species.index('H')] = -0.25
        # atm.alpha[species.index('H2')] = -0.25
    
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
        For all the species in vulcan_cfg.condesne_sp, pre-calculating the  
        saturation varpor pressure (in dyne/cm2) and storing in atm.sat_p. 
        '''
        # the list that the data has been coded 
        sat_sp_list = ["H2O",'NH3']
        
        for sp in vulcan_cfg.condesne_sp:
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
                
                
                # ice from Ackerman&Marley (2001) in cgs unit
                #c0 = 6111.5; c1 = 23.036; c2 = -333.7; c3 = 279.82
                
                # water from Yaws (1999)
                #c0 = 2.98605e1  c1 = -3.1522e3  c2 = -7.30370e0  c3 = 2.4247e-9  c4 = 1.8090e-6
                
                # ice in Moses 2005
                moses_ice = 10**( 12.537 - 2663.5/(T+273.) ) *10
                
                #saturate_p = (T<0)*( c0 * np.exp( (c1*T + T**2/c2)/(T + c3) ) ) +\
                #(T>0)*(a_water[0] + a_water[1]*T + a_water[2]*T**2 + a_water[3]*T**3 + a_water[4]*T**4 + a_water[5]*T**5 + a_water[6]*T**6)*1.e3
                
                # ice from Ackerman&Marley (2001)
                c0 = 6111.5; c1 = 23.036; c2 = -333.7; c3 = 279.82
                # water from Ackerman&Marley (2001)
                w0 = 6112.1; c1 = 18.729; c2 = -227.3; c3 = 257.87
    
                saturate_p = (T<0)*( c0 * np.exp( (c1*T + T**2/c2)/(T + c3) ) ) + (T>0)*(c0 * np.exp( (c1*T + T**2/c2)/(T + c3) ) )
                
                atm.sat_p[sp] = saturate_p
                #atm.sat_p[sp] = c0 * np.exp( (c1*T + T**2/c2)/(T + c3) )
                #atm.sat_p[sp] = moses_ice
            
            if sp == "NH3":
                # from Weast (1971) in bar
                c0 = 10.53; c1 = -2161.0;  c2 = -86596.0
                saturate_p = np.exp(c0 + c1/T + c2/T**2)
                atm.sat_p[sp] = saturate_p * 1.e6
        

if __name__ == "__main__":
    print("This module stores classes for constructing atmospheric structure \
    and initializing its chemical composition from the desinated elemental abudance.")