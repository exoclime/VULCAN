import numpy as np
from numpy import polynomial
import scipy
from scipy import interpolate
import scipy.optimize as sop
import vulcan_cfg
from phy_const import kb, Navo
from vulcan_cfg import nz
import chem_funs
from chem_funs import ni, nr  # number of species and reactions in the network
species = chem_funs.spec_list

###
compo = np.genfromtxt(vulcan_cfg.com_file,names=True,dtype=None)
compo_row = list(compo['species'])
###


class InitialAbun(object):
    """
    Calculating the appropriate initial mixing ratios with the assigned elemental abundance
    """
    
    def __init__(self):
        self.ini_m = [0.9,0.1,0.,0.] # initial guess
        self.EQ_ini_file = vulcan_cfg.EQ_ini_file
        
        self.atom_list = vulcan_cfg.atom_list

    def abun_lowT(self, x):
        """
        calculating the initial mixing ratios of the following 4 molecules (with CH4) 
        satisfying the assigned elemental abundance
        x1:H2 x2:H2O x3:CH4 x4:He
        """
        O_H, C_H, He_H = vulcan_cfg.O_H, vulcan_cfg.C_H, vulcan_cfg.He_H
        x1,x2,x3,x4 = x
        f1 = x1+x2+x3+x4-1.
        f2 = x2 - (2*x1+2*x2+4*x3)*O_H
        f3 = x3 - (2*x1+2*x2+4*x3)*C_H
        f4 = x4 - (2*x1+2*x2+4*x3)*He_H
        return f1,f2,f3,f4
        
    def abund_highT(self, x):
        """
        calculating the initial mixing ratios of the following 4 molecules (with CO) 
        satisfying the assigned elemental abundance
        x1:H2 x2:H2O x3:CO x4:He
        """
        x1,x2,x3,x4 = x
        f1 = x1+x2+x3+x4-1.
        f2 = x2+x3 - (2*x1+2*x2)*O_H
        f3 = x3 - (2*x1+2*x2)*C_H
        f4 = x4 - (2*x1+2*x2)*He_H
        return f1,f2,f3,f4
    
    def abund_eq9(self, sp_ini, tp, pp):
        '''
        Using the 6-molecule function to solve 5th order polynomial equation for CO mixing ratio (Heng & Tsai 2016).
        Besides the 6 molecules, several additional molecules are non self-consistently "post-processed."
        '''
        
        from chem_funs import gibbs_sp
        
        n_eq, mix_eq = {}, {}
        n_o, n_c, n_n = vulcan_cfg.O_H, vulcan_cfg.C_H, 0.
               
        kk1 = lambda T, pbar: np.exp( -( gibbs_sp('CO',T) + 3*gibbs_sp('H2',T) -1*gibbs_sp('CH4',T)-1*gibbs_sp('H2O',T) ) ) * (pbar)**-2
        kk2 = lambda T: np.exp( -( gibbs_sp('CO',T) + gibbs_sp('H2O',T) -1*gibbs_sp('CO2',T) -1*gibbs_sp('H2',T) ) )
        kk3 = lambda T, pbar: np.exp( -( gibbs_sp('C2H2',T) + 3*gibbs_sp('H2',T) -2*gibbs_sp('CH4',T) ) ) * (pbar)**-2
        kk4 = lambda T, pbar: np.exp( -( gibbs_sp('C2H2',T) + gibbs_sp('H2',T) -gibbs_sp('C2H4',T) ) ) * (pbar)**-1
        kk5 = lambda T, pbar: np.exp( -( gibbs_sp('N2',T) + 3*gibbs_sp('H2',T) -2*gibbs_sp('NH3',T) ) ) * (pbar)**-2
        kk6 = lambda T, pbar: np.exp( -( gibbs_sp('HCN',T) + 3*gibbs_sp('H2',T) -gibbs_sp('NH3',T) -gibbs_sp('CH4',T) ) ) * (pbar)**-2
     
        k1 = kk1(tp,pp)
        k5 = kk5(tp,pp)
        k6 = kk6(tp,pp)
        d2 = 1.0 + 2.0*k1*(n_o+n_c)
        a0 = 256.0*(k1**3)*k5*(n_o**3)*n_c*n_c
        a1 = 32.0*((k1*n_o)**2)*n_c*( k6 - 4.0*k5*( d2 + k1*n_c ) )
        a2 = 16.0*k1*k5*n_o*( 8.0*k1*k1*n_o*n_c + d2*d2 + 4.0*k1*d2*n_c ) + 8.0*k1*k6*n_o*( 2.0*k6*(n_c-n_n) - 2.0*k1*n_c - d2 )
        a3 = -8.0*k1*k5*( 4.0*k1*d2*n_o + 8.0*k1*k1*n_o*n_c + d2*d2 ) + 4.0*k1*k6*( 2.0*k1*n_o + d2 ) + 4.0*k6*k6*( 2.0*k1*n_n - d2 )
        a4 = 16.0*k1*k1*k5*( k1*n_o + d2 ) + 4.0*k1*k6*(k6-k1)
        a5 = -8.0*(k1**3)*k5
        result = polynomial.polynomial.polyroots([a0,a1,a2,a3,a4,a5])
        result = result[result.real > 0.0]
        result = result[result.real < 2.0*n_o]

        try: # multiple positive realroots 
            result = result[0]#.real
        except: # for only single positive real root
            pass
        n_eq['CO'] = result.real
        
        # for no root found, meaning the mixing ratio of CO is extreamly low
        if not n_eq['CO']: 
            n_eq['CO']= 0.
            print ('No equilibrium solution found for CO! Setting to zero.')
              
        k2 = kk2(tp)
        k3 = kk3(tp,pp)
        k4 = kk4(tp,pp)
        
        # H2O
        c2 = 1.0/kk2(tp)
        n_eq['H2O'] = (2.0*n_o - n_eq['CO'])/(1.0 + 2.0*c2*n_eq['CO'])
        
        # CH4
        n_eq['CH4'] = n_eq['CO']/k1/n_eq['H2O']
        
        # CO2       
        n_eq['CO2'] = n_eq['CO']*n_eq['H2O']/k2

        # C2H2
        n_eq['C2H2'] = k3*n_eq['CH4']**2
        
        # C2H4
        n_eq['C2H4'] = n_eq['C2H2']/k4
        
        # non self-consistent "post-processing" for other species
        
        # kk7: C2H2 + 0.5 * H2 -> C2H3
        kk7 = lambda T, pbar: np.exp( -( gibbs_sp('C2H3',T) -gibbs_sp('C2H2',T) -0.5*gibbs_sp('H2',T) ) ) * (pbar)**0.5
        # kk8: C2H2 + 1.5 * H2 -> C2H5
        kk8 = lambda T, pbar: np.exp( -( gibbs_sp('C2H5',T) -gibbs_sp('C2H2',T) -1.5*gibbs_sp('H2',T) ) ) * (pbar)**1.5
        # kk9: C2H2 + 2 * H2 -> C2H6
        kk9 = lambda T, pbar: np.exp( -( gibbs_sp('C2H6',T) -gibbs_sp('C2H2',T) -2*gibbs_sp('H2',T) ) ) * (pbar)**2
        # kk10: CO + 1.5 * H2 -> CH2OH
        kk10 = lambda T, pbar: np.exp( -( gibbs_sp('CH2OH',T) -gibbs_sp('CO',T) -1.5*gibbs_sp('H2',T) ) ) * (pbar)**1.5
        # kk11: CO + 2 * H2 -> CH3OH
        kk11 = lambda T, pbar: np.exp( -( gibbs_sp('CH3OH',T) -gibbs_sp('CO',T) -2*gibbs_sp('H2',T) ) ) * (pbar)**2
        # kk12: CO + 0.5 * H2 -> HCO
        kk12 = lambda T, pbar: np.exp( -( gibbs_sp('HCO',T) -gibbs_sp('CO',T) -0.5*gibbs_sp('H2',T) ) ) * (pbar)**0.5
        # kk13: CO + H2 -> H2CO
        kk13 = lambda T, pbar: np.exp( -( gibbs_sp('H2CO',T) -gibbs_sp('CO',T) -gibbs_sp('H2',T) ) ) * (pbar)**1
        
        k7, k8, k9, k10, k11, k12, k13 = kk7(tp,pp), kk8(tp,pp), kk9(tp,pp), kk10(tp,pp), kk11(tp,pp), kk12(tp,pp), kk13(tp,pp)
        
        n_eq['C2H3'] = n_eq['C2H2']*k7
        n_eq['C2H5'] = n_eq['C2H2']*k8
        n_eq['C2H6'] = n_eq['C2H2']*k9
        n_eq['CH2OH'] = n_eq['CO']*k10
        n_eq['CH3OH'] = n_eq['CO']*k11
        n_eq['HCO'] = n_eq['CO']*k12
        n_eq['H2CO'] = n_eq['CO']*k13
              
        norm = sum([n_eq[sp] for sp in sp_ini]) + 1. + vulcan_cfg.He_H # 1. is from H2
        for sp in sp_ini: # normalizing so that the total mixing ratio equals 1
            mix_eq[sp] = n_eq[sp]/norm 
        mix_eq['H2'] = 1./norm
        mix_eq['He'] = vulcan_cfg.He_H/norm
    
        return mix_eq
    
    def ini_mol(self):
        return np.array(sop.fsolve(self.abun_lowT, self.ini_m))
       
    def ini_y(self, data_var, data_atm): 
        # initial mixing ratios of the molecules
        
        ini_mol = self.ini_mol() 
           
        ini = np.zeros(ni)
        y_ini = data_var.y
        gas_tot = data_atm.M
    
        for i in range(nz):
            
            if vulcan_cfg.ini_mix == 'EQ':
                nine_sp = ['CH4', 'CO', 'CO2', 'H2O', 'C2H2', 'C2H4', 'NH3', 'N2', 'HCN', 'C2H4', 'C2H5', 'C2H6'\
            , 'CH2OH', 'CH3OH', 'HCO', 'H2CO']
                # sp_ini are the species actually selected by the user
                sp_ini = [sp for sp in species if sp in nine_sp]               
                ini_eq = self.abund_eq9(sp_ini, data_atm.Tco[i], data_atm.pco[i]/1.E6)
                y_ini[i,:] = ini               
                for sp in sp_ini + ['H2', 'He']:
                    y_ini[i,species.index(sp)] = ini_eq[sp]*gas_tot[i]
            
            elif vulcan_cfg.ini_mix == 'CH4': 
                y_ini[i,:] = ini
                y_ini[i,species.index('H2')] = ini_mol[0]*gas_tot[i]; y_ini[i,species.index('H2O')] = ini_mol[1]*gas_tot[i]; y_ini[i,species.index('CH4')] = ini_mol[2]*gas_tot[i]
                # assign rest of the particles to He
                y_ini[i,species.index('He')] = gas_tot[i] - np.sum(y_ini[i,:])
            
            elif vulcan_cfg.ini_mix == 'CO':
                y_ini[i,:] = ini
                y_ini[i,species.index('H2')] = ini_mol[0]*gas_tot[i]; y_ini[i,species.index('H2O')] = ini_mol[1]*gas_tot[i]; y_ini[i,species.index('CO')] = ini_mol[2]*gas_tot[i]
                # assign rest of the particles to He
                y_ini[i,species.index('He')] = gas_tot[i] - np.sum(y_ini[i,:])
        
        ysum = np.sum(y_ini, axis=1).reshape((-1,1))
        # storing ymix
        data_var.ymix = y_ini/ysum
        
        return data_var
        


    def ele_sum(self, data_var):
        
        for atom in self.atom_list:
            data_var.atom_ini[atom] = np.sum([compo[compo_row.index(species[i])][atom] * data_var.y[:,i] for i in range(ni)])
            data_var.atom_loss[atom] = 0.

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
    
        

    def load_TPK(self, data_atm, output):
        
        PTK_fun = {}
        
        if self.type == 'isothermal': 
            data_atm.Tco = np.repeat(vulcan_cfg.Tiso,nz)
            data_atm.Kzz = np.repeat(vulcan_cfg.const_Kzz,nz-1)
            
        elif self.type == 'analytical': 
            
            # plotting T-P on the fly                               
            para_atm = vulcan_cfg.para_anaTP 
            
            # return the P-T function
            PTK_fun['pT'] = lambda pressure: self.TP_H14(pressure, *para_atm)        
            data_atm.Tco = PTK_fun['pT'](data_atm.pco)
            data_atm.Kzz = np.repeat(vulcan_cfg.const_Kzz,nz-1)
            
        elif self.type == 'file':
            
            if self.Kzz_prof == 'const':     
                atm_table = np.genfromtxt(vulcan_cfg.atm_file, names=True, dtype=None, skip_header=1)
                p_file, T_file = atm_table['Pressure'], atm_table['Temp']
            
            elif self.Kzz_prof == 'file':
                atm_table = np.genfromtxt(vulcan_cfg.atm_file, names=True, dtype=None, skip_header=1)
                p_file, T_file, Kzz_file = atm_table['Pressure'], atm_table['Temp'], atm_table['Kzz']
 
            else: raise IOError ('\n"Kzz_prof" (the type of Kzz profile) cannot be recongized.\nPlease assign it as "file" or "const" in vulcan_cfg.')

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
               
        # calculating and storing M(the third body)
        data_atm.M = data_atm.pco/(kb*data_atm.Tco)
        data_atm.n_0 = data_atm.M.copy()
        
        # plot T-P profile
        if vulcan_cfg.plot_TP == True: output.plot_TP(data_atm)
        
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
        """ calculating the molar mass of each species?"""
        return compo['mass'][compo_row.index(sp)]

    def mean_mass(self, var, atm, ni):
        mu = np.zeros(nz)
        for i in range(ni):
            mu += self.mol_mass(species[i]) * var.ymix[:,i]
        atm.mu = mu
        return atm        
        
    def f_mu_dz(self, data_var, data_atm):  
            
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

        dzi = 0.5*(dz + np.roll(dz,1))
        dzi = dzi[1:]
        # for the j grid, dzi[j] from the grid above and dz[j-1] from the grid below
        
        # updating and storing dz and dzi
        data_atm.dz = dz
        data_atm.dzi = dzi
        
        return data_atm

if __name__ == "__main__":
    print("This module stores classes for constructing atmospheric structure \
    and initializing its chemical composition from the desinated elemental abudance.")