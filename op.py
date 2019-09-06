# ==============================================================================
# Module includes all the numerical functions VULCAN. 
# Copyright (C) 2016 Shang-Min Tsai (Shami)                    
# ==============================================================================
# ReadRate() reads in the chemical network and construct the rate constants based
# on the T-P sturcture.
# Integration() is the backbone of integrating for one time step
# ODESolver() contains the functions for solving system of ODEs (e.g. dy/dt, Jacobian, etc.)
# ==============================================================================

import numpy as np
import scipy
from scipy import sparse
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.legend as lg
import time, os, pickle
import csv, ast

#from builtins import input
#from collections import defaultdict
# TODO :test the TODO buldle

import vulcan_cfg
try: from PIL import Image
except ImportError: 
    try: import Image
    except: vulcan_cfg.use_PIL = False

import build_atm
import chem_funs
from chem_funs import ni, nr  # number of species and reactions in the network

from phy_const import kb, Navo, hc, ag0 # hc is used to convert to the actinic flux

from vulcan_cfg import nz

# imported functions 
chemdf = chem_funs.chemdf
neg_achemjac = chem_funs.neg_symjac
compo = build_atm.compo
compo_row = build_atm.compo_row

species = chem_funs.spec_list


class ReadRate(object):
    
    """
    to read in rate constants from the network file and compute the reaction rates for the corresponding Tco and pco 
    """
    
    def __init__(self):
        
        self.i = 1
        # flag of trimolecular reaction
        self.re_tri, self.re_tri_k0 = False, False
        self.list_tri = []

        
    def read_rate(self, var, atm):
        
        Rf, Rindx, a, n, E, a_inf, n_inf, E_inf, k, k_fun, k_inf, kinf_fun, k_fun_new, pho_rate_index = \
        var.Rf, var.Rindx, var.a, var.n, var.E, var.a_inf, var.n_inf, var.E_inf, var.k, var.k_fun, var.k_inf, var.kinf_fun, var.k_fun_new,\
         var.pho_rate_index
        
        i = self.i
        re_tri, re_tri_k0 = self.re_tri, self.re_tri_k0
        list_tri = self.list_tri
        
        Tco = atm.Tco.copy()
        M = atm.M.copy()
        
        special_re = False
        conden_re = False
        photo_re = False
        end_re = False
        br_read  = False
        
        photo_sp = [] 
               
        with open(vulcan_cfg.network) as f:
            all_lines = f.readlines()
            for line_indx, line in enumerate(all_lines):
                
                # switch to 3-body and dissociation reations 
                if line.startswith("# 3-body"): 
                    re_tri = True
                    
                if line.startswith("# 3-body reactions without high-pressure rates"):
                    re_tri_k0 = True
                    
                elif line.startswith("# special"): 
                    special_re = True # switch to reactions with special forms (hard coded)  
                
                elif line.startswith("# condensation"): 
                    special_re = False 
                    conden_re = True
                    var.conden_indx = i
                    
                elif line.startswith("# photo"):
                    conden_re = False
                    special_re = False # turn off reading in the special form
                    photo_re = True
                    var.photo_indx = i  
                
                elif line.startswith("# re_end"):
                    end_re = True
                    
                elif line.startswith("# braching info start"):
                    br_read = True
       
                elif line.startswith("# braching info end"):
                    br_read = False
                    
                # skip common lines and blank lines
                # ========================================================================================
                if not line.startswith("#") and line.strip() and special_re == False and conden_re == False and photo_re == False and end_re == False: # if not starts
                    
                    Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                    li = line.partition(']')[-1].strip()
                    columns = li.split()
                    Rindx[i] = int(line.partition('[')[0].strip())
                              
                    a[i] = float(columns[0])
                    n[i] = float(columns[1])
                    E[i] = float(columns[2])
                
                    # switching to trimolecular reactions (len(columns) > 3 for those with high-P limit rates)   
                    if re_tri == True and re_tri_k0 == False:
                        a_inf[i] = float(columns[3])
                        n_inf[i] = float(columns[4])
                        E_inf[i] = float(columns[5])
                        list_tri.append(i) 
                    
                    if columns[-1].strip() == 'He': re_He = i
                    elif columns[-1].strip() == 'ex1': re_CH3OH = i
                
                    # Note: make the defaut i=i
                    k_fun[i] = lambda temp, mm, i=i: a[i] *temp**n[i] * np.exp(-E[i]/temp)
                
                
                    if re_tri == False:
                        k[i] = k_fun[i](Tco, M)
                    
                    # for 3-body reactions, also calculating k_inf
                    elif re_tri == True and len(columns)>=6:
        
        
                        kinf_fun[i] = lambda temp, i=i: a_inf[i] *temp**n_inf[i] * np.exp(-E_inf[i]/temp)
                        k_fun_new[i] = lambda temp, mm, i=i: (a[i] *temp**n[i] * np.exp(-E[i]/temp))/(1 + (a[i] *temp**n[i] * np.exp(-E[i]/temp))*mm/(a_inf[i] *temp**n_inf[i] * np.exp(-E_inf[i]/temp)) ) 
        
                        #k[i] = k_fun_new[i](Tco, M)
                        k_inf = a_inf[i] *Tco**n_inf[i] * np.exp(-E_inf[i]/Tco)
                        k[i] = k_fun[i](Tco, M)
                        k[i] = k[i]/(1 + k[i]*M/k_inf )
        
        
                    else: # for 3-body reactions without high-pressure rates
                        k[i] = k_fun[i](Tco, M)
                           
                    ### TEST CAPPING
                    # k[i] = np.minimum(k[i],1.E-11)
                    ###
    
                    i += 2
                    # end if not 
                 # ========================================================================================    
                elif special_re == True and line.strip() and not line.startswith("#") and end_re == False:

                    Rindx[i] = int(line.partition('[')[0].strip())
                    Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                
                    if Rf[i] == 'OH + CH3 + M -> CH3OH + M':
                        print ('Using special form for the reaction: ' + Rf[i])
                    
                        k[i] = 1.932E3*Tco**-9.88 *np.exp(-7544./Tco) + 5.109E-11*Tco**-6.25 *np.exp(-1433./Tco)
                        k_inf = 1.031E-10 * Tco**-0.018 *np.exp(16.74/Tco)
                        k[i] = k[i]/(1 + k[i]*M/k_inf )
                    
                        k_fun[i] = lambda temp, mm, i=i: 1.932E3 *temp**-9.88 *np.exp(-7544./temp) + 5.109E-11*temp**-6.25 *np.exp(-1433./temp)
                        kinf_fun[i] = lambda temp, mm, i=i: 1.031E-10 * temp**-0.018 *np.exp(16.74/temp)
                        k_fun_new[i] = lambda temp, mm, i=i: (1.932E3 *temp**-9.88 *np.exp(-7544./temp) + 5.109E-11*temp**-6.25 *np.exp(-1433./temp))/\
                        (1 + (1.932E3 *temp**-9.88 *np.exp(-7544./temp) + 5.109E-11*temp**-6.25 *np.exp(-1433./temp)) * mm / (1.031E-10 * temp**-0.018 *np.exp(16.74/temp)) )
                    
                    i += 2
                
                # Testing condensation
                elif conden_re == True and line.strip() and not line.startswith("#") and end_re == False:
                    Rindx[i] = int(line.partition('[')[0].strip())
                    Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                    
                    var.conden_re_list.append(i)
                    k[i] = np.zeros(nz)
                    k[i+1] = np.zeros(nz)
                           
                    i += 2
                    
                # setting photo dissociation reactions to zeros
                elif photo_re == True and line.strip() and not line.startswith("#") and end_re == False:
                    
                    k[i] = np.zeros(nz)
                    Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                    
                    # adding the photo species
                    photo_sp.append(Rf[i].split()[0])
                    
                    li = line.partition(']')[-1].strip()
                    columns = li.split()
                    Rindx[i] = int(line.partition('[')[0].strip())
                    # columns[0]: the species being dissocited; branch index: columns[1]
                    pho_rate_index[(columns[0],int(columns[1]))] = Rindx[i]
                    
                    # store the number of branches
                    var.n_branch[columns[0]] = int(columns[1])
                    
                    i += 2
                
                # end_re == True
                elif br_read == True and not line.startswith("#"):
                    # read in the quantum yields of photolysis reactions
                    sp_list = line.partition(':')
                    sp = sp_list[0]
                    lists = sp_list[-1]
                    wavelen_yield = lists.partition(';')
                    # wavelen_yield is tuple of string in wavelength seitch, ;, Q yield e.g. ('[165.]', ';', '[(1.,0),(0,1.)]')
                    var.wavelen[sp] = ast.literal_eval(wavelen_yield[0].strip())
                    var.br_ratio[sp] = ast.literal_eval(wavelen_yield[-1].strip())
                
        k_fun.update(k_fun_new)
    
        # store k into data_var
        # remeber k_fun has not removed reactions from remove_list
        var.k = k
        var.k_fun = k_fun
        var.kinf_fun = kinf_fun
        
        var.photo_sp = set(photo_sp)

        return var
        
    def rev_rate(self, var, atm):
        # ToDo
        # Re-organising these
        if vulcan_cfg.use_photo == False and vulcan_cfg.use_condense == False:
            rev_list = range(2,nr+1,2)
        elif vulcan_cfg.use_photo == False and vulcan_cfg.use_condense == True:
            rev_list = range(2,var.conden_indx,2)
            # setting the reversal zeros
            for i in range(var.conden_indx+1, nr+1,2):
                var.k[i] = np.zeros(nz)
        elif vulcan_cfg.use_photo == True and vulcan_cfg.use_condense == False:
            try: rev_list = range(2,var.photo_indx,2)
            except: print ('No photodissociation reaction provided while "use_photo" set True in vulcan_cfg.\n'); raise
            # setting the reversal zeros
            for i in range(var.photo_indx+1, nr+1,2):
                var.k[i] = np.zeros(nz)
        else: #vulcan_cfg.use_photo == True and vulcan_cfg.use_condense == True
            try: rev_list = range(2,var.conden_indx,2)
            except: print ('No condensation reaction provided while "use_photo" set True in vulcan_cfg.\n'); raise            
            # setting the reversal of both condensation and photolisys rates to zero
            for i in range(var.photo_indx+1, nr+1,2):
                var.k[i] = np.zeros(nz)
                
        Tco = atm.Tco.copy()
        
        # reversing rates and storing into data_var
        print ('Rates greater than 1e-6:')
        for i in rev_list: 
            var.k_fun[i] = lambda temp, mm, i=i: var.k_fun[i-1](temp, mm)/chem_funs.Gibbs(i-1,temp)
            var.k[i] = var.k[i-1]/chem_funs.Gibbs(i-1,Tco)
            
            if np.any(var.k[i] > 1.e-6): print ('R' + str(i) + ':  ' + str(np.amax(var.k[i])) )
            if np.any(var.k[i-1] > 1.e-6): print ('R' + str(i-1) + ':  ' + str(np.amax(var.k[i-1])) )
        return var
        
    
    def remove_rate(self, var):
        
        for i in vulcan_cfg.remove_list:
            var.k[i] = np.repeat(0.,nz)
            var.k_fun[i] = lambda temp, mm, i=i: np.repeat(0.,nz)
            
        return var
        
    def read_rateFun(self, var):
        '''
        Reading in the reaction network and returning only the functions (k_fun)
        Used for pathway analysis
        '''
        
        Rf, Rindx, a, n, E, a_inf, n_inf, E_inf, k_fun, kinf_fun, k_fun_new = \
        var.Rf, var.Rindx, var.a, var.n, var.E, var.a_inf, var.n_inf, var.E_inf, var.k_fun, var.kinf_fun,  var.k_fun_new
        
        i = self.i
        re_tri, re_tri_k0 = self.re_tri, self.re_tri_k0
        list_tri = self.list_tri
        
        special_re = False
        conden_re = False
        photo_re = False
        end_re = False
        br_read  = False
        
        photo_sp = []
               
        with open(vulcan_cfg.network) as f:
            for line in f.readlines():
                
                # switch to 3-body and dissociation reations 
                if line.startswith("# 3-body"): 
                    re_tri = True
                    
                if line.startswith("# 3-body reactions without high-pressure rates"):
                    re_tri_k0 = True
                    
                elif line.startswith("# special"): 
                    special_re = True # switch to reactions with special forms (hard coded)                   
                
                elif line.startswith("# photo"):
                    special_re = False # turn off reading in the special form
                    photo_re = True
                    var.photo_indx = i  
                
                elif line.startswith("# re_end"):
                    end_re = True
                    
                elif line.startswith("# braching info start"):
                    br_read = True
       
                elif line.startswith("# braching info end"):
                    br_read = False
                       
                # skip common lines and blank lines
                # ========================================================================================
                if not line.startswith("#") and line.strip() and special_re == False and photo_re == False and end_re == False: # if not starts
                    
                    Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                    li = line.partition(']')[-1].strip()
                    columns = li.split()
                    Rindx[i] = int(line.partition('[')[0].strip())
                              
                    a[i] = float(columns[0])
                    n[i] = float(columns[1])
                    E[i] = float(columns[2])
                
                    # switching to trimolecular reactions (len(columns) > 3 for those with high-P limit rates)   
                    if re_tri == True and re_tri_k0 == False:
                        a_inf[i] = float(columns[3])
                        n_inf[i] = float(columns[4])
                        E_inf[i] = float(columns[5])
                        list_tri.append(i) 
                    
                    if columns[-1].strip() == 'He': re_He = i
                    elif columns[-1].strip() == 'ex1': re_CH3OH = i
                
                    # Note: make the defaut i=i
                    k_fun[i] = lambda temp, mm, i=i: a[i] *temp**n[i] * np.exp(-E[i]/temp)
                    
                    
                    # for 3-body reactions, also calculating k_inf
                    if re_tri == True and len(columns)>=6:
          
                        kinf_fun[i] = lambda temp, i=i: a_inf[i] *temp**n_inf[i] * np.exp(-E_inf[i]/temp)
                        k_fun_new[i] = lambda temp, mm, i=i: (a[i] *temp**n[i] * np.exp(-E[i]/temp))/(1 + (a[i] *temp**n[i] * np.exp(-E[i]/temp))*mm/(a_inf[i] *temp**n_inf[i] * np.exp(-E_inf[i]/temp)) ) 
                        
                    i += 2
                    # end if not 
                 # ========================================================================================    
                elif special_re == True and line.strip() and not line.startswith("#") and end_re == False:

                    Rindx[i] = int(line.partition('[')[0].strip())
                    Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                
                    if Rf[i] == 'OH + CH3 + M -> CH3OH + M':
                        print ('Using special form for the reaction: ' + Rf[i])
                        k_fun[i] = lambda temp, mm, i=i: 1.932E3 *temp**-9.88 *np.exp(-7544./temp) + 5.109E-11*temp**-6.25 *np.exp(-1433./temp)
                        kinf_fun[i] = lambda temp, mm, i=i: 1.031E-10 * temp**-0.018 *np.exp(16.74/temp)
                        k_fun_new[i] = lambda temp, mm, i=i: (1.932E3 *temp**-9.88 *np.exp(-7544./temp) + 5.109E-11*temp**-6.25 *np.exp(-1433./temp))/\
                        (1 + (1.932E3 *temp**-9.88 *np.exp(-7544./temp) + 5.109E-11*temp**-6.25 *np.exp(-1433./temp)) * mm / (1.031E-10 * temp**-0.018 *np.exp(16.74/temp)) )
                    
                    i += 2
                
                # setting photo dissociation reactions to zeros
                elif photo_re == True and line.strip() and not line.startswith("#") and end_re == False:
                    
                    #k[i] = np.zeros(nz)
                    Rf[i] = line.partition('[')[-1].rpartition(']')[0].strip()
                    
                    # adding the photo species
                    photo_sp.append(Rf[i].split()[0])
                    
                    li = line.partition(']')[-1].strip()
                    columns = li.split()
                    Rindx[i] = int(line.partition('[')[0].strip())
                    # columns[0]: the species being dissocited; branch index: columns[1]
                    pho_rate_index[(columns[0],int(columns[1]))] = Rindx[i]
                    
                    # store the number of branches
                    var.n_branch[columns[0]] = int(columns[1])
                    
                    i += 2
                
                # end_re == True
                elif br_read == True and not line.startswith("#"):
                    # read in the quantum yields of photolysis reactions
                    sp_list = line.partition(':')
                    sp = sp_list[0]
                    lists = sp_list[-1]
                    wavelen_yield = lists.partition(';')
                    # wavelen_yield is tuple of string in wavelength seitch, ;, Q yield e.g. ('[165.]', ';', '[(1.,0),(0,1.)]')
                    var.wavelen[sp] = ast.literal_eval(wavelen_yield[0].strip())
                    var.br_ratio[sp] = ast.literal_eval(wavelen_yield[-1].strip())
                    
        k_fun.update(k_fun_new)
    
        # store k_fun into data_var
        var.k_fun = k_fun
        var.kinf_fun = kinf_fun

        return var
        
    def rev_rateFun(self, var):
        '''
        Revarsing only the functions of forward rates (k_fun) and the T, P values (at the examined level)
        Used for pathway analysis
        '''
        
        rev_list = range(2,nr+1,2)
                
        # reversing rates and storing into data_var
        for i in rev_list: 
            var.k_fun[i] = lambda temp, mm, i=i: var.k_fun[i-1](temp, mm)/chem_funs.Gibbs(i-1,temp)
            
        return var
    
    def make_bins_read_cross(self,var):
        '''
        determining the bin range and only use the min and max wavelength that the molecules absorb
        to avoid photons with w0=1 (pure scatteing) in certain wavelengths
        '''
        photo_sp = list(var.photo_sp)
        sp0 = photo_sp[0]
        
        cross_raw, scat_raw = {}, {}
        # reading in cross sections into dictionary
        for n, sp in enumerate(photo_sp):
            try: cross_raw[sp] = np.genfromtxt(vulcan_cfg.cross_folder+sp+'_cross.csv',dtype=float,delimiter=',',skip_header=1, names = ['lambda','cross','disso'])
            except: print ('\nMissing the cross section from ' + sp)
            if cross_raw[sp]['cross'][0] == 0 or cross_raw[sp]['cross'][-1] ==0:
                raise IOError ('\n Please remove the zeros in the cross file of ' + sp)
            
            if n==0: # the first species
                bin_min = cross_raw[sp]['lambda'][0]
                bin_max = cross_raw[sp]['lambda'][-1]
            else:
                sp_min, sp_max = cross_raw[sp]['lambda'][0], cross_raw[sp]['lambda'][-1]
                if sp_min < bin_min: bin_min = sp_min
                if sp_max > bin_max: bin_max = sp_max
        
        # constraining the bin_min and bin_max by the default values defined in store.py
        bin_min = max(bin_min, var.def_bin_min)
        bin_max = min(bin_max, var.def_bin_max)
        
        var.dbin = vulcan_cfg.dbin
        bins = np.arange(bin_min,bin_max, var.dbin)
        bins[-1] -= 1.E-10  # for solaving the floating pb in binary when doing read_cross 
        var.bins = bins
        var.nbin = len(bins)
        
        # all variables that depend on nbins
        # the direct beam (staggered)
        var.sflux = np.zeros( (nz+1, var.nbin) )
        # the diffusive flux (staggered)
        var.dflux_u, var.dflux_d = np.zeros( (nz+1, var.nbin) ), np.zeros( (nz+1, var.nbin) )
        # the total actinic flux (non-staggered)
        var.aflux = np.zeros( (nz, var.nbin) )
        # the total actinic flux from the previous calculation 
        prev_aflux = np.zeros( (nz, var.nbin) )
        
        # staggered
        var.tau = np.zeros( (nz+1, var.nbin) )
        # the stellar flux at TOA
        var.sflux_top = np.zeros(var.nbin)
        
        
        # read_cross
        # creat a dict of cross section with key=sp and values=bins in data_var
        var.cross = dict([(sp, np.zeros(var.nbin)) for sp in var.photo_sp])
        
        #read cross of disscoiation
        var.cross_J = dict([(sp, np.zeros(var.nbin)) for sp in var.photo_sp])
        
        var.cross_scat = dict([(sp, np.zeros(var.nbin)) for sp in vulcan_cfg.scat_sp])
        
        for sp in var.photo_sp:
            # for values outside the boundary => fill_value = 0
            inter_cross = interpolate.interp1d(cross_raw[sp]['lambda'], cross_raw[sp]['cross'], bounds_error=False, fill_value=0)
            inter_cross_J = interpolate.interp1d(cross_raw[sp]['lambda'], cross_raw[sp]['disso'], bounds_error=False, fill_value=0)
            
            for n, ld in enumerate(bins):
                var.cross[sp][n] = inter_cross(ld)
                var.cross_J[sp][n] = inter_cross_J(ld)
                
        # reading in cross sections of Rayleigh Scattering
        for sp in vulcan_cfg.scat_sp:
            scat_raw[sp] = np.genfromtxt(vulcan_cfg.cross_folder + 'rayleigh/' + sp+'_scat.txt',dtype=float,\
            skip_header=1, names = ['lambda','cross'])

            # for values outside the boundary => fill_value = 0
            inter_scat = interpolate.interp1d(scat_raw[sp]['lambda'], scat_raw[sp]['cross'], bounds_error=False, fill_value=0)
            
            for n, ld in enumerate(bins):
                var.cross_scat[sp][n] = inter_scat(ld)
        
    # def read_cross(self, var):
    #     # flux stored in a 2D array: lambda, flux
    #     # using flux[np.where(flux[:,0]==lambda),1]
    #
    #     # creat a dict of cross section with key=sp and values=bins in data_var
    #     var.cross = dict([(sp, np.zeros(var.nbin)) for sp in var.photo_sp])
    #     var.cross_scat = dict([(sp, np.zeros(var.nbin)) for sp in vulcan_cfg.scat_sp])
    #
    #     for sp in var.photo_sp:
    #         # for values outside the boundary => fill_value = 0
    #         inter_cross = interpolate.interp1d(cross_raw[sp]['lambda'], cross_raw[sp]['cross'], bounds_error=False, fill_value=0)
    #         for n, ld in enumerate(bins):
    #             var.cross[sp][n] = inter_cross(ld)
        
        
        # cross_raw, scat_raw = {}, {}
        # # creat a dict of cross section with key=sp and values=bins in data_var
        # var.cross = dict([(sp, np.zeros(var.nbin)) for sp in var.photo_sp])
        # var.cross_scat = dict([(sp, np.zeros(var.nbin)) for sp in vulcan_cfg.scat_sp])
        # bins = var.bins
        #
        # # reading in cross sections into dictionary
        # for sp in var.photo_sp:
        #     cross_raw[sp] = np.genfromtxt(vulcan_cfg.cross_folder+sp+'_cross.csv',dtype=float,delimiter=',',skip_header=1, names = ['lambda','cross'])
        #
        #     # for values outside the boundary => fill_value = 0
        #     inter_cross = interpolate.interp1d(cross_raw[sp]['lambda'], cross_raw[sp]['cross'], bounds_error=False, fill_value=0)
        #     for n, ld in enumerate(bins):
        #         var.cross[sp][n] = inter_cross(ld)
            
        # reading in cross sections of Rayleigh Scattering
        for sp in vulcan_cfg.scat_sp:
            scat_raw[sp] = np.genfromtxt(vulcan_cfg.cross_folder + 'rayleigh/' + sp+'_scat.txt',dtype=float,\
            skip_header=1, names = ['lambda','cross'])

            # for values outside the boundary => fill_value = 0
            inter_scat = interpolate.interp1d(scat_raw[sp]['lambda'], scat_raw[sp]['cross'], bounds_error=False, fill_value=0)
            
            for n, ld in enumerate(bins):
                var.cross_scat[sp][n] = inter_scat(ld)

class Integration(object):
    """
    time-stepping until the stopping criteria (steady-state) is satisfied
    #all the operators required in the continuity equation: dn/dt + dphi/dz = P - L
    #or class incorporating the esential numerical operations?
    """
    
    def __init__(self, odesolver, output):

        self.mtol = vulcan_cfg.mtol
        self.atol = vulcan_cfg.atol
        self.output = output
         
        self.odesolver = odesolver
        self.non_gas_sp = vulcan_cfg.non_gas_sp
        self.use_settling = vulcan_cfg.use_settling 
        self.g = vulcan_cfg.g
        
        # TEST
        if vulcan_cfg.use_photo == True: self.update_photo_frq = vulcan_cfg.ini_update_photo_frq
        
        
    def __call__(self, var, atm, para, make_atm):
        
        use_print_prog, use_live_plot = vulcan_cfg.use_print_prog, vulcan_cfg.use_live_plot
        
        while not self.stop(var, para, atm): # Looping until the stop condition is satisfied
            
            var = self.backup(var)
            
            # updating tau, flux, and the photolosys rate
            # TEST
            if vulcan_cfg.use_photo == True and var.longdy < vulcan_cfg.yconv_min*10. and var.longdydt < 1.e-6:  
                self.update_photo_frq = vulcan_cfg.final_update_photo_frq
                if para.switch_final_photo_frq == False:
                    print ('update_photo_frq changed to ' + str(vulcan_cfg.final_update_photo_frq) +'\n')
                    para.switch_final_photo_frq = True
            
            if vulcan_cfg.use_photo == True and para.count % self.update_photo_frq == 0:
                self.odesolver.compute_tau(var, atm)
                self.odesolver.compute_flux(var, atm)
                self.odesolver.compute_J(var, atm)
            
            # Testing condensation
            # updating the condensation rates
            if vulcan_cfg.use_condense == True and var.t > vulcan_cfg.start_conden_time:
                var = self.conden(var,atm)
            
                
            # integrating one step
            var, para = self.odesolver.one_step(var, atm, para)

            if para.count % vulcan_cfg.update_frq == 0: # updating mu and dz (dzi) due to diffusion
                atm = self.update_mu_dz(var, atm, make_atm)
            
            # Test?
            # MAINTAINING HYDROSTATIC BALANCE
            if vulcan_cfg.use_condense == True:
                var.y[:,atm.exc_conden] = np.vstack(atm.n_0)*var.ymix[:,atm.exc_conden]
            else:
                var.y = np.vstack(atm.n_0)*var.ymix
            
            # calculating the changing of y
            var = self.f_dy(var, para)
            
            # save values of the current step
            var, para = self.save_step(var, para)
            
            # adjusting the step-size
            var = self.odesolver.step_size(var, para)
            
            if use_print_prog == True and para.count % vulcan_cfg.print_prog_num==0:
                self.output.print_prog(var)
                
            if vulcan_cfg.use_live_flux == True and vulcan_cfg.use_photo == True and para.count % vulcan_cfg.live_plot_frq ==0:
                #plt.figure('flux')
                self.output.plot_flux_update(var, atm, para)
                
            if use_live_plot == True and para.count % vulcan_cfg.live_plot_frq ==0:
                #plt.figure('mix')
                self.output.plot_update(var, atm, para)
            
            
        
    def backup(self, var):
        var.y_prev = np.copy(var.y)
        var.dy_prev = np.copy(var.dy)
        var.atom_loss_prev = var.atom_loss.copy()
        return var
        
    def update_mu_dz(self, var, data_atm, make_atm): #y, ni, spec, Tco, pco
        
        # gravity
        g = make_atm.g
        
        # calculating mu (mean molecular weight)
        data_atm = make_atm.mean_mass(var, data_atm, ni)
        data_atm.Hp = kb*data_atm.Tco/(data_atm.mu/Navo*g)
        
        # calculating the pressure at interface
        data_atm = make_atm.f_pico(data_atm)
        
        # updating dz, zco, and dzi
        dz, dzi, pico = data_atm.dz, data_atm.dzi, data_atm.pico
        dz = data_atm.Hp* np.log(pico[:-1]/np.roll(pico,-1)[:-1])
        data_atm.zco = np.insert(np.cumsum(dz),0, 0.) # cumulative sum of dz
        
        # for the j grid, dzi[j] counts the distance from the grid above and dz[j-1] from the grid below
        dzi = 0.5*(dz + np.roll(dz,1))
        dzi = dzi[1:]
        
        data_atm.dz, data_atm.dzi, data_atm.pico = dz, dzi, pico
        
        # Calculating the settling velocity
        if self.use_settling == True:
            # Using Gao 2018 (50)
            r_p = 1.e-4
            R_uni = kb*Navo # the universal gas const
        
            Ti = data_atm.Ti
            mui = 0.5*(data_atm.mu + np.roll(data_atm.mu,-1))
            mui = mui[:-1]
            pi = data_atm.pico[1:-1]
            
            # TEST
            # mass density of the air
            rho_a = pi/Ti/(R_uni/mui)
            data_atm.rho_a = rho_a
            for sp in self.non_gas_sp:
                if sp == 'H2O_l_s':
                    rho_p = 0.9
                    # negative for downward
                    data_atm.vs[:,species.index(sp)] = -0.5*rho_p*g*r_p / rho_a *(np.pi*mui /(2*R_uni*Ti))**0.5
                    
        return data_atm
    

    # function calculating the change of y
    def f_dy(self, var, para): # y, y_prev, ymix, yconv, count, dt
        if para.count == 0: 
            var.dy, var.dydt = 1., 1.
            return var
        y, ymix, y_prev = var.y, var.ymix, var.y_prev    
        dy =  np.abs(y - y_prev)
        dy[ymix < vulcan_cfg.mtol] = 0   
        dy[y < vulcan_cfg.atol] = 0 
        dy = np.amax( dy[y>0]/y[y>0] )
        
        var.dy, var.dydt = dy, dy/var.dt
        
        return var
    
    
    def conv(self, var, para, atm, out=False, print_freq=100):
        '''
        funtion returns TRUE if the convergence condition is satisfied
        '''
        st_factor, mtol_conv, atol, yconv_cri, slope_cri, yconv_min =\
         vulcan_cfg.st_factor, vulcan_cfg.mtol_conv, vulcan_cfg.atol, vulcan_cfg.yconv_cri, vulcan_cfg.slope_cri, vulcan_cfg.yconv_min
        y, ymix, y_time, t_time = var.y.copy(), var.ymix.copy(), var.y_time, var.t_time
        count = para.count
        
        #slope_min = min( np.amin(atm.Kzz)/np.amax(0.1*atm.Hp)**2 , 1.e-8)
        slope_min = min( np.amin(atm.Kzz/(0.1*atm.Hp[:-1])**2) , 1.e-8)
        slope_min = max(slope_min, 1.e-12)
        #print ('slope_min= ' + "{:.2e}".format(slope_min))
        
        # if t < trun_min: indx = -100
        indx = np.abs(t_time-var.t*st_factor).argmin()   
     
        if indx == para.count-1: indx-=1  #Important!! For dt larger than half of the runtime (count-1 is the last one) 
        
        longdy = np.abs((y_time[count-1] - y_time[indx])/np.vstack(atm.n_0))
        longdy[ymix < mtol_conv] = 0
        longdy[y < atol] = 0
        
        if vulcan_cfg.use_condense == True: 
            for sp in self.non_gas_sp:
                longdy[:,species.index(sp)] = 0
         
        #indx_max = np.nanargmax(longdy[ymix>0]/ymix[ymix>0]) # this index is effectively from species[ymix>0]
        longdy = np.amax( longdy[ymix>0]/ymix[ymix>0] )
        longdydt = longdy/(t_time[-1]-t_time[indx])
        # store longdy and longdydt
        var.longdy, var.longdydt = longdy, longdydt
        
        #if para.count % print_freq==0: print ('from... ' + str(int(indx_max/ni)) + ' , ' + species[indx_max%ni])
        
        if (longdy < yconv_cri and longdydt < slope_cri or longdy < yconv_min and longdydt < slope_min) and var.aflux_change<vulcan_cfg.flux_cri: 
            return True
        return False
    
    def stop(self, var, para, atm):
        '''
        To check the convergence criteria and stop the integration 
        '''
        if var.t > vulcan_cfg.trun_min and para.count > vulcan_cfg.count_min and self.conv(var, para, atm):
            print ('Integration successful with ' + str(para.count) + ' steps and long dy, long dydt = ' + str(var.longdy) + ' ,' + str(var.longdydt) + '\nActinic flux change: ' + '{:.2E}'.format(var.aflux_change)) 
            self.output.print_end_msg(var, para)
            para.end_case = 1
            return True
        elif var.t > vulcan_cfg.runtime:
            print ("After ------- %s seconds -------" % ( time.time()- para.start_time ) + ' s CPU time')
            print ('Integration not completed...\nMaximal allowed runtime exceeded ('+ \
            str (vulcan_cfg.runtime) + ' sec)!')
            para.end_case = 2
            return True
        elif para.count > vulcan_cfg.count_max:
            print ("After ------- %s seconds -------" % ( time.time()- para.start_time ) + ' s CPU time')
            print ('Integration not completed...\nMaximal allowed steps exceeded (' + \
            str (vulcan_cfg.count_max) + ')!')
            para.end_case = 3
            return True
    
    def save_step(self, var, para):
        '''
        save current values of y and add 1 to the counter
        '''
        var.t += var.dt   
        para.count += 1
        
        # tmp = list(var.y)
        #if para.count % self.y_time_freq ==0:
        var.y_time.append(var.y)
        #var.ymix_time.append(var.ymix.copy())
        var.t_time.append(var.t)
    
        # only used in PI_control
        var.dy_time.append(var.y)
        var.dydt_time.append(var.dydt)
        var.atom_loss_time.append(list(var.atom_loss.values()) )
        
        return var, para
        
    
    # TESTing condensation
    def conden(self, var, atm):
        '''
        Updating the condensation reactions according to the new number density
        using the timescale given in (15) Hu 2012
        '''
        for re in var.conden_re_list:
            if var.Rf[re] == 'H2O -> H2O_l_s':
                m = 18./Navo
                rho_p = 0.9
                r_p = 1.e-4 # assuming 1 micron from the stratospheric aerosols in Table 2. [Toon & Farlow]
            
                rate_c = m/(4*rho_p)*(8*kb*atm.Tco/np.pi/m)**0.5 *(var.y[:,species.index('H2O')]-atm.sat_p['H2O']/kb/atm.Tco)/r_p
                
                #var.rate_c = rate_c
                
                var.k[re] = rate_c
                var.k[re+1] = rate_c
                
                # positive: condensation
                var.k[re] = np.maximum(var.k[re], 0)
                # negative: evaporation
                var.k[re+1] = np.minimum(var.k[re+1], 0)
                var.k[re+1] = np.abs(var.k[re+1])
                
                # TEST capping
                rate_mix = np.amax(atm.Kzz/(atm.Hpi)**2) # the max rate (shortest tau_dyn)
                var.k[re] = np.minimum(var.k[re], rate_mix)
                var.k[re+1] = 0
                
                # var.k[re] *= 1e-6
#                 var.k[re+1] *= 1e-20
                
                # /kb/data_atm.Tco is to convert p_sat into n_sat
                #print ('water conden')
                #print (var.k[re])
            
        
            elif var.Rf[re] == 'NH3 -> NH3_l':
                m = 17./Navo
                rho_p = 0.7
                r_p = 1.e-4 # assuming 1 micron
            
                rate_c = m/(4*rho_p)*(8*kb*atm.Tco/np.pi/m)**0.5 *(var.y[:,species.index('NH3')]-atm.sat_p['NH3']/kb/atm.Tco)/r_p
        
                var.k[re] = rate_c
                var.k[re+1] = rate_c
                
                # positive: condensation
                var.k[re] = np.maximum(var.k[re], 0)
                # negative: evaporation
                var.k[re+1] = np.minimum(var.k[re+1], 0)
                var.k[re+1] = np.abs(var.k[re+1])
                
                # TEST capping
                #var.k[re] = np.minimum(var.k[re], 1.e-8)
                #var.k[re+1] = np.minimum(var.k[re+1], 1.e-30)
                
                var.k[re] *= 1e-6
                var.k[re+1] *= 1e-20

        # for sp in vulcan_cfg.condesne_sp:
#             atm.sat_mix[sp] = atm.sat_p[sp]/atm.pco
#             pre_conden = np.copy(var.y[:,species.index(sp)])
#             var.y[:,species.index(sp)] = np.minimum(atm.n_0 * atm.sat_mix[sp], var.y[:,species.index(sp)])
#             # storing the removed species
#             var.y_conden[:,species.index(sp)] += np.abs(pre_conden - var.y[:,species.index(sp)])
        
        # re-scaling every species
        #var.ymix = var.y/np.vstack(np.sum(var.y,axis=1))
    
        return var

    
class ODESolver(object):
    
    def __init__(self): # do I always need to update var, atm, para ?
        
        self.mtol = vulcan_cfg.mtol
        self.atol = vulcan_cfg.atol
        self.non_gas_sp = vulcan_cfg.non_gas_sp       
        
    def diffdf_no_mol(self, y, atm): 
        """
        function of eddy diffusion without molecular diffusion, with zero-flux boundary conditions and non-uniform grids (dzi)
        in the form of Aj*y_j + Bj+1*y_j+1 + Cj-1*y_j-1
        """
        y = y.copy()
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            ysum = np.sum(y[:,atm.exc_conden], axis=1)
        else: ysum = np.sum(y, axis=1)
        # TEST condensation excluding non-gaseous species
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
        vz = atm.vz.copy()
        
        A, B, C = np.zeros(nz), np.zeros(nz), np.zeros(nz)

        A[0] = -1./(dzi[0])*(Kzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[0]     
        B[0] = 1./(dzi[0])*(Kzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[1] 
        C[0] = 0 
        A[nz-1] = -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-1] 
        B[nz-1] = 0 
        C[nz-1] = 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-2] 
        
        # vertical adection with zero-flux B.C. 
        A[0] += -( (vz[0]>0)*vz[0] )/dzi[0]
        B[0] += -( (vz[0]<0)*vz[0] )/dzi[0]
        A[-1] += ( (vz[-1]<0)*vz[-1] )/dzi[-1]
        C[-1] += ( (vz[-1]>0)*vz[-1] )/dzi[-1]
        # vertical adection
        
        for j in range(1,nz-1):  
            dz_ave = 0.5*(dzi[j-1] + dzi[j])
            A[j] = -2./(dzi[j-1] + dzi[j])* ( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j]+ysum[j-1])/2. ) /ysum[j]  
            B[j] = 2./(dzi[j-1] + dzi[j])*Kzz[j]/dzi[j] *(ysum[j+1]+ysum[j])/2. /ysum[j+1]
            C[j] = 2./(dzi[j-1] + dzi[j])*Kzz[j-1]/dzi[j-1] *(ysum[j]+ysum[j-1])/2. /ysum[j-1]
            
            # vertical adection
            A[j] += -( (vz[j]>0)*vz[j] - (vz[j-1]<0)*vz[j-1] )/dz_ave
            B[j] += -( (vz[j]<0)*vz[j] )/dz_ave
            C[j] += ( (vz[j-1]>0)*vz[j-1] )/dz_ave
            # vertical adection
            
        tmp0 = A[0]*y[0] + B[0]*y[1]
        tmp1 = np.ndarray.flatten( (np.vstack(A[1:nz-1])*y[1:(nz-1)] + np.vstack(B[1:nz-1])*y[1+1:(nz-1)+1] + np.vstack(C[1:nz-1])*y[1-1:(nz-1)-1]) )
        tmp2 = (A[nz-1]*y[nz-1] +C[nz-1]*y[nz-2]) 
        diff = np.append(np.append(tmp0, tmp1), tmp2)
        diff = diff.reshape(nz,ni)
        
        if vulcan_cfg.use_topflux == True:
            # Don't forget dz!!! -d phi/ dz
            ### the const flux has no contribution to the jacobian ### 
            diff[-1] += atm.top_flux / atm.dz[-1]
        if vulcan_cfg.use_botflux == True:
            ### the deposition term needs to be included in the jacobian!!!   
            diff[0] += (atm.bot_flux - y[0]*atm.bot_vdep) / atm.dz[0]
            
        return diff
    
    def diffdf(self, y, atm): 
        """
        function of eddy diffusion including molecular diffusion, with zero-flux boundary conditions and non-uniform grids (dzi)
        in the form of Aj*y_j + Bj+1*y_j+1 + Cj-1*y_j-1
        """
        
        y = y.copy()
        
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            ysum = np.sum(y[:,atm.exc_conden], axis=1)
        else: ysum = np.sum(y, axis=1)
        # TEST condensation excluding non-gaseous species
    
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
        vz = atm.vz.copy()
        Dzz = atm.Dzz.copy()
        alpha = atm.alpha.copy()
        Tco = atm.Tco.copy()
        ms = atm.ms.copy()
        Hp = atm.Hp.copy()
        g = vulcan_cfg.g
        
        # define T_1/2 for the molecular diffusion
        Ti = 0.5*(Tco + np.roll(Tco,-1))
        Ti = Ti[:-1] 
        Hpi = 0.5*(Hp + np.roll(Hp,-1))
        Hpi = Hpi[:-1]
        # store Ti and Hpi
        atm.Ti = Ti
        atm.Hpi = Hpi
        
        A, B, C = np.zeros(nz), np.zeros(nz), np.zeros(nz)
        Ai, Bi, Ci = [ np.zeros((nz,ni)) for i in range(3)]
        
        A[0] = -1./(dzi[0])*(Kzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[0]     
        B[0] = 1./(dzi[0])*(Kzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[1] 
        C[0] = 0 
        A[nz-1] = -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-1] 
        B[nz-1] = 0 
        C[nz-1] = 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-2] 
        
        # vertical adection (with closed B.C.) 
        A[0] += -( (vz[0]>0)*vz[0] )/dzi[0]
        B[0] += -( (vz[0]<0)*vz[0] )/dzi[0]
        A[-1] += ( (vz[-1]<0)*vz[-1] )/dzi[-1]
        C[-1] += ( (vz[-1]>0)*vz[-1] )/dzi[-1]
        # vertical adection
        
        # shape of ni-long 1D array
        Ai[0] = -1./(dzi[0])*(Dzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[0] +\
        1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] )  
        Bi[0] = 1./(dzi[0])*(Dzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[1] +\
        1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] )
        Ci[0] = 0 
        Ai[nz-1] = -1./(dzi[-1])*(Dzz[nz-2]/dzi[-1]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-1] \
        -1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] )
        Bi[nz-1] = 0
        Ci[nz-1] = 1./(dzi[-1])*(Dzz[nz-2]/dzi[-1]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-2] \
        -1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] )
        
        for j in range(1,nz-1):
            dz_ave = 0.5*(dzi[j-1] + dzi[j])
            A[j] = -1./dz_ave * ( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j]+ysum[j-1])/2. ) /ysum[j]  
            B[j] = 1./dz_ave * Kzz[j]/dzi[j] *(ysum[j+1]+ysum[j])/2. /ysum[j+1]
            C[j] = 1./dz_ave * Kzz[j-1]/dzi[j-1] *(ysum[j]+ysum[j-1])/2. /ysum[j-1]
            
            # vertical adection
            A[j] += -( (vz[j]>0)*vz[j] - (vz[j-1]<0)*vz[j-1] )/dz_ave
            B[j] += -( (vz[j]<0)*vz[j] )/dz_ave
            C[j] += ( (vz[j-1]>0)*vz[j-1] )/dz_ave
            # vertical adection
            
            # Ai in the shape of nz*ni and Ai[j] in the shape of ni 
            Ai[j] = -1./dz_ave * ( Dzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Dzz[j-1]/dzi[j-1]*(ysum[j]+ysum[j-1])/2. ) /ysum[j]  
            Bi[j] = 1./dz_ave * Dzz[j]/dzi[j] *(ysum[j+1]+ysum[j])/2. /ysum[j+1]
            Ci[j] = 1./dz_ave * Dzz[j-1]/dzi[j-1] *(ysum[j]+ysum[j-1])/2. /ysum[j-1]
            
            Ai[j] += 1./(2.*dz_ave)*( Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] ) \
            - Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+ alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] ) ) #/ysum[j]
            Bi[j] += 1./(2.*dz_ave)* Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] )
            Ci[j] += -1./(2.*dz_ave)* Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] )
 
        tmp0 = (A[0] + Ai[0])*y[0] + (B[0] + Bi[0])*y[1] # shape of ni-long 1D array  
        tmp1 = np.ndarray.flatten( (np.vstack(A[1:nz-1])*y[1:(nz-1)] + np.vstack(B[1:nz-1])*y[1+1:(nz-1)+1] + np.vstack(C[1:nz-1])*y[1-1:(nz-1)-1]) ) 
        tmp1 += np.ndarray.flatten( Ai[1:nz-1]*y[1:(nz-1)] + Bi[1:nz-1]*y[1+1:(nz-1)+1] + Ci[1:nz-1]*y[1-1:(nz-1)-1] ) # shape of (nz-2,ni)
        tmp2 = (A[nz-1] + Ai[nz-1])*y[nz-1] + (C[nz-1] + Ci[nz-1])*y[nz-2]
        diff = np.append(np.append(tmp0, tmp1), tmp2)
        diff = diff.reshape(nz,ni)

        if vulcan_cfg.use_topflux == True:
            # Don't forget dz!!! -d phi/ dz
            ### the const flux has no contribution to the jacobian ### 
            diff[-1] += atm.top_flux / atm.dz[-1]
        if vulcan_cfg.use_botflux == True:
            ### the deposition term needs to be included in the jacobian!!!   
            diff[0] += (atm.bot_flux - y[0]*atm.bot_vdep) / atm.dz[0]
        
        return diff
    
    
    def diffdf_settling(self, y, atm): 
        """
        function of eddy diffusion including molecular diffusion and the settling velocity for particles, with zero-flux boundary conditions and non-uniform grids (dzi)
        in the form of Aj*y_j + Bj+1*y_j+1 + Cj-1*y_j-1
        """
        
        y = y.copy()
        
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            ysum = np.sum(y[:,atm.exc_conden], axis=1)
        else: ysum = np.sum(y, axis=1)
        # TEST condensation excluding non-gaseous species
    
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
        vz = atm.vz.copy()
        Dzz = atm.Dzz.copy()
        vs = atm.vs.copy()
        alpha = atm.alpha.copy()
        Tco = atm.Tco.copy()
        ms = atm.ms.copy()
        Hp = atm.Hp.copy()
        g = vulcan_cfg.g
        
        # define T_1/2 for the molecular diffusion
        Ti = 0.5*(Tco + np.roll(Tco,-1))
        Ti = Ti[:-1] 
        Hpi = 0.5*(Hp + np.roll(Hp,-1))
        Hpi = Hpi[:-1]
        # store Ti and Hpi
        atm.Ti = Ti
        atm.Hpi = Hpi
        
        A, B, C = np.zeros(nz), np.zeros(nz), np.zeros(nz)
        Ai, Bi, Ci = [ np.zeros((nz,ni)) for i in range(3)]
        
        A[0] = -1./(dzi[0])*(Kzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[0]     
        B[0] = 1./(dzi[0])*(Kzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[1] 
        C[0] = 0 
        A[nz-1] = -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-1] 
        B[nz-1] = 0 
        C[nz-1] = 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-2] 
        
        # vertical adection (with closed B.C.) 
        A[0] += -( (vz[0]>0)*vz[0] )/dzi[0]
        B[0] += -( (vz[0]<0)*vz[0] )/dzi[0]
        A[-1] += ( (vz[-1]<0)*vz[-1] )/dzi[-1]
        C[-1] += ( (vz[-1]>0)*vz[-1] )/dzi[-1]
        # vertical adection
        
        # shape of ni-long 1D array
        # Including the settling velocity of the particles
        Ai[0] = -1./(dzi[0])*(Dzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[0] +\
        1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] )  -( (vs[0]>0)*vs[0] )/dzi[0]  
        Bi[0] = 1./(dzi[0])*(Dzz[0]/dzi[0]) *(ysum[1]+ysum[0])/2. /ysum[1] +\
        1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] )  -( (vs[0]<0)*vs[0] )/dzi[0]
        Ci[0] = 0 
        Ai[nz-1] = -1./(dzi[-1])*(Dzz[nz-2]/dzi[-1]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-1] \
        -1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] )  +( (vs[-1]<0)*vs[-1] )/dzi[-1]
        Bi[nz-1] = 0
        Ci[nz-1] = 1./(dzi[-1])*(Dzz[nz-2]/dzi[-1]) *(ysum[nz-1]+ysum[nz-2])/2. /ysum[nz-2] \
        -1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] )  +( (vs[-1]>0)*vs[-1] )/dzi[-1]
        
        for j in range(1,nz-1):
            dz_ave = 0.5*(dzi[j-1] + dzi[j])
            A[j] = -1./dz_ave * ( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j]+ysum[j-1])/2. ) /ysum[j]  
            B[j] = 1./dz_ave * Kzz[j]/dzi[j] *(ysum[j+1]+ysum[j])/2. /ysum[j+1]
            C[j] = 1./dz_ave * Kzz[j-1]/dzi[j-1] *(ysum[j]+ysum[j-1])/2. /ysum[j-1]
            
            # vertical adection
            A[j] += -( (vz[j]>0)*vz[j] - (vz[j-1]<0)*vz[j-1] )/dz_ave
            B[j] += -( (vz[j]<0)*vz[j] )/dz_ave
            C[j] += ( (vz[j-1]>0)*vz[j-1] )/dz_ave
            # vertical adection
            
            # Ai in the shape of nz*ni and Ai[j] in the shape of ni 
            # Including the settling velocity of the particles
            Ai[j] = -1./dz_ave * ( Dzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Dzz[j-1]/dzi[j-1]*(ysum[j]+ysum[j-1])/2. ) /ysum[j] -( (vs[j]>0)*vs[j] + (vs[j-1]<0)*vs[j-1] )/dz_ave
            Bi[j] = 1./dz_ave * Dzz[j]/dzi[j] *(ysum[j+1]+ysum[j])/2. /ysum[j+1]  -( (vs[j]<0)*vs[j] )/dz_ave
            Ci[j] = 1./dz_ave * Dzz[j-1]/dzi[j-1] *(ysum[j]+ysum[j-1])/2. /ysum[j-1]  +( (vs[j-1]>0)*vs[j-1] )/dz_ave
            
            Ai[j] += 1./(2.*dz_ave)*( Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] ) \
            - Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+ alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] ) ) #/ysum[j]
            Bi[j] += 1./(2.*dz_ave)* Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] )
            Ci[j] += -1./(2.*dz_ave)* Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] )
 
        tmp0 = (A[0] + Ai[0])*y[0] + (B[0] + Bi[0])*y[1] # shape of ni-long 1D array  
        tmp1 = np.ndarray.flatten( (np.vstack(A[1:nz-1])*y[1:(nz-1)] + np.vstack(B[1:nz-1])*y[1+1:(nz-1)+1] + np.vstack(C[1:nz-1])*y[1-1:(nz-1)-1]) ) 
        tmp1 += np.ndarray.flatten( Ai[1:nz-1]*y[1:(nz-1)] + Bi[1:nz-1]*y[1+1:(nz-1)+1] + Ci[1:nz-1]*y[1-1:(nz-1)-1] ) # shape of (nz-2,ni)
        tmp2 = (A[nz-1] + Ai[nz-1])*y[nz-1] + (C[nz-1] + Ci[nz-1])*y[nz-2]
        diff = np.append(np.append(tmp0, tmp1), tmp2)
        diff = diff.reshape(nz,ni)

        if vulcan_cfg.use_topflux == True:
            # Don't forget dz!!! -d phi/ dz
            ### the const flux has no contribution to the jacobian ### 
            diff[-1] += atm.top_flux / atm.dz[-1]
        if vulcan_cfg.use_botflux == True:
            ### the deposition term needs to be included in the jacobian!!!   
            diff[0] += (atm.bot_flux - y[0]*atm.bot_vdep) / atm.dz[0]
        
        return diff
        
        
    def jac_tot(self, var, atm):
        """
        jacobian matrix for dn/dt + dphi/dz = P - L (including molecular diffusion)
        zero-flux BC:  1st derivitive of y is zero
        """
        
        y = var.y.copy()
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            ysum = np.sum(y[:,atm.exc_conden], axis=1)
        else: ysum = np.sum(y, axis=1)
        # TEST condensation excluding non-gaseous species
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
        Dzz = atm.Dzz.copy()
        vz = atm.vz.copy()
        alpha = atm.alpha.copy()
        Tco = atm.Tco.copy()
        mu, ms = atm.mu.copy(),  atm.ms.copy()
        g = vulcan_cfg.g
        
        # define T_1/2 for the molecular diffusion
        #Ti = 0.5*(Tco + np.roll(Tco,-1))
        #Ti = Ti[:-1]
        
        Ti = atm.Ti.copy()
        Hpi = atm.Hpi.copy()
        
        
        dfdy = achemjac(y, atm.M, var.k)
        j_indx = []
        
        for j in range(nz):
            j_indx.append( np.arange(j*ni,j*ni+ni) )

        for j in range(1,nz-1): 
            # excluding the buttom and the top cell
            # at j level consists of ni species
            dz_ave = 0.5*(dzi[j-1] + dzi[j]) 
            dfdy[j_indx[j], j_indx[j]] +=  -1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j] -( (vz[j]>0)*vz[j] - (vz[j-1]<0)*vz[j-1] )/dz_ave  
            dfdy[j_indx[j], j_indx[j+1]] += 1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) ) -( (vz[j]<0)*vz[j] )/dz_ave
            dfdy[j_indx[j], j_indx[j-1]] += 1./dz_ave*( Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) ) +( (vz[j-1]>0)*vz[j-1] )/dz_ave
            
            # [j_indx[j], j_indx[j]] has size ni*ni
            dfdy[j_indx[j], j_indx[j]] +=  -1./dz_ave*( Dzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Dzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j]\
            +1./(2.*dz_ave)*( Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] ) \
            - Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] ) )
            dfdy[j_indx[j], j_indx[j+1]] += 1./dz_ave*( Dzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) ) \
            +1./(2.*dz_ave)* Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] ) 
            dfdy[j_indx[j], j_indx[j-1]] += 1./dz_ave*( Dzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) ) \
            -1./(2.*dz_ave)* Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] )

              
        dfdy[j_indx[0], j_indx[0]] += -1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[0]) -( (vz[0]>0)*vz[0] )/dzi[0]
        dfdy[j_indx[0], j_indx[0]] += -1./(dzi[0])*(Dzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[0]) \
        +1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] ) 
        # deposition velocity
        if vulcan_cfg.use_botflux == True: dfdy[j_indx[0], j_indx[0]] += -1.*atm.bot_vdep / atm.dz[0]
        
        dfdy[j_indx[0], j_indx[1]] += 1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[1]) -( (vz[0]<0)*vz[0] )/dzi[0] 
        dfdy[j_indx[0], j_indx[1]] += 1./(dzi[0])*(Dzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[1]) \
        +1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] )

        dfdy[j_indx[nz-1], j_indx[nz-1]] += -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[nz-1]) +( (vz[-1]<0)*vz[-1] )/dzi[-1] 
        dfdy[j_indx[nz-1], j_indx[nz-1]] += -1./(dzi[nz-2])*(Dzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/(2.*ysum[nz-1]) \
        - 1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] )
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] += 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2])* (ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[(nz-1)-1]) +( (vz[-1]>0)*vz[-1] )/dzi[-1]  
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] += 1./(dzi[nz-2])*(Dzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/(2.*ysum[(nz-1)-1]) \
                -1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] )

        return dfdy
    
    def lhs_jac_tot(self, var, atm):      
        """
        directly constructing lhs = 1./(r*h)*sparse.identity(ni*nz) - dfdy
        jacobian matrix for dn/dt + dphi/dz = P - L (including molecular diffusion)
        zero-flux BC:  1st derivitive of y is zero
        """
        y = var.y.copy()
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            ysum = np.sum(y[:,atm.exc_conden], axis=1)
        else: ysum = np.sum(y, axis=1)
        # TEST condensation excluding non-gaseous species
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
        Dzz = atm.Dzz.copy()
        vz = atm.vz.copy()
        alpha = atm.alpha.copy()
        Tco = atm.Tco.copy()
        mu, ms = atm.mu.copy(),  atm.ms.copy()
        g = vulcan_cfg.g

        Ti = atm.Ti.copy()
        Hpi = atm.Hpi.copy()

        # c0 = 1./(r*h) where r = 1. + 1./2.**0.5
        r = 1. + 1./2.**0.5
        c0 = 1./(r*var.dt)
        dfdy = neg_achemjac(y, atm.M, var.k)
        np.fill_diagonal(dfdy, c0 + np.diag(dfdy)) 
        j_indx = []
        
        for j in range(nz):
            j_indx.append( np.arange(j*ni,j*ni+ni) )

        for j in range(1,nz-1):
            # excluding the buttom and the top cell
            # at j level consists of ni species
            dz_ave = 0.5*(dzi[j-1] + dzi[j])
            dfdy[j_indx[j], j_indx[j]] -=  -1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j] -( (vz[j]>0)*vz[j] - (vz[j-1]<0)*vz[j-1] )/dz_ave
            dfdy[j_indx[j], j_indx[j+1]] -= 1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) ) -( (vz[j]<0)*vz[j] )/dz_ave
            dfdy[j_indx[j], j_indx[j-1]] -= 1./dz_ave*( Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) ) +( (vz[j-1]>0)*vz[j-1] )/dz_ave

            # [j_indx[j], j_indx[j]] has size ni*ni
            dfdy[j_indx[j], j_indx[j]] -=  -1./dz_ave*( Dzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Dzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j]\
            +1./(2.*dz_ave)*( Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] ) \
            - Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] ) )
            dfdy[j_indx[j], j_indx[j+1]] -= 1./dz_ave*( Dzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) ) \
            +1./(2.*dz_ave)* Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] )
            dfdy[j_indx[j], j_indx[j-1]] -= 1./dz_ave*( Dzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) ) \
            -1./(2.*dz_ave)* Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] )
    
        dfdy[j_indx[0], j_indx[0]] -= -1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[0]) -( (vz[0]>0)*vz[0] )/dzi[0]
        dfdy[j_indx[0], j_indx[0]] -= -1./(dzi[0])*(Dzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[0]) \
        +1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] ) 
        # deposition velocity
        if vulcan_cfg.use_botflux == True: dfdy[j_indx[0], j_indx[0]] -= -1.*atm.bot_vdep / atm.dz[0]
        
        dfdy[j_indx[0], j_indx[1]] -= 1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[1]) -( (vz[0]<0)*vz[0] )/dzi[0]
        dfdy[j_indx[0], j_indx[1]] -= 1./(dzi[0])*(Dzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[1]) \
        +1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] )

        dfdy[j_indx[nz-1], j_indx[nz-1]] -= -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[nz-1]) +( (vz[-1]<0)*vz[-1] )/dzi[-1]  
        dfdy[j_indx[nz-1], j_indx[nz-1]] -= -1./(dzi[nz-2])*(Dzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/(2.*ysum[nz-1]) \
        - 1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] )
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] -= 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2])* (ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[(nz-1)-1]) +( (vz[-1]>0)*vz[-1] )/dzi[-1]  
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] -= 1./(dzi[nz-2])*(Dzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/(2.*ysum[(nz-1)-1]) \
                -1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] )

        return dfdy
    
    
    def lhs_jac_no_mol(self, var, atm):      
        """
        directly constructing lhs = 1./(r*h)*sparse.identity(ni*nz) - dfdy 
        jacobian matrix for dn/dt + dphi/dz = P - L (WITHOUT molecular diffusion)
        zero-flux BC:  1st derivitive of y is zero
        """
        y = var.y.copy()
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            ysum = np.sum(y[:,atm.exc_conden], axis=1)
        else: ysum = np.sum(y, axis=1)
        # TEST condensation excluding non-gaseous species
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
        vz = atm.vz.copy()
        Tco = atm.Tco.copy()
        mu, ms = atm.mu.copy(),  atm.ms.copy()
        g = vulcan_cfg.g

        r = 1. + 1./2.**0.5
        c0 = 1./(r*var.dt)
        dfdy = neg_achemjac(y, atm.M, var.k)
        np.fill_diagonal(dfdy, c0 + np.diag(dfdy)) 
        j_indx = []
        
        for j in range(nz):
            j_indx.append( np.arange(j*ni,j*ni+ni) )

        for j in range(1,nz-1):
            # excluding the buttom and the top cell
            # at j level consists of ni species
            dz_ave = 0.5*(dzi[j-1] + dzi[j])
            dfdy[j_indx[j], j_indx[j]] -=  -1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j] -( (vz[j]>0)*vz[j] - (vz[j-1]<0)*vz[j-1] )/dz_ave
            dfdy[j_indx[j], j_indx[j+1]] -= 1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) ) -( (vz[j]<0)*vz[j] )/dz_ave
            dfdy[j_indx[j], j_indx[j-1]] -= 1./dz_ave*( Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) ) +( (vz[j-1]>0)*vz[j-1] )/dz_ave
    
        dfdy[j_indx[0], j_indx[0]] -= -1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[0]) -( (vz[0]>0)*vz[0] )/dzi[0]
        # deposition velocity
        if vulcan_cfg.use_botflux == True: dfdy[j_indx[0], j_indx[0]] -= -1.*atm.bot_vdep / atm.dz[0]
        
        dfdy[j_indx[0], j_indx[1]] -= 1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[1]) -( (vz[0]<0)*vz[0] )/dzi[0]

        dfdy[j_indx[nz-1], j_indx[nz-1]] -= -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[nz-1]) +( (vz[-1]<0)*vz[-1] )/dzi[-1] 
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] -= 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2])* (ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[(nz-1)-1]) +( (vz[-1]>0)*vz[-1] )/dzi[-1]  

        return dfdy
    
    def lhs_jac_fix_all_bot(self, var, atm):
        """
        directly constructing lhs = 1./(r*h)*sparse.identity(ni*nz) - dfdy
        jacobian matrix for dn/dt + dphi/dz = P - L (including molecular diffusion)
        Fixed all species BC: all species at bottom (y[0]) remains fixed
        """
        y = var.y.copy()
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            ysum = np.sum(y[:,atm.exc_conden], axis=1)
        else: ysum = np.sum(y, axis=1)
        # TEST condensation excluding non-gaseous species
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
        Dzz = atm.Dzz.copy()
        vz = atm.vz.copy()
        alpha = atm.alpha.copy()
        Tco = atm.Tco.copy()
        mu, ms = atm.mu.copy(),  atm.ms.copy()
        g = vulcan_cfg.g

        Ti = atm.Ti.copy()
        Hpi = atm.Hpi.copy()

        r = 1. + 1./2.**0.5
        c0 = 1./(r*var.dt)
        dfdy = neg_achemjac(y, atm.M, var.k)
        np.fill_diagonal(dfdy, c0 + np.diag(dfdy)) 
        j_indx = []
        
        for j in range(nz):
            j_indx.append( np.arange(j*ni,j*ni+ni) )

        for j in range(1,nz-1):
            # excluding the buttom and the top cell
            # at j level consists of ni species
            dz_ave = 0.5*(dzi[j-1] + dzi[j])
            dfdy[j_indx[j], j_indx[j]] -=  -1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j] -( (vz[j]>0)*vz[j] - (vz[j-1]<0)*vz[j-1] )/dz_ave
            dfdy[j_indx[j], j_indx[j+1]] -= 1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) ) -( (vz[j]<0)*vz[j] )/dz_ave
            dfdy[j_indx[j], j_indx[j-1]] -= 1./dz_ave*( Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) ) +( (vz[j-1]>0)*vz[j-1] )/dz_ave

            # [j_indx[j], j_indx[j]] has size ni*ni
            dfdy[j_indx[j], j_indx[j]] -=  -1./dz_ave*( Dzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Dzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j]\
            +1./(2.*dz_ave)*( Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] ) \
            - Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] ) )
            dfdy[j_indx[j], j_indx[j+1]] -= 1./dz_ave*( Dzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) ) \
            +1./(2.*dz_ave)* Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] )
            dfdy[j_indx[j], j_indx[j-1]] -= 1./dz_ave*( Dzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) ) \
            -1./(2.*dz_ave)* Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] )
    
        # deposition velocity (off with fixed all BC)
        # if vulcan_cfg.use_botflux == True: dfdy[j_indx[0], j_indx[0]] -= -1.*atm.bot_vdep / atm.dz[0]
        
        # Fix bottom BC
        #print (dfdy[:, j_indx[0]])
        dfdy[:, j_indx[0]] = 0.
        
        dfdy[j_indx[0], j_indx[1]] -= 1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[1]) -( (vz[0]<0)*vz[0] )/dzi[0]
        dfdy[j_indx[0], j_indx[1]] -= 1./(dzi[0])*(Dzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[1]) \
        +1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] )

        dfdy[j_indx[nz-1], j_indx[nz-1]] -= -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[nz-1]) +( (vz[-1]<0)*vz[-1] )/dzi[-1] 
        dfdy[j_indx[nz-1], j_indx[nz-1]] -= -1./(dzi[nz-2])*(Dzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/(2.*ysum[nz-1]) \
        - 1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] )
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] -= 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2])* (ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[(nz-1)-1]) +( (vz[-1]>0)*vz[-1] )/dzi[-1]  
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] -= 1./(dzi[nz-2])*(Dzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/(2.*ysum[(nz-1)-1]) \
                -1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] )

        return dfdy
        
    def lhs_jac_no_mol_fix_all_bot(self, var, atm):      
        """
        directly constructing lhs = 1./(r*h)*sparse.identity(ni*nz) - dfdy 
        jacobian matrix for dn/dt + dphi/dz = P - L (WITHOUT molecular diffusion)
        Fixed all species BC: all species at bottom (y[0]) remains fixed
        """
        y = var.y.copy()
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            ysum = np.sum(y[:,atm.exc_conden], axis=1)
        else: ysum = np.sum(y, axis=1)
        # TEST condensation excluding non-gaseous species
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
        vz = atm.vz.copy()
        Tco = atm.Tco.copy()
        mu, ms = atm.mu.copy(),  atm.ms.copy()
        g = vulcan_cfg.g

        r = 1. + 1./2.**0.5
        c0 = 1./(r*var.dt)
        dfdy = neg_achemjac(y, atm.M, var.k)
        np.fill_diagonal(dfdy, c0 + np.diag(dfdy)) 
        j_indx = []
        
        for j in range(nz):
            j_indx.append( np.arange(j*ni,j*ni+ni) )

        for j in range(1,nz-1):
            # excluding the buttom and the top cell
            # at j level consists of ni species
            dz_ave = 0.5*(dzi[j-1] + dzi[j])
            dfdy[j_indx[j], j_indx[j]] -=  -1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j] -( (vz[j]>0)*vz[j] - (vz[j-1]<0)*vz[j-1] )/dz_ave
            dfdy[j_indx[j], j_indx[j+1]] -= 1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) ) -( (vz[j]<0)*vz[j] )/dz_ave
            dfdy[j_indx[j], j_indx[j-1]] -= 1./dz_ave*( Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) ) +( (vz[j-1]>0)*vz[j-1] )/dz_ave
    
        #dfdy[j_indx[0], j_indx[0]] -= -1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[0]) -( (vz[0]>0)*vz[0] )/dzi[0]
        # deposition velocity (off with fixed all BC)
        # if vulcan_cfg.use_botflux == True: dfdy[j_indx[0], j_indx[0]] -= -1.*atm.bot_vdep / atm.dz[0]
        
        # Fix bottom BC
        dfdy[:, j_indx[0]] = 0.
        
        dfdy[j_indx[0], j_indx[1]] -= 1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[1]) -( (vz[0]<0)*vz[0] )/dzi[0]

        dfdy[j_indx[nz-1], j_indx[nz-1]] -= -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[nz-1]) +( (vz[-1]<0)*vz[-1] )/dzi[-1] 
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] -= 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2])* (ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[(nz-1)-1]) +( (vz[-1]>0)*vz[-1] )/dzi[-1]  

        return dfdy
        
    def lhs_jac_settling(self, var, atm):      
        """
        directly constructing lhs = 1./(r*h)*sparse.identity(ni*nz) - dfdy
        jacobian matrix for dn/dt + dphi/dz = P - L (including molecular diffusion and gravitation settling for particles)
        zero-flux BC:  1st derivitive of y is zero
        """
        y = var.y.copy()
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            ysum = np.sum(y[:,atm.exc_conden], axis=1)
        else: ysum = np.sum(y, axis=1)
        # TEST condensation excluding non-gaseous species
        dzi = atm.dzi.copy()
        Kzz = atm.Kzz.copy()
        Dzz = atm.Dzz.copy()
        vz = atm.vz.copy()
        vs = atm.vs.copy()
        alpha = atm.alpha.copy()
        Tco = atm.Tco.copy()
        mu, ms = atm.mu.copy(),  atm.ms.copy()
        g = vulcan_cfg.g

        Ti = atm.Ti.copy()
        Hpi = atm.Hpi.copy()

        # c0 = 1./(r*h) where r = 1. + 1./2.**0.5
        r = 1. + 1./2.**0.5
        c0 = 1./(r*var.dt)
        dfdy = neg_achemjac(y, atm.M, var.k)
        np.fill_diagonal(dfdy, c0 + np.diag(dfdy)) 
        j_indx = []
        
        for j in range(nz):
            j_indx.append( np.arange(j*ni,j*ni+ni) )

        for j in range(1,nz-1):
            # excluding the buttom and the top cell
            # at j level consists of ni species
            dz_ave = 0.5*(dzi[j-1] + dzi[j])
            dfdy[j_indx[j], j_indx[j]] -=  -1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j] -( (vz[j]>0)*vz[j] - (vz[j-1]<0)*vz[j-1] )/dz_ave
            dfdy[j_indx[j], j_indx[j+1]] -= 1./dz_ave*( Kzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) ) -( (vz[j]<0)*vz[j] )/dz_ave
            dfdy[j_indx[j], j_indx[j-1]] -= 1./dz_ave*( Kzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) ) +( (vz[j-1]>0)*vz[j-1] )/dz_ave

            # [j_indx[j], j_indx[j]] has size ni*ni
            dfdy[j_indx[j], j_indx[j]] -=  -1./dz_ave*( Dzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/2. + Dzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/2. ) /ysum[j]\
            +1./(2.*dz_ave)*( Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] ) \
            - Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] ) )  -( (vs[j]>0)*vs[j] - (vs[j-1]<0)*vs[j-1] )/dz_ave
            dfdy[j_indx[j], j_indx[j+1]] -= 1./dz_ave*( Dzz[j]/dzi[j]*(ysum[j+1]+ysum[j])/(2.*ysum[j+1]) ) \
            +1./(2.*dz_ave)* Dzz[j]*(-1./Hpi[j]+ms*g/(Navo*kb*Ti[j])+alpha/Ti[j]*(Tco[j+1]-Tco[j])/dzi[j] )  -( (vs[j]<0)*vs[j] )/dz_ave
            dfdy[j_indx[j], j_indx[j-1]] -= 1./dz_ave*( Dzz[j-1]/dzi[j-1]*(ysum[j-1]+ysum[j])/(2.*ysum[j-1]) ) \
            -1./(2.*dz_ave)* Dzz[j-1]*(-1./Hpi[j-1]+ms*g/(Navo*kb*Ti[j-1])+alpha/Ti[j-1]*(Tco[j]-Tco[j-1])/dzi[j-1] )  +( (vs[j-1]>0)*vs[j-1] )/dz_ave
    
        dfdy[j_indx[0], j_indx[0]] -= -1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[0]) -( (vz[0]>0)*vz[0] )/dzi[0]
        dfdy[j_indx[0], j_indx[0]] -= -1./(dzi[0])*(Dzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[0]) \
        +1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] )  -( (vs[0]>0)*vs[0] )/dzi[0]
        # deposition velocity
        if vulcan_cfg.use_botflux == True: dfdy[j_indx[0], j_indx[0]] -= -1.*atm.bot_vdep / atm.dz[0]
        
        dfdy[j_indx[0], j_indx[1]] -= 1./(dzi[0])*(Kzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[1]) -( (vz[0]<0)*vz[0] )/dzi[0] 
        dfdy[j_indx[0], j_indx[1]] -= 1./(dzi[0])*(Dzz[0]/dzi[0]) * (ysum[1]+ysum[0])/(2.*ysum[1]) \
        +1./(dzi[0])* Dzz[0]/2.*(-1./Hpi[0]+ms*g/(Navo*kb*Ti[0])+alpha/Ti[0]*(Tco[1]-Tco[0])/dzi[0] ) -( (vs[0]<0)*vs[0] )/dzi[0]

        dfdy[j_indx[nz-1], j_indx[nz-1]] -= -1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2]) *(ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[nz-1]) +( (vz[-1]<0)*vz[-1] )/dzi[-1]  
        dfdy[j_indx[nz-1], j_indx[nz-1]] -= -1./(dzi[nz-2])*(Dzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/(2.*ysum[nz-1]) \
        - 1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] ) +( (vs[-1]<0)*vs[-1] )/dzi[-1]
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] -= 1./(dzi[nz-2])*(Kzz[nz-2]/dzi[nz-2])* (ysum[(nz-1)-1]+ysum[nz-1])/(2.*ysum[(nz-1)-1]) +( (vz[-1]>0)*vz[-1] )/dzi[-1]   
        dfdy[j_indx[nz-1], j_indx[(nz-1)-1]] -= 1./(dzi[nz-2])*(Dzz[nz-2]/dzi[nz-2]) *(ysum[nz-1]+ysum[nz-2])/(2.*ysum[(nz-1)-1]) \
                -1./(dzi[-1])* Dzz[-1]/2.*(-1./Hpi[-1]+ms*g/(Navo*kb*Ti[-1])+alpha/Ti[-1]*(Tco[-1]-Tco[-2])/dzi[-1] ) +( (vs[-1]>0)*vs[-1] )/dzi[-1]

        return dfdy
            
        
        
    def clip(self, var, para, atm, pos_cut = vulcan_cfg.pos_cut, nega_cut = vulcan_cfg.nega_cut):
        '''
        function to clip samll and negative values
        and to calculate the particle loss
        '''
        y, ymix = var.y, var.ymix.copy()
         
        para.small_y += np.abs(np.sum(y[np.logical_and(y<pos_cut, y>=0)]))
        para.nega_y += np.abs(np.sum(y[np.logical_and(y>nega_cut, y<=0)]))
        y[np.logical_and(y<pos_cut, y>=nega_cut)] = 0.
        
        # Also setting y=0 when ymix<mtol
        y[np.logical_and(ymix<self.mtol, y<0)] = 0.
        
        var = self.loss(var)
        
        # store y and ymix
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            var.y, var.ymix = y, var.y/np.vstack(np.sum(var.y[:,atm.exc_conden],axis=1)) 
        else: var.y, var.ymix = y, y/np.vstack(np.sum(y,axis=1))
        # TEST condensation excluding non-gaseous species
        
        return var , para
        
    def loss(self, data_var): 
        
        y = data_var.y
        atom_list = vulcan_cfg.atom_list
        
        # changed atom_tot to dictionary atom_sum
        atom_sum = data_var.atom_sum
        
        for atom in atom_list:
            data_var.atom_sum[atom] = np.sum([compo[compo_row.index(species[i])][atom] * data_var.y[:,i] for i in range(ni)])
            data_var.atom_loss[atom] = (data_var.atom_sum[atom] - data_var.atom_ini[atom])/data_var.atom_ini[atom]

        return data_var
        
    def step_ok(self, var, para, loss_eps = vulcan_cfg.loss_eps, rtol = vulcan_cfg.rtol):
        if np.all(var.y>=0) and np.amax( np.abs( np.fromiter(var.atom_loss.values(),float) - np.fromiter(var.atom_loss_prev.values(),float) ) )<loss_eps and para.delta<=rtol:
            return True
        else:
            return False
            
    def step_reject(self, var, para, loss_eps = vulcan_cfg.loss_eps, rtol = vulcan_cfg.rtol):
        
        if para.delta > rtol: # truncation error larger than the tolerence value
            para.delta_count += 1
            
        elif np.any(var.y < 0):             
            para.nega_count += 1
            if vulcan_cfg.use_print_prog == True:
                self.print_nega(var,para) # print the info for the negative solutions (where y < 0)
            # print input: y, t, count, dt
            

        else: # meaning np.amax( np.abs( np.abs(y_loss) - np.abs(loss_prev) ) )<loss_eps
            para.loss_count +=1
            if vulcan_cfg.use_print_prog == True:
                self.print_lossBig(para)
        
        
        var = self.reset_y(var) # reset y and dt to the values at previous step
        
        if var.dt < vulcan_cfg.dt_min:
            var.dt = vulcan_cfg.dt_min
            var.y[var.y<0] = 0. # clipping of negative values
            print ('Keep producing negative values! Clipping negative solutions and moving on!')
            return True
        
        return False
            
    def reset_y(self, var, dt_reduc = vulcan_cfg.dt_var_min):
        '''
        reset y and reduce dt by dt_reduc
        '''
        
        # reset and store y and dt
        var.y = var.y_prev
        var.dt *= dt_reduc
        # var.dt = np.maximum(var.dt, vulcan_cfg.dt_min)

        return var
        
    def print_nega(self, data_var, data_para): 
        
        nega_i = np.where(data_var.y<0)
        print ('Negative y at time ' + str("{:.2e}".format(data_var.t)) + ' and step: ' + str(data_para.count) )
        print ('Negative values:' + str(data_var.y[data_var.y<0]) )
        print ('from levels: ' + str(nega_i[0]) )
        print ('species: ' + str([species[_] for _ in nega_i[1]]) )
        print ('dt= ' + str(data_var.dt))
        print ('...reset dt to dt*0.2...')
        print ('------------------------------------------------------------------')
    
    def print_lossBig(self, para):
        
        print ('partical conservation is violated too large (by numerical errors)')
        print ('at step: ' + str(para.count))
        print ('------------------------------------------------------------------')
        
    def thomas_vec(a, b, c, d): 
        '''
        Thomas vectorized solver, a b c d refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
        d is a matrix
        not used in this current version
        '''
        # number of equations
        nf = len(a) 
        aa, bb, cc, dd = map(np.copy, (a, b, c, d))  
        # d needs to reshape
        dd = dd.reshape(nf,-1)
        #C' and D'
        cp = [cc[0]/bb[0]]; dp = [dd[0]/bb[0]]  
        x = np.zeros((nf, np.shape(dd)[1]))
  
        for i in range(1, nf-1):
            cp.append( cc[i]/(bb[i] - aa[i]*cp[i-1]) ) 
            dp.append( (dd[i] - aa[i]*dp[i-1])/(bb[i] - aa[i]*cp[i-1]) )  
   
        dp.append( (dd[(nf-1)] - aa[(nf-1)]*dp[(nf-1)-1])/(bb[(nf-1)] - aa[(nf-1)]*cp[(nf-1)-1]) ) # nf-1 is the last element
        x[nf-1] = dp[nf-1]/1
        for i in range(nf-2, -1, -1):
            x[i] = dp[i] - cp[i]*x[i+1]
        
        return x
    
    ### photo-calculation starts from here
    
    # def tot_cross(self, var):
#     ''' compute the total cross section from all species '''
#         for sp in vulcan_cfg.photo_sp:
#             cross[sp]
    

    def compute_tau(self, var, atm):
        ''' compute the optical depth '''
        
        # reset to zero
        var.tau.fill(0)
        
        for j in range(nz-1,-1,-1):
            
            for sp in var.photo_sp:
            # summing over all photo species
                var.tau[j] += var.y[j,species.index(sp)] * atm.dz[j] * var.cross[sp] # only the j-th laye
            
            for sp in vulcan_cfg.scat_sp: # scat_sp are not necessary photo_sp, e.g. He
                var.tau[j] += var.y[j,species.index(sp)] * atm.dz[j] * var.cross_scat[sp]
            # adding the layer above at the end of species loop   
            var.tau[j] += var.tau[j+1]
               
    # Lines like chi = zeta_m**2*tran**2 - zeta_p**2 doing large np 2D array multiplication
    # can be sped up with cython           
    def compute_flux(self, var, atm): # Vectorise this loop!  
        # change it to stagerred grids
        # top: stellar flux
        # bottom BC: zero upcoming flux
        
        # Note!!! Matej's mu is defined in the outgoing hemisphere so his mu<0
        # My cos[sl_angle] is always 0<=mu<=1
        # Converting my mu to Matej's mu (e.g. 45 deg -> 135 deg)
      
        mu_ang = -1.*np.cos(vulcan_cfg.sl_angle)
        edd = vulcan_cfg.edd
        tau = var.tau
        
        # delta_tau (length nz) is used in the transmission function
        delta_tau = tau - np.roll(tau,-1,axis=0) # np.roll(tau,-1,axis=0) are the upper layers
        delta_tau = delta_tau[:-1]
        
        
        # single-scattering albedo
        nbins = len(var.bins)
        tot_abs, tot_scat = np.zeros((nz, nbins)), np.zeros((nz, nbins))
        for sp in var.photo_sp: 
            tot_abs += np.vstack(var.ymix[:,species.index(sp)])*var.cross[sp] # nz * nbins
        for sp in vulcan_cfg.scat_sp: 
            tot_scat += np.vstack(var.ymix[:,species.index(sp)])*var.cross_scat[sp]

        total = tot_abs + tot_scat
        
        w0 = tot_scat  / (tot_abs + tot_scat) # 2D: nz * nbins
        # tot_abs + tot_scat can be zero when certain gas (e.g. H2) does not exist
        
        # Replace nan with zero and inf with very large numbers
        w0 = np.nan_to_num(w0)
        
        # to avoit w0=1
        w0 = np.minimum(w0,1.-1.E-8)

        # sflux: the direct beam; dflux: diffusive flux
        ''' Beer's law for the intensity'''
        var.sflux = var.sflux_top *  np.exp(-1.*tau/np.cos(vulcan_cfg.sl_angle) ) 
        # converting the intensity to flux for the raditive transfer calculation
        dir_flux = var.sflux*np.cos(vulcan_cfg.sl_angle) # multiplied by the zenith angle for calculating the diffuse flux
        
        # scattering
        # the transmission function (length nz)
        if ag0 == 0: # to save memory
            tran = np.exp( -1./edd *(1.- w0)**0.5 * delta_tau ) # 2D: nz * nbins
            zeta_p = 0.5*( 1. + (1.-w0)**0.5 )
            zeta_m = 0.5*( 1. - (1.-w0)**0.5 )
            ll = -1.*w0/( 1./mu_ang**2 -1./edd**2 *(1.-w0) )
            g_p = 0.5*( ll*(1./edd+1./mu_ang) )
            g_m = 0.5*( ll*(1./edd-1./mu_ang) ) 

        else:
            tran = np.exp( -1./edd *( (1.- w0*ag0)*(1.- w0) )**0.5 * delta_tau )
            zeta_p = 0.5*( 1. + ((1.-w0)/(1-w0*ag0))**0.5 )
            zeta_m = 0.5*( 1. - ((1.-w0)/(1-w0*ag0))**0.5 ) 
            ll = ( (1.-w0)*(1-w0*ag0) - 1.)/( 1./mu_ang**2 -1./edd**2 *(1.-w0)*(1-w0*ag0) )
            g_p = 0.5*( ll*(1./edd+1/(mu_ang*(1.-w0*ag0))) + w0*ag0*mu_ang/(1.-w0*ag0)  )
            g_m = 0.5*( ll*(1./edd-1/(mu_ang*(1.-w0*ag0))) - w0*ag0*mu_ang/(1.-w0*ag0)  )
        
        
        # to avoit zero denominator
        ll = np.minimum(ll, 1.e10)
        ll = np.maximum(ll, -1.e10)
        

        # 2D: nz * nbins
        chi = zeta_m**2*tran**2 - zeta_p**2
        xi = zeta_p*zeta_m*(1.-tran**2)
        phi = (zeta_m**2-zeta_p**2)*tran
        
        # 2D: nz * nbins
        i_u = phi*g_p*dir_flux[:-1] - (xi*g_m+chi*g_p)*dir_flux[1:]
        i_d = phi*g_m*dir_flux[1:] - (chi*g_m+xi*g_p)*dir_flux[:-1]
        # sflux[1:] are all the layers above and sflux[:-1] are all the layers abelow
        
        var.zeta_m = zeta_m
        var.zeta_p = zeta_p
        var.tran = tran
        
        
        
        #starting recording time
        #start_time = timeit.default_timer()

        
        # propagating downward layer by layer and then upward
        # var.dflux_d and var.dflux_p are defined at the interfaces (staggerred)
        # the rest is defined in the center of the layer
        for j in range(nz-1,-1,-1): # dflux_d goes from the second top interface (nz+1 interfaces) 
            var.dflux_d[j] = 1./chi[j]*(phi[j]*var.dflux_d[j+1] - xi[j]*var.dflux_u[j] + i_d[j]/mu_ang )
        for j in range(1,nz+1):        
            var.dflux_u[j] = 1./chi[j-1]*(phi[j-1]*var.dflux_u[j-1] - xi[j-1]*var.dflux_d[j] + i_u[j-1]/mu_ang )
        

        #print ("time passed...")
        #print (timeit.default_timer() - start_time)
        
        
        # old
        # # the average intensity (not flux!) of the direct beam
#         ave_int = 0.5*( var.sflux[:-1] + var.sflux[1:])
#         tot_int = (ave_int + 0.5*(var.dflux_u[:-1] + var.dflux_u[1:] + var.dflux_d[1:] + var.dflux_d[:-1]) )/edd
#         # devided by the Eddington coefficient to recover the intensity
        
        
        # the average flux from the direct beam
        # !!! WITHOUT multiplied by the cos zenith angle (flux per unit area perpendicular to the direction of propagationat) !!! 
        ave_dir_flux = 0.5*( var.sflux[:-1] + var.sflux[1:]) 
        # devided by the Eddington coefficient to recover the intensity then multiplied by 4pi to get the integrated flux
        tot_flux = ave_dir_flux + 0.5*(var.dflux_u[:-1] + var.dflux_u[1:] + var.dflux_d[1:] + var.dflux_d[:-1])/edd 
        
    
        # ### Debug
        
        #var.ave_int = ave_int
        
        # var.ll = ll
        # var.chi=chi
        # var.phi=phi
        # var.xi = xi
        #
        # var.i_u = i_u
        # var.i_d = i_d
        # var.w0 = w0
        # var.tot_abs = tot_abs
        # var.tot_scat = tot_scat
        # var.tran = tran
        # var.delta_tau = delta_tau
        
        ### Debug
        if np.any(tot_flux< -1.e-20):
            print (tot_flux[tot_flux<-1.e-20])
            raise IOError ('\nNegative diffusive flux! ')
        
         
        # store the previous actinic flux into prev_aflux
        var.prev_aflux = np.copy(var.aflux)
        # converting to the actinic flux and storing the current flux
        var.aflux = tot_flux / (hc/var.bins)
        # the change of the actinic flux
        var.aflux_change = np.nanmax( np.abs(var.aflux-var.prev_aflux)[var.aflux>vulcan_cfg.flux_atol]/var.aflux[var.aflux>vulcan_cfg.flux_atol] )
        
        #print ('aflux change: ' + '{:.4E}'.format(var.aflux_change) )
        
        
    def compute_J(self, var, atm): # the vectorised version
        # change it to trapezoid integral!!!
        flux = var.aflux
        
        #cross = var.cross
        cross = var.cross_J
        
        bins = var.bins
        n_branch = var.n_branch
        wavelen = var.wavelen
        br_ratio = var.br_ratio

        # reset to zeros every time
        var.J_sp = dict([( (sp,bn) , np.zeros(nz)) for sp in var.photo_sp for bn in range(n_branch[sp]+1) ])

        for sp in var.photo_sp:
            # shape: flux (nz,nbin) cross (nbin)


            # I want to parallelize this bit
            # for n in range(var.nbin):
            # I want to parallelize this bit

            # the number of wavelength regions
            wl_num = len(wavelen[sp])

            # convert to actinic flux *1/(hc/ld)
            if wl_num == 0:
                for nbr in range(1, n_branch[sp]+1):
                    # axis=1 is to sum over all wavelength 
                    var.J_sp[(sp, nbr)] = np.sum( br_ratio[sp][0][nbr-1] * flux * cross[sp] * var.dbin, axis=1)
                    var.J_sp[(sp, nbr)] -= 0.5* br_ratio[sp][0][nbr-1] * flux[:,0] * cross[sp][0] * var.dbin +\
                    0.5* br_ratio[sp][0][nbr-1] * flux[:,-1] * cross[sp][-1] * var.dbin

            elif wl_num == 1:
                for nbr in range(1, n_branch[sp]+1):
                    var.J_sp[(sp, nbr)] = np.sum( br_ratio[sp][0][nbr-1] * (bins<wavelen[sp]) * flux * cross[sp] * var.dbin + br_ratio[sp][1][nbr-1] * (bins>=wavelen[sp]) * flux * cross[sp] * var.dbin, axis=1)
                    var.J_sp[(sp, nbr)] -= 0.5* ( br_ratio[sp][0][nbr-1] * (bins[0]<wavelen[sp]) * flux[:,0] * cross[sp][0] * var.dbin +\
                        br_ratio[sp][1][nbr-1] * (bins[0]>=wavelen[sp]) * flux[:,0] * cross[sp][0] * var.dbin + br_ratio[sp][0][nbr-1] * (bins[-1]<wavelen[sp]) * flux[:,-1] * cross[sp][-1] * var.dbin +\
                        br_ratio[sp][1][nbr-1] * (bins[-1]>=wavelen[sp]) * flux[:,-1] * cross[sp][-1] * var.dbin )
 
            elif wl_num >= 2:
                for nbr in range(1, n_branch[sp]+1):
                    # the first wavelength region
                    var.J_sp[(sp, nbr)] = np.sum( br_ratio[sp][0][nbr-1] * (bins<wavelen[sp][0]) * flux * cross[sp] * var.dbin ,axis=1)
                    # the last wavelength region
                    var.J_sp[(sp, nbr)] += np.sum( br_ratio[sp][-1][nbr-1] * (bins>=wavelen[sp][-1]) * flux * cross[sp] * var.dbin, axis=1)
                    
                    # trapezoid integral
                    var.J_sp[(sp, nbr)] -= 0.5*(br_ratio[sp][0][nbr-1] * (bins[0]<wavelen[sp][0]) * flux[:,0] * cross[sp][0] * var.dbin + br_ratio[sp][0][nbr-1] * (bins[0]<wavelen[sp][0]) * flux[:,-1] * cross[sp][-1] * var.dbin)
                    var.J_sp[(sp, nbr)] -= 0.5*(br_ratio[sp][-1][nbr-1] * (bins[-1]>=wavelen[sp][-1]) * flux[:,0] * cross[sp][0] * var.dbin + br_ratio[sp][-1][nbr-1] * (bins[-1]>=wavelen[sp][-1]) * flux[:,-1] * cross[sp][-1] * var.dbin)
                    
                    for region_num in range(1,wl_num):
                        
                        var.J_sp[(sp, nbr)] += np.sum(br_ratio[sp][region_num][nbr-1] * (np.array(bins>=wavelen[sp][region_num-1]) & np.array(bins<wavelen[sp][region_num])) * flux * cross[sp] * var.dbin, axis=1)
                        # trapezoid integral
                        var.J_sp[(sp, nbr)] -= 0.5*(br_ratio[sp][region_num][nbr-1] * (bins[0]>=wavelen[sp][region_num-1] and bins[0]<wavelen[sp][region_num]) * flux[:,0] * cross[sp][0] * var.dbin +\
                        br_ratio[sp][region_num][nbr-1] * (bins[-1]>=wavelen[sp][region_num-1] and bins[-1]<wavelen[sp][region_num]) * flux[:,-1] * cross[sp][-1] * var.dbin)

            else:
                raise IOError ('\n Too many wavelength switches! Check the photo-network file.')

            # end of the loop: for sp in var.photo_sp:
            # if sp in lya_ratio.keys():
                # define the range of Ly-alpha


            # 0 is the total dissociation rate
            # summing all branches
            for nbr in range(1, n_branch[sp]+1):
                var.J_sp[(sp, 0)] += var.J_sp[(sp, nbr)]
                # incoperating J into rate coefficients
                if var.pho_rate_index[(sp, nbr)] not in vulcan_cfg.remove_list:
                    var.k[ var.pho_rate_index[(sp, nbr)]  ] = var.J_sp[(sp, nbr)]
          

class Ros2(ODESolver):
    '''
    class inheritance from ODEsolver for 2nd order Rosenbrock solver 
    '''
    def __init__(self):
        ODESolver.__init__(self)
         
    def store_bandM(self, a, nb, nn):
        """
        store block-tridiagonal matrix(bandwidth=1) into diagonal ordered form 
        (http://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_banded.html) 
        a : square block-tridiagonal matirx
        nb: size of the block matrix (number of species)
        nn: number of the block matrices (number of layers)
        """
    
        # band width (treat block-banded as banded matrix)
        bw = 2*nb-1 
        ab = np.zeros((2*bw+1,nb*nn))

        # first 2 columns
        for i in range(0,2*nb):
            ab[-(2*nb+i):,i] = a[0:2*nb+i,i]
    
        # middle
        for i in range(2*nb, nn*nb-2*nb):
            ab[:,i] = a[(i-2*nb+1):(i-2*nb+1)+(2*bw+1),i] 
    
        # last 2 columns
        for ne,i in enumerate(range(nn*nb-2*nb,nn*nb)):
            ab[:(2*bw+1 -ne),i] = a[-(2*bw+1 -ne):,i]
            
        return (ab, bw)
    
    def solver(self, var, atm, para):
        """
        2nd order Rosenbrock [Verwer et al. 1997] with banded-matrix solver
        with switches to include the molecular diffusion or not
        """
        
        y, ymix, h, k = var.y, var.ymix, var.dt, var.k
        M, dzi, Kzz = atm.M, atm.dzi, atm.Kzz
            
        if vulcan_cfg.use_moldiff == True and vulcan_cfg.use_settling == False:
            diffdf = self.diffdf
            jac_tot = self.lhs_jac_tot
        elif vulcan_cfg.use_moldiff == True and vulcan_cfg.use_settling == True:
            diffdf = self.diffdf_settling
            jac_tot = self.lhs_jac_settling
        else:
            diffdf = self.diffdf_no_mol
            jac_tot = self.lhs_jac_no_mol
    
        r = 1. + 1./2.**0.5

        df = chemdf(y,M,k).flatten() + diffdf(y, atm).flatten()
        lhs = jac_tot(var, atm)
        
        lhs_b, bw = self.store_bandM(lhs,ni,nz)
        k1_flat = scipy.linalg.solve_banded((bw,bw),lhs_b,df)
        
        k1 = k1_flat.reshape(y.shape)
        
        yk2 = y + k1/r
        df = chemdf(yk2,M,k).flatten() + diffdf(yk2, atm).flatten()

        rhs = df - 2./(r*h)*k1_flat
        k2 = scipy.linalg.solve_banded((bw,bw),lhs_b,rhs)
        k2 = k2.reshape(y.shape)
        
        sol = y + 3./(2.*r)*k1 + 1/(2.*r)*k2
        if vulcan_cfg.use_fix_sp_bot: # if use_fix_sp_bot = {} (empty), it returns false
            for sp in vulcan_cfg.use_fix_sp_bot.keys():
                sol[0,species.index(sp)] = vulcan_cfg.use_fix_sp_bot[sp]*atm.n_0[0]
        # surface sink (and top BC for the particles?)
        if vulcan_cfg.use_condense == True:
            #pass
            for sp in self.non_gas_sp:
                sol[0,species.index(sp)] = 0
                # TEST top sink
                # sol[-1,species.index(sp)] = 0
                
        delta = np.abs(sol-yk2)
        delta[ymix < self.mtol] = 0
        delta[sol < self.atol] = 0

        # TEST condensation
        if vulcan_cfg.use_condense == True:
            for sp in vulcan_cfg.non_gas_sp:
                delta[:,species.index(sp)] = 0
        # TEST condensation
        delta = np.amax( delta[sol>0]/sol[sol>0] )

        var.y = sol
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            var.ymix = var.y/np.vstack(np.sum(var.y[:,atm.exc_conden],axis=1)) 
        else: var.ymix = var.y/np.vstack(np.sum(var.y,axis=1))
        # TEST condensation excluding non-gaseous species
        para.delta = delta    
        
        return var, para
        
    def solver_fix_all_bot(self, var, atm, para):
        """
        2nd order Rosenbrock [Verwer et al. 1997] with banded-matrix solver
        with switches to include the molecular diffusion or not
        """
        
        y, ymix, h, k = var.y, var.ymix, var.dt, var.k
        M, dzi, Kzz = atm.M, atm.dzi, atm.Kzz
        
        # store the fixed bottom level
        bottom = np.copy(ymix[0])
        
        if vulcan_cfg.use_moldiff == True:
            diffdf = self.diffdf
            jac_tot = self.lhs_jac_fix_all_bot
        else:
            diffdf = self.diffdf_no_mol
            jac_tot = self.lhs_jac_no_mol_fix_all_bot
    
        r = 1. + 1./2.**0.5

        df = chemdf(y,M,k).flatten() + diffdf(y, atm).flatten()
        lhs = jac_tot(var, atm)
        
        lhs_b, bw = self.store_bandM(lhs,ni,nz)
        k1_flat = scipy.linalg.solve_banded((bw,bw),lhs_b,df)
        
        k1 = k1_flat.reshape(y.shape)
        
        yk2 = y + k1/r
        df = chemdf(yk2,M,k).flatten() + diffdf(yk2, atm).flatten()
        
        rhs = df - 2./(r*h)*k1_flat
        k2 = scipy.linalg.solve_banded((bw,bw),lhs_b,rhs)
        k2 = k2.reshape(y.shape)
        
        sol = y + 3./(2.*r)*k1 + 1/(2.*r)*k2  
        # for sp in :
#             sol[0,species.index(sp)] = bottom*atm.n_0[0,species.index(sp)]
        sol[0] = bottom*atm.n_0[0] 
        
        delta = np.abs(sol-yk2)
        delta[ymix < self.mtol] = 0
        delta[sol < self.atol] = 0

        # TEST condensation
        # if vulcan_cfg.use_condense == True:
#             for sp in vulcan_cfg.non_gas_sp:
#                 delta[:,species.index(sp)] = 0
        # TEST condensation
        delta = np.amax( delta[sol>0]/sol[sol>0] )

        var.y = sol
        
        # TEST condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            var.ymix = var.y/np.vstack(np.sum(var.y[:,atm.exc_conden],axis=1)) 
        else: var.ymix = var.y/np.vstack(np.sum(var.y,axis=1))
        # TEST condensation excluding non-gaseous species

        para.delta = delta    
        
        return var, para   
      
    
    def naming_solver(self, para):
        
        # if vulcan_cfg.use_fix_all_bot == True:
        #     if vulcan_cfg.use_moldiff == True: print ('Use fixed bottom BC and molecular diffusion.')
        #     else: print ('Use fixed bottom BC and No molecular diffusion.')
        #     para.solver_str = 'solver_fix_all_bot'
            
        #else:
        if vulcan_cfg.use_moldiff == True: print ('Include molecular diffusion.')
        else: print ('No molecular diffusion.')
        para.solver_str = 'solver'
        
        
    def one_step(self, var, atm, para):

        while True:
            
           var, para =  getattr(self, para.solver_str)(var, atm, para)
           
           # clipping small negative values and also calculating atomic loss (atom_loss)  
           var , para = self.clip(var, para, atm) 
            
           if self.step_ok(var, para): break
           elif self.step_reject(var, para): break # giving up and moving on
               
        return var, para                    
        
    def step_size(self, var, para, dt_var_min = vulcan_cfg.dt_var_min, dt_var_max = vulcan_cfg.dt_var_max, dt_min = vulcan_cfg.dt_min, dt_max = vulcan_cfg.dt_max):  
        """
        step-size control by delta(truncation error) for the Rosenbrock method
        """
        h = var.dt
        delta = para.delta
        rtol = vulcan_cfg.rtol
               
        if delta==0: delta = 0.01*rtol
        h_factor = 0.9*(rtol/delta)**0.5 # 0.9 is simply a safety factor
        h_factor = np.maximum(h_factor, dt_var_min)    
        h_factor = np.minimum(h_factor, dt_var_max)    
        
        h *= h_factor
        h = np.maximum(h, dt_min)
        h = np.minimum(h, dt_max)
        
        # store the adopted dt
        var.dt = h
        
        return var
            
    
class Output(object):
    
    def __init__(self):
        
        output_dir, out_name, plot_dir = vulcan_cfg.output_dir, vulcan_cfg.out_name, vulcan_cfg.plot_dir

        if not os.path.exists(output_dir): os.makedirs(output_dir)
        if not os.path.exists(plot_dir): os.makedirs(plot_dir)
        if vulcan_cfg.use_save_movie == True:
            if not os.path.exists(vulcan_cfg.movie_dir): os.makedirs(vulcan_cfg.movie_dir)
        
        if os.path.isfile(output_dir+out_name):
            # Fix Python 3.x and 2.x.
            # try: input = raw_input
            # except NameError: pass
            # input("  The output file: " + str(out_name) + " already exists.\n"
            #           "  Press enter to overwrite the existing file,\n"
            #           "  or Ctrl-Z and Return to leave and choose a different out_name in vulcan_cfg.")
            
            print ('Warning... the output file: ' + str(out_name) + ' already exists.\n')
        
    def print_prog(self, var):
        print ('time: ' +str("{:.2e}".format(var.t)) + '  and dt= ' + str(var.dt) )
        print ('longdy=' + str(var.longdy) + ' and longdy/dt= ' + str(var.longdydt) )
        
        # TEST
        #print ('dn/dt =')
        #print ('dn/dt =')
        
        print ('------------------------------------------------------' )

    def print_end_msg(self, var, para ): 
        print ("After ------- %s seconds -------" % ( time.time()- para.start_time ) + ' s CPU time') 
        print ('VULCAN has sucesfully run to steady-state with ' + str(para.count) + ' steps and ' + str("{:.2e}".format(var.t)) + ' s' )
        print ('long dy = ' + str(var.longdy) + ' and long dy/dt = ' + str(var.longdydt) )
        
        print ('total atom loss:')
        for atom in vulcan_cfg.atom_list: print (atom + ': ' + str(var.atom_loss[atom]) + ' ')
      
        print ('negative solution counter:')
        print (para.nega_count)
        print ('loss rejected counter:')
        print (para.loss_count)
        print ('delta rejected counter:')
        print (para.delta_count)
        if vulcan_cfg.use_shark == True: print ("It's a long journey to this shark planet. Don't stop bleeding.")
        print ('------ Live long and prosper \V/ ------') 
        
        
        
    def save_cfg(self, dname):
        output_dir, out_name = vulcan_cfg.output_dir, vulcan_cfg.out_name
        if not os.path.exists(output_dir):
            print ('The output directory assigned in vulcan_cfg.py does not exist.')
            print( 'Directory ' , output_dir,  " created.")
            os.mkdir(output_dir)

        # copy the vulcan_cfg.py file
        with open('vulcan_cfg.py' ,'r') as f:
            cfg_str = f.read()
        with open(dname + '/' + output_dir + "cfg_" + out_name[:-3] + "txt", 'w') as f: f.write(cfg_str)
    
    def save_out(self, var, atm, para, dname): 
        output_dir, out_name = vulcan_cfg.output_dir, vulcan_cfg.out_name
        output_file = dname + '/' + output_dir + out_name
        
        if not os.path.exists(output_dir):
            print ('The output directory assigned in vulcan_cfg.py does not exist.')
            print( 'Directory ' , output_dir,  " created.")
            os.mkdir(output_dir)
            
        # convert lists into numpy arrays
        for key in var.var_evol_save:
            as_nparray = np.array(getattr(var, key))
            setattr(var, key, as_nparray)
        
        # plotting
        if vulcan_cfg.use_plot_evo == True: 
            self.plot_evo(var, atm)
        if vulcan_cfg.use_plot_end == True:
            self.plot_end(var, atm, para)
        else: plt.close()
        
        # making the save dict
        var_save = {'species':species, 'nr':nr}
        
        for key in var.var_save:
            var_save[key] = getattr(var, key)
        if vulcan_cfg.save_evolution == True:
            # slicing time-sequential data to reduce ouput filesize
            fq = vulcan_cfg.save_evo_frq
            for key in var.var_evol_save:
                as_nparray = getattr(var, key)[::fq]
                setattr(var, key, as_nparray)
                var_save[key] = getattr(var, key)

        with open(output_file, 'wb') as outfile:
            if vulcan_cfg.output_humanread == True: # human-readable form, less efficient 
                outfile.write(str({'variable': var_save, 'atm': vars(atm), 'parameter': vars(para)}))
            else:
                # the protocol must be <= 2 for python 2.X
                pickle.dump( {'variable': var_save, 'atm': vars(atm), 'parameter': vars(para) }, outfile, protocol=2)
                # how to add  'config': vars(vulcan_cfg) ?
        
            
    def plot_update(self, var, atm, para):
        
        images = []
        colors = ['b','g','r','c','m','y','k','orange','pink', 'grey',\
        'darkred','darkblue','salmon','chocolate','mediumspringgreen','steelblue','plum','hotpink']
        
        tex_labels = {'H':'H','H2':'H$_2$','O':'O','OH':'OH','H2O':'H$_2$O','CH':'CH','C':'C','CH2':'CH$_2$','CH3':'CH$_3$','CH4':'CH$_4$','HCO':'HCO','H2CO':'H$_2$CO', 'C4H2':'C$_4$H$_2$',\
        'C2':'C$_2$','C2H2':'C$_2$H$_2$','C2H3':'C$_2$H$_3$','C2H':'C$_2$H','CO':'CO','CO2':'CO$_2$','He':'He','O2':'O$_2$','CH3OH':'CH$_3$OH','C2H4':'C$_2$H$_4$','C2H5':'C$_2$H$_5$','C2H6':'C$_2$H$_6$','CH3O': 'CH$_3$O'\
        ,'CH2OH':'CH$_2$OH', 'NH3':'NH$_3$'}
        
        plt.figure('live mixing ratios')
        plt.ion()
        color_index = 0
        for color_index, sp in enumerate(vulcan_cfg.plot_spec):
            if sp in tex_labels: sp_lab = tex_labels[sp]
            else: sp_lab = sp
            if color_index == len(para.tableau20): # when running out of colors
                para.tableau20.append(tuple(np.random.rand(3)))
            if vulcan_cfg.plot_height == False:
                line, = plt.plot(var.ymix[:,species.index(sp)], atm.pco/1.e6, color = para.tableau20[color_index], label=sp_lab)
                plt.gca().set_yscale('log')
                plt.gca().invert_yaxis()
                plt.ylabel("Pressure (bar)")
                plt.ylim((vulcan_cfg.P_b/1.E6,vulcan_cfg.P_t/1.E6))
            else: # plotting with height
                line, = plt.plot(var.ymix[:,species.index(sp)], atm.zmco/1.e5, color = para.tableau20[color_index], label=sp_lab)
                plt.ylim((atm.zco[0]/1e5,80))
                plt.ylabel("Height (km)")
                
            images.append((line,))
        
        plt.title(str(para.count)+' steps and ' + str("{:.2e}".format(var.t)) + ' s' )
        plt.gca().set_xscale('log')         
        plt.xlim(1.E-16, 1.)
        plt.legend(frameon=0, prop={'size':14}, loc=3)
        plt.xlabel("Mixing Ratios")
        plt.show(block=0)
        plt.pause(0.001)
        if vulcan_cfg.use_save_movie == True: 
            plt.savefig( vulcan_cfg.movie_dir+str(int(para.count/vulcan_cfg.save_movie_rate))+'.jpg', dpi=200)
        plt.clf()
    
    def plot_flux_update(self, var, atm, para):
        
        images = []
        plt.ion()
        
        # fig.add_subplot(121) fig.add_subplot(122)
        
        line1, = plt.plot(np.sum(var.dflux_u,axis=1), atm.pico/1.e6, label='up flux')
        line2, = plt.plot(np.sum(var.dflux_d,axis=1), atm.pico/1.e6, label='down flux', ls='--', lw=1.2)
        line3, = plt.plot(np.sum(var.sflux,axis=1), atm.pico/1.e6, label='stellar flux', ls=':', lw=1.5)
            
        images.append((line1,line2))        
        
        plt.title(str(para.count)+' steps and ' + str("{:.2e}".format(var.t)) + ' s' )
        plt.gca().set_xscale('log')       
        plt.gca().set_yscale('log') 
        plt.gca().invert_yaxis() 
        plt.xlim(xmin=1.E-8)
        plt.ylim((atm.pico[0]/1.e6,atm.pico[-1]/1.e6))
        plt.legend(frameon=0, prop={'size':14}, loc=3)
        plt.xlabel("Diffusive flux")
        plt.ylabel("Pressure (bar)")
        plt.show(block=0)
        plt.pause(0.1)
        if vulcan_cfg.use_flux_movie == True: plt.savefig( 'plot/movie/flux-'+str(para.count)+'.jpg')
        
        plt.clf()
        
    def plot_end(self, var, atm, para):
        
        plot_dir = vulcan_cfg.plot_dir
        colors = ['b','g','r','c','m','y','k','orange','pink', 'grey',\
        'darkred','darkblue','salmon','chocolate','mediumspringgreen','steelblue','plum','hotpink']
        
        plt.figure('live mixing ratios')
        color_index = 0
        for sp in vulcan_cfg.plot_spec:
            if vulcan_cfg.plot_height == False:
                line, = plt.plot(var.ymix[:,species.index(sp)], atm.pco/1.e6, color = colors[color_index], label=sp)
                plt.gca().set_yscale('log')
                plt.gca().invert_yaxis() 
                plt.ylabel("Pressure (bar)")
                plt.ylim((vulcan_cfg.P_b/1.E6,vulcan_cfg.P_t/1.E6))
            else: # plotting with height
                line, = plt.plot(var.ymix[:,species.index(sp)], atm.zmco/1.e5, color = colors[color_index], label=sp)
                plt.ylim((atm.zco[0]/1e5,80))
                plt.ylabel("Height (km)")
            color_index +=1
                  
        plt.title(str(para.count)+' steps and ' + str("{:.2e}".format(var.t)) + ' s' )
        plt.gca().set_xscale('log')
        plt.xlim(1.E-20, 1.)
        plt.legend(frameon=0, prop={'size':14}, loc=3)
        plt.xlabel("Mixing Ratios")
        plt.savefig(plot_dir + 'mix.png')       
        if vulcan_cfg.use_live_plot == True:
            # plotting in the same window of real-time plotting
            plt.draw()
        elif vulcan_cfg.use_PIL == True: # plotting in a new window with PIL package            
            plot = Image.open(plot_dir + 'mix.png')
            plot.show()
            plt.close()
            
    def plot_evo(self, var, atm, plot_j=-1, dn=1):
        
        plot_spec = vulcan_cfg.plot_spec
        plot_dir = vulcan_cfg.plot_dir
        plt.figure('evolution')
        
        ymix_time = np.array(var.y_time/atm.n_0[:,np.newaxis])
        
        for i,sp in enumerate(vulcan_cfg.plot_spec):
            plt.plot(var.t_time[::dn], ymix_time[::dn,plot_j,species.index(sp)],c = plt.cm.rainbow(float(i)/len(plot_spec)),label=sp)

        plt.gca().set_xscale('log')       
        plt.gca().set_yscale('log') 
        plt.xlabel('time')
        plt.ylabel('mixing ratios')
        plt.ylim((1.E-30,1.))
        plt.legend(frameon=0, prop={'size':14}, loc='best')
        plt.savefig(plot_dir + 'evo.png')
        if vulcan_cfg.use_PIL == True:
            plot = Image.open(plot_dir + 'evo.png')
            plot.show()
            plt.close()
        # else: plt.show(block = False)
    
    def plot_evo_inter(self, var, atm, plot_j=-1, dn=1):
        '''
        plot the evolution when the code is interrupted
        '''
        var.t_time = np.array(var.t_time)
        ymix_time = np.array(var.y_time/atm.n_0[:,np.newaxis])
        
        plot_spec = vulcan_cfg.plot_spec
        plot_dir = vulcan_cfg.plot_dir
        plt.figure('evolution')
    
        for i,sp in enumerate(vulcan_cfg.plot_spec):
            plt.plot(var.t_time[::dn], ymix_time[::dn,plot_j,species.index(sp)],c = plt.cm.rainbow(float(i)/len(plot_spec)),label=sp)

        plt.gca().set_xscale('log')       
        plt.gca().set_yscale('log') 
        plt.xlabel('time')
        plt.ylabel('mixing ratios')
        plt.ylim((1.E-30,1.))
        plt.legend(frameon=0, prop={'size':14}, loc='best')
        plt.savefig(plot_dir + 'evo.png')
        if vulcan_cfg.use_PIL == True:
            plot = Image.open(plot_dir + 'evo.png')
            plot.show()
            plt.close()
    
    def plot_TP(self, atm):
        plot_dir = vulcan_cfg.plot_dir
        plt.figure('TP')
        if vulcan_cfg.plot_height == False:
            plt.semilogy( atm.Tco, atm.pco/1.e6, c='black')
            plt.gca().invert_yaxis()
            plt.ylim((vulcan_cfg.P_b/1.E6,vulcan_cfg.P_t/1.E6))
            plt.ylabel("Pressure (bar)")
        else: # plotting with height
            
            line, = plt.plot(atm.Tco, atm.zmco/1.e5, c='black')            
            plt.ylabel("Height (km)")
        
        plt.xlabel("Temperature (K)")
        plot_name = plot_dir + 'TP.png'
        plt.savefig(plot_name)
        if vulcan_cfg.use_PIL == True:        
            plot = Image.open(plot_name)
            plot.show()
            # close the matplotlib window
            plt.close()
        else: plt.show(block = False)
            
class SemiEU(ODESolver):
    '''
    class inheritance from ODEsolver for semi-implicit Euler solver 
    '''
    def __init__(self):
        ODESolver.__init__(self)
           
    def solver(self, var, atm):
        """
        semi-implicit Euler solver (1st order) 
        """
        y, ymix, h, k = var.y, var.ymix, var.dt, var.k
        M, dzi, Kzz = atm.M, atm.dzi, atm.Kzz

        diffdf = self.diffdf
        jac_tot = self.jac_tot
        
        df = chemdf(y,M,k).flatten() + diffdf(var, atm).flatten()
        dfdy = jac_tot(var, atm)        
        aa = np.identity(ni*nz) - h*dfdy
        aa = scipy.linalg.solve(aa,df)
        aa = aa.reshape(y.shape)
        y = y + aa*h
        
        var.y = y
        var.ymix = var.y/np.vstack(np.sum(var.y,axis=1))
        
        return var

    def step_ok(self, var, para, loss_eps = vulcan_cfg.loss_eps):
        if np.all(var.y>=0) and np.amax( np.abs( np.fromiter(var.atom_loss.values(),float) - np.fromiter(var.atom_loss_prev.values(),float) ) )<loss_eps and para.delta<=rtol:
            return True
        else:
            return False    
    
    def one_step(self, var, atm, para):
        
        while True:
           var = self.solver(var, atm)
           
           # clipping small negative values and also calculating atomic loss (atom_loss)  
           var , para = self.clip(var, para, atm) 
            
           if self.step_ok(var, para): break
           elif self.step_reject(var, para): break # giving up and moving on
               
        return var, para                    
    
    def step_size(self, var, para):
        '''
        PID control required for all semi-Euler like methods
        '''
        dt_var_min, dt_var_max, dt_min, dt_max = vulcan_cfg.dt_var_min, vulcan_cfg.dt_var_max, vulcan_cfg.dt_min, vulcan_cfg.dt_max
        PItol = vulcan_cfg.PItol
        dy, dy_prev, h = var.dy, var.dy_prev, var.dt
        
        if dy == 0 or dy_prev == 0: 
            var.dt = np.minimum(h*2.,dt_max)
            return var
            
        if para.count > 0:
            
            h_factor = (dy_prev/dy)**0.075 * (PItol/dy)**0.175
            h_factor = np.maximum(h_factor, dt_var_min)    
            h_factor = np.minimum(h_factor, dt_var_max)
            h *= h_factor
            h = np.maximum(h, dt_min)
            h = np.minimum(h, dt_max)
        
        # store the adopted dt
        var.dt = h

        return var


class SparSemiEU(SemiEU):
    '''
    class inheritance from SemiEU.
    It is the same semi-implicit Euler solver except for utilizing sparse-matrix solvers
    '''
    def __init__(self):
        SemiEU.__init__(self)
    
    # override solver
    def solver(self, var, atm):
        """
        sparse-matrix semi-implicit Euler solver (1st order) 
        """
        y, ymix, h, k = var.y, var.ymix, var.dt, var.k
        M, dzi, Kzz = atm.M, atm.dzi, atm.Kzz

        diffdf = self.diffdf
        jac_tot = self.jac_tot
        
        df = chemdf(y,M,k).flatten() + diffdf(var, atm).flatten()
        dfdy = jac_tot(var, atm)        
                
        aa = sparse.csc_matrix( np.identity(ni*nz) - h*dfdy )
        aa = sparse.linalg.spsolve(aa,df)
        aa = aa.reshape(y.shape)
        y = y + aa*h
        
        var.y = y
        var.ymix = var.y/np.vstack(np.sum(var.y,axis=1))

        return var 
    

### back-up methods: extrapolation semi_implicit Euler ###