# ==============================================================================
# Module contains classes that store all the variables used in VULCAN. 
# Copyright (C) 2016 Shang-Min Tsai (Shami)                    
# ==============================================================================
# There are three classes: Variables, AtmData, and Parameters, which store the main
# physical variables, atmospheric variables, and numerical parameters, respectively.
# ==============================================================================

import numpy as np
import scipy
import vulcan_cfg
from vulcan_cfg import nz
from chem_funs import ni, nr, spec_list  # number of species and reactions in the network

class Variables(object):
    """
    store the essential variables for calculation  
    """
    def __init__(self): # self means the created object instance
        self.k = {}  # rate coefficients: k[1] is the rate constant of R1 reaction at every level (same shape as Tco and pco)
        self.y = np.zeros((nz, ni)) # current number density in the shape of (number of vertical levels, number of species)
        self.y_prev = np.zeros((nz, ni)) # number density at the previous step
        self.ymix = np.zeros((nz, ni)) # current mixing ratios
        self.y_ini = np.zeros((nz, ni)) # the initial number density 
        self.t = 0 # integration time
        self.dt = vulcan_cfg.dttry # time step
        self.dy = 1. # the max change of y from the previous step to current step 
        self.dy_prev = 1. # the max change of y from i-2 step to i-1 step 
        self.dydt = 1. # the max change of dy/dt
        self.longdy = 1. # the max change of y over a period of time 
        # for checking the steady state (dn in (10) in Tsai et al 2017)
        self.longdydt = 1. # the max change of dy/dt over a period of time 
        # (dn/dt in (11) in Tsai et al 2017)
        #self.slope_min = 0 # 
        
        self.dy_time = [] # storing dy at each step 
        self.dydt_time = []
        self.atim_loss_time = []
        self.ymix_time = [] # storing the mixing ratio at each step 
        self.y_time = [] # storing the number density at each step 
        self.t_time = [] # storing the time  at each step 
        self.dt_time = [] # storing the time step  at each step 
        self.atom_loss_time = [] # storing the loss of atoms at each step
        
        self.atom_ini = {}
        self.atom_sum = {}
        self.atom_loss = {}
        self.atom_loss_prev = {}
        self.atom_conden = {}  
        
        self.Rf = {}  # reaction equation
        self.Rindx = {}
        self.a, self.n, self.E, self.a_inf, self.n_inf, self.E_inf,= [{} for i in range(6)]
        self.k, self.k_fun, self.k_inf = [{} for i in range(3)] 
        self.photo_sp = set()  
        self.pho_rate_index, self.n_branch, self.wavelen, self.br_ratio = {}, {}, {}, {}
        
        self.kinf_fun = {}
        self.k_fun_new = {}
        
        # the max change of the actinic flux (for convergence)
        # if photochemistry is off, the value remaines 0 for checking convergence
        self.aflux_change = 0.
        
        self.def_bin_min = 2. 
        self.def_bin_max = 800.1

        # Define what variables to save in the output file!
        self.var_save = ['k','y','ymix','y_ini','t','dt','longdy','longdydt',\
        'atom_ini','atom_sum','atom_loss','atom_conden','aflux_change','Rf'] 
        if vulcan_cfg.use_photo == True: self.var_save.extend(['nbin','bins','dbin','tau','sflux','aflux','cross','cross_scat','cross_J', 'J_sp','wavelen','n_branch','br_ratio'])
        self.var_evol_save = ['y_time','t_time']
        self.conden_re_list = []        
        
        ### ### ### ### ### ### ### ### ### ### ###
        # List the names variables defined in op here!
        ### ### ### ### ### ### ### ### ### ### ###
        
        # Old stuff
        # self.yconv = 1.
        # self.yconv_prev = 1.
        # self.y_conden = np.zeros((nz, ni))  # to store the species removed by condensation      
        # photo data initiated in read_cross in op.py

class AtmData(object):
    """
    store the data of atmospheric structure  
    """
    def __init__(self):
        self.pco = np.logspace(np.log10(vulcan_cfg.P_b),np.log10(vulcan_cfg.P_t),nz) # pressure grids
        self.pico = np.empty(nz+1) # pressure grids at the interface
        self.dz = np.empty(nz) # height grids
        self.dzi = np.empty(nz-1) # height grids at the interface
        self.zco = np.empty(nz+1) # not used in calculation
        self.zmco = np.empty(nz)  
        self.Tco = np.empty(nz)  # temperature grids
        self.Kzz = np.zeros(nz-1) 
        self.vz = np.zeros(nz-1) # vertical velocity
        self.M = np.empty(nz) # the number density of the third body
        self.n_0 = np.empty(nz) # the total number density
        self.Hp = np.empty(nz) # the scale height
        self.mu = np.empty(nz) # mean molecular weight
        self.ms = np.empty(ni) # molecular weight for every species
        self.Dzz = np.zeros((nz-1,ni)) # molecular diffusion (nz,ni)
        self.vs = np.zeros((nz-1,ni)) # the settling velocity
        self.alpha = -0.25*np.ones(ni) # thermal diffusion factor = -0.25 for every species except for H and H2 (defined in mol_diff() in build_atm.py)
        
        self.top_flux = np.zeros(ni) # the assigned flux at the top boundary
        self.bot_flux = np.zeros(ni) # the assigned flux at the bottom boundary
        self.bot_vdep = np.zeros(ni) # the deposition velocity at the bottom boundary (surface sink)
        self.bot_fix_sp = np.zeros(ni) # the fixed mixing ratios for certain species at the bottom boundary
        
        self.sat_p = {}
        self.sat_mix = {}
        
        # condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            self.exc_conden = [_ for _ in range(ni) if spec_list[_] not in vulcan_cfg.non_gas_sp]
        # TEST condensation excluding non-gaseous species
        

class Parameters(object):
    """
    store the overall parameters for numerical method and counters  
    """
    def __init__(self):
        self.nega_y = 0
        self.small_y = 0
        self.delta = 0
        self.count = 0  # number of steps taken
        self.nega_count =0
        self.loss_count = 0
        self.delta_count= 0
        self.end_case = 0
        self.solver_str = '' # for assigning the name of solver
        self.switch_final_photo_frq = False
        
        # These are the "Tableau 20" colors as RGB.    
        self.tableau20 = [(31, 119, 180),(255, 127, 14),(44, 160, 44),(214, 39, 40),(148, 103, 189),(140, 86, 75), (227, 119, 194),(127, 127, 127),(188, 189, 34),(23, 190, 207),\
        (174, 199, 232),(255, 187, 120),(152, 223, 138),(255, 152, 150),(197, 176, 213),(196, 156, 148),(247, 182, 210),(199, 199, 199),(219, 219, 141),(158, 218, 229)]     
  
        # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
        for i in range(len(self.tableau20)):    
            r, g, b = self.tableau20[i]    
            self.tableau20[i] = (r / 255., g / 255., b / 255.)