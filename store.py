import numpy as np
import scipy
import vulcan_cfg
from vulcan_cfg import nz
from chem_funs import ni, nr, spec_list  # number of species and reactions in the network

class Variables(object):
    """
    store the essential variables for calculation  
    """
    def __init__(self): # self means the created object
        #self.nz = vulcan_cfg.nz
        self.k = {}  # k[] follows the same shape as Tco in read_rate() in op.py
        self.y = np.zeros((nz, ni))
        self.y_prev = np.zeros((nz, ni))
        self.ymix = np.zeros((nz, ni))
        self.y_ini = np.zeros((nz, ni))
        self.y_conden = np.zeros((nz, ni))  # to store the species removed by condensation
        
        #self.yconv = 1.
        #self.yconv_prev = 1.
        
        # y_loss has been changed to atom_loss
        #self.y_loss = 0
        
        self.t = 0
        self.dt = vulcan_cfg.dttry
        self.dy = 1.
        self.dy_prev = 1.
        self.dydt = 1.
        self.longdy = 1.
        self.longdydt = 1.
        
        self.dy_time = []
        self.dydt_time = []
        self.atim_loss_time = []
        self.ymix_time = []
        self.y_time = []
        self.t_time = []
        self.dt_time = []
        self.atom_loss_time = []
        
        # self.atom_tot = [0] * na
        self.atom_ini = {}
        self.atom_sum = {}
        self.atom_loss = {}
        self.atom_loss_prev = {}
        #  self.mean_mass = np.zeros(nz) 
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
        
        self.def_bin_min = 2. #20.
        self.def_bin_max = 800.1
        
        
        ### ### ### ### ### ### ### ### ### ### ###
        # List the names variables defined in op here!
        ### ### ### ### ### ### ### ### ### ### ###
        
        # User define what to save!
        self.var_save = ['k','y','ymix','y_ini','y_conden','t','dt','longdy','longdydt',\
        'atom_ini','atom_sum','atom_loss','atom_conden','aflux_change','Rf'] 
        if vulcan_cfg.use_photo == True: self.var_save.extend(['nbin','bins','dbin','tau','sflux','aflux','cross','cross_scat','cross_J', 'J_sp','wavelen','n_branch','br_ratio'])
        self.var_evol_save = ['y_time','t_time'] #'atom_loss_time'
        self.conden_re_list = []        
        
        # photo data initiated in read_cross in op.py
        # self.cross = {}
        # self.cross_scat = {}
        
        # self.dbin = vulcan_cfg.dbin #0.2 # d-lambda
#         self.bins = np.arange(20.,400.1, self.dbin)
#         self.bins[-1] -= 1.E-10   # for solaving the floating pb in binary when doing read_cross
#         self.nbin = len(self.bins)

        # # the direct beam (staggered)
#         self.sflux = np.empty( (nz+1, self.nbin) )
#         # the diffusive flux (staggered)
#         self.dflux_u, self.dflux_d = np.zeros( (nz+1, self.nbin) ), np.zeros( (nz+1, self.nbin) )
#         # the total actinic flux (non-staggered)
#         self.aflux = np.empty( (nz, self.nbin) )
#         # staggered
#         self.tau = np.zeros( (nz+1, self.nbin) )
#
#         self.sflux_top = np.empty(self.nbin)
        # the total actinic flux (center, non-staggered) from the previous calculation 
        # self.prev_aflux = np.empty( (nz, self.nbin) )
        

class AtmData(object):
    """
    store the data of atmospheric structure  
    """
    def __init__(self):
        self.pco = np.logspace(np.log10(vulcan_cfg.P_b),np.log10(vulcan_cfg.P_t),nz)
        self.pico = np.empty(nz+1) # need the bottom and top layers to calculate hydrostatic relation of dz
        self.dz = np.empty(nz)
        self.dzi = np.empty(nz-1)
        self.zco = np.empty(nz+1) # not used in calculation
        self.zmco = np.empty(nz)
        self.Tco = np.empty(nz)
        self.Kzz = np.zeros(nz-1)
        self.vz = np.zeros(nz-1)
        self.M = np.empty(nz)
        self.n_0 = np.empty(nz)
        self.Hp = np.empty(nz)
        self.mu = np.empty(nz)
        self.ms = np.empty(ni) # molecular weight for every species
        self.Dzz = np.zeros((nz-1,ni))
        self.vs = np.zeros((nz-1,ni)) # the settling velocity
        self.alpha = -0.25*np.ones(ni) # thermal diffusion factor = -0.25 for every species except for H and H2 (defined in mol_diff() in build_atm.py)
        
        self.top_flux = np.zeros(ni)
        self.bot_flux = np.zeros(ni)
        self.bot_vdep = np.zeros(ni)
        self.bot_fix_sp = np.zeros(ni)
        
        self.sat_p = {}
        self.sat_mix = {}
        
        # condensation excluding non-gaseous species
        if vulcan_cfg.use_condense == True:
            self.exc_conden = [_ for _ in range(ni) if spec_list[_] not in vulcan_cfg.non_gas_sp]
        # TEST condensation excluding non-gaseous species
        
        # # put in build atm ?
        # if vulcan_cfg.atm_base == 'H2'
        #     self.Dzz = dict([( (sp, ) , np.zeros(nz)) for sp in var.photo_sp for bn in range(n_branch[sp]+1) ])
        #

class Parameters(object):
    """
    store the overall parameters for numerical method and counters  
    """
    def __init__(self):
        self.nega_y = 0
        self.small_y = 0
        self.delta = 0
        self.count = 0
        self.nega_count =0
        self.loss_count = 0
        self.delta_count= 0
        self.end_case = 0
        self.solver_str = '' # for assigning the name of solver
        self.switch_final_photo_frq = False