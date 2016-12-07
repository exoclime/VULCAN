import numpy as np
import scipy
#from vulcan_cfg import *
import vulcan_cfg
#from phy_const import *
#import HOC, HOC_jac
from vulcan_cfg import nz
from vulcan_cfg import na # the number of atoms
from chem_funs import ni, nr  # number of species and reactions in the network


class Variables(object):
    """
    store the essential variables for calculation  
    """
    def __init__(self): # self means the created object
        #self.nz = vulcan_cfg.nz
        self.k = {}
        self.y = np.zeros((nz, ni))
        self.y_prev = np.zeros((nz, ni))
        self.ymix = np.zeros((nz, ni))
        self.ysum = 0
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
        self.mean_mass = np.zeros(nz) 
        
        
        self.Rf = {}  # reaction equation
        self.Rindx = {}
        self.a, self.n, self.E, self.a_inf, self.n_inf, self.E_inf,= [{} for i in range(6)]
        self.k, self.k_fun, self.k_inf = [{} for i in range(3)]   
        
        
        ### ???
        self.kinf_fun = {}
        self.k_fun_new = {}
    

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
        self.Tco = np.empty(nz)
        self.Kzz = np.empty(nz-1)
        self.M = np.empty(nz)
        self.n_0 = np.empty(nz)
        self.Hp = np.empty(nz)
        self.mu = np.empty(nz)

                


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
        
