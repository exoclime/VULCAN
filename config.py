import os
VULCAN_DIR = os.path.dirname(os.path.abspath(__file__)) + "/"

class Config:
    """
    Configuration class for VULCAN model.
    """

    def __init__(self):
        self.atom_list               = ['H', 'O', 'C']
        self.network                 = VULCAN_DIR+'thermo/CHO_photo_network.txt'
        self.use_lowT_limit_rates    = False

        self.atm_base                = 'H2'
        self.rocky                   = True           # for the surface gravity
        self.nz                      = 61   # number of vertical layers
        self.P_b                     = 320906663.80470425  # pressure at the bottom (dyne/cm^2)
        self.P_t                     = 10.0  # pressure at the top (dyne/cm^2)
        self.atm_type                = 'file'
        self.atm_file                = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem/profile.dat'

        self.sflux_file              = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem/star.dat'
        self.top_BC_flux_file        = VULCAN_DIR+'atm/BC_top.txt' # the file for the top boundary conditions
        self.bot_BC_flux_file        = VULCAN_DIR+'atm/BC_bot.txt' # the file for the lower boundary conditions

        self.output_dir              = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem/'
        self.plot_dir                = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem/'
        self.movie_dir               = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem//frames/'
        self.out_name                = 'vulcan.pkl'

        # ====== Setting up the elemental abundance ======
        self.ini_mix = 'table'
        self.const_mix = { 'H2O':1.34606878e-02, 'CO2':5.60647590e-03, 'H2':2.74728526e-01, 'CH4':4.88058267e-04, 'CO':7.01674285e-01, 'N2':4.86994497e-04, 'NH3':3.47522764e-03, 'S2':2.34261857e-09, 'SO2':6.09690088e-08, 'H2S':7.96813790e-05 }
        self.vul_ini = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem/vmrs.dat'


        # ====== Setting up photochemistry ======
        self.use_ion         = False
        self.use_photo       = True
        self.r_star          = 0.2768848325427627     # stellar radius (R_sun)
        self.Rp              = 817738393.0      # Planetary radius (cm)
        self.orbit_radius    = 0.048833377211965893    # planet-star distance in A.U.
        self.gs              = 1275.58438      # surface gravity (cm/s^2)  (HD189:2140  HD209:936)
        self.sl_angle        = 0.955393232541696   # the zenith angle
        self.f_diurnal       = 0.25
        self.scat_sp         = ['H2', 'O2', 'CO2']
        self.T_cross_sp      = []

        self.edd             = 0.5 # the Eddington coefficient
        self.dbin1           = 0.1  # the uniform bin width < dbin_12trans (nm)
        self.dbin2           = 2.   # the uniform bin width > dbin_12trans (nm)
        self.dbin_12trans    = 240. # the wavelength switching from dbin1 to dbin2 (nm)

        # the frequency to update the actinic flux and optical depth
        self.ini_update_photo_frq    = 100
        self.final_update_photo_frq  = 5

        # ====== Mixing processes ======
        self.use_moldiff = True
        self.use_vz      = True
        self.vz_prof     = 'const'  # Options: 'const' or 'file'
        self.const_vz    = 0.0 # (cm/s)
        self.use_Kzz     = True
        self.Kzz_prof    = 'file' # Options: 'const','file'
        self.const_Kzz   = 100000.0 # Only reads when Kzz_prof = 'const'
        self.K_max       = 1e5        # for Kzz_prof = 'Pfunc'
        self.K_p_lev     = 0.1      # for Kzz_prof = 'Pfunc'
        self.update_frq  = 50    # frequency for updating dz and dzi due to change of mu

        # ====== Setting up the boundary conditions ======
        self.use_topflux     = False
        self.use_botflux     = False
        self.use_fix_sp_bot  = {  } # fixed mixing ratios at the lower boundary
        self.diff_esc        = [] # species for diffusion-limit escape at TOA
        self.max_flux        = 1e13  # upper limit for the diffusion-limit fluxes

        # ====== Reactions to be switched off  ======
        self.remove_list = [] # in pairs e.g. [1,2]

        # == Condensation ======
        self.use_condense        = False
        self.use_settling        = False
        self.start_conden_time   = 1e10
        self.condense_sp         = []
        self.non_gas_sp          = []
        self.fix_species         = []      # fixed the condensable species after condensation-evapoation EQ has reached
        self.fix_species_time    = 0  # after this time to fix the condensable species

        # ====== steady state check ======
        self.st_factor = 0.5
        self.conv_step = 100

        # ====== Setting up numerical parameters for the ODE solver ======
        self.ode_solver      = 'Ros2' # case sensitive
        self.trun_min        = 1e2
        self.runtime         = 1e22
        self.use_print_prog  = True
        self.use_print_delta = False
        self.print_prog_num  = 20  # print the progress every x steps
        self.dttry           = 1e-6
        self.dt_min          = 1e-8
        self.dt_max          = 1e18
        self.dt_var_max      = 2.
        self.dt_var_min      = 0.5

        self.count_min       = 120
        self.count_max       = int(3E4)
        self.atol            = 5.E-2 # Try decreasing this if the solutions are not stable
        self.mtol            = 1.E-22
        self.mtol_conv       = 1.E-20
        self.pos_cut         = 0
        self.nega_cut        = -1.
        self.loss_eps        = 1e-1
        self.yconv_cri       = 0.05  # for checking steady-state
        self.slope_cri       = 0.0001  # for checking steady-state
        self.yconv_min       = 0.5
        self.flux_cri        = 0.1
        self.flux_atol       = 1. # the tol for actinc flux (# photons cm-2 s-1 nm-1)
        self.conver_ignore   = [] # added 2023. to get rid off non-convergent species, e.g. HC3N without sinks

        # ====== Setting up numerical parameters for Ros2 ODE solver ======
        self.rtol             = 0.7 # relative tolerence for adjusting the stepsize
        self.post_conden_rtol = 0.1 # switched to this value after fix_species_time

        # ====== Setting up for output and plotting ======
        self.plot_TP         = False
        self.use_live_plot   = False
        self.use_live_flux   = False
        self.use_plot_end    = False
        self.use_plot_evo    = False
        self.use_save_movie  = False
        self.use_flux_movie  = False
        self.plot_height     = False
        self.use_PIL         = False
        self.live_plot_frq   = 50
        self.save_movie_rate = 50
        self.y_time_freq     = 1  #  storing data for every 'y_time_freq' step
        self.plot_spec       = ['H2',  'H', 'H2O', 'CH4', 'CO', 'CO2', 'C2H2']

        # output:
        self.output_humanread = False
        self.use_shark        = False
        self.save_evolution   = False   # save the evolution of chemistry (y_time and t_time) for every save_evo_frq step
        self.save_evo_frq     = 10
