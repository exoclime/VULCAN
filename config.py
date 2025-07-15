from paths import VULCAN_DIR
import os

class Config:
    """
    Configuration class for VULCAN model.
    """

    def __init__(self):
        self.atom_list               = ['H', 'O', 'C', 'N']
        self.network                 = VULCAN_DIR+'thermo/NCHO_photo_network.txt'
        self.use_lowT_limit_rates    = False

        self.atm_base                = 'H2'
        self.rocky                   = False           # for the surface gravity
        self.nz                      = 100   # number of vertical layers
        self.P_b                     = 1e9 # pressure at the bottom (dyne/cm^2)
        self.P_t                     = 1e-1  # pressure at the top (dyne/cm^2)

        # Set T(p) from file
        self.atm_type                = 'file'
        self.atm_file                = VULCAN_DIR+'atm/atm_HD189_Kzz.txt'

        # Analytical T(p) parameters
        #  T_int, T_irr, ka_L, ka_S, beta_S, beta_L (details see Heng et al. 2014)
        self.para_anaTP              = [120., 1500., 0.1, 0.02, 1., 1.]

        # Isothermal T(p)
        self.Tiso                    = 1000.0
            
        self.sflux_file              = VULCAN_DIR+'atm/stellar_flux/sflux-HD189_Moses11.txt'
        self.top_BC_flux_file        = VULCAN_DIR+'atm/BC_top.txt' # the file for the top boundary conditions
        self.bot_BC_flux_file        = VULCAN_DIR+'atm/BC_bot.txt' # the file for the lower boundary conditions

        self.output_dir              = VULCAN_DIR + "output/"
        self.plot_dir                = self.output_dir + "plot/"
        self.out_name                = 'example.pkl'


        # ====== Setting up the elemental abundance ======
        self.ini_mix = 'const_mix'
        self.const_mix = { 'H2O':3e-2, 'CO2':5e-3, 'H2':0.9, 'CH4':3e-2, 'CO':1e-3, 'N2':4e-3, 'NH3':3e-2}
        self.vul_ini = '_unset_'


        # ====== Setting up photochemistry ======
        self.use_ion         = False
        self.use_photo       = True
        self.r_star          = 0.805     # stellar radius (R_sun)
        self.Rp              = 1.138*7.1492E9      # Planetary radius (cm)
        self.orbit_radius    = 0.03142    # planet-star distance in A.U.
        self.gs              = 2140.       # surface gravity (cm/s^2)  (HD189:2140  HD209:936)
        self.sl_angle        = 58 /180.*3.14159  # the zenith angle
        self.f_diurnal       = 1.0
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
        self.Kzz_prof    = 'Pfunc' # Options: 'const','file'
        self.const_Kzz   = 1e10 # Only reads when Kzz_prof = 'const'
        self.K_max       = 1e5        # for Kzz_prof = 'Pfunc'
        self.K_p_lev     = 0.1      # for Kzz_prof = 'Pfunc'
        self.update_frq  = 50    # frequency for updating dz and dzi due to change of mu

        # ====== Setting up the boundary conditions ======
        self.use_topflux     = False
        self.use_botflux     = False
        self.use_fix_sp_bot  = {  } # fixed mixing ratios at the lower boundary
        self.diff_esc        = ['H'] # species for diffusion-limit escape at TOA
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

        # ====== coupling to AGNI ======
        self.agni_call_frq  = 10
        self.spectral_file  = "res/spectral_files/Dayspring/48/Dayspring.sf"
        self.use_rayleigh   = True
        self.surf_albedo    = 0.0
        self.Tsurf_guess    = 2000.0
        self.agni_atol      = 1e-1
        self.agni_rtol      = 1e-2
        self.agni_nlev      = 35
        self.solve_rce      = True


        # ====== steady state check ======
        self.st_factor      = 0.5
        self.conv_step      = 100

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
        self.plot_dpi        = 130
        self.plot_TP         = True
        self.use_live_plot   = True
        self.use_live_flux   = False
        self.use_plot_end    = False
        self.use_plot_evo    = False
        self.use_flux_movie  = False
        self.plot_height     = False
        self.live_plot_frq   = 50
        self.y_time_freq     = 1  #  storing data for every 'y_time_freq' step
        self.plot_spec       = ['H2',  'H', 'H2O', 'CH4', 'CO', 'CO2', 'C2H2']

        # output:
        self.output_humanread = False
        self.save_evolution   = False   # save the evolution of chemistry (y_time and t_time) for every save_evo_frq step
        self.save_evo_frq     = 10

    def write_file(self):
        fpath = os.path.join(self.output_dir,"vulcan_cfg.txt")
        with open(fpath,'w') as hdl:
            hdl.write("# Variables set in vulcan_cfg object \n")
            for attr, value in self.__dict__.items():
                hdl.write("%-22s  %-s \n"%(attr,str(value)))
