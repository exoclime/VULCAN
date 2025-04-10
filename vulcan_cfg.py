# VULCAN CONFIGURATION FILE
# CREATED AUTOMATICALLY BY PROTEUS

atom_list               = ['H', 'O', 'C']
network                 = 'thermo/CHO_photo_network.txt'
use_lowT_limit_rates    = False
gibbs_text              = 'thermo/gibbs_text.txt' # (all the nasa9 files must be placed in the folder: thermo/NASA9/)
cross_folder            = 'thermo/photo_cross/'
com_file                = 'thermo/all_compose.txt'

atm_base                = 'H2'
rocky                   = True           # for the surface gravity
nz                      = 61   # number of vertical layers
P_b                     = 320906663.80470425  # pressure at the bottom (dyne/cm^2)
P_t                     = 10.0  # pressure at the top (dyne/cm^2)
atm_type                = 'file'
atm_file                = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem/profile.dat'

sflux_file              = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem/star.dat'
top_BC_flux_file        = 'atm/BC_top.txt' # the file for the top boundary conditions
bot_BC_flux_file        = 'atm/BC_bot.txt' # the file for the lower boundary conditions

output_dir              = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem/'
plot_dir                = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem/'
movie_dir               = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem//frames/'
out_name                = 'vulcan.pkl'

# ====== Setting up the elemental abundance ======
ini_mix = 'table'
const_mix = { 'H2O':1.34606878e-02, 'CO2':5.60647590e-03, 'H2':2.74728526e-01, 'CH4':4.88058267e-04, 'CO':7.01674285e-01, 'N2':4.86994497e-04, 'NH3':3.47522764e-03, 'S2':2.34261857e-09, 'SO2':6.09690088e-08, 'H2S':7.96813790e-05 }
vul_ini = '/Users/nichollsh/Projects/PROTEUS/output/physical_agni/offchem/vmrs.dat'


# ====== Setting up photochemistry ======
use_ion         = False
use_photo       = True
r_star          = 0.2768848325427627     # stellar radius (R_sun)
Rp              = 817738393.0      # Planetary radius (cm)
orbit_radius    = 0.048833377211965893    # planet-star distance in A.U.
gs              = 1275.58438      # surface gravity (cm/s^2)  (HD189:2140  HD209:936)
sl_angle        = 0.955393232541696   # the zenith angle
f_diurnal       = 0.25
scat_sp         = ['H2', 'O2', 'CO2']
T_cross_sp      = []

edd             = 0.5 # the Eddington coefficient
dbin1           = 0.1  # the uniform bin width < dbin_12trans (nm)
dbin2           = 2.   # the uniform bin width > dbin_12trans (nm)
dbin_12trans    = 240. # the wavelength switching from dbin1 to dbin2 (nm)

# the frequency to update the actinic flux and optical depth
ini_update_photo_frq    = 100
final_update_photo_frq  = 5

# ====== Mixing processes ======
use_moldiff = True

use_vz      = True
vz_prof     = 'const'  # Options: 'const' or 'file'
const_vz    = 0.0 # (cm/s)

use_Kzz     = True
Kzz_prof    = 'file' # Options: 'const','file'
const_Kzz   = 100000.0 # Only reads when Kzz_prof = 'const'
K_max       = 1e5        # for Kzz_prof = 'Pfunc'
K_p_lev     = 0.1      # for Kzz_prof = 'Pfunc'

update_frq  = 50    # frequency for updating dz and dzi due to change of mu

# ====== Setting up the boundary conditions ======
use_topflux     = False
use_botflux     = False
use_fix_sp_bot  = {  } # fixed mixing ratios at the lower boundary
diff_esc        = [] # species for diffusion-limit escape at TOA
max_flux        = 1e13  # upper limit for the diffusion-limit fluxes

# ====== Reactions to be switched off  ======
remove_list = [] # in pairs e.g. [1,2]

# == Condensation ======
use_condense        = False
use_settling        = False
start_conden_time   = 1e10
condense_sp         = []
non_gas_sp          = []
fix_species         = []      # fixed the condensable species after condensation-evapoation EQ has reached
fix_species_time    = 0  # after this time to fix the condensable species

# ====== steady state check ======
st_factor = 0.5
conv_step = 100

# ====== Setting up numerical parameters for the ODE solver ======
ode_solver      = 'Ros2' # case sensitive
trun_min        = 1e2
runtime         = 1.E22
use_print_prog  = True
use_print_delta = False
print_prog_num  = 20  # print the progress every x steps
dttry           = 1.E-6
dt_min          = 1.E-8
dt_max          = runtime*1e-4
dt_var_max      = 2.
dt_var_min      = 0.5

count_min       = 120
count_max       = int(3E4)
atol            = 5.E-2 # Try decreasing this if the solutions are not stable
mtol            = 1.E-22
mtol_conv       = 1.E-20
pos_cut         = 0
nega_cut        = -1.
loss_eps        = 1e-1
yconv_cri       = 0.05  # for checking steady-state
slope_cri       = 0.0001  # for checking steady-state
yconv_min       = 0.5
flux_cri        = 0.1
flux_atol       = 1. # the tol for actinc flux (# photons cm-2 s-1 nm-1)
conver_ignore   = [] # added 2023. to get rid off non-convergent species, e.g. HC3N without sinks

# ====== Setting up numerical parameters for Ros2 ODE solver ======
rtol             = 0.7 # relative tolerence for adjusting the stepsize
post_conden_rtol = 0.1 # switched to this value after fix_species_time

# ====== Setting up for output and plotting ======
plot_TP         = False
use_live_plot   = False
use_live_flux   = False
use_plot_end    = False
use_plot_evo    = False
use_save_movie  = False
use_flux_movie  = False
plot_height     = False
use_PIL         = False
live_plot_frq   = 50
save_movie_rate = live_plot_frq
y_time_freq     = 1  #  storing data for every 'y_time_freq' step
plot_spec       = ['H2',  'H', 'H2O', 'CH4', 'CO', 'CO2', 'C2H2']
# output:
output_humanread = False
use_shark        = False
save_evolution   = False   # save the evolution of chemistry (y_time and t_time) for every save_evo_frq step
save_evo_frq     = 10
