# ============================================================================= 
# Configuration file of VULCAN:  
# ============================================================================= 

# ====== Setting up the elements included in the network ======
atom_list = ['H', 'O', 'C', 'He', 'N']
# ====== Setting up paths and filenames for the input and output files  ======
# input:
network = 'thermo/NCHO_earth_photo_network.txt'
gibbs_text = 'thermo/gibbs_text.txt' # (all the nasa9 files must be placed in the folder: thermo/NASA9/)
cross_folder = 'thermo/photo_cross/'
com_file = 'thermo/all_compose.txt'
atm_file = 'atm/atm_Earth_Jan_Kzz.txt'
sflux_file = 'atm/stellar_flux/VPL_solar_100nmcut.txt' # This is the flux density at the stellar surface
top_BC_flux_file = 'atm/BC_top.txt'
bot_BC_flux_file = 'atm/BC_bot_Earth.txt'
vul_ini = 'output/No-condense-Earth.vul'
# output:
output_dir = 'output/'
plot_dir = 'plot/'
movie_dir = 'plot/movie/earth/'
out_name =  'rtol005-2nd-Earth.vul'

# ====== Setting up the elemental abundance ======
use_solar = False # True: using the solar abundance from Table 10. K.Lodders 2009; False: using the customized elemental abundance. 
# customized elemental abundance (only reads when use_solar = False)
O_H = 6.0618E-4 *(0.793)  
C_H = 2.7761E-4  
N_H = 8.1853E-5
He_H = 0.09691
ini_mix = 'vulcan_ini' # Options: 'EQ', 'const_mix', 'vulcan_ini' (for 'vulcan_ini, the T-P grids have to be exactly the same)
# Initialsing uniform (constant with pressure) mixing ratios (only reads when ini_mix = const_mix)
const_mix = {'CH4':2.7761E-4*2, 'O2':4.807e-4, 'He':0.09691, 'N2':8.1853E-5, 'H2':1. -2.7761E-4*2*4/2} 

# ====== Setting up photochemistry ======
use_photo = True
# astronomy input
r_star = 1 # stellar radius in solar radius
orbit_radius = 1 # planet-star distance in A.U.
sl_angle = 48 /180.*3.14159 # the zenith angle of the star in degree
# radiation parameters 
excit_sp = ['O_1', 'CH2_1'] # N_D not included due to lack of NASA9 Gibbs energy
scat_sp = ['N2', 'O2'] # the bulk compositions that contribute to Rayleigh scattering
edd = 0.669 #(cos(48 deg) ) # the Eddington coefficient 
dbin = 0.2  # the uniform bin width
# frequency to update the flux and optical depth
ini_update_photo_frq = 100
final_update_photo_frq = 10

# ====== Setting up parameters for the atmosphere ======
atm_base = 'N2' #Options: 'H2', 'N2', 'O2', 'CO2 -- the bulk gas of the atmosphere: affects molecular diffsion
nz = 100   # number of vertical layers
P_b = 1.E6 # pressure at the bottom (dyne/cm^2)
P_t = 1.e0 # pressure at the top (dyne/cm^2)
use_Kzz = True
use_moldiff = True
use_vz = False
atm_type = 'file' # Options: 'isothermal', 'analytical', or 'file'
Kzz_prof = 'file' # Options: 'const' or 'file'
vz_prof = 'const' # Options: 'const' or 'file'
g = 1000.         # gravity (cm/s^2)  (HD189:2140  HD209:936)
Tiso = 3000. # only reads when atm_type = 'isothermal'
# setting the parameters for the analytical T-P from (126)in Heng et al. 2014. Only reads when atm_type = 'analytical' 
# T_int, T_irr, ka_L, ka_S, beta_S, beta_L
para_warm = [120., 1500., 0.1, 0.02, 1., 1.]
para_anaTP = para_warm
const_Kzz = 1.E6 # (cm^2/s) Only reads when use_Kzz = True and Kzz_prof = 'const'
const_vz = 0 # (cm/s) Only reads when use_vz = True and vz_prof = 'const'

# frequency for updating dz and dzi due to change of mu
update_frq = 100 

# ====== Setting up the boundary conditions ======
# Boundary Conditions:
use_topflux = False
use_botflux = True
#use_fix_all_bot = True
use_fix_sp_bot = {'H2O':0.01}

# ====== Reactions to be switched off  ======
remove_list = [] # in pairs e.g. [1,2]

# == Condensation (Ongoing testing!)  ======
use_condense = True
use_settling = False
start_conden_time = 1e6
condesne_sp = ["H2O"]    # , 'NH3'
non_gas_sp = ["H2O_l_s"]

# ====== steady state check ======
st_factor = 0.95

# ====== Setting up numerical parameters for the ODE solver ====== 
ode_solver = 'Ros2' # case sensitive
use_print_prog = True
print_prog_num = 500  # every x steps to print progress
dttry = 1.E-10
trun_min = 1e2
runtime = 1.E22
dt_min = 1.E-14
dt_max = runtime*1e-5
dt_var_max = 2.
dt_var_min = 0.5
count_min = 120
count_max = int(2E5)
atol = 1.E-2 # Try decreasing this if the solutions are not stable
mtol = 1.E-22
mtol_conv = 1.E-16
pos_cut = 0
nega_cut = -1.
loss_eps = 1e-1
yconv_cri = 0.01 # for checking steady-state
slope_cri = 1.e-4
yconv_min = 0.1
flux_cri = 5.e-2 
flux_atol = 1. # the tol for actinc flux (# photons cm-2 s-1 nm-1)
# ====== Setting up numerical parameters for Ros2 ODE solver ====== 
rtol = 0.05
# ====== Setting up numerical parameters for SemiEu/SparSemiEU ODE solver (Not used) ====== 
PItol = 0.1

# ====== Setting up for output and plotting ======
# plotting:
plot_TP = True
use_live_plot = True
use_live_flux = False
use_plot_end = True
use_plot_evo = True
use_save_movie = False
use_flux_movie = False
plot_height = False
use_PIL = True 
live_plot_frq = 10
save_movie_rate = live_plot_frq
y_time_freq = 1  #  storing data for every 'y_time_freq' step
plot_spec = ['H2O', 'H2O_l_s', 'CO2', 'CH4', 'NO', 'NO2', 'O2', 'O3','N2O', 'NH3']  
# output:
output_humanread = False
save_evolution = False
save_evo_frq = 10
