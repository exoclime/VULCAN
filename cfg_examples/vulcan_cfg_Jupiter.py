# ============================================================================= 
# Configuration file of VULCAN:  
# ============================================================================= 

# ====== Setting up the elements included in the network ======
atom_list = ['H', 'O', 'C', 'N']
use_lowT_limit_rates = True

# ====== Setting up paths and filenames for the input and output files  ======
# input:
network = 'thermo/NCHO_photo_network_lowT_Jupiter.txt' # NCHO_photo_network_lowT_Jupiter.txt
gibbs_text = 'thermo/gibbs_text.txt' # (all the nasa9 files must be placed in the folder: thermo/NASA9/)
cross_folder = 'thermo/photo_cross/'
com_file = 'thermo/all_compose.txt'
atm_file = 'atm/Jupiter_deep_top.txt' # TP and Kzz (optional) file
sflux_file = 'atm/stellar_flux/Gueymard_solar.txt' # This is the flux density at the stellar surface
top_BC_flux_file = 'atm/BC_top_Jupiter.txt' # the file for the top boundary conditions
bot_BC_flux_file = 'atm/BC_bot.txt' # the file for the lower boundary conditions
vul_ini = 'output/' # the file to initialize the abundances for ini_mix = 'vulcan_ini'
# output:
output_dir = 'output/'
plot_dir = 'plot/'
movie_dir = 'plot/movie/'
out_name =  'Jupiter.vul' # output file name

# ====== Setting up the elemental abundance ======
use_solar = False # True: using the solar abundance from Table 10. K.Lodders 2009; False: using the customized elemental abundance. 
# customized elemental abundance (only read when use_solar = False)
O_H =  6.0618E-4  *0.5
C_H =  1.19e-3 # X4.28 solar (Atreya 2020 review)
N_H =  8.1853E-5 *2.78 # Li 2017 Juno 
# S_H = 4.45E-5   # (Atreya 2020 review)
He_H = 0.0785   # about 0.8 solar (Atreya 2020 review)
ini_mix = 'EQ' # Options: 'EQ', 'const_mix', 'vulcan_ini', 'table' (for 'vulcan_ini, the T-P grids have to be exactly the same)
fastchem_met_scale = 1. # scaling factor for other elements in fastchem (e.g., if fastchem_met_scale = 0.1, other elements such as Si and Mg will take 0.1 solar values)


# Initialsing uniform (constant with pressure) mixing ratios (only reads when ini_mix = const_mix)
const_mix = {} 

# ====== Setting up photochemistry ======
use_photo = True
# astronomy input
r_star = 1 # stellar radius in solar radius
Rp = 1.*7.1492E9 # Planetary radius (cm) (for computing gravity)
orbit_radius = 5.2 # planet-star distance in A.U.
sl_angle = 60. /180.*3.14159 # the zenith angle of the star in degree (usually 58 deg for the dayside average)
f_diurnal = 0.5 # to account for the diurnal average of solar flux (i.e. 0.5 for Earth; 1 for tidally-locked planets) 
scat_sp = ['H2', 'He'] # the bulk gases that contribute to Rayleigh scattering
T_cross_sp = [] # warning: slower start! available atm: 'CO2','H2O','NH3', 'SH','H2S','SO2', 'S2', 'COS', 'CS2'

edd = 0.5 # the Eddington coefficient 
dbin1 = 0.1  # the uniform bin width < dbin_12trans (nm)
dbin2 = 2.   # the uniform bin width > dbin_12trans (nm)
dbin_12trans = 240. # the wavelength switching from dbin1 to dbin2 (nm)

# the frequency to update the actinic flux and optical depth
ini_update_photo_frq = 100
final_update_photo_frq = 5

# ====== Setting up ionchemistry ======
use_ion = False
if use_photo == False and use_ion == True:
    print ('Warning: use_ion = True but use_photo = False')
# photoionization needs to run together with photochemistry


# ====== Setting up parameters for the atmosphere ======
atm_base = 'H2' #Options: 'H2', 'N2', 'O2', 'CO2 -- the bulk gas of the atmosphere: changes the molecular diffsion, thermal diffusion factor, and settling velocity
rocky = False # for the surface gravity
nz = 150   # number of vertical layers
P_b = 5e9  # pressure at the bottom (dyne/cm^2)
P_t = 1e-2 # pressure at the top (dyne/cm^2)
use_Kzz = True
use_moldiff = True
use_vm_mol = False # use upwind scheme for molecular diffusion -- under testing
use_vz = False
atm_type = 'file'  # Options: 'isothermal', 'analytical', 'file', or 'vulcan_ini' 'table'
Kzz_prof = 'file' # Options: 'const','file' or 'Pfunc' (Kzz increased with P^-0.4)
K_max = 1e5        # for Kzz_prof = 'Pfunc'
K_p_lev = 0.1      # for Kzz_prof = 'Pfunc'
vz_prof = 'const'  # Options: 'const' or 'file'
gs = 2479.         # surface gravity (cm/s^2)  (HD189:2140  HD209:936)
Tiso = 1000 # only read when atm_type = 'isothermal'
# setting the parameters for the analytical T-P from (126)in Heng et al. 2014. Only reads when atm_type = 'analytical' 
# T_int, T_irr, ka_L, ka_S, beta_S, beta_L
para_warm = [120., 1500., 0.1, 0.02, 1., 1.]
para_anaTP = para_warm
const_Kzz = 1.E10 # (cm^2/s) Only reads when use_Kzz = True and Kzz_prof = 'const'
const_vz = 0 # (cm/s) Only reads when use_vz = True and vz_prof = 'const'

# frequency for updating dz and dzi due to change of mu
update_frq = 100 

# ====== Setting up the boundary conditions ======
# Boundary Conditions:
use_topflux = True
use_botflux = False
use_fix_sp_bot = {} # fixed mixing ratios at the lower boundary
diff_esc = [] # species for diffusion-limit escape at TOA
max_flux = 1e13  # upper limit for the diffusion-limit fluxes

# ====== Reactions to be switched off  ======
remove_list = [] # in pairs e.g. [1,2]

# == Condensation ======
use_condense = True
use_relax = ["H2O" , "NH3"]
humidity = 1.
use_settling = True

condense_sp = ["H2O" , "NH3"]      
non_gas_sp = ['H2O_l_s', 'NH3_l_s']
r_p = {'H2O_l_s': 5e-5, 'NH3_l_s': 5e-5}
rho_p = {'H2O_l_s': 0.9, 'NH3_l_s': 0.7}
humidity = 1. # for water on Earth only
fix_species = ["H2O" , "NH3", 'H2O_l_s', 'NH3_l_s']      # fixed the condensable species after condensation-evapoation EQ has reached  
start_conden_time = 1e6
stop_conden_time = 1e8  
fix_species_time = stop_conden_time # after this time to fix the condensable species
fix_species_from_coldtrap_lev = True

use_ini_cold_trap = True  # cold trap from the get go
use_sat_surfaceH2O = False # for terrestriall atmospheres where H2O is controlled by surface T

# ====== steady state check ======
st_factor = 0.5
conv_step = 500

# ====== Setting up numerical parameters for the ODE solver ====== 
ode_solver = 'Ros2' # case sensitive
use_print_prog = True
use_print_delta = False
print_prog_num = 500  # print the progress every x steps 
dttry = 1.E-10
trun_min = 1e2
runtime = 1.E22
dt_min = 1.E-14
dt_max = runtime*1e-5
dt_var_max = 2.
dt_var_min = 0.5
count_min = 120
count_max = int(1E4)
atol = 1.E-1 # Try decreasing this if the solutions are not stable
mtol = 1.E-18
mtol_conv = 1.E-16
pos_cut = 0
nega_cut = -1.
loss_eps = 1e12 # for using BC
yconv_cri = 0.01 # for checking steady-state
slope_cri = 1.e-4
yconv_min = 0.1
flux_cri = 0.1
flux_atol = 1. # the tol for actinc flux (# photons cm-2 s-1 nm-1)

# ====== Setting up numerical parameters for Ros2 ODE solver ====== 
rtol = 0.05             # relative tolerence for adjusting the stepsize 
post_conden_rtol = 0.15 # switched to this value after fix_species_time

# ====== Setting up for ouwtput and plotting ======
# plotting:
plot_TP = False
use_live_plot = True
use_live_flux = False
use_plot_end = False
use_plot_evo = False
use_save_movie = False
use_flux_movie = False
plot_height = False
use_PIL = True 
live_plot_frq = 10
save_movie_rate = live_plot_frq
y_time_freq = 1  #  storing data for every 'y_time_freq' step
plot_spec = ['H2O', 'H2O_l_s', 'CH4', 'CO', 'C2H2', 'NH3', 'NH3_l_s' ]
# output:
output_humanread = False
use_shark = False
save_evolution = False   # save the evolution of chemistry (y_time and t_time) for every save_evo_frq step
save_evo_frq = 10
