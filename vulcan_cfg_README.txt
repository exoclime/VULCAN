# ============================================================================= 
# Configuration file of VULCAN:  
# ============================================================================= 

# ====== Setting up the elements included in the network ======
atom_list = ['H', 'O', 'C', 'N'] # only for the purpose of checking element conservation
# ====== Setting up paths and filenames for the input and output files  ======
# input:
network = 'thermo/NCHO_photo_network.txt' # the path to the chemical network file
gibbs_text = 'thermo/gibbs_text.txt' # (all the nasa9 files must be placed in the folder: thermo/NASA9/)
cross_folder = 'thermo/photo_cross/' # the path to the photolysis cross sections
com_file = 'thermo/all_compose.txt'  # the file for basic chemistry peroperties 
atm_file = 'atm/atm_HD189_Kzz.txt' # TP and Kzz (optional) file
sflux_file = 'atm/stellar_flux/sflux-HD189_Moses11.txt' # the file for stellar flux (the flux density is defined at the stellar surface)
top_BC_flux_file = 'atm/BC_top.txt' # the file for the top boundary conditions
bot_BC_flux_file = 'atm/BC_bot.txt' # the file for the lower boundary conditions
vul_ini = 'output/HD189-nominal.vul' # the file to initialize the abundances for ini_mix = 'vulcan_ini'
# output:
output_dir = 'output/'    # output path
plot_dir = 'plot/'        # output plot path
movie_dir = 'plot/movie/' # output path for live–plotting for making movies 
out_name =  'HD189.vul'   # output file name

# ====== Setting up the elemental abundance ======
use_solar = False # True: using the solar abundance from Table 10. K.Lodders 2009; False: using the customized elemental abundance. 
# Below are customized elemental abundances (only read when use_solar = False)
# X_H is the X to H ratio
O_H = 6.0618E-4 *(0.85) #
C_H = 2.7761E-4  
N_H = 8.1853E-5
S_H = 1.3183E-5
He_H = 0.09692

ini_mix = 'EQ' # The initial abundances (Options: 'EQ', 'const_mix', 'vulcan_ini', 'table' ...for 'vulcan_ini, the T-P grids have to be exactly the same)
# only for ini_mix = 'const_mix': initialsing uniform (constant with pressure) mixing ratios
const_mix = {'CH4':2.7761E-4*2, 'O2':4.807e-4, 'He':0.09691, 'N2':8.1853E-5, 'H2':0.9} 

# ====== Setting up photochemistry ======
use_photo = True  # include photochemistry or not
# Below are astronomy input
r_star = 0.805 # stellar radius (solar radius)
Rp = 1.138*7.1492E9 # planetary radius (cm) (only for computing g(z))
orbit_radius = 0.03142 # planet-star distance (AU)
sl_angle = 58 /180.*3.14159 # the zenith angle of the star in degree (58 deg for the dayside average)
f_diurnal = 1. # the factor mutiplying too all photolysis rates (to account for the diurnal average of solar flux, i.e. 0.5 for Earth; 1 for tidally-locked planets) 
scat_sp = ['H2', 'He'] # the bulk gases that contribute to Rayleigh scattering
T_cross_sp = [] # T-dependent cross section (warning: slower start! available atm: 'CO2','H2O','NH3', 'SH','H2S','SO2', 'S2', 'COS', 'CS2')

edd = 0.5 # the Eddington coefficient 
dbin1 = 0.1  # the uniform bin width (nm) in the VUV ( < dbin_12trans )
dbin2 = 2.   # the uniform bin width (nm) in the ~MUV-NUV ( > dbin_12trans )
dbin_12trans = 240. # the wavelength (nm) switching from dbin1 to dbin2 

# the frequency to update the actinic flux and optical depth
ini_update_photo_frq = 100 
final_update_photo_frq = 5

# ====== Setting up ionchemistry ======
use_ion = False  # include ionchemistry 
if use_photo == False and use_ion == True:
    print ('Warning: use_ion = True but use_photo = False')
# photoionization usually needs to run together with photochemistry


# ====== Setting up parameters for the atmosphere ======
atm_base = 'H2' # Options: 'H2', 'N2', 'O2', 'CO2 -- the bulk gas of the atmosphere: changes the molecular diffsion, thermal diffusion factor, and settling velocity
rocky = False # only for setting up the surface gravity
nz = 120   # number of vertical layers
P_b = 1e9  # pressure at the bottom (dyne/cm^2)
P_t = 1e-2 # pressure at the top (dyne/cm^2)
use_Kzz = True # include eddy diffusion 
use_moldiff = True # include molecular diffusion
use_vz = False     # include vertical advection
atm_type = 'file'  # the way to prescribe P-T (options: 'isothermal', 'analytical', 'file', 'vulcan_ini', or 'table')
Kzz_prof = 'Pfunc' # the types of eddy-diffusion profile (options: 'const', 'file' or 'Pfunc' -- Kzz increased with P^-0.4 above "K_p_lev" level)
K_max = 1e5        # for Kzz_prof = 'Pfunc': Kzz of the constant bit (cm^2/s)
K_p_lev = 0.1      # for Kzz_prof = 'Pfunc': the transition level from constant Kzz to P^-0.4 Kzz; usually the boundary of convective and radiative zone (bar)
vz_prof = 'const'  # the types verical advection (options: 'const' or 'file')
gs = 2140.         # surface gravity (cm/s^2) 
Tiso = 1000 # only read when atm_type = 'isothermal'

# setting the parameters for the analytical T-P from (126)in Heng et al. 2014. Only reads when atm_type = 'analytical' 
# T_int, T_irr, ka_L, ka_S, beta_S, beta_L (details see Heng et al. 2014)
para_warm = [120., 1500., 0.1, 0.02, 1., 1.]
para_anaTP = para_warm
const_Kzz = 1.E10 # Kzz (cm^2/s) for constant eddy diffusion (only reads when use_Kzz = True and Kzz_prof = 'const')
const_vz = 0 # verical wind (cm/s) (only reads when use_vz = True and vz_prof = 'const')

# frequency (every X steps) for updating dz and dzi due to change of mean molecular weight
update_frq = 100 

# ====== Setting up the boundary conditions ======
use_topflux = False # apply the upper boundary conditions from "top_BC_flux_file"
use_botflux = False # apply the lower boundary conditions from "top_BC_flux_file"
use_fix_sp_bot = {} # fixed mixing ratios at the lower boundary (e.g. use_fix_sp_bot = {'H2O':0.01} prescribes [H2O] = 0.01 at the surface)
diff_esc = ['H'] # species for diffusion-limited escape at TOA
max_flux = 1e13  # (cm^2/s) upper limit for the diffusion-limit fluxes

# ====== Reactions to be switched off  ======
remove_list = [] # removing the reactions (usually in forward/reverse pairs  e.g. [1,2])

# == Condensation ======
use_condense = False    # include condensation reactions
use_settling = False    # include gravitational setteling of particles
start_conden_time = 0   # time to start condensation
condense_sp = []        # condensable species
non_gas_sp = []         # non-gaseous species (e.g.: H2O_l_s)
fix_species = []        # the condensable species to be fixed after condensation-evapoation EQ has reached  
fix_species_time = 0    # after this time to fix the condensable species

# ====== steady state check ======
st_factor = 0.5         # checking steady-state for st_factor of the integration time
conv_step = 500         # capping the step to check for steady-state (e.g. conv_step = 500  means using the difference between N-500 and N at most) 
 
# ====== Setting up numerical parameters for the ODE solver ====== 
ode_solver = 'Ros2'     # default: the 2nd-order Rosenberg solver
use_print_prog = True   # option to print some integration info 
use_print_delta = False # option to print delta (truncation error)
print_prog_num = 500    # print the progress every X steps 
dttry = 1.E-10          # starting timestep (s)
trun_min = 1e2          # mininum of total model time (s)
runtime = 1.E22         # maximum of total model time (s)
dt_min = 1.E-14         # mininum of timestep (s)
dt_max = runtime*1e-3   # maximum of timestep (s)
dt_var_max = 2.         # constraints of varying the timestep
dt_var_min = 0.5        # constraints of varying the timestep
count_min = 120         # mininum of total number of steps
count_max = int(3E4)    # maximum of total number of steps
# the model will run for at least count_min steps and trun_min seconds and will stop when either time exceeds runtime or steps have reached count_max 

atol = 1.E-1            # absolute tolerance  (try decreasing this if the solutions are not stable)
mtol = 1.E-22           # relative tolerance for adjusting stepsize 
mtol_conv = 1.E-20      # relative tolerance for checking steady-state
# the solver will ignore the abundance below atol (molecule/cm^2) and mixing ratio below mtol 

pos_cut = 0             # clipping small numbers (normally no need) 
nega_cut = -1.          # clipping small negative numbers
loss_eps = 1e-1         # criteria for element conservation
yconv_cri = 0.01        # dy for checking steady-state
slope_cri = 1.e-4       # dy/dt for checking steady-state
yconv_min = 0.1         # alternative condition for checking steay-state
flux_atol = 1.          # the absolute tolerence for actinc flux (# photons cm-2 s-1 nm-1)        
flux_cri = 0.1          # the thershold of change of actinic flux for checking steady-state
# in conv(), the steady-state is considered reached when (dy < yconv_cri and dy/dt < slope_cri or dy < yconv_min and dydt < slope_min) and aflux_change < flux_cri 


# ====== Setting up numerical parameters for Ros2 ODE solver ====== 
rtol = 0.5             # relative tolerence for adjusting the stepsize 
post_conden_rtol = 0.1 # switched to this value after fix_species_time

# ====== Setting up for ouwtput and plotting ======
# plotting:
plot_TP = False        # plotting TP at the start    
use_live_plot = True   # turning on live-time plot
use_live_flux = False  # turning on live-time plot for the actinic flux
use_plot_end = False    
use_plot_evo = False   # plot the abundance evolution at j level at the end
use_save_movie = False # saving the live-time plotted files in the folder asigned in movie_dir/
use_flux_movie = False
plot_height = False    # use height coordinate instead of pressure
use_PIL = True         
live_plot_frq = 10     # every N step to make the live-time plot 
save_movie_rate = live_plot_frq
y_time_freq = 1        #  storing data for every 'y_time_freq' step
plot_spec = ['H2O', 'H', 'CH4', 'CO', 'CO2', 'HCN', 'NH3' ] # the list of species for the live-time plot
# output:
output_humanread = False
use_shark = False
save_evolution = False   # save the evolution of chemistry (y_time and t_time) for every save_evo_frq step
save_evo_frq = 10        # every N step to save when save_evolution = True
