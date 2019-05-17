# ============================================================================= 
# Summary for the configuration file of VULCAN (vulcan_cfg.py)
# More detailed instructions can be found in README on GitHub.
# The variable types are indicated in the brackets:
# int: integer  str: string  float: float  bool: Boolean   
# ============================================================================= 

# ====== Setting up the elements included in the network ======
atom_list = ['H', 'O', 'C', 'He', 'N'] : included elements (list) --for the purpose of checking atom conservation

# ====== Set up paths and filenames for the input and output files  ======
network = 'CHO_network.txt' : the filename of the desinated chemical network (str)
gibbs_text = 'thermo/gibbs_text.txt' : the input text file for prepine.py (str)

cross_folder = 'thermo/leiden_cross/' : the path of the folder for the UV cross sections (str)
com_file = 'thermo/HOC_compose.txt' : the file name of the basic property of chemical species (str)
atm_file = 'atm/atm_HD189_Kzz.txt' : the file name of the input T-P(-Kzz) profile (str)

sflux_file = 'atm/flux-HD189_Moses11.txt' : The stellar flux at the stellar surface (str)
top_BC_mix_file = 'atm/BC_top.txt' : the file name of the top boundary for constant flux/mixing boundary conditions (str)
bot_BC_mix_file = 'atm/BC_bot.txt' : the file name of the bottom boundary for constant flux/mixing boundary conditions (str)
vul_ini = 'output/moses_HD189.vul' : the file name of the vulcan output as the intial mixing-ratio profiles, for ini_mix = 'vulcan_ini' (str)

output_dir = 'output/' : the name of the output directory (str)
plot_dir = 'plot/' the name of the plot directory (str)
movie_dir = 'plot/movie/new-HD189/' : the name of the movie directory, required for use_save_movie = True (str)
out_name = 'HD189_K9.vul' : the name of the output file (str)
out_y_time_freq = 10 : the frequency (for every X steps) to store the temporal data (int)

# ====== Setting up the photochemistry ======
use_photo = True 
excit_sp = ['O_1', 'CH2_1'] # N_D to avoid in the initial abundances by fc
scat_sp = ['H2', 'He'] # # the molecules that contribute to Rayleigh scattering
r_star = 0.752 #0.752 HD209: 1.118
orbit_radius = 0.03142 #0.03142 # planet-star distance in A.U.
sl_angle = 48 /180.*3.14159 # the zenith angle of the star  
edd = 0.669 #(cos(48 deg) )  # the Eddington coefficient 
dbin = 0.2

# ====== Setting up the elemental abundance ======
# default: solar abundance (from K.Lodders 2009)
O_H = 6.0618E-4 : O/H ratio (float)
C_H = 2.7761E-4 : C/H ratio (float)
N_H = 8.1853E-5 : N/H ratio (float)
He_H = 0.09691  : He/H ratio (float)
ini_mix = 'fc'  : defines the way to set the initial composition, can be 'fc', 'const_mix', 'const_lowT', 'fc_precal', or 'vulcan_ini'
				  'fc': Using chemical equilibrium by running fastchem
				  'const_mix' : assigning constant-with-altitude mixing ratios
const_mix = {'CO2':1.81E-3, 'He':0.136,'N2':0.8, 'H2':1.-1.81E-3-0.136-0.8 } : for ini_mix = 'const_mix'


# ====== Reactions to be switched off  ======
remove_list = [] : the reactions (by the index in the chemical network) to be set to zero
				   should usually be a pair for the forward and reverse reactions, e.g. [1,2] (list)

# ====== Setting up parameters for the atmosphere ======
atm_base = 'H2' : the bulk gas of the atmosphere: affects molecular diffsion and some 3-body reactions (only 'H2' for the moment)
nz = 160 : number of vertical levels (int)
use_Kzz = True : the switch of eddy diffusion (bool)
use_vz = 0 : the switch of vertical advection (bool)
use_moldiff = True : the switch of molecular diffusion (bool)

# Boundary Conditions:
use_topflux = False
use_botflux = False
use_fix_bot = False 

atm_type = 'file' : can be 'isothermal', 'analytical', or 'file' (str)
                   'isothermal': constant T with the assigned value Tiso 
				   'analytical': analytical T-P profile from (126) in Heng, Mendonca & Lee (2014)
				   'file': read in temperature, pressure (and eddy diffusion) from the input T-P(-Kzz) file 

Kzz_prof = 'file': can be 'const' or 'file' (str) 
				   'const': constant eddy diffusion
				   'file': Kzz read in from the input T-P(-Kzz) file (only when atm_type = 'file') 

vz_prof = 'const' : can be 'const' or 'file' (str) 

g = 2140. : gravitational acceleration (cm/s^2)  (float)	
		   (2140 for HD 189733b and 936 for HD 209458b)			   			   
Tiso = 1000. : the constant value of temperature (K) for atm_type = 'isothermal' (float)
para_anaTP = [0., 1500., 0.01, 0.001, 1., 1.] : T_int, T_irr, ka_L, ka_S, beta_S, beta_L 
			 , the six parameters of the analytical T-P profile
const_Kzz = 1.E9 :  the constant value of Kzz (cm^2/s) for Kzz_prof = 'const'
const_vz = 0 : the constant value of vz (cm/s) for vz_prof = 'const'

P_b = 1.E9 : pressure (dyne/cm^2) at the bottom (float)
P_t = 1.E2 : pressure (dyne/cm^2) at the top (float)


# == Condensation still testing...
use_condense = False : the switch of condensation (bool)
#start_conden_time = 1e10
condesne_sp = ["H2O"] : the list of condense species
non_gas_sp = ["H2O_l_s"] : the list of non-gaseous (condensed species)

# ====== Setting up general parameters for the ODE solver ====== 
ode_solver = 'Ros2' : can be 'Ros2', 'SemiEU', or 'SparSemiEU'  (str)
					 For Rosenbrock, semi-implicit Euler, or semi-implicit Euler with sparse matrix respectivey. 
					 'Ros2' is highly recommended. 	
use_print_prog = False : print out ODE time-stepping information (bool)
use_height = False : live plotting in hight instead of pressure (bool)
print_prog_num = 200 : the frequency (for every X steps) to print out the ODE info (int)
use_live_plot = True : real-time plotting (bool)
use_save_movie = False
use_flux_movie = False
live_plot_frq = 10 : the frequency (for every X steps) of the real-time plotting (int)

save_movie_rate = live_plot_frq


use_plot_end = True : plotting the steady-state mixing ratios in the end (bool)
use_plot_evo = False : plotting the temporal evolution at certain level in the end (bool)
plot_TP = True : plotting the input T-P profile (bool)
output_humanread = False : 
plot_spec = ['H', 'H2', 'CH3', 'CH4', 'C2H2', 'CO', 'CH3OH', 'CH2OH', 'He'] : the species to be plotted for the temporal evolution
live_plot_spec = ['H', 'H2', 'H2O', 'CH4', 'CO', 'CO2', 'C2H2', 'C2H4', 'C2H6', 'CH3OH'] : the species to be plotted for the real-time and steady-state plotting


# ====== steady state check ======
st_factor = 0.5  : the steady-state factor i.e. t - st_factor * t will be used to check for steady-state 
				  May try larger st_factor when it is not converged after long time, especially when T < 1000K (float)
count_min = 120 : the minimum steps to stop (int)
 
# ====== Setting up numerical parameters for the ODE solver ====== 
dttry = 1.E-8 : the initial value of the stepsize (sec)  (float)
dt_std = 1. : the final value of the stepsize to decrease to for the steady-state (float)
trun_min = 1.E3 : the minimum runtime (float)
runtime = 1.E24 : the maximum runtime (float)
dt_min = 1.E-14 : the minimum stepsize (float)
dt_max = runtime*0.01 : the maximum stepsize (float)
dt_var_max = 2. : the maximum increment of the stepsize (float)
dt_var_min = 0.2 : the maximum decrease of the stepsize (float)
atol = 1.E-3 : the lower limit of the absolute number density to compute determine accept/reject the solution after each step (float)
			  i.e. ignore the trace species with number density below this value when calculating the variation of the solution
mtol = 1.E-20 : the lower limit of the mixing ratio to compute determine accept/reject the solution after each step (float) 
mtol_conv = 1.E-26 : the lower limit of the mixing ratio to check for convergence (float)
pos_cut = 0 : the positive value for clipping (float)
nega_cut = -1. : the negative value for clipping (float)
loss_eps = 1e-4 : the tolerance for partical conservation (float)
yconv_cri = 0.05 : the relative variation for checking the steady-state (float)
slope_cri = 1.e-4 : the relative slope (devided by time) for checking the steady-state (float)

yconv_min = 0.1
flux_cri = 5.e-2 : the critaria for the actinic flux to converge
flux_atol = 1. # the tol for actinc flux (# photons cm-2 s-1 nm-1)

count_max = int(1E5)
count_max = int(2E4) : the maximum steps (int)
update_frq = 100 : the frequency (for every X steps) to update the layer width, dz and dzi, due to the change of mean molecular weight (int)

# ====== Setting up numerical parameters for Ros2 ODE solver ====== 
rtol = 0.05 : the relative tolerance for the truncation error (try smaller value if numerical unstable)

# ====== Setting up numerical parameters for SemiEu/SparSemiEU ODE solver ====== 
PItol = 0.1 : the relative tolerance for the PID control

use_PIL = True : use PIL package (for plotting in the native picture viewer)

Note:
1. All the thermodynamics files are placed in the folder ./thermo and the NASA9 polynomials are in ./thermo/NASA9/
2. The atmospheric T-P(-Kzz) files are placed in ./atm

