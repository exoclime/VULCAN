# ============================================================================= 
# Configuration file of VULCAN:  
# ============================================================================= 

# ====== Set up paths and filenames for the input and output files  ======
network = 'CHO_network.txt'
gibbs_text = 'thermo/gibbs_text.txt'
# all the nasa9 files must be placed in the folder: thermo/NASA9/
com_file = 'thermo/HOC_compose.txt'
atm_file = 'atm/atm_HD189_Kzz.txt'
output_dir = 'output/'
plot_dir = 'plot/'
out_name = 'test.vul'
# storing data for every 'out_y_time_freq' step  
out_y_time_freq = 10 
EQ_ini_file = ''

# ====== Setting up the elemental abundance ======
na = 4 # na: The number of elements. Default is 4: H,O,C,He
atom_list = ['H', 'O', 'C', 'He']
# default: solar abundance (from K.Lodders 2009)
O_H = 6.0618E-4
C_H = 2.7761E-4 
He_H = 0.09691
ini_mix = 'EQ' # 'EQ', 'CH4', or 'CO'

# ====== Reactions to be switched off  ======
remove_list = []

# ====== Setting up parameters for the atmosphere ======
nz = 100
use_Kzz = True
atm_type = 'file' # 'isothermal', 'analytical', or 'file'
Kzz_prof = 'const' # 'const' or 'file'
Tiso = 1000.
# T_int, T_irr, ka_L, ka_S, beta_S, beta_L
para_anaTP = [0., 1500., 0.01, 0.001, 1., 1.]
const_Kzz = 1.E9 # (cm^2/s)
g = 2140 # (cm/s^2)
P_b = 1.E9 #(dyne/cm^2)
P_t = 1.E2  

# ====== Setting up general parameters for the ODE solver ====== 
ode_solver = 'Ros2' # case sensitive
use_print_prog = False
print_prog_num = 200
use_live_plot = True
live_plot_frq = 10
use_plot_end = True
use_plot_evo = True
plot_TP = True
output_humanread = False
plot_spec = ['H', 'H2', 'CH3', 'CH4', 'C2H2', 'CO', 'CH3OH', 'CH2OH', 'He']
live_plot_spec = ['H', 'H2', 'H2O', 'CH4', 'CO', 'CO2', 'C2H2', 'C2H4', 'C2H6', 'CH3OH']

# ====== steady state check ======
st_factor = 0.2  
# Try larger st_factor when T < 1000K
count_min = 100

# ====== Setting up numerical parameters for the ODE solver ====== 
dttry = 1.E-8
dt_std = 1.
trun_min = 1e2
runtime = 1.E24
dt_min = 1.E-14
dt_max = runtime*0.01
dt_var_max = 2.
dt_var_min = 0.2
atol = 1.E-3
mtol = 1.E-20
mtol_conv = 1.E-26
pos_cut = 0
nega_cut = -0.1
loss_eps = 1e-4
yconv_cri = 0.05 # for checking steady-state
slope_cri = 1.e-4
count_max = int(1E4)
update_frq = 100 # for updating dz and dzi due to change of mu

# ====== Setting up numerical parameters for Ros2 ODE solver ====== 
rtol = 0.05

# ====== Setting up numerical parameters for SemiEu/SparSemiEU ODE solver ====== 
PItol = 0.1

use_PIL = True