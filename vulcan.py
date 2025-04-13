#!/usr/bin/env python3

# ==============================================================================
# This is the main file of VULCAN: the chemical kinetics code.
# Copyright (C) 2016 Shang-Min Tsai (Shami)
# ==============================================================================

# Import system modules
import importlib
import time
import sys

# Import some VULCAN modules
from paths import COM_FILE
from make_chem_funs import make_all
from config import Config

# import the configuration inputs
def main(vulcan_cfg:Config):

    # Reload chem_funs
    # importlib.reload(chem_funs)

    # Import modules for running VULCAN
    import store, build_atm, op

    ### read in the basic chemistry data
    with open(COM_FILE, 'r') as f:
        columns = f.readline() # reading in the first line
        num_ele = len(columns.split())-2 # number of elements (-2 for removing "species" and "mass")
    type_list = ['int' for i in range(num_ele)]
    type_list.insert(0,'U20'); type_list.append('float')

    ### create the instances for storing the variables and parameters
    data_var = store.Variables(vulcan_cfg)
    data_atm = store.AtmData(vulcan_cfg)
    data_para = store.Parameters(vulcan_cfg)

    # record starting CPU time
    data_para.start_time = time.time()

    # create atmosphere object
    make_atm = build_atm.Atm(vulcan_cfg)

    # for plotting and printing
    output = op.Output(vulcan_cfg)

    # save a copy of the config file
    output.save_cfg()

    # construct pico
    data_atm = make_atm.f_pico(data_atm)

    # construct Tco and Kzz
    data_atm =  make_atm.load_TPK(data_atm)

    # construct Dzz (molecular diffusion)
    # Only setting up ms (the species molecular weight) if vulcan_cfg.use_moldiff == False
    make_atm.mol_diff(data_atm)

    # calculating the saturation pressure
    if vulcan_cfg.use_condense:
        make_atm.sp_sat(data_atm)

    # for reading rates
    rate = op.ReadRate(vulcan_cfg)

    # read-in network and calculating forward rates
    data_var = rate.read_rate(data_var, data_atm)

    # for low-T rates e.g. Jupiter
    if vulcan_cfg.use_lowT_limit_rates:
        data_var = rate.lim_lowT_rates(data_var, data_atm)

    # reversing rates
    data_var = rate.rev_rate(data_var, data_atm)

    # removing rates
    data_var = rate.remove_rate(data_var)

    # initialing y and ymix (the number density and the mixing ratio of every species)
    ini_abun = build_atm.InitialAbun(vulcan_cfg)
    data_var = ini_abun.ini_y(data_var, data_atm)

    # storing the initial total number of atmos
    data_var = ini_abun.ele_sum(data_var)

    # calculating mean molecular weight, dz, and dzi and plotting TP
    data_atm = make_atm.f_mu_dz(data_var, data_atm, output)

    # specify the BC
    make_atm.BC_flux(data_atm)


    # ============== Execute VULCAN  ==============
    # time-steping in the while loop until conv() returns True or count > count_max

    # setting the numerical solver to the desinated one in vulcan_cfg
    solver = op.Ros2(vulcan_cfg)

    # Setting up for photo chemistry
    if vulcan_cfg.use_photo:
        rate.make_bins_read_cross(data_var, data_atm)
        #rate.read_cross(data_var)
        make_atm.read_sflux(data_var, data_atm)

        # computing the optical depth (tau), flux, and the photolisys rates (J) for the first time
        solver.compute_tau(data_var, data_atm)
        solver.compute_flux(data_var, data_atm)
        solver.compute_J(data_var, data_atm)
        # they will be updated in op.Integration by the assigned frequence

        # removing rates
        data_var = rate.remove_rate(data_var)

    # Assgining the specific solver corresponding to different B.C.s
    integ = op.Integration(solver, output, vulcan_cfg)
    solver.naming_solver(data_para)

    # Running the integration loop
    integ(data_var, data_atm, data_para, make_atm)

    # Save result to disk
    output.save_out(data_var, data_atm, data_para)


if __name__ == "__main__":
    # Entry point for the script when run directly

    # Make config
    vulcan_cfg = Config()

    # no arguments or not setting '-n' (no re-making chem_funs.py) option
    if len(sys.argv) < 2 or sys.argv[1] != '-n':
        # running prepipe to construch chem_funs.py
        print ('Making chem_funs.py ...')
        make_all(vulcan_cfg)

    # Run model
    main(vulcan_cfg)
