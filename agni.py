# Wrapper functions for coupling VULCAN to AGNI

# Import system modules
import logging
import os
import numpy as np
from juliacall import Main as jl
from scipy.interpolate import PchipInterpolator
from scipy.integrate import trapezoid

# initialise logger for vulcan
log = logging.getLogger("fwl."+__name__)

# Import some VULCAN modules
import paths
from config import Config
from store import AtmData, Variables
from phy_const import r_sun, au
from chem_funs import spec_list as gas_list

# constants
AGNI_LOGFILE_NAME = "agni_recent.log"

def activate_julia(vulcan_cfg:Config):

    log.info("Activating Julia environment")
    jl.seval("using Pkg")
    jl.Pkg.activate(paths.AGNI_DIR)

    # Plotting configuration
    jl.seval('ENV["GKSwstype"] = "100"')
    jl.seval("using Plots")
    jl.seval('default(label=nothing, dpi=250)')

    # Import AGNI
    jl.seval("import AGNI")

    # Setup logging from AGNI
    #    This handle will be kept open throughout, so the file
    #    should not be deleted at runtime. However, it will be emptied when appropriate.
    verbosity = 1
    logpath = os.path.join(vulcan_cfg.output_dir, AGNI_LOGFILE_NAME)
    jl.AGNI.setup_logging(logpath, verbosity)

    log.debug("AGNI will log to '%s'"%logpath)

def init_agni_atmos(vulcan_cfg:Config, atm:AtmData, var:Variables):
    log.debug("New AGNI atmosphere")

    atmos = jl.AGNI.atmosphere.Atmos_t()

    # Stellar spectrum path (scaled to stellar surface)
    sflux_path  = vulcan_cfg.sflux_file

    # Calculate instellation by scaling spectrum to TOA
    sflux_data = np.loadtxt(sflux_path).T # erg/s/cm^2/nm
    sflux_integ = trapezoid(sflux_data[1], x=sflux_data[0]) # erg/s/cm^2
    sflux_integ *= (vulcan_cfg.r_star*r_sun)**2  / (vulcan_cfg.orbit_radius*au)**2
    sflux_integ *= 0.001 # W/m^2

    # Spectral file path
    try_spfile = os.path.join(vulcan_cfg.output_dir , "runtime.sf")
    if os.path.exists(try_spfile):
        # exists => don't modify it
        input_sf =      try_spfile
        input_star =    ""
    else:
        # doesn't exist => AGNI will copy it + modify as required
        input_sf =      os.path.join(paths.AGNI_DIR,vulcan_cfg.spectral_file)
        input_star =    sflux_path

    # composition set initially well-mixed
    vol_dict = {}
    for g in gas_list:
        if "_" in g:
            continue
        vol_dict[g] = np.median(var.ymix[:,gas_list.index(g)])

    # set condensation
    condensates = []

    # Boundary pressures
    p_surf = atm.pico[0] * 1e-6 # convert to bar
    p_top  = atm.pico[-1]* 1e-6 # convert to bar

    # Setup struct
    jl.AGNI.atmosphere.setup_b(atmos,
                        paths.AGNI_DIR, vulcan_cfg.output_dir, input_sf,

                        sflux_integ,
                        vulcan_cfg.f_diurnal,
                        0.0,
                        vulcan_cfg.sl_angle*np.pi/180.0, # convert to degrees

                        vulcan_cfg.Tsurf_guess,
                        vulcan_cfg.gs*0.01, # convert to SI, m/s^2
                        vulcan_cfg.Rp*0.01, # convert to SI, m

                        vulcan_cfg.agni_nlev,
                        p_surf,
                        p_top,

                        vol_dict, "",

                        flag_rayleigh  = vulcan_cfg.use_rayleigh,
                        flag_cloud     = False,
                        overlap_method = 'ee',

                        albedo_s=vulcan_cfg.surf_albedo,
                        surface_material='greybody',
                        condensates=condensates,
                        use_all_gases=False,
                        fastchem_work = '',
                        real_gas = False,

                        skin_d=0.01,
                        skin_k=2.0,
                        tmp_magma=vulcan_cfg.Tsurf_guess, 
                        tmp_floor=50.0
                        )

    # Allocate arrays
    jl.AGNI.atmosphere.allocate_b(atmos,input_star)

    # Set temperature profile to initial guess
    tmp_top = 500.0
    tmp_top = min(tmp_top, vulcan_cfg.Tsurf_guess)
    log.debug("Initialised log-linear (top = %.2f K)"%tmp_top)
    jl.AGNI.setpt.loglinear_b(atmos, -0.5 * vulcan_cfg.Tsurf_guess)
    jl.AGNI.setpt.stratosphere_b(atmos, tmp_top)

    return atmos


def deallocate_atmos(atmos):
    """
    Deallocate atmosphere struct
    """
    jl.AGNI.atmosphere.deallocate_b(atmos)


def _solve_energy(atmos, vulcan_cfg:Config,):
    """Use AGNI to solve for energy-conserving solution. """

    # atmosphere solver plotting frequency
    modplot = 0
    if vulcan_cfg.use_live_plot:
        modplot = 1

    # tracking
    agni_success = False  # success?
    attempts = 0          # number of attempts so far

    # make attempts
    while not agni_success:
        attempts += 1
        log.info("Attempt %d" % attempts)

        # default parameters
        linesearch  = 2
        easy_start  = True
        dx_max      = 300.0
        ls_increase = 1.04
        perturb_all = True
        max_steps   = 70

        # try different solver parameters if struggling
        if attempts == 2:
            linesearch  = 1
            dx_max     *= 2.0
            ls_increase = 1.1

        log.debug("Solver parameters:")
        log.debug("    ls_method=%d, easy_start=%s, dx_max=%.1f, ls_increase=%.2f"%(
            linesearch, str(easy_start), dx_max, ls_increase
        ))

        # Try solving temperature profile
        agni_success = jl.AGNI.solver.solve_energy_b(atmos,
                            sol_type  = int(3),
                            method    = int(1),
                            chem_type = int(0),

                            conduct=False, convect=True, sens_heat=True,
                            latent=False, rainout=True,

                            max_steps=int(max_steps), max_runtime=900.0,
                            conv_atol=float(vulcan_cfg.agni_atol),
                            conv_rtol=float(vulcan_cfg.agni_rtol),

                            ls_increase=float(ls_increase), ls_method=int(linesearch),
                            dx_max=float(dx_max), easy_start=easy_start,
                            perturb_all=perturb_all,

                            save_frames=False, modplot=int(modplot),
                            plot_jacobian=False
                            )

        # Model status check
        if agni_success:
            # success
            log.info("Attempt %d succeeded" % attempts)
            break
        else:
            # failure
            log.warning("Attempt %d failed" % attempts)

            # Max attempts
            if attempts >= 2:
                log.error("Maximum attempts when executing AGNI")
                break
    return atmos

def _solve_once(atmos):
    """Use AGNI to solve radiative transfer with prescribed T(p) profile"""

    #    dry convection
    jl.AGNI.setpt.dry_adiabat_b(atmos)

    #    temperature floor in stratosphere
    jl.AGNI.setpt.stratosphere_b(atmos, 0.5)

    # solve fluxes
    jl.AGNI.energy.calc_fluxes_b(atmos, False, True, False, False, calc_cf=True)

    # fill kzz values
    jl.AGNI.energy.fill_Kzz_b(atmos)

    return atmos

def run_agni(atmos, vulcan_cfg:Config, atm:AtmData, var:Variables):
    """Run AGNI atmosphere model.  """

    # Inform
    log.info("Running climate calculation...")

    # update gas compositions (interpolate from VULCAN to AGNI)
    for g in atmos.gas_names:
        p_vul = atm.pco[::-1]*0.1 # convert to Pa
        g_vul = np.clip(var.ymix[::-1,gas_list.index(g)], 1e-30, 1.0)
        g_itp = PchipInterpolator(p_vul, np.log10(g_vul))

        atmos.gas_vmr[g][:]  = 10**g_itp(atmos.p)[:]
        atmos.gas_ovmr[g][:] = atmos.gas_vmr[g][:]

    jl.AGNI.plotting.plot_vmr(atmos,os.path.join(vulcan_cfg.plot_dir,"_agni_vmr.png"))

    # Run model
    if vulcan_cfg.solve_rce:
        log.info("Using nonlinear solver to conserve fluxes")
        atmos = _solve_energy(atmos, vulcan_cfg)
    else:
        log.info("Using prescribed temperature profile")
        atmos = _solve_once(atmos)

    # Make plots
    jl.AGNI.plotting.plot_pt(atmos,os.path.join(vulcan_cfg.plot_dir, "_agni_tp.png"))

    # Write output data
    # ncdf_path = os.path.join(dirs["output"],"data",time_str+"_atm.nc")
    # jl.AGNI.save.write_ncdf(atmos, ncdf_path)

    # Parse result (interpolate from AGNI to VULCAN)
    t_itp = PchipInterpolator(atmos.p, atmos.tmp)
    atm.Tco[:] = t_itp(atm.pco*0.1)[:] 

    log.info('------------------------------------------------------------------------')
    return atmos
