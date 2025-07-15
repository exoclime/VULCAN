import glob
import logging
import os
from typing import TYPE_CHECKING

import numpy as np
from juliacall import Main as jl
from scipy.interpolate import PchipInterpolator

if TYPE_CHECKING:
    from proteus.config import Config

log = logging.getLogger("fwl."+__name__)


def sync_log_files(outdir:str):
    # Logfile paths
    agni_logpath = os.path.join(outdir, AGNI_LOGFILE_NAME)
    logpath = GetLogfilePath(outdir, GetCurrentLogfileIndex(outdir))

    # Copy logfile content
    with open(agni_logpath, "r") as infile:
        inlines = infile.readlines()

        with open(logpath, "a") as outfile:
            for i,line in enumerate(inlines):
                # First line of agni logfile has NULL chars at the start, for some reason
                if i == 0:
                    line = "[" + line.split("[", 1)[1]
                # copy the line
                outfile.write(line)

    # Remove logfile content
    with open(agni_logpath, "w") as hdl:
        hdl.write("")

def activate_julia(dirs:dict):

    log.info("Activating Julia environment")
    jl.seval("using Pkg")
    jl.Pkg.activate(dirs["agni"])

    # Plotting configuration
    jl.seval('ENV["GKSwstype"] = "100"')
    jl.seval("using Plots")
    jl.seval('default(label=nothing, dpi=250)')

    # Import AGNI
    jl.seval("import AGNI")

    # Setup logging from AGNI
    #    This handle will be kept open throughout the PROTEUS simulation, so the file
    #    should not be deleted at runtime. However, it will be emptied when appropriate.
    verbosity = 1
    logpath = os.path.join(dirs["output"], AGNI_LOGFILE_NAME)
    jl.AGNI.setup_logging(logpath, verbosity)

    log.debug("AGNI will log to '%s'"%logpath)


def _construct_voldict(hf_row:dict, dirs:dict):

    # get from hf_row
    vol_dict = {}
    vol_sum = 0.0
    for vol in gas_list:
        vol_dict[vol] = hf_row[vol+"_vmr"]
        vol_sum += vol_dict[vol]

    # Check that the total VMR is not zero
    if vol_sum < 1e-4:
        UpdateStatusfile(dirs, 20)
        raise ValueError("All volatiles have a volume mixing ratio of zero")

    return vol_dict


def init_agni_atmos(dirs:dict, config:Config, hf_row:dict):
    """Initialise atmosphere struct for use by AGNI.

    Does not set the temperature profile.

    Parameters
    ----------
        dirs : dict
            Dictionary containing paths to directories
        config : Config
            Configuration options and other variables
        hf_row : dict
            Dictionary containing simulation variables for current iteration

    Returns
    ----------
        atmos : atmosphere.Atmos_t
            Atmosphere struct

    """

    log.debug("New AGNI atmosphere")

    atmos = jl.AGNI.atmosphere.Atmos_t()

    # Stellar spectrum path
    sflux_files = glob.glob(os.path.join(dirs["output"], "data", "*.sflux"))
    sflux_times = [ int(s.split("/")[-1].split(".")[0]) for s in sflux_files]
    sflux_path  = os.path.join(dirs["output"],
                                "data", "%d.sflux"%int(sorted(sflux_times)[-1]))

    # Spectral file path
    try_spfile = os.path.join(dirs["output"] , "runtime.sf")
    if os.path.exists(try_spfile):
        # exists => don't modify it
        input_sf =      try_spfile
        input_star =    ""
    else:
        # doesn't exist => AGNI will copy it + modify as required
        input_sf =      get_spfile_path(dirs["fwl"], config)
        input_star =    sflux_path

    # composition
    vol_dict = _construct_voldict(hf_row, dirs)

    # set condensation
    condensates = []
    if config.atmos_clim.agni.condensation:
        if len(vol_dict) == 1:
            # single-gas case
            condensates = list(vol_dict.keys())
        else:
            # get sorted gases (in order of decreasing VMR at the surface)
            vol_sorted = sorted(vol_dict.items(), key=lambda item: item[1])[::-1]

            # Set gases as condensates...
            condensates = ['H2O'] # always prefer H2O
            for v in vol_sorted:
                # add gas if it has non-zero abundance
                if (v[1] > 1e-30) and (v[0] not in condensates):
                    condensates.append(v[0])

            # Remove the least abundant gas from the list, so that we have something to
            #     fill the background with if everything else condenses at a given level
            condensates = condensates[:-1]

    # Chemistry
    chem_type = config.atmos_clim.agni.chemistry
    include_all = False
    fc_dir = "_unset"
    if chem_type == 'eq':
        # equilibrium
        include_all = True
        condensates = []

        # working folder for fastchem coupling
        fc_dir = create_tmp_folder()
        log.debug("Fastchem work folder: '%s'"%fc_dir)

    # Surface single-scattering albedo
    surface_material = config.atmos_clim.agni.surf_material
    if "greybody" in str(surface_material).lower():
        # Grey value
        surface_material = "greybody"
        log.debug("Using grey single-scattering surface properties")

    else:
        # Empirical values
        log.debug(f"Using '{surface_material}' single-scattering surface properties")
        surface_material = os.path.join(dirs["fwl"], surface_material)
        if not os.path.isfile(surface_material):
            raise FileNotFoundError(surface_material)

    # Boundary pressures
    p_surf = hf_row["P_surf"]
    p_top  = config.atmos_clim.agni.p_top
    p_surf = max(p_surf, p_top * 1.1) # this will happen if the atmosphere is stripped

    # Setup struct
    jl.AGNI.atmosphere.setup_b(atmos,
                        dirs["agni"], dirs["output"], input_sf,

                        hf_row["F_ins"],
                        config.orbit.s0_factor,
                        config.atmos_clim.albedo_pl,
                        config.orbit.zenith_angle,

                        hf_row["T_surf"],
                        hf_row["gravity"], hf_row["R_int"],

                        int(config.atmos_clim.agni.num_levels),
                        p_surf,
                        p_top,

                        vol_dict, "",

                        flag_rayleigh = config.atmos_clim.rayleigh,
                        flag_cloud    = config.atmos_clim.cloud_enabled,
                        overlap_method  = config.atmos_clim.agni.overlap_method,

                        albedo_s=config.atmos_clim.surf_greyalbedo,
                        surface_material=surface_material,
                        condensates=condensates,
                        use_all_gases=include_all,
                        fastchem_work = fc_dir,
                        real_gas = config.atmos_clim.agni.real_gas,

                        skin_d=config.atmos_clim.surface_d,
                        skin_k=config.atmos_clim.surface_k,
                        tmp_magma=hf_row["T_surf"], tmp_floor=config.atmos_clim.tmp_minimum
                        )

    # Allocate arrays
    jl.AGNI.atmosphere.allocate_b(atmos,input_star)

    # Set temperature profile from old NetCDF if it exists
    nc_files = glob.glob(os.path.join(dirs["output"],"data","*_atm.nc"))
    if len(nc_files) > 0:
        log.debug("Load NetCDF profile")

        nc_times = [ int(s.split("/")[-1].split("_")[0]) for s in nc_files]
        nc_path  = os.path.join(dirs["output"],
                                "data", "%.0f_atm.nc"%int(sorted(nc_times)[-1]))
        jl.AGNI.setpt.fromncdf_b(atmos, nc_path)

    # Otherwise, set to initial guess
    else:
        tmp_top = 400.0
        tmp_top = min(tmp_top, hf_row["T_surf"])
        log.debug("Initialised log-linear (top = %.2f K)"%tmp_top)
        jl.AGNI.setpt.loglinear_b(atmos, -0.5 * hf_row["T_surf"])
        jl.AGNI.setpt.stratosphere_b(atmos, tmp_top)

    # Logging
    sync_log_files(dirs["output"])

    return atmos


def deallocate_atmos(atmos):
    """
    Deallocate atmosphere struct
    """
    jl.AGNI.atmosphere.deallocate_b(atmos)
    safe_rm(str(atmos.fastchem_work))


def update_agni_atmos(atmos, hf_row:dict, dirs:dict, transparent:bool):
    """Update atmosphere struct.

    Sets the new boundary conditions and composition.

    Parameters
    ----------
        atmos : AGNI.atmosphere.Atmos_t
            AGNI atmosphere struct
        hf_row : dict
            Dictionary containing simulation variables for current iteration
        dirs : dict
            Directories dictionary
        transparent : bool
            If True, the atmosphere is configured to be transparent to radiation

    Returns
    ----------
        atmos : AGNI.atmosphere.Atmos_t
            Atmosphere struct
    """

    # ---------------------
    # Update instellation flux
    atmos.instellation = float(hf_row["F_ins"])

    # ---------------------
    # Update compositions
    vol_dict = _construct_voldict(hf_row, dirs)
    for g in vol_dict.keys():
        atmos.gas_vmr[g][:]  = vol_dict[g]
        atmos.gas_ovmr[g][:] = vol_dict[g]

    # ---------------------
    # Update surface temperature(s)
    atmos.tmp_surf  = float(hf_row["T_surf"] )
    atmos.tmp_magma = float(hf_row["T_magma"])

    # ---------------------
    # Transparent mode?
    if transparent:
        jl.AGNI.atmosphere.make_transparent_b(atmos)
        atmos.tmp[:]  = float(atmos.tmp_surf)
        atmos.tmpl[:] = float(atmos.tmp_surf)
        return atmos

    # ---------------------
    # Store old/current log-pressure vs temperature arrays
    p_old = list(atmos.p)
    t_old = list(atmos.tmp)
    nlev_c = len(p_old)

    #    extend to lower pressures
    p_old = [p_old[0]/10] + p_old
    t_old = [t_old[0]]    + t_old

    #    extend to higher pressures
    p_old = p_old + [p_old[-1]*10]
    t_old = t_old + [t_old[-1]]

    #    create interpolator
    itp = PchipInterpolator(np.log10(p_old), t_old)

    # ---------------------
    # Update surface pressure [Pa] and generate new grid
    atmos.p_boa = 1.0e5 * float(hf_row["P_surf"])
    jl.AGNI.atmosphere.generate_pgrid_b(atmos)

    # ---------------------
    # Set temperatures at all levels
    for i in range(nlev_c):
        atmos.tmp[i]  = float( itp(np.log10(atmos.p[i]))  )
        atmos.tmpl[i] = float( itp(np.log10(atmos.pl[i])) )
    atmos.tmpl[-1]    = float( itp(np.log10(atmos.pl[-1])))

    return atmos


def _solve_energy(atmos, loops_total:int, dirs:dict, config:Config):
    """Use AGNI to solve for energy-conserving solution.

    Parameters
    ----------
        atmos : AGNI.atmosphere.Atmos_t
            Atmosphere struct
        loops_total : int
            Model total loops counter.
        dirs : dict
            Dictionary containing paths to directories
        config : Config
            Configuration options and other variables

    Returns
    ----------
        atmos : AGNI.atmosphere.Atmos_t
            Atmosphere struct
    """

    # atmosphere solver plotting frequency
    modplot = 0
    plot_jacobian = False
    if config.params.out.logging == "DEBUG":
        modplot = 1
        plot_jacobian = True

    # tracking
    agni_success = False  # success?
    attempts = 0          # number of attempts so far

    # make attempts
    while not agni_success:
        attempts += 1
        log.info("Attempt %d" % attempts)

        # default parameters
        linesearch  = 2
        easy_start  = False
        dx_max      = config.interior.spider.tsurf_atol*2+5.0
        ls_increase = 1.01
        perturb_all = True
        max_steps   = 70
        chem_type   = int(config.atmos_clim.agni.chemistry_int)

        # first few iterations
        if loops_total < 3:
            dx_max = 200.0
            ls_increase = 1.1
            max_steps   = 200

        # very first iteration parameters
        if loops_total == 0:
            easy_start  = True
            dx_max      = 300.0

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
                            sol_type  = int(config.atmos_clim.surf_state_int),
                            method    = int(1),
                            chem_type = chem_type,

                            conduct=False, convect=True, sens_heat=True,
                            latent=config.atmos_clim.agni.condensation, rainout=True,

                            max_steps=int(max_steps), max_runtime=900.0,
                            conv_atol=float(config.atmos_clim.agni.solution_atol),
                            conv_rtol=float(config.atmos_clim.agni.solution_rtol),

                            ls_increase=float(ls_increase), ls_method=int(linesearch),
                            dx_max=float(dx_max), easy_start=easy_start,
                            perturb_all=perturb_all,

                            save_frames=False, modplot=int(modplot),
                            plot_jacobian=plot_jacobian
                            )

        # Move AGNI logfile content into PROTEUS logfile
        sync_log_files(dirs["output"])

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

def _solve_once(atmos, config:Config):
    """Use AGNI to solve radiative transfer with prescribed T(p) profile

    Parameters
    ----------
        atmos : AGNI.atmosphere.Atmos_t
            Atmosphere struct
        config : Config
            PROTEUS config object

    Returns
    ----------
        atmos : AGNI.atmosphere.Atmos_t
            Atmosphere struct
    """

    # set temperature profile
    #    rainout volatiles
    rained = jl.AGNI.setpt.prevent_surfsupersat_b(atmos)
    rained = bool(rained)
    if rained:
        log.info("    gases are condensing at the surface")
    #    dry convection
    jl.AGNI.setpt.dry_adiabat_b(atmos)
    #    condensation above
    if config.atmos_clim.agni.condensation:
        for gas in gas_list:
            jl.AGNI.setpt.saturation_b(atmos, str(gas))
    #    temperature floor in stratosphere
    jl.AGNI.setpt.stratosphere_b(atmos, 0.5)

    # do chemistry
    chem_int = config.atmos_clim.agni.chemistry_int
    if chem_int > 0:
        jl.AGNI.chemistry.fastchem_eqm_b(atmos, chem_int, True)

    # solve fluxes
    jl.AGNI.energy.calc_fluxes_b(atmos, False, True, False, False, calc_cf=True)

    # fill kzz values
    jl.AGNI.energy.fill_Kzz_b(atmos)

    return atmos

def _solve_transparent(atmos, config:Config):
    """
    Use AGNI to solve for the surface temperature under a transparent atmosphere

    Parameters
    ----------
        atmos : AGNI.atmosphere.Atmos_t
            Atmosphere struct
        config : Config
            PROTEUS config object

    Returns
    ----------
        atmos : AGNI.atmosphere.Atmos_t
            Atmosphere struct
    """

    atol = float(config.atmos_clim.agni.solution_atol)
    rtol = float(config.atmos_clim.agni.solution_rtol)
    max_steps = 120

    jl.AGNI.solver.solve_transparent_b(atmos,
                                        sol_type=int(config.atmos_clim.surf_state_int),
                                        conv_atol=atol, conv_rtol=rtol,
                                        max_steps=int(max_steps))
    return atmos

def run_agni(atmos, loops_total:int, dirs:dict, config:Config,
                hf_row:dict, transparent:bool):
    """Run AGNI atmosphere model.

    Calculates the temperature structure of the atmosphere and the fluxes, etc.
    Stores the new flux boundary condition to be provided to SPIDER.

    Parameters
    ----------
        atmos : AGNI.atmosphere.Atmos_t
            Atmosphere struct
        loops_total : int
            Model total loops counter.
        dirs : dict
            Dictionary containing paths to directories
        config : Config
            Configuration options and other variables
        hf_row : dict
            Dictionary containing simulation variables for current iteration
        transparent : bool
            Find solution assuming a transparent atmosphere

    Returns
    ----------
        atmos : AGNI.atmosphere.Atmos_t
            Atmosphere struct
        output : dict
            Output variables, as a dictionary
    """

    # Inform
    log.debug("Running AGNI...")
    time_str = "%d"%hf_row["Time"]

    # Solve atmosphere
    if bool(atmos.transparent):
        # no opacity
        log.info("Using transparent solver")
        atmos = _solve_transparent(atmos, config)

    else:
        # has opacity
        if config.atmos_clim.agni.solve_energy:
            log.info("Using nonlinear solver to conserve fluxes")
            atmos = _solve_energy(atmos, loops_total, dirs, config)
        else:
            log.info("Using prescribed temperature profile")
            atmos = _solve_once(atmos, config)

    # Write output data
    ncdf_path = os.path.join(dirs["output"],"data",time_str+"_atm.nc")
    jl.AGNI.save.write_ncdf(atmos, ncdf_path)

    # Make plots
    if multiple(loops_total, config.params.out.plot_mod):
        cff = os.path.join(dirs["output/plots"], f"plot_cff.{config.params.out.plot_fmt}")
        jl.AGNI.plotting.plot_contfunc1(atmos, cff)

    # ---------------------------
    # Calculate observables
    # ---------------------------

    # observed height and derived bulk density
    jl.AGNI.atmosphere.calc_observed_rho_b(atmos)
    rho_obs = float(atmos.transspec_rho)
    p_obs   = float(atmos.transspec_p) # set by peak of contribution function
    r_obs   = float(atmos.transspec_r)

    # ---------------------------
    # Parse results
    # ---------------------------

    log.debug("Parse results")
    tot_flux =      np.array(atmos.flux_tot)
    LW_flux_up =    np.array(atmos.flux_u_lw)
    SW_flux_up =    np.array(atmos.flux_u_sw)
    SW_flux_down =  np.array(atmos.flux_d_sw)
    T_surf =        float(atmos.tmp_surf)

    # New flux from SOCRATES
    F_atm_new = tot_flux[0]

    # Enforce positive limit on F_atm, if enabled
    if config.atmos_clim.prevent_warming:
        F_atm_lim = max( 1e-8 , F_atm_new )
    else:
        F_atm_lim = F_atm_new
    if not np.isclose(F_atm_lim , F_atm_new ):
        log.warning("Change in F_atm [W m-2] limited in this step!")
        log.warning("    %g  ->  %g" % (F_atm_new , F_atm_lim))

    log.info("    T_surf = %.3f K"%float(atmos.tmp_surf))
    log.info("    R_obs  = %.3f km"%float(r_obs/1e3))
    log.info("    F_top  = %.2e W m-2"%float(tot_flux[0]))
    log.info("    F_bot  = %.2e W m-2"%float(tot_flux[-1]))

    # XUV height in atm
    if config.escape.module == 'zephyrus':
        # escape level set by zephyrus config
        p_xuv = config.escape.zephyrus.Pxuv # [bar]
    else:
        # escape level set to surface
        p_xuv = hf_row["P_surf"] # [bar]
    p_xuv, r_xuv = get_radius_from_pressure(atmos.p, atmos.r, p_xuv*1e5) # [Pa], [m]

    # final things to store
    output = {}
    output["F_atm"]  = F_atm_lim
    output["F_olr"]  = LW_flux_up[0]
    output["F_sct"]  = SW_flux_up[0]
    output["T_surf"] = T_surf
    output["p_obs"]  = p_obs/1e5 # convert [Pa] to [bar]
    output["R_obs"]  = r_obs
    output["rho_obs"]= rho_obs
    output["albedo"] = SW_flux_up[0]/SW_flux_down[0]
    output["p_xuv"]  = p_xuv/1e5        # Closest pressure from Pxuv    [bars]
    output["R_xuv"]  = r_xuv            # Radius at Pxuv                [m]

    return atmos, output
