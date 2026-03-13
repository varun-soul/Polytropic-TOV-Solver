import numpy as np
from scipy.integrate import solve_ivp
import eos
import events
import equations as eqns
import units

'''
    This is the actual integrator. It will integrate the TOV
    equations from r~0.

    > PARAMETERS <

    1. rho_c = Central Pressure
    2. K = Polytropic Constant
    3. Gamma = Polytropic Index
    4. use_log = Flag for using log-pressure
    5. P_floor = Minimum pressure for termination
    6. r_max = Maximum Radius for integration

    The module uses rho_c, which gets converted in P_c internally.

    > RETURNS <

    1. Radius of NS (R)
    2. Mass of NS (M)
    3. Full sol object
'''

################
## INTEGRATOR ##
################

r_start = 1e-6
r_max = 50.0

def solver(rho_c, K, Gamma, method = 'direct', P_floor = 1e-10,
           r_max = r_max, r_start = r_start):


    r_start = 1e-6 # Set to << 1 to avoid r = 0 singularity.

    ########################
    ## INITIAL CONDITIONS ##
    ########################

    P_c = eos.pressure_from_density(rho_c, K, Gamma)
    eps_c = eos.energy_density(P_c, K, Gamma)

    # Near r = 0, we use Taylor's Expansion for initial mass.
    m_start = (4.0/3.0)*(np.pi)*(r_start**3)*eps_c

    ###############
    ## INTEGRATE ##
    ###############
      
    if method == 'direct':

        y0 = [P_c, m_start]
        RHS = eqns.TOV_RHS_direct
        event = events.surface_event_direct(P_floor)

        sol = solve_ivp(
            fun      = RHS,
            t_span   = (r_start, r_max),
            y0       = y0,
            method   = 'DOP853',
            events   = event,
            args     = (K, Gamma),
            rtol     = 1e-12,
            atol     = 1e-14,
            max_step = 1e-4,        # prevents stepping over the surface
            dense_output = True,
        )

    elif method == 'log':

        y0 = [np.log(P_c), m_start]
        RHS = eqns.TOV_RHS_log
        event = events.surface_event_log(P_floor)

        sol = solve_ivp(
            fun      = RHS,
            t_span   = (r_start, r_max),
            y0       = y0,
            method   = 'DOP853',
            events   = event,
            args     = (K, Gamma),
            rtol     = 1e-12,
            atol     = 1e-14,
            max_step = 1e-4,        # prevents stepping over the surface
            dense_output = True,
        )
    
    elif method == 'enthalpy':

        h_c = eos.enthalpy_from_pressure(P_c, K, Gamma)
        y0 = [h_c, m_start]
        RHS = eqns.TOV_RHS_enthalpy
        event = events.surface_event_enthalpy()

        # print(f"  h_c     = {h_c:.6f}")
        # print(f"  P_c     = {P_c:.6e}")
        # print(f"  m_start = {m_start:.6e}")

        sol = solve_ivp(
            fun      = RHS,
            t_span   = (r_start, r_max),
            y0       = y0,
            method   = 'DOP853',
            events   = event,
            args     = (K, Gamma),
            rtol     = 1e-12,
            atol     = 1e-14,
            max_step = 1e-4,        # prevents stepping over the surface
            dense_output = True,
        )


    ################
    ## EXTRACTION ##
    ################

    if sol.t_events[0].size > 0:
        # Terminated cleanly at surface event
        R = sol.t_events[0][0]
        M = sol.y_events[0][0][1]

    else:
        # Hit r_max without finding surface — star may be too large
        R = sol.t[-1]
        M = sol.y[1, -1]
    
        print(f"  Warning: surface not found before r_max={r_max}")

    return R, M, sol   
