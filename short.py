import numpy as np
import eos
import units

#######################################
## CONSTRUCTING SHORT FILES FOR GR1D ##
#######################################

import directory
import os 
import solver
import numpy as np

'''
This module contains the code for construction of initial stellar
profiles which will be fed to the GR1D code.
'''
r0 = solver.r_start

def create(sol, K, Gamma, name, flag):

    filename = os.path.join(directory.prof_dir, name)

    zones = int(input("  Enter number of required zones: "))

    r_max = sol.t_events[0][0]
    r_uniform = np.linspace(r0, r_max, zones)
    y_uniform = sol.sol(r_uniform)
    
    if flag == 1:
        r_array = r_uniform*units.LENGTH_CGS
        M_array = y_uniform[1]*units.M_CGS
        P_array = y_uniform[0]
        rho_array = eos.density_from_pressure(P_array, K, Gamma)*units.RHO_CGS

    elif flag == 2:
        r_array = r_uniform*units.LENGTH_CGS
        M_array = y_uniform[1]*units.M_CGS
        P_array = y_uniform[0]
        P_array = np.exp(P_array)
        rho_array = eos.density_from_pressure(P_array, K, Gamma)*units.RHO_CGS  

    elif flag == 3:
        r_array = r_uniform*units.LENGTH_CGS
        M_array = y_uniform[1]*units.M_CGS
        P_array = eos.pressure_from_enthalpy(y_uniform[0], K, Gamma)
        rho_array = eos.density_from_pressure(P_array, K, Gamma)*units.RHO_CGS 

    padding = str(input("  Do you want Atmosphere Padding? (y/n)"))

    if padding == "y":

        atm_size = float(input("  Extend atmosphere up till (CGS): "))
        dr = int(input("  Enter the atmospheric grid spacing: "))
        atm_rho = float(input("  Enter atmospheric density: "))

        atm_array = np.arange(r_array[-1] + dr, atm_size, dr)
        atm_rho_array = np.full(len(atm_array), atm_rho)

        M_surface = M_array[-1]
        atm_mass_array = M_surface + (4.0*np.pi/3.0)*(atm_rho)*(atm_array - r_array[-1])

        r_array = np.concatenate([r_array, atm_array])
        M_array = np.concatenate([M_array, atm_mass_array])
        rho_array = np.concatenate([rho_array, atm_rho_array])

    N = len(r_array)

    with open(filename, 'w') as f:
        f.write(f"{N-1}\n")

        for i in range(N-1):
            line = f"{i + 1:<6}    {M_array[i]:<24.18e}    {r_array[i]:<24.18e}    {0.0:<8.1f}   {rho_array[i]:<24.18e}  {0.0:<8.1f}   {0.0:<8.1f}   {0.0:<8.1f}\n"
            f.write(line)






    

    






