import numpy as np
import eos

#############################
## SETTING UP THE TOV EQNS ##
#############################

'''
We will be working with two formulations of the TOV equations.
(1) In the first, we shall be considering a direct implementation of
    the integration with P. 

(2) In the second, we shall be considering evolving with log(P) in rder 
    smoothly integrate near the surface. This will prevent our solver
    from taking large steps which might have strange/unphysical results.
'''

# (An enthalpy formulation may also be considered)

def TOV_RHS_direct(r, y, K, Gamma):

    # Here, y = [P, m]
    # Function returns: dP/dr, dm/dr

    P, m = y

    if P <= 0.0:
        return [0.0, 0.0]
    
    eps = eos.energy_density(P, K, Gamma)

    denom = r*(r-2*m)
    dPdr = -(eps + P)*(m + (4*np.pi*r**3)*P)/denom
    dmdr = (4*np.pi*r**2)*eps

    return [dPdr, dmdr]

def TOV_RHS_log(r, y, K, Gamma):

    # Here, y = [q, m], where q = log(P)
    # Function returns [dq/dr, dm/dr]

    q, m = y

    P = np.exp(q)
    eps = eos.energy_density(P, K, Gamma)

    denom = r*(r-2*m)
    dPdr = -(eps + P)*(m + (4*np.pi*r**3)*P)/denom
    dqdr = dPdr/P
    dmdr = (4*np.pi*r**2)*eps

    return [dqdr, dmdr]

def TOV_RHS_enthalpy(r, y, K, Gamma):

    # Here, y = [h, m], where h = Enthalpy
    # Function returns [dhdr, dmdr]

    h, m = y

    if h <= 1.0:
        return [0.0, 0.0]

    P = eos.pressure_from_enthalpy(h, K, Gamma)
    eps = eos.energy_density(P, K, Gamma)

    denom = r*(r-2*m)
    dhdr = -h*(m + 4*np.pi*(r**3)*P)/(denom)
    dmdr = 4*(np.pi)*(r**2)*(eps)

    return [dhdr, dmdr]
