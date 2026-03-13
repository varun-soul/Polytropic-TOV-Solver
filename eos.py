import numpy as np

###############################
## SETTING UP POLYTROPIC EOS ##
###############################

'''
Polytropic Equation of State:
P = (K)*(rho)^(Gamma)
'''

def pressure_from_density(rho, K, Gamma):
    return (K)*(rho)**(Gamma)

def density_from_pressure(P, K, Gamma):
    return (P/K)**(1.0/Gamma)

def energy_density(P, K, Gamma):
    rho = density_from_pressure(P, K, Gamma)
    return rho + P/(Gamma - 1.0)

def enthalpy_from_pressure(P, K, Gamma):
    rho = density_from_pressure(P, K, Gamma)
    return 1 + ((Gamma)/(Gamma - 1))*(P/rho)

def pressure_from_enthalpy(h, K, Gamma):
    return ((Gamma - 1.0) / Gamma * (h - 1.0)) ** (Gamma / (Gamma - 1.0)) \
           / K ** (1.0 / (Gamma - 1.0))    
