from scipy import constants

'''
This file contains the true values for constants in CGS units.
It also contains the conversion factors required for converting
from GEOMETRIC UNITS to CGS UNITS.
'''

##############################
### SETTING UP UNITS (CGS) ###
##############################

c = constants.c*100 # cgs
G = constants.G*1000 # cgs
M_sun = (1.988416)*(10**33)

####################################
### CONVERSIONS FROM GEOM TO CGS ###
####################################

LENGTH_CGS = (G*M_sun)/(c**2)
RHO_CGS = (M_sun)/(LENGTH_CGS**3)
PRESSURE_CGS = (RHO_CGS)*(c**2)
M_CGS = M_sun
K_CGS = (c**2)/(RHO_CGS)

print(1.455e5/K_CGS)
