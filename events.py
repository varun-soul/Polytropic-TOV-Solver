import numpy as np

#############################
## DEFINING SURFACE EVENTS ##
#############################

'''
This module defines the surface events for implementations
of the TOV Solver (direct & log)
'''

def surface_event_direct(P_floor):
    def event(r,y, *args):
        return y[0]-P_floor
    
    event.terminal = True
    event.direction = -1
    return event

def surface_event_log(P_floor):
    def event(r, y, *args):
        return y[0] - np.log(P_floor)
    
    event.terminal = True
    event.direction = -1
    return event

def surface_event_enthalpy():
    def event(r, y, *args):
        return y[0] - 1.0      # fires when h == 1
    
    event.terminal  = True
    event.direction = -1
    return event

