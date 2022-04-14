'''
Project based on Bounded Martian satellite relative motion by Marcus and Gurfil

Focuses on the study of satellite relative motion around Mars
'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

'''
ideas:
-similar to orbit propagation, a,e,i,da,de,di as function into odeint and plot
'''
# INPUT PARAMETERS
# mu = 0.042828e6          # universal gravitational parameter, km^3/s^2
# r_eq = 3396.2            # equatorial radius of Mars, km

# 2 MARS GRAVITATIONAL HARMONICS
# Table 1: Coefficients of the gravitational harmonics of Mars, normalized to J2
J2 = 1
J3 = 1.612e-2
J4 = -7.889e-3

# Table 2: Coefficients of the gravitational harmonics of Earth, normalized to J2 (use?)

# 3 THE ORBITAL EVOLUTION
# 3.1 Model Formulation
def func(a, e, i, O, w, f):

    # parameters + r eqn
    mu = 0.042828e6          # universal gravitational parameter, km^3/s^2
    r_eq = 3396.2            # equatorial radius of Mars, km
    r = a*(1 - e**2) / (1 + e*np.cos(f))
    # eqn 3
    R = - mu*r_eq**2 * (4*J2*r**2 * (3*np.sin(i)**2*np.sin(f+w)**2 - 1)  
        + 4*J3*r*r_eq*(5*np.sin(i)**2*np.sin(f+w)**2 - 3) *np.sin(i)*np.sin(f+w) 
        + J4*r_eq**2 *(35*np.sin(i)**4*np.sin(f+w)**4 - 30*np.sin(i)**2 *np.sin(f+w)**2 +3)) / 8*r**5
    # eqn 4
    a_dot = 0
    # eqn 5
    e_dot = r_eq**3*np.sin(i)*np.cos(w) * np.sqrt((1-e*e)*mu/a**(11)) * (6*a*J3*(5*np.sin(i)**2 - 4) + 
            15*e*r_eq*J4*np.sin(i)*np.sin(w)*(7*np.sin(i)**2 - 6)) / 16
    # eqn 6
    i_dot = -e*mu*r_eq**3*np.cos(i)*np.cos(w) * (6*a*J3*(1 - 5*np.cos(i)**2 + 15*e*r_eq*J4*np.sin(i)
            *np.sin(w)*(7*np.sin(i)**2 - 6))) / (16*np.sqrt((1-e*e)*mu*a**(11)))
    # eqn 7
    O_dot = -mu*r_eq**2 * (12*a**2*np.sin(i)*(3*e**2 + 2)*J2 + 6*a*e*r_eq*np.sin(w)*(15*np.sin(i)**2 - 4)*J3 +
            15*r_eq**2*np.sin(i)*(14*e**2*np.sin(i)**2*np.sin(w)**2 - 4 - 17*e**2 + 7*np.sin(i)**2 
            - 6*e**2*np.sin(w)**2 + 28*e**2*np.sin(i)**2)*J4) / (16*np.tan(i)*np.sqrt((1-e*e)*mu*a**(11)))
    # eqn 8
    w_dot = mu*r_eq**2 * (12*(e*e - 5*np.sin(i**2)+4)*a**2*e*J2 + 
            6*a*r_eq*np.sin(w)*(4*np.sin(i) - 10*e*e*np.sin(i)**3 - 4*e*e*np.csc(i) - 5*np.sin(i)**3 + 15*e*e*np.sin(i))*J3 +
            15*e*r_eq**2*(28*e**2*np.sin(i)**2 -6*e**2*np.sin(w)**2 - 14*e**2*np.sin(i)**4 - 8 - 
            7*e**2*np.sin(i)**4*np.sin(w)**2 + 14*e**2*np.sin(i)**2*np.sin(w)**2 + 28*np.sin(i)**2 - 
            7*np.sin(i)**4*np.sin(w)**2 - 21*np.sin(i)**4 + 6*np.sin(i)**2*np.sin(w)**2 - 13*e*e)*J4) / (16*np.tan(i)*np.sqrt((1-e*e)*mu*a**(11)))
    # eqn 9
    M_dot = (1/(32*e))*np.sqrt(mu/a**(11))*(32*a**4*e + 24*a**2*e*r_eq**2*(2 + 8*e**2 - 12*e**2*np.sin(i)**2 - 3*np.sin(i)**2)*J2 +
             12*a*r_eq**3*np.sin(i)*np.sin(w)*(36*e**2 -45*e**2*np.sin(i)**2 + 5*np.sin(i)**2 - 4)*J3 +
             15*e*r_eq**4*(72*e**2*np.sin(i)**2*np.sin(w)**2 - 168*e**2*np.sin(i)**4 + 204*e**2*np.sin(i)**2 -
             84*e**2*np.sin(i)**4*np.sin(w)**2 + 14*np.sin(i)**4*np.sin(w)**2 - 7*np.sin(i)**4 - 48*e*e - 
             12*np.sin(i)**2*np.sin(w)**2 + 6*np.sin(i)**2)*J4)

    return R, a_dot, e_dot, i_dot, O_dot, w_dot, M_dot
