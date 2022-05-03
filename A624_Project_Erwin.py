'''
Project based on Bounded Martian satellite relative motion by Marcus and Gurfil

Focuses on the study of satellite relative motion around Mars
'''

import numpy as np
from scipy.integrate import odeint
from mpl_toolkits import mplot3d 
import matplotlib.pyplot as plt

'''
ideas:
-similar to orbit propagation, a,e,i,da,de,di as function into odeint and plot
'''
# INPUT PARAMETERS
mu = 0.042828e6          # universal gravitational parameter, km^3/s^2
# r_eq = 3396.2            # equatorial radius of Mars, km

# 2 MARS GRAVITATIONAL HARMONICS
# Table 1: Coefficients of the gravitational harmonics of Mars, normalized to J2
J2 = 1
J3 = 1.612e-2
J4 = -7.889e-3

# Table 2: Coefficients of the gravitational harmonics of Earth, normalized to J2 (use?)

# 3 THE ORBITAL EVOLUTION
# 3.1 Model Formulation
def mean_orbital_element_rates(y, t, mu):
	
    # grab inital/updated variables
    a = y[0]; e = y[1]; i = y[2]; O = y[3]; w = y[4]; M = y[5]

	# givens
    r_eq = 3396.2            	# equatorial radius of Mars, km
    # n = np.sqrt(mu/(a**3))      # 1/s, mean motion
    # M = M0 + n*(t-t0)
    # f = M + (2*e - (e**3)/4)*np.sin(M) + (5*(e**2) * np.sin(2*M))/4 + (13*(e**3) * np.sin(3*M))/12
    # r = a*(1 - e**2) / (1 + e*np.cos(f))

    # eqn 3
    '''R = - mu*r_eq**2 * (4*J2*r**2 * (3*np.sin(i)**2*np.sin(f+w)**2 - 1)  
        + 4*J3*r*r_eq*(5*np.sin(i)**2*np.sin(f+w)**2 - 3) *np.sin(i)*np.sin(f+w) 
        + J4*r_eq**2 *(35*np.sin(i)**4*np.sin(f+w)**4 - 30*np.sin(i)**2 *np.sin(f+w)**2 +3)) / 8*r**5'''
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
            6*a*r_eq*np.sin(w)*(4*np.sin(i) - 10*e*e*np.sin(i)**3 - 4*e*e*(1/np.sin(i)) - 5*np.sin(i)**3 + 15*e*e*np.sin(i))*J3 +
            15*e*r_eq**2*(28*e**2*np.sin(i)**2 -6*e**2*np.sin(w)**2 - 14*e**2*np.sin(i)**4 - 8 - 
            7*e**2*np.sin(i)**4*np.sin(w)**2 + 14*e**2*np.sin(i)**2*np.sin(w)**2 + 28*np.sin(i)**2 - 
            7*np.sin(i)**4*np.sin(w)**2 - 21*np.sin(i)**4 + 6*np.sin(i)**2*np.sin(w)**2 - 13*e*e)*J4) / (16*np.tan(i)*np.sqrt((1-e*e)*mu*a**(11)))

    # eqn 12 -- not sure which w_dot to use...
    '''w_dot = 3*mu*r_eq**2 * (2*J2*a**2 * e*(4-5*(np.sin(i))**2) - a*J3*r_eq*np.sin(i)*(4-5*(np.sin(i))**2) + 
                5*J4*e*r_eq**2 * (17*(np.sin(i))**2 - 14*(np.sin(i))**4 - 4) ) / (8*e*np.sqrt((1-e*e)*mu*a**(11)))'''
    # eqn 9
    M_dot = (1/(32*e))*np.sqrt(mu/a**(11))*(32*a**4*e + 24*a**2*e*r_eq**2*(2 + 8*e**2 - 12*e**2*np.sin(i)**2 - 3*np.sin(i)**2)*J2 +
             12*a*r_eq**3*np.sin(i)*np.sin(w)*(36*e**2 -45*e**2*np.sin(i)**2 + 5*np.sin(i)**2 - 4)*J3 +
             15*e*r_eq**4*(72*e**2*np.sin(i)**2*np.sin(w)**2 - 168*e**2*np.sin(i)**4 + 204*e**2*np.sin(i)**2 -
             84*e**2*np.sin(i)**4*np.sin(w)**2 + 14*np.sin(i)**4*np.sin(w)**2 - 7*np.sin(i)**4 - 48*e*e - 
             12*np.sin(i)**2*np.sin(w)**2 + 6*np.sin(i)**2)*J4)

    yDot = np.array([a_dot, e_dot, i_dot, O_dot, w_dot, M_dot])

    return yDot

# find position and velocity at each point
def position_velocity(r, theta, RAAN, ArgPeri, i, mu, h, e):

    '''takes r: position norm (km), theta: arument of perigee + true anomaly (rad), RAAN: (rad), i: inclination (rad), 
       mu: gravitational constant (km^3/s^2), h: angular momentum norm (km^2/s), and e: eccentricity. 
       Returns position and velocity '''

    # now, we can find radius and velocity (Battin Problem 3-21)
    rx = r*(np.cos(RAAN)*np.cos(theta) - np.sin(RAAN)*np.sin(theta)*np.cos(i))
    ry = r*(np.sin(RAAN)*np.cos(theta) + np.cos(RAAN)*np.sin(theta)*np.cos(i))
    rz = r*(np.sin(theta)*np.sin(i))

    vx = -(mu/h)*(np.cos(RAAN)*(np.sin(theta) + e*np.sin(ArgPeri)) + np.sin(RAAN)*(np.cos(theta) + e*np.cos(ArgPeri))*np.cos(i))
    vy = -(mu/h)*(np.sin(RAAN)*(np.sin(theta) + e*np.sin(ArgPeri)) - np.cos(RAAN)*(np.cos(theta) + e*np.cos(ArgPeri))*np.cos(i))
    vz =  (mu/h)*(np.cos(theta) + e*np.cos(ArgPeri)*np.sin(i))

    rvVector = np.array([rx, ry, rz, vx, vy, vz]) # r in km & v in km/s, respectively

    return rvVector

# propogate orbit using Keplerian parameters
def keplerian_propagation(t0, tf, mu, a, e, i, O, w, updated_elems, sol_array):

    '''takes t0: initial time (s), tf: final time (s), mu: gravitational constant (km^3/s^2), a: semi-major axis (km),
       e: eccentricity, i: inclination (rad), w: argument of perigee (rad), RAAN: (rad), sol_array: empty array. 
       Returns keplerian propagation in sol_array.'''

    for t in np.arange(t0, tf):

        a, e, i, O, w, M = updated_elems[t]   # grab orbital elements at current time

        # https://en.wikipedia.org/wiki/True_anomaly
        f = M + (2*e - (e**3)/4)*np.sin(M) + (5*(e**2) * np.sin(2*M))/4 + (13*(e**3) * np.sin(3*M))/12

        # now we can find theta and radius
        theta = w + f                 # angle, rad
        p = a*(1 - e**2)              # semi-latus rectum, km
        r = p/(1 + e*np.cos(f))       # radius, km
        h = np.sqrt(p*mu)             # angular momentum, km^2/s

        # find updated position and velocity, add to solution array
        temp = position_velocity(r, theta, O, w, i, mu, h, e)
        sol_array.append(temp)
            
    return sol_array

# given initial mean orbital elements 
a0_mean = 3897e3				# m, initial semi-major axis
e0_mean = 0.00546				# initial eccentricity
i0_mean = np.deg2rad(50)		# rad, initial inclination
# given initial oscullating elements
a0_osc = 3891815				# m, initial semi-major axis
e0_osc = 0.00462				# initial eccentricity
i0_osc = np.deg2rad(49.968)		# rad, initial inclination
# given initial elements (same for mean and osc.)
O0 = np.deg2rad(120)			# rad, initial RAAN
w0 = np.deg2rad(270)			# rad, initial argument of perigee
f0 = np.deg2rad(0)				# rad, initial true anomaly
# initial mean anomalies
M0_mean = np.arctan2((np.sqrt(1-e0_mean*e0_mean*np.sin(f0))/(1+e0_mean*np.cos(f0))), (e0_mean+np.cos(f0)/(1+e0_mean*np.cos(f0)))) - e0_mean*(np.sqrt(1-e0_mean*e0_mean*np.sin(f0))/(1+e0_mean*np.cos(f0)))
M0_osc = np.arctan2((np.sqrt(1-e0_osc*e0_osc*np.sin(f0))/(1+e0_osc*np.cos(f0))), (e0_osc+np.cos(f0)/(1+e0_osc*np.cos(f0)))) - e0_osc*(np.sqrt(1-e0_osc*e0_osc*np.sin(f0))/(1+e0_osc*np.cos(f0)))

# update elements over time
t0 = 0                                       # s, initial time
tf = 365*24*3600                             # s, final time: one year
t = np.linspace(t0, tf, tf)                  # s, integrate across 8 orbital periods
# mean initial elements
y0_mean = np.array([a0_mean, e0_mean, i0_mean, O0, w0, M0_mean])            # initial orital elements, e unchanged
y0_mean_div5 = np.array([a0_mean, (e0_mean/5), i0_mean, O0, w0, M0_mean])   # initial orital elements, e/5
y0_mean_x5 = np.array([a0_mean, (5*e0_mean), i0_mean, O0, w0, M0_mean])     # initial orital elements, 5*e
# oscullating initial elements
y0_osc = np.array([a0_osc, e0_osc, i0_osc, O0, w0, M0_osc])                 # initial orital elements, e unchanged
y0_osc_div5 = np.array([a0_osc, (e0_osc/5), i0_osc, O0, w0, M0_osc])        # initial orital elements, e/5
y0_osc_x5 = np.array([a0_osc, (5*e0_osc), i0_osc, O0, w0, M0_osc])          # initial orital elements, 5*e 
# updated mean elements
updated_elems_mean = odeint(mean_orbital_element_rates, y0_mean, t, args=(mu, ), rtol=1e-6, atol=1e-6)
updated_elems_mean_div5 = odeint(mean_orbital_element_rates, y0_mean_div5, t, args=(mu, ), rtol=1e-6, atol=1e-6)
updated_elems_mean_x5 = odeint(mean_orbital_element_rates, y0_mean_x5, t, args=(mu, ), rtol=1e-6, atol=1e-6)
# updated oscullating elements
updated_elems_osc = odeint(mean_orbital_element_rates, y0_osc, t, args=(mu, ), rtol=1e-6, atol=1e-6)
updated_elems_osc_div5 = odeint(mean_orbital_element_rates, y0_osc_div5, t, args=(mu, ), rtol=1e-6, atol=1e-6)
updated_elems_osc_x5 = odeint(mean_orbital_element_rates, y0_osc_x5, t, args=(mu, ), rtol=1e-6, atol=1e-6)

# plot Figure 1 (pg. 28): Argument of Periapsis
fig = plt.figure('Figure 1')
plt.plot(t, np.rad2deg(updated_elems_mean[:,4]), label='$w_{mean} \;when\; e_{mean,t=0} = e_{frozen}$') 
plt.plot(t, np.rad2deg(updated_elems_mean_div5[:,4]), label='$w_{mean} \;when\; e_{mean,t=0} = e_{frozen}/5$') 
plt.plot(t, np.rad2deg(updated_elems_mean_x5[:,4]), label='$w_{mean} \;when\; e_{mean,t=0} = e_{frozen}*5$') 
plt.plot(t, np.rad2deg(updated_elems_osc[:,4]), label='$w_{osc} \;when\; e_{mean,t=0} = e_{frozen}$') 
plt.plot(t, np.rad2deg(updated_elems_osc_div5[:,4]), label='$w_{osc} \;when\; e_{mean,t=0} = e_{frozen}/5$') 
plt.plot(t, np.rad2deg(updated_elems_osc_x5[:,4]), label='$w_{osc} \;when\; e_{mean,t=0} = e_{frozen}*5$') 
plt.xlabel('Time (seconds)')
plt.ylabel('Argument of Periapsis (deg)')
plt.legend(loc='lower center')
plt.show()

# plot Figure 2 (pg. 28): Eccentricity
fig = plt.figure('Figure 2')
plt.plot(t, updated_elems_mean[:,1], label='$w_{mean} \; when\; e_{mean,t=0} = e_{frozen}$') 
plt.plot(t, updated_elems_mean_div5[:,1], label='$w_{mean} \; when\; e_{mean,t=0} = e_{frozen}/5$') 
plt.plot(t, updated_elems_mean_x5[:,1], label='$w_{mean} \;when\; e_{mean,t=0} = e_{frozen}*5$') 
plt.plot(t, updated_elems_osc[:,1], label='$w_{osc} \;when\; e_{mean,t=0} = e_{frozen}$') 
plt.plot(t, updated_elems_osc_div5[:,1], label='$w_{osc} \;when\; e_{mean,t=0} = e_{frozen}/5$') 
plt.plot(t, updated_elems_osc_x5[:,1], label='$w_{osc} \;when\; e_{mean,t=0} = e_{frozen}*5$') 
plt.xlabel('Time (seconds)')
plt.ylabel('Eccentricity')
plt.legend(loc='lower center')
plt.show()


'''
# propogate orbits
emp_sol_mean = [] # empty array for solution
kep_sol_mean = keplerian_propagation(t0, tf, mu, a0_mean, e0_mean, i0_mean, O0, w0, updated_elems_mean, emp_sol_mean)
kep_sol_mean = np.reshape(kep_sol_mean, (len(kep_sol_mean), 6))

# plot orbital position
fig = plt.figure('Orbital Propagation')
ax = plt.axes(projection='3d')
ax.plot3D(kep_sol_mean[:, 0], kep_sol_mean[:,1], kep_sol_mean[:,2], label='Mean')

ax.set_xlabel('rx (m)')
ax.set_ylabel('ry (m)')
ax.set_zlabel('rz (m)')
ax.set_title('Orbital Propagation')
plt.legend()
plt.show()'''



# 4.1 ANALYSIS


# FROM APPENDIX
ef = 'frozen ecentricity'
A_RAAN = 84*(3*ef**2 + 2)*a*a*np.sin(i)*J2 - 54*(15*(np.sin(i))**2 + 4)*a*ef*r_eq*J3 + 165*(42*ef*ef*(np.sin(i))**2 - 23*ef*ef + 7*(np.sin(i))**2 - 4)*r_eq**2 * np.sin(i)*J4
E_RAAN = 12*(3*ef**2 - 8)*a*a*ef*np.sin(i)*J2 + 6*(15*(np.sin(i))**2 - 4)*a*r_eq*J3 + (636*ef**2*(np.sin(i))**2 - 345*ef**2 - 16*(np.sin(i))**2 + 750)*ef*r_eq**2*np.sin(i)*J4
