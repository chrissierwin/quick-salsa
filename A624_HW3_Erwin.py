# A624 HW 3, due April 21, 2022

import numpy as np
from scipy.integrate import odeint
from mpl_toolkits import mplot3d 
import matplotlib.pyplot as plt

# Consider the perturbed two body problem, where the rates of the osculating orbital
# elements are governed by the following Gauss' variational equations, given by:
def Gauss_Vari(y, t, t0):

    # grab inital/updated variables
    O = y[0]; i = y[1]; a = y[2]; w = y[3]; e = y[4]; f = y[5]

    # constants / givens
    mu = 3.986e14                         # m^3/s^2, Earth's gravitational parameter
    p = a*(1 - e*e)                       # m, semi-latus rectum
    r = p/(1 + e*np.cos(f))               # m, radius
    n = np.sqrt(mu/(a**3))                # s^-1, mean motion
    eta = np.sqrt(1 - e*e)                # unitless parameter
    b = a*eta                             # m, semi-minor axis
    theta = w + f                         # rad
    r_eq = 6378.1370e3                    # m, e mean equatorial radius of the Earth
    J2 = 1082.63e-6                       # variation in Earth's Oblateness

    # Gauss' variational equations
    dOdt = -3*J2*n*(a*a/(b*r))*(r_eq/r)**2 * np.sin(theta)**2 * np.cos(i)
    didt = -3*J2/4 *n*(a*a/(b*r))*(r_eq/r)**2 * np.sin(2*theta) * np.sin(2*i)
    dwdt =  3*J2/2 *n*(p/(r*r*e*eta**3))*(r_eq/r)**2 * (2*r*e*np.cos(i)**2*np.sin(theta)**2 - (p+r)*np.sin(f)*np.sin(i)**2*np.sin(2*theta) + p*np.cos(f)*(1-3*np.sin(i)**2*np.sin(theta)**2))
    dadt = -3*J2*n*(a**4/(b*r*r))*(r_eq/r)**2 * (e*np.sin(f)*(1-3*np.sin(theta)**2*np.sin(i)) + (p/r)*np.sin(2*theta)*np.sin(i)**2)
    dedt = -3*J2/2 *n*(a*a/(b*r))*(r_eq/r)**2 * ((p/r)*np.sin(f)*(1-3*np.sin(theta)**2*np.sin(i)**2) + (e+np.cos(f)*(2+e*np.cos(f)))*np.sin(2*theta)*np.sin(i)**2)
    dM0dt = 3*J2/2 *n*(p/(r*r*e*eta**2))*(r_eq/r)**2 * ((p+r)*np.sin(f)*np.sin(i)**2*np.sin(2*theta) + (2*r*e - p*np.cos(f))*(1-3*np.sin(i)**2*np.sin(theta)**2)) + (1.5*n/a)*(t-t0)*dadt

    yDot = np.append(dOdt, didt, dwdt, dadt, dedt, dM0dt)

    return yDot

# find mean and eccentric anomalies
def meanAnomaly(t0, t, n):
    '''takes t0: time at periapse (s), t: current time (s), and n: mean motion (s). returns M: mean anomaly (rad)'''
    return n*(t - t0) 

def func(M, E, e): 
    '''takes M: mean anomaly (rad), E: eccentric anomaly (rad), e: eccentricity'''	
    return (E - e*np.sin(E) - M)

def derivFunc(E, e): 
    '''takes E: eccentric anomaly (rad), e: eccentricity'''
    return (1 - e*np.cos(E))

# function to find the root 
def newtonRaphson(M, E, e): 
    '''uses the root finding method, Newton Raphson, to return the eccentric anomaly (rad).
       takes M: mean anomaly (rad), E: initial guess of eccentric anomaly (rad), and e: eccentricity'''
    
    h = func(M, E, e) / derivFunc(E, e) 
    
	# convergence criteria 
    while abs(h) >= 0.000001: 
        h = func(M, E, e)/derivFunc(E, e) 
		
		# E(k+1) = E(k) - f(E) / f'(E) 
        E = E - h 

    return E

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

    rVector = np.array([rx, ry, rz]) # km
    vVector = np.array([vx, vy, vz]) # km/s
    rvVector = np.array([rx, ry, rz, vx, vy, vz])

    return rvVector

# prpogate orbit using Keplerian parameters
def keplerian_propagation(t0, tf, mu, a, e, i, w, RAAN, sol_array):

    '''takes t0: initial time (s), tf: final time (s), mu: gravitational constant (km^3/s^2), a: semi-major axis (km),
       e: eccentricity, i: inclination (rad), w: argument of perigee (rad), RAAN: (rad), sol_array: empty array. 
       Returns keplerian propagation in sol_array.'''

    n = np.sqrt(mu/(a**3))                # sec
    E0 = meanAnomaly(t0, t0, n) + 0.25*e  # initial guess, rad

    for t in np.arange(t0, tf):
        M = meanAnomaly(t0, t, n)                     # mean anomaly, rad
        E = newtonRaphson(M, E0, e)                   # eccentric anomaly, rad

        if E >=0 and E <= (2*np.pi):

            # v = np.arccos((np.cos(E) - e) / (1 - e*np.cos(E))) # true anomaly, rad
            # https://en.wikipedia.org/wiki/True_anomaly
            v = M + (2*e - (e**3)/4)*np.sin(M) + (5*(e**2) * np.sin(2*M))/4 + (13*(e**3) * np.sin(3*M))/12
            #v_all.append(v)

            # now we can find theta and radius

            p = a*(1 - e**2)           # semi-latus rectum, km
            h = np.sqrt(p*mu)          # angular momentum, km^2/s

            # now we can find theta and radius
            theta = w + v            # angle, rad
            r = p/(1 + e*np.cos(v))  # radius, km

            # re-run with new E
            E0 = E      

            temp = position_velocity(r, theta, RAAN, w, i, mu, h, e)
            sol_array.append(temp)
            
    return sol_array


# PROBLEM 1 =======================================
'''Propagate an orbit using the Gauss' equations for the J2 perturbed two body problem, 
whose initial osculating elements are given below for 8 orbits starting from the 
initial true anomaly. Plot the time history of the osculating elements. Calculate a 
straight line fit of the orbital elements as a function of time and plot the straight 
line. For each straight line fit, assess the std deviation of the fit and report its values.'''

# orbital parameters 
mu = 3.986004415e14          # universal gravitational parameter, m^3/s^2
rEarth = 6378.135e3          # radius of Earth, m

# givens
O0 = np.deg2rad(45)  # rad, initial RAAN
i0 = np.deg2rad(30)  # rad, initial inclination
a0 = 7378e3          # m, initial semi-major axis
w0 = np.deg2rad(270) # rad, initial arg. of perigee
e0 = 0.01            # initial eccentricity
f0 = np.pi/10        # rad, initial true anomaly

# update elements over time
P = 2*np.pi*np.sqrt(a0**3/3.986e14)
y0 = np.append(O0, i0, a0, w0, e0, f0)      # initial orital elements
t0 = 0.                                     # s, initial time
t = np.linspace(t0, 8*P, 8*P)               # s, integrate across 8 orbital periods
updated_elems = odeint(Gauss_Vari, y0, t, args = (t0,))

# propogate orbits
emp_sol = [] # empty array for solution
kep_sol = keplerian_propagation(0, P, mu, a0, e0, i0, w0, O0, emp_sol)
kep_sol = np.reshape(kep_sol, (len(kep_sol), 6))

# plot orbital position
fig = plt.figure('Orbital Propagation')
ax = plt.axes(projection='3d')
ax.plot3D(kep_sol[:, 0], kep_sol[:,1], kep_sol[:,2], label='Orbit')

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = rEarth*np.cos(u)*np.sin(v)
y = rEarth*np.sin(u)*np.sin(v)
z = rEarth*np.cos(v)
ax.plot_wireframe(x, y, z, color="g", label='Earth')

ax.set_xlabel('rx (m)')
ax.set_ylabel('ry (m)')
ax.set_zlabel('rz (m)')
ax.set_title('Orbital Propagation')
plt.legend()
plt.show()


plt.plot(t, t)
