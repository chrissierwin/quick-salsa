# -*- coding: utf-8 -*-
"""
Created Mar 2022

@author: Chrissi Erwin

Exam 1, Part B: Machine Problem
"""

import numpy as np
from scipy.integrate import odeint
from mpl_toolkits import mplot3d 
import matplotlib.pyplot as plt

# orbital parameters/givens =======
mu = 398600.4415          # universal gravitational parameter, km^3/s^2
rEarth = 6378.135         # radius of Earth, km

r0_vec = np.array([7800, 5400, 1800])    # initial orbit, km
v0_vec = np.array([9.2, 7.0107, 2.3042]) # initial velocity, km/s

# generate useful functions ======
def calc_elems(r_vec, v_vec):

    '''Take position and velocity vectors and returns norms/orbital elements '''

    # begin w/ vector norms
    r = np.linalg.norm(r_vec) # km
    v = np.linalg.norm(v_vec) # km/s

    # use the vis-viva eqn to find semi-major axis
    a = (2/r - (v**2 / mu))**(-1) # km

    # calculate the angular momentum
    h_vec = np.cross(r_vec,v_vec) # km^2 / s
    h = np.linalg.norm(h_vec)     # km^2 / s

    # use eqn for semi-latus rectum to find eccentricity
    e = np.sqrt(1 - (h**2)/(mu*a)) # unit-less

    return r, v, h, a, e

def orbit_shape(e):
    '''Takes eccentricity and describes orbital shape'''

    if e > 0 and e < 1:
        print('The orbit is elliptical')
    elif e == 1:
        print('The orbit is parabolic')
    elif e > 1:
        print('The orbit is hyperbolic')


def orbital_angles(r_vec, v_vec, a, mu):

    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)

    # inclination calc
    i = np.arccos(h_vec[2]/h)
    
    # RAAN calcs
    n_vec = np.array([-h_vec[1], h_vec[0], 0])
    n = np.linalg.norm(n_vec)
    if n_vec[1] >= 0:
        Omega = np.arccos(n_vec[0] / n)
    elif n_vec[1] < 0:
        Omega = 2*np.pi - np.arccos(n_vec[0] / n)

    # Arg. Peri. calcs
    e_vec = (np.cross(v_vec,h_vec) - mu*r_vec/r) / mu 
    e = np.linalg.norm(e_vec) # check for matching values: correct
    if e_vec[2] >= 0:
        w = np.arccos( np.dot(n_vec, e_vec) / (n * e) ) 
    elif e_vec[2] < 0:
        w = 2*np.pi - np.arccos( np.dot(n_vec, e_vec) / (n * e) ) 

    # true anomaly calcs
    rv = np.dot(r_vec, v_vec)
    if rv >= 0:
        nu = np.arccos( np.dot(e_vec, r_vec) / (e * r) ) 
    elif rv < 0:
        nu = 2*np.pi - np.arccos( np.dot(e_vec, r_vec) / (e * r) )  

    # period calc
    P = 2*np.pi*np.sqrt((a**3)/mu)

    # print orbital elements
    print('Semi-major axis (a, km): \t', a)
    print('Eccentricity (e): \t\t', e)
    print('Inclination (i, deg): \t\t', np.rad2deg(i))
    print('RAAN (Omega, deg): \t\t', np.rad2deg(Omega))
    print('Argument of perigee (w, deg): \t', np.rad2deg(w))
    print('True anomaly (v, rad): \t\t', np.rad2deg(nu))
    print('Orbital period (P, s): \t\t', P, '\n') # period

    return a, e, i, Omega, w, nu, P

def twoBody(y, t, mu):

    '''return 3D vector of velocity and accel. for a two body problem'''

    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)
    
    yDot = np.empty((6,))
    
    yDot[0] = y[3]
    yDot[1] = y[4]
    yDot[2] = y[5]
    yDot[3] = (-mu/(r**3))*y[0]
    yDot[4] = (-mu/(r**3))*y[1]
    yDot[5] = (-mu/(r**3))*y[2]
    
    return yDot


# part 1 ==========================
# calculate eccentricity to find orbit shape using orbital elements

# orbital elements, shape, and print values
r0, v0, h0, a0, e0 = calc_elems(r0_vec,v0_vec)
orbit_shape(e0)
a0, e0, i0, Omega0, w0, nu0, P0 = orbital_angles(r0_vec, v0_vec, a0, mu)
print('Since hyperbolic orbits are used for escape trajectories (and that they do '
      'not close) the period does not exist')

# part 2 ==========================
# pick an orbit that is close to the initial conditions to create an elliptical orbit
print('\nNow, reduce the initial velocity in order to create an elliptical orbit.\n')

# new vectors
rn_vec = np.array([15000, 14000, 1800])    # increased position, km 
vn_vec = np.array([0.01, 4.30299270602785441, 0.03]) # reduced velocity, km/s

# new orbital elements, shape, and print values
rn, vn, hn, an, en = calc_elems(rn_vec,vn_vec)
orbit_shape(en)
an, en, inn, Omegan, wn, nun, Pn = orbital_angles(rn_vec, vn_vec, an, mu)

print('Assumptions: the velocity can be reduced, and position can be increased in order to create '
      'an elliptical orbit that does not pass through the Earth. Spherical Earth is also assumed. ')

# part 3 ==========================
# propagate the orbit you found in item 2 for a period of 5 orbital time periods

tol = 1e-6
print('\nODE propagation tolerance: \t', tol)

# integration parameters and ODE integration (initial)
y0 = np.append(r0_vec, v0_vec)
t = np.linspace(0., (.5*3600), 1000) # propagate .5 hours
sol0 = odeint(twoBody, y0, t, args=(mu, ), atol=tol)

# integration parameters and ODE integration (new)
y0 = np.append(rn_vec, vn_vec)
t = np.linspace(0., 5*Pn, 27612) # propagate 5 orbits
sol = odeint(twoBody, y0, t, args=(mu, ), atol=tol)

# cartsian plot
fig = plt.figure(1)
ax = plt.axes(projection='3d')
ax.plot3D(sol0[:, 0], sol0[:,1], sol0[:,2], label='Initial')
ax.plot3D(sol[:, 0], sol[:,1], sol[:,2], label='New')

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = rEarth*np.cos(u)*np.sin(v)
y = rEarth*np.sin(u)*np.sin(v)
z = rEarth*np.cos(v)
ax.plot_wireframe(x, y, z, color="g", label='Earth')

ax.set_xlabel('rx (km)')
ax.set_ylabel('ry (km)')
ax.set_zlabel('rz (km)')
ax.set_title('ODE Propagation')
plt.legend()
plt.show()


# ========================================================================================
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



# part 4 ==========================


def twoBody(y, t, mu):

    '''return 3D vector of velocity and accel. for a two body problem'''

    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)
    
    yDot = np.empty((6,))
    
    yDot[0] = y[3]
    yDot[1] = y[4]
    yDot[2] = y[5]
    yDot[3] = (-mu/(r**3))*y[0]
    yDot[4] = (-mu/(r**3))*y[1]
    yDot[5] = (-mu/(r**3))*y[2]
    
    return yDot



# position vector in terms of k1 and k2
def k_position(a, e, d1, d2):

    # common eqns
    p = a*(1 - e*e)
    k_sqrt = np.sqrt(1 - d1*d1) * np.sqrt(1 - d2*d2)

    # build position vectors
    r_const = a*( np.sqrt(1 - d1*d2 + k_sqrt) ) / (np.sqrt(2) * (2 + d1 + d2))
    rx = r_const * ( np.sqrt(1 - d1*d1) + np.sqrt(1 - d2*d2) )
    ry = r_const * (d1 - d2)
    rz = 0

    # build velocity vectors
    v_const = np.sqrt(mu/(2*p))
    vx = v_const * -np.sqrt(1 - d1*d2 - k_sqrt)
    vy = v_const * ( np.sqrt(1 + d1*d2 - k_sqrt) + np.sqrt(1 + d1*d2 + k_sqrt) )
    vz = 0

    rvVector = np.array([rx, ry, rz, vx, vy, vz]) # km

    return rvVector


# position vector in terms of k1 and k2
def k_position(a, k1, k2):

    k_sqrt = np.sqrt(1 - k1*k1) * np.sqrt(1 - k2*k2)
    rx = a*(np.sqrt(1 + k1*k2 + k_sqrt) + np.sqrt(1 + k1*k2 - k_sqrt)) / np.sqrt(2)
    ry = a*(k1 - k2) / 2

    rVector = np.array([rx, ry]) # km

    return rVector

'''# Kepler's equation in terms of k1 and k2
def k_Kepler(k1, k2):

    kep = np.arctan((k1 - k2) / (np.sqrt(1 - k1*k1) + np.sqrt(1 - k2*k2))) - (np.sqrt(1 - k1*k1) - np.sqrt(1 - k2*k2))/2

    return kep'''

# propagate orbit in terms of k
def k_propagation(t0, tf, mu, a, e, sol_array):

    n = np.sqrt(mu/(a**3))                  # sec
    E0 = meanAnomaly(t0, t0, n) + 0.25*e    # initial guess, rad
    phi = np.arcsin(e)                      # rad

    for t in np.arange(t0, tf):
        M = meanAnomaly(t0, t, n)                     # mean anomaly, rad (note: this is eqaul to the exam equation in k-terms)
        E = newtonRaphson(M, E0, e)                   # eccentric anomaly, rad

        if E >=0 and E <= (2*np.pi):

            # now that eccentric anomaly is known, k1 and k2 can be solved
            k1 = np.sin(phi + E)#; print(k1)
            k2 = np.sin(phi - E)#; print(k2)

            # re-run with new E
            E0 = E      

            temp = k_position(a, k1, k2)
            sol_array.append(temp)
            
    return sol_array

'''
empa_sol = [] ; phi_all = []
k_sol = k_propagation(0, 2*Pn, mu, an, en, empa_sol)      # CHANGE to 5 OrbitS !!!!!!
k_sol = np.reshape(k_sol, (len(k_sol), 2))

# plot orbital position
fig = plt.figure(2)
ax = plt.axes(projection='3d')
ax.plot3D(k_sol[:, 0], k_sol[:,1], label='K-Sol')#, k_sol[:,2], label='K-Sol')

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = rEarth*np.cos(u)*np.sin(v)
y = rEarth*np.sin(u)*np.sin(v)
z = rEarth*np.cos(v)
ax.plot_wireframe(x, y, z, color="g", label='Earth')

ax.set_xlabel('rx (km)')
ax.set_ylabel('ry (km)')
ax.set_zlabel('rz (km)')
ax.set_title('K-Propagation (From Part A)')
plt.legend()
plt.show()
'''

# ========================================================================================





# =====================================


# part 5 ==========================

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

'''
def orbital_elements(rvvector, mu):
    
    rVec = rvvector[0:3]
    vVec = rvvector[3:6]
    hVec = np.cross(rVec,vVec)
    
    r = np.linalg.norm(rVec)
    v = np.linalg.norm(vVec)
    h = np.linalg.norm(hVec)
    
    p = h**2 / mu
    a = 1/(2/r - (v**2 / mu))
    e = np.sqrt(1 - (p/a))
    
    return a,e
'''

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

            # now we can find more parameters
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



emp_sol = [] #; v_all = []
kep_sol = keplerian_propagation(0, Pn, mu, an, en, inn, wn, Omegan, emp_sol)
kep_sol = np.reshape(kep_sol, (len(kep_sol), 6))

# plot orbital position
fig = plt.figure(3)
ax = plt.axes(projection='3d')
ax.plot3D(kep_sol[:, 0], kep_sol[:,1], kep_sol[:,2], label='New')

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = rEarth*np.cos(u)*np.sin(v)
y = rEarth*np.sin(u)*np.sin(v)
z = rEarth*np.cos(v)
ax.plot_wireframe(x, y, z, color="g", label='Earth')

ax.set_xlabel('rx (km)')
ax.set_ylabel('ry (km)')
ax.set_zlabel('rz (km)')
ax.set_title('Keplerian Propagation')
plt.legend()
plt.show()



# part 6 ==========================
# compare the accuracies of the orbit propagation

print('\nTo compare the accuracies of each solver, it is assumed that the Keplerian orbit '
      'propagation is correct (the accepted value) and the ODE orbital propagation is the '
      'value being tested (experimental value). By comparing the two propagations, the '
      'error can be calculated. ')


def accuracy_checker(accepted_value, experimental_value, error_list):

    '''takes the accepted value, experimental value, and returns the error list'''

    # ensure that accepted_value and experimental_value are of the same length
    # note: here t is used since the length is based on orbital period
    for t in np.arange(len(experimental_value)):    

        error = (accepted_value[t] - experimental_value[t]) / accepted_value[t]

        error_list.append(error)

    return error_list

# set empty arrays to be plotted
rx_error = []; ry_error = []; rz_error = []
vx_error = []; vy_error = []; vz_error = []

# check accuracy of position and velocity
rx_error = accuracy_checker(kep_sol[:, 0], sol[:, 0], rx_error)
ry_error = accuracy_checker(kep_sol[:, 1], sol[:, 1], ry_error)
rz_error = accuracy_checker(kep_sol[:, 2], sol[:, 2], rz_error)
vx_error = accuracy_checker(kep_sol[:, 3], sol[:, 3], vx_error)
vy_error = accuracy_checker(kep_sol[:, 4], sol[:, 4], vy_error)
vz_error = accuracy_checker(kep_sol[:, 5], sol[:, 5], vz_error)

# plot accuracy
fig3 = plt.figure(4)
plt.subplot(2,3,1)
plt.plot(t, rx_error)
plt.xlabel('time (s)')
plt.ylabel('x-position error')

plt.subplot(2,3,2)
plt.plot(t, ry_error)
plt.xlabel('time (s)')
plt.ylabel('y-position error')

plt.subplot(2,3,3)
plt.plot(t, rz_error)
plt.xlabel('time (s)')
plt.ylabel('z-position error')

plt.subplot(2,3,4)
plt.plot(t, vx_error)
plt.xlabel('time (s)')
plt.ylabel('x-velocity error')

plt.subplot(2,3,5)
plt.plot(t, vy_error)
plt.xlabel('time (s)')
plt.ylabel('y-velocity error')

plt.subplot(2,3,6)
plt.plot(t, vz_error)
plt.xlabel('time (s)')
plt.ylabel('z-velocity error')

plt.suptitle('Propagation Accuracy (5 Orbital Periods)')
plt.tight_layout()
plt.show()

# end