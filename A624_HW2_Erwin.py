# -*- coding: utf-8 -*-
"""
Created on Feb 2022

@author: Chrissi Erwin

Simulate the deployment of each of the first 4 vehicles by determining each of their orbits
and propagating them for 1 orbital time period (for their corresponding orbital time period)
"""

import numpy as np
from scipy.integrate import odeint
from mpl_toolkits import mplot3d 
import matplotlib.pyplot as plt

# orbital parameters =============
# 8 satellites
mu = 398600.4415          # universal gravitational parameter, km
rEarth = 6378.135         # radius of Earth, km

e = 0.01                   # eccentricity
rp = 7100 + rEarth         # radius at perigee, km
a = rp/(1-e)               # semi-major axis, km
i = np.deg2rad(10.0)       # inclination, deg2rad
RAAN = np.deg2rad(5.0)     # right ascension of the ascending node, deg2rad
ArgPeri = np.deg2rad(15.0) # arg. of perigee, deg2rad
#phi = np.deg2rad(0.0)      # true anomaly, deg2rad

p = a*(1 - e**2)           # semi-latus rectum, km
h = np.sqrt(p*mu)          # angular momentum, km^2/s



# =====================================

def func(E): 
	return (E - e*np.sin(E) - M)

def derivFunc(E): 
	return (1 - e*np.cos(E))

# function to find the root 
def newtonRaphson(E, e): 
    h = func(E) / derivFunc(E) 
    
	# convergence criteria is 10e-4
    while abs(h) >= 0.0001: 
        h = func(E)/derivFunc(E) 
		
		# E(k+1) = E(k) - f(E) / f'(E) 
        E = E - h 

    return E




# Part 1 =========================


# NOTE: CHANGE TO r,theta COORDINATES !!!!!!
def position_velocity(M, RAAN, ArgPeri, i, mu, h, e):

    # initial value to find eccentric anomaly 
    E0 = M + 0.25*e
    E = newtonRaphson(E0, e) 

    # find true anomaly 
    phi = np.arccos((np.cos(E) - e) / (1 - e*np.cos(E))) # true anomaly, rad
    E = np.arccos((e + np.cos(phi))/(1 + e*np.cos(phi))) # eccentric anomaly, rad

    # now we can find theta and radius
    theta = ArgPeri + phi      # angle, rad
    r = p/(1 + e*np.cos(phi))  # radius, km

    # now, we can find radius and velocity
    rx = r*(np.cos(RAAN)*np.cos(theta) - np.sin(RAAN)*np.sin(theta)*np.cos(i))
    ry = r*(np.sin(RAAN)*np.cos(theta) + np.cos(RAAN)*np.sin(theta)*np.cos(i))
    rz = r*(np.sin(theta)*np.sin(i))

    vx = -(mu/h)*(np.cos(RAAN)*(np.sin(theta) + e*np.sin(ArgPeri)) + np.sin(RAAN)*(np.cos(theta) + e*np.cos(ArgPeri))*np.cos(i))
    vy = -(mu/h)*(np.sin(RAAN)*(np.sin(theta) + e*np.sin(ArgPeri)) - np.cos(RAAN)*(np.cos(theta) + e*np.cos(ArgPeri))*np.cos(i))
    vz = +(mu/h)*(np.cos(theta) + e*np.sin(ArgPeri)*np.sin(i))

    rVector = np.array([rx, ry, rz]) # km
    vVector = np.array([vx, vy, vz]) # km/s

    return np.append(rVector, vVector)


# return 3D vector of velocity and accel.
# for a two body problem
def twoBody(y, t, mu):
    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)
    
    yDot = np.empty((6,))
    
    yDot[0] = y[3]
    yDot[1] = y[4]
    yDot[2] = y[5]
    yDot[3] = (-mu/(r**3))*y[0]
    yDot[4] = (-mu/(r**3))*y[1]
    yDot[5] = (-mu/(r**3))*y[2]
    
    return yDot
    
# integration parameters and ODE integration
for t in np.arange(4):
    
    n = np.sqrt(mu/a**3) # sec
    t0 = 0		         # sec at perigee
    t = 300 * 60         # min to sec
    M = n*(t-t0)         # mean anomaly
    
    
y0 = position_velocity(RAAN, ArgPeri, i, mu, h, e)
tf = 4. * 86400                 # days to sec
t = np.linspace(0., tf, 1000)   # days

sol = odeint(twoBody, y0, t, args=(mu, ))

# plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(sol[:, 0], sol[:,1], sol[:,2])

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = rEarth*np.cos(u)*np.sin(v)
y = rEarth*np.sin(u)*np.sin(v)
z = rEarth*np.cos(v)
ax.plot_wireframe(x, y, z, color="g")

plt.show()


# Part 2 =========================




