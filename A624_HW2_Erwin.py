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
mu = 398600.4415          # universal gravitational parameter, km^3/s^2
rEarth = 6378.135         # radius of Earth, km

e = 0.01                   # eccentricity
rp = 7100                  # radius at perigee, km
a = rp/(1-e)               # semi-major axis, km
i = np.deg2rad(10.0)       # inclination, deg2rad
RAAN = np.deg2rad(5.0)     # right ascension of the ascending node, deg2rad
ArgPeri = np.deg2rad(15.0) # arg. of perigee, deg2rad
#phi = np.deg2rad(0.0)     # true anomaly, deg2rad

p = a*(1 - e**2)           # semi-latus rectum, km
h = np.sqrt(p*mu)          # angular momentum, km^2/s

n = np.sqrt(mu/(a**3))     # sec
t0 = 0		               # sec at perigee
tf = np.pi*1/n
#tf = 8. * 3600             # hours to sec         CHECK TIMES AND FIX ORBITAL PERIOD; CHECK RESULTS OF Y ARRAY

# =====================================

# find mean and eccentric anomalies
def meanAnomaly(t, n):
    return n*(t) #* 2*np.pi

def func(E, e): 
	return (E - e*np.sin(E) - M)

def derivFunc(E, e): 
	return (1 - e*np.cos(E))

# function to find the root 
def newtonRaphson(E, e): 
    h = func(E, e) / derivFunc(E, e) 
    
	# convergence criteria 
    while abs(h) >= 0.000001: 
        h = func(E, e)/derivFunc(E, e) 
		
		# E(k+1) = E(k) - f(E) / f'(E) 
        E = E - h 

    return E


def position_velocity(r, theta, RAAN, ArgPeri, i, mu, h, e):

    # # initial value to find eccentric anomaly 
    # E0 = M + 0.25*e
    # E = newtonRaphson(E0, e) 

    # # find true anomaly 
    # phi = np.arccos((np.cos(E) - e) / (1 - e*np.cos(E))) # true anomaly, rad
    # E = np.arccos((e + np.cos(phi))/(1 + e*np.cos(phi))) # eccentric anomaly, rad

    # # now we can find theta and radius
    # theta = ArgPeri + phi      # angle, rad
    # r = p/(1 + e*np.cos(phi))  # radius, km

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

def position_velocity2(r, mu, e, h, phi):

    rVector = np.array([r*np.cos(phi), r*np.sin(phi), 0])
    vVector = np.array([(-mu/h)*np.sin(phi), (mu/h)*(e + np.cos(phi)), 0])
    
    return [rVector, vVector]

def twoBody(y, t, mu):

    # return 3D vector of velocity and accel.
    # for a two body problem

    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)
    
    yDot = np.empty((6,))
    
    yDot[0] = y[3]
    yDot[1] = y[4]
    yDot[2] = y[5]
    yDot[3] = (-mu/(r**3))*y[0]
    yDot[4] = (-mu/(r**3))*y[1]
    yDot[5] = (-mu/(r**3))*y[2]
    
    return yDot


# y = tf*np.zeros((6,1))
# orbitalElements = tf*np.zeros((6,1))
y = []; phi_all = []; E_all = []
E0 = meanAnomaly(t0, n) + 0.25*e
for t in np.arange(t0, tf):
    M = meanAnomaly(t, n)                      # mean anomaly, rad
    print(M)                                    # FIX MEAN ANOMALY TO REACH 2pi !!!!!!
    E = newtonRaphson(E0, e)                   # eccentric anomaly, rad

    if E >=0 and E <= (2*np.pi):

        E_all.append(E)
        phi = np.arccos((np.cos(E) - e) / (1 - e*np.cos(E))) # true anomaly, rad
        phi_all.append(phi)

        # now we can find theta and radius
        theta = ArgPeri + phi      # angle, rad
        r = p/(1 + e*np.cos(phi))  # radius, km

        # re-run with new E
        E0 = E      

        y.append(position_velocity(r, theta, RAAN, ArgPeri, i, mu, h, e))
        #y.append(position_velocity2(r, mu, e, h, phi))
    
print(len(phi_all))
print(len(E_all))
y = np.reshape(y, (len(E_all), 6))




    
    
# integration parameters and ODE integration
#for t in np.arange(4):  
    
#y0 = position_velocity(RAAN, ArgPeri, i, mu, h, e)
#tf = 4. * 86400                 # days to sec
#tf = 2*np.pi*np.sqrt(a**3/mu)
#t = np.linspace(0., tf, 1000)   # days

#sol = odeint(twoBody, y, t, args=(mu, ))

# plot orbital position
fig = plt.figure()
ax = plt.axes(projection='3d')
# ax.plot3D(sol[:, 0], sol[:,1], sol[:,2])
ax.plot3D(y[:, 0], y[:,1], y[:,2])

# # draw sphere
# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
# x = rEarth*np.cos(u)*np.sin(v)
# y = rEarth*np.sin(u)*np.sin(v)
# z = rEarth*np.cos(v)
# ax.plot_wireframe(x, y, z, color="g")

plt.show()








print('Deployment maifest:')
print('True anomolies:')
print('Orbital elements:')
print('The post deployment orbital elements were arrived at through the '
      'use of the mean anomaly to ')
#print('phi\n',phi_all)

# Part 2 =========================




