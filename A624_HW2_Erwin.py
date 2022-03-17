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
satellite_num = 8         # 8 satellites
mu = 398600.4415          # universal gravitational parameter, km^3/s^2
rEarth = 6378.135         # radius of Earth, km

e = 0.01                   # eccentricity
rp = 7100                  # radius at perigee, km
a = rp/(1-e)               # semi-major axis, km
i = np.deg2rad(10.0)       # inclination, deg2rad
RAAN = np.deg2rad(5.0)     # right ascension of the ascending node, deg2rad
ArgPeri = np.deg2rad(15.0) # arg. of perigee, deg2rad
deltaV = 0.05              # km/s

p = a*(1 - e**2)           # semi-latus rectum, km
h = np.sqrt(p*mu)          # angular momentum, km^2/s

n = np.sqrt(mu/(a**3))     # sec
t0 = 0		               # sec at perigee
tf = 2*np.pi*np.sqrt((a**3)/mu)
tq = tf/satellite_num

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

# empty arrays for later use
LV = []; sat1 = []; sat2 = []; sat3 = []; sat4 = []
phi_all = []; E_all = []

E0 = meanAnomaly(t0, n) + 0.25*e
for t in np.arange(t0, tf):
    M = meanAnomaly(t, n)                      # mean anomaly, rad
    E = newtonRaphson(E0, e)                   # eccentric anomaly, rad

    if E >=0 and E <= (2*np.pi):

        E_all.append(E)
        # phi = np.arccos((np.cos(E) - e) / (1 - e*np.cos(E))) # true anomaly, rad
        # https://en.wikipedia.org/wiki/True_anomaly
        phi = M + (2*e - (e**3)/4)*np.sin(M) + (5*(e**2) * np.sin(2*M))/4 + (13*(e**3) * np.sin(3*M))/12
        phi_all.append(phi)

        # now we can find theta and radius
        theta = ArgPeri + phi      # angle, rad
        r = p/(1 + e*np.cos(phi))  # radius, km

        # re-run with new E
        E0 = E      

        LV.append(position_velocity(r, theta, RAAN, ArgPeri, i, mu, h, e))

#print(phi_all[3000:3200],'\n',E_all[3000:3200])

# data for satellite 1
E1 = []; phi1 = []
for t in np.arange(1*tq, tf):
    M = meanAnomaly(t, n)                      # mean anomaly, rad
    E = newtonRaphson(E0, e)                   # eccentric anomaly, rad
    E1.append(E)

    if E >=0 and E <= (2*np.pi):

        phi = M + (2*e - (e**3)/4)*np.sin(M) + (5*(e**2) * np.sin(2*M))/4 + (13*(e**3) * np.sin(3*M))/12
        phi1.append(phi)

        # now we can find theta and radius
        theta = ArgPeri + phi      # angle, rad
        r = p/(1 + e*np.cos(phi))  # radius, km

        # re-run with new E
        E0 = E      

        temp = position_velocity(r, theta, RAAN, ArgPeri, i, mu, h, e)
        temp[1] = temp[1] + deltaV
        sat1.append(temp)

# data for satellite 2
E2 = []; phi2 = []
for t in np.arange(2*tq, tf):
    M = meanAnomaly(t, n)                      # mean anomaly, rad
    E = newtonRaphson(E0, e)                   # eccentric anomaly, rad
    E2.append(E)

    if E >=0 and E <= (2*np.pi):

        phi = M + (2*e - (e**3)/4)*np.sin(M) + (5*(e**2) * np.sin(2*M))/4 + (13*(e**3) * np.sin(3*M))/12
        phi2.append(phi)

        # now we can find theta and radius
        theta = ArgPeri + phi      # angle, rad
        r = p/(1 + e*np.cos(phi))  # radius, km

        # re-run with new E
        E0 = E      

        temp = position_velocity(r, theta, RAAN, ArgPeri, i, mu, h, e)
        temp[1] = temp[1] + deltaV
        sat2.append(temp)

# data for satellite 3
E3 = []; phi3 = []
for t in np.arange(3*tq, tf):
    M = meanAnomaly(t, n)                      # mean anomaly, rad
    E = newtonRaphson(E0, e)                   # eccentric anomaly, rad
    E3.append(E)

    if E >=0 and E <= (2*np.pi):

        phi = M + (2*e - (e**3)/4)*np.sin(M) + (5*(e**2) * np.sin(2*M))/4 + (13*(e**3) * np.sin(3*M))/12
        phi3.append(phi)

        # now we can find theta and radius
        theta = ArgPeri + phi      # angle, rad
        r = p/(1 + e*np.cos(phi))  # radius, km

        # re-run with new E
        E0 = E      

        temp = position_velocity(r, theta, RAAN, ArgPeri, i, mu, h, e)
        temp[1] = temp[1] + deltaV
        sat3.append(temp)

# data for satellite 4
E4 = []; phi4 = []
for t in np.arange(4*tq, tf):
    M = meanAnomaly(t, n)                      # mean anomaly, rad
    E = newtonRaphson(E0, e)                   # eccentric anomaly, rad
    E4.append(E)

    if E >=0 and E <= (2*np.pi):

        phi = M + (2*e - (e**3)/4)*np.sin(M) + (5*(e**2) * np.sin(2*M))/4 + (13*(e**3) * np.sin(3*M))/12
        phi4.append(phi)

        # now we can find theta and radius
        theta = ArgPeri + phi      # angle, rad
        r = p/(1 + e*np.cos(phi))  # radius, km

        # re-run with new E
        E0 = E      

        temp = position_velocity(r, theta, RAAN, ArgPeri, i, mu, h, e)
        temp[1] = temp[1] + deltaV
        sat4.append(temp)

# reshape solutions
LV = np.reshape(LV, (len(E_all), 6))
sat1 = np.reshape(sat1, (len(E1), 6))
sat2 = np.reshape(sat2, (len(E2), 6))
sat3 = np.reshape(sat3, (len(E3), 6))
sat4 = np.reshape(sat4, (len(E4), 6))


# plot orbital position
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(LV[:, 0], LV[:,1], LV[:,2], label='LV Orbit')
ax.plot3D(sat1[:, 0], sat1[:,1], sat1[:,2])
ax.plot3D(sat2[:, 0], sat2[:,1], sat2[:,2])
ax.plot3D(sat3[:, 0], sat3[:,1], sat3[:,2])
ax.plot3D(sat4[:, 0], sat4[:,1], sat4[:,2])

ax.scatter(LV[0,0], LV[1,1], LV[2,2], label='Periapse')
ax.scatter(sat1[0, 0], sat1[1,1], sat1[2,2], label = 'Sat1')
ax.scatter(sat2[0, 0], sat2[1,1], sat2[2,2], label = 'Sat2')
ax.scatter(sat3[0, 0], sat3[1,1], sat3[2,2], label = 'Sat3')
ax.scatter(sat4[0, 0], sat4[1,1], sat4[2,2], label = 'Sat4')

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = rEarth*np.cos(u)*np.sin(v)
y = rEarth*np.sin(u)*np.sin(v)
z = rEarth*np.cos(v)
ax.plot_wireframe(x, y, z, color="yellowgreen", label='Earth')

ax.set_xlabel('rx')
ax.set_ylabel('ry')
ax.set_zlabel('rz')
plt.legend()
plt.show()


# results
t = np.arange(t0, tf)

print('Deployment maifest:')
print('----------------------')
print('True anomalies:')
print('time (s):','\t',t[1000], '\t\t',t[3072], '\t\t', t[5300])
print(u'\u03BD (rad):','\t', phi_all[1000],'\t', phi_all[3072],'\t', phi_all[5300]) # true anomaly, nu
print('----------------------')
print('Orbital elements:')
print('a (km):','\t', a)
print('e:','\t\t', e)
print('i (rad):','\t', i)
print(u'\u038F (rad):','\t', RAAN) # RAAN, Omega
print(u'\u03C9 (rad):','\t', ArgPeri) # Argument of perigee, omega
print('\nThe post deployment orbital elements were arrived at through the '
      'use of time (t) and time at periapse (t0) to find the mean anomaly (M). From M '
      'the eccentric (E) and true anomalies (v) to be found, from which '
      'the rest of the orbital elements (a, e, i, RAAN, w, v) are known. '
      'Note that the only orbital elements changing with time are the mean, '
      'eccentric, and true anomalies. Since these are all related only the true '
      'anomaly is shown.')


# Part 2 =========================




