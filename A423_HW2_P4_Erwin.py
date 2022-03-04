# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 20:55:24 2021

@author: Chrissi Erwin
"""

import numpy as np
from scipy.integrate import odeint
from mpl_toolkits import mplot3d 
import matplotlib.pyplot as plt

# orbital parameters =============
a = 17000                 # km
e = 0.35
i = np.deg2rad(21.5)      # deg2rad
RAAN = np.deg2rad(40)     # deg2rad
ArgPeri = np.deg2rad(70)  # deg2rad
phi = np.deg2rad(60)      # deg2rad

mu = 398600.4415          # km
rEarth = 6378.135         # km


# part a =========================
p = a*(1 - e**2)

r = p/(1 + e*np.cos(phi))
phiDot = np.sqrt(mu*p)/r**2
rDot = phiDot * (p*e*np.sin(phi))/((1 + e*np.cos(phi))**2)

rVector = np.array([r*np.cos(phi), r*np.sin(phi), 0])
vVector = np.array([rDot*np.cos(phi) - r*phiDot*np.sin(phi), rDot*np.sin(phi) + r*phiDot*np.cos(phi), 0])

print('rVec: ', rVector) # km
print('vVec: ', vVector) # km/s


# part b =========================

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
y0 = np.append(rVector, vVector)
tf = 8. * 3600               # hours to sec
t = np.linspace(0., tf, 100) # hours

sol = odeint(twoBody, y0, t, args=(mu, ))

# cartsian plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(sol[:, 0], sol[:,1], sol[:,2])
plt.show()
'''
plt.plot(sol[:, 0], sol[:,1], sol[:,2])
plt.title('Earth Orbiting Satellite')
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.grid('True')
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.show()


# polar graphing
nStep = len(sol)
rNorm = np.zeros((nStep,))
theta = np.zeros((nStep,))

for i in range(nStep):
    x = sol[i, 0]
    y = sol[i, 1]
    rNorm[i] = np.linalg.norm(sol[i, 0:2])
    theta[i] = np.arctan2(sol[i, 4], sol[i, 3])


ax = plt.subplot(111, projection = 'polar')
ax.set_title('Earth Orbiting Satellite')
#ax.plot(theta[0:25], rNorm[0:25], 'b', theta[25:50], rNorm[25:50], 'g', theta[50:75], rNorm[50:75], 'r', theta[75:100], rNorm[75:100], 'y')
ax.plot(theta, rNorm, label='8 hour orbit [km]')
'''



