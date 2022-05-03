# A 624 Final Exam, due May 2, 2022
# Problems are taken from Chapter 10 of Schaub and Junkins

import numpy as np
from scipy.integrate import odeint
from mpl_toolkits import mplot3d 
import matplotlib.pyplot as plt


# COMMON FUNCTIONS

def twoBody(y, t, mu):

    '''return 3D vector of velocity and accel. for cylindrical coordinates'''

    r = np.sqrt(y[0]**2 + y[1]**2 + y[2]**2)
    
    yDot = np.empty((6,))
    
    yDot[0] = y[3]
    yDot[1] = y[4]
    yDot[2] = y[5]
    yDot[3] = (-mu/(r**3))*y[0]
    yDot[4] = (-mu/(r**3))*y[1]
    yDot[5] = (-mu/(r**3))*y[2]
    
    return yDot



# eqn 10.9, generic time varying function
def f(t, rij, rijo):
    return rij / rijo

# initial position and velocity vectors function
def init_position_velocity_vec(t, rij, rijo):

    # updated relative position vectors
    r12 = rij[0]; r13 = rij[1]; r23 = rij[2]
    r12o = rijo[0]; r13o = rijo[1]; r23o = rijo[2]

    # initial radial distances, eqn. 10.11
    r1o = np.sqrt(m2**2 * r12o**2 + m3**2 * r13o**2 + 2*m2*m3*r12o*r13o*np.cos(alpha)) / M
    r2o = np.sqrt(m1**2 * r12o**2 + m3**2 * r23o**2 - 2*m1*m3*r12o*r23o*np.cos(alpha)) / M
    r3o = np.sqrt(m1**2 * r13o**2 + m2**2 * r23o**2 + 2*m1*m2*r13o*r23o*np.cos(alpha)) / M

    # transfer to radial distances (m)
    r1_t = r1o * f(t, r12, r12o)        # r1(t), eqn. 10.12
    r2_t = r2o * f(t, r13, r13o)        # r2(t), eqn. 10.13
    r3_t = r3o * f(t, r23, r23o)        # r3(t), eqn. 10.14
    #print(r1_t); print(r2_t); print(r3_t)
    #print(r1o); print(r2o); print(r3o)

    #w = r1_dot_mag*np.sin(gamma) / r1_t       # 1/s, angular acceleration (v/r)
    w = np.sqrt(G*M / r12**3)

    # eqn. 10.35, rearranged (m/s)
    r1_dot_t = np.sqrt(( r1_dot_mag**2 - r1_t**2 * w**2 )) 
    r2_dot_t = np.sqrt(( r2_dot_mag**2 - r2_t**2 * w**2 )) 
    r3_dot_t = np.sqrt(( r3_dot_mag**2 - r3_t**2 * w**2 )) 

    # print(( r1_dot_mag**2 - r1_t**2 * w**2 ))
    # print(( r2_dot_mag**2 - r2_t**2 * w**2 ))
    # print(( r3_dot_mag**2 - r3_t**2 * w**2 ))
    # print(r1_dot_t); print(r2_dot_t); print(r3_dot_t)

    '''
    # eqn. 10.36, rearranged (m/s)
    r1_dot = r1_t * w / np.tan(gamma)
    r2_dot = r2_t * w / np.tan(gamma)
    r3_dot = r3_t * w / np.tan(gamma)'''
    print(r3_t*w)

    # build position and velocity vectors
    r1_vec = np.array([ r1_t, 0, 0 ])               # eqn. 10.17
    r2_vec = np.array([ r2_t, 0, 0 ]) 
    r3_vec = np.array([ r3_t, 0, 0 ]) 
    v1_vec = np.array([ r1_dot_t, r1_t*w, 0 ])      # eqn. 10.18
    v2_vec = np.array([ r2_dot_t, r2_t*w, 0 ])
    v3_vec = np.array([ r3_dot_t, r3_t*w, 0 ])  
    '''
    rv1_vec = np.append(r1_vec, v1_vec)
    rv2_vec = np.append(r2_vec, v2_vec)
    rv3_vec = np.append(r3_vec, v3_vec)

    rv_vec = np.append(rv1_vec, rv2_vec, rv3_vec)'''

    return r1_vec, v1_vec, r2_vec, v2_vec, r3_vec, v3_vec

# returns position and velocity vectors for problem 1
def pos_vel_vecs(G, masses, rijo, ri_dot_mag, angles):

    # redefine input vectors
    M = masses[0]; m1 = masses[1]; m2 = masses[2]; m3 = masses[3]
    r12o = rijo[0]; r13o = rijo[1]; r23o = rijo[2]
    r1_dot_mag = ri_dot_mag[0]; r2_dot_mag = ri_dot_mag[1]; r3_dot_mag = ri_dot_mag[2]
    gamma = angles[0]; alpha1 = angles[1]; alpha2 = angles[2]; alpha3 = angles[3]

    # angular velocity (eqn. 10.66)
    w = np.sqrt(G*(m1+m2) / r12o**3)

    # initial position vectors, r1=r1o, r2=r2o, r3=r3o (eqns. 10.11-10.14)
    r1 = np.sqrt(m2**2 * r12o**2 + m3**2 * r13o**2 + 2*m2*m3*r12o*r13o*np.cos(alpha1)) / M
    r2 = np.sqrt(m1**2 * r12o**2 + m3**2 * r23o**2 - 2*m1*m3*r12o*r23o*np.cos(alpha2)) / M
    r3 = np.sqrt(m1**2 * r13o**2 + m2**2 * r23o**2 + 2*m1*m2*r13o*r23o*np.cos(alpha3)) / M

    # initial velocity vectors (eqn. 10.35)
    try:
        r1_dot = np.sqrt(r1_dot_mag**2 - r1*r1*w*w)
        r2_dot = np.sqrt(r2_dot_mag**2 - r2*r2*w*w)
        r3_dot = np.sqrt(r3_dot_mag**2 - r3*r3*w*w)
    finally:
        r1_dot = np.sqrt(np.abs(r1_dot_mag**2 - r1*r1*w*w))
        r2_dot = np.sqrt(np.abs(r2_dot_mag**2 - r2*r2*w*w))
        r3_dot = np.sqrt(np.abs(r3_dot_mag**2 - r3*r3*w*w))

    # CYLINDRICAL 2 CARTESIAN COORDINATES
    # init empty pos/vel vector + knowns
    rvVector1 = np.zeros((6,))
    rvVector2 = np.zeros((6,))
    rvVector3 = np.zeros((6,))
    theta = gamma
    thetaDot = w

    # mass 1
    rvVector1[0] = r1*np.cos(theta)
    rvVector1[1] = r1*np.sin(theta)
    # rz will stay zero
    rvVector1[3] = r1_dot*np.cos(theta) - r1*thetaDot*np.sin(theta)
    rvVector1[4] = r1_dot*np.sin(theta) + r1*thetaDot*np.cos(theta)
    # vz will stay zero

    # mass 2
    rvVector2[0] = r2*np.cos(theta)
    rvVector2[1] = r2*np.sin(theta)
    # rz will stay zero
    rvVector2[3] = r2_dot*np.cos(theta) - r2*thetaDot*np.sin(theta)
    rvVector2[4] = r2_dot*np.sin(theta) + r2*thetaDot*np.cos(theta)
    # vz will stay zero

    # mass 3
    rvVector3[0] = r3*np.cos(theta)
    rvVector3[1] = r3*np.sin(theta)
    # rz will stay zero
    rvVector3[3] = r3_dot*np.cos(theta) - r3*thetaDot*np.sin(theta)
    rvVector3[4] = r3_dot*np.sin(theta) + r3*thetaDot*np.cos(theta)
    # vz will stay zero

    return rvVector1, rvVector2, rvVector3


# PROBLEM 1: EXAMPLE 10.1 (pg. 584) ==============================
'''To illustrate the gen. equilateral solution of the three-body problem, the motion 
of the following three-body system is numerically solved using eq. 10.39.'''

# givens
m1 = 5.967e23           # kg (1/10 of Earth's mass)
m2 = 7.35e22            # kg (the moon's mass)
m3 = 3.657e22           # kg (half of the moon's mass)

l0 = 1e9                # m, initial length where l0 = r12 = r23 = r13
gamma = np.deg2rad(40)  # rad, angle of initial velocity to respective radial position vector
alpha = np.deg2rad(60)  # rad, flight path angle
r1_dot_mag = 29.8659    # m/s, magnitude of inital velocity vector 1
r2_dot_mag = 189.181    # m/s, magnitude of inital velocity vector 2
r3_dot_mag = 195.552    # m/s, magnitude of inital velocity vector 3

# mass eqns (kg)
M = m1 + m2 + m3                            # eqn. 10.2
M1 = (m2**2 + m3**2 + m2*m3)**(3/2) / M**2  # eqn 10.40
M2 = (m1**2 + m3**2 + m1*m3)**(3/2) / M**2  # eqn 10.41
M3 = (m2**2 + m1**2 + m1*m2)**(3/2) / M**2  # eqn 10.42

# gravitational parameters
G = 6.67430e-11     # m^3/kg*s^2, gravitational constant
mu1 = G*M1          # m^3/s^2, gravitational parameter for mass 1
mu2 = G*M2          # m^3/s^2, gravitational parameter for mass 2
mu3 = G*M3          # m^3/s^2, gravitational parameter for mass 3


# input parameters for all orbits
masses = np.array([M, m1, m2, m3])  # kg, masses
# rij  = np.array([l0, l0, l0])       # m, relative distance, initial update
rijo = np.array([l0, l0, l0])       # m, initial relative distance
ri_dot_mag = np.array([r1_dot_mag, r2_dot_mag, r3_dot_mag]) # initial velocity vector magnitudes
angles = np.array([gamma, np.deg2rad(60), np.deg2rad(60), np.deg2rad(60)])   # initial angles
tol = 1e-6                          # ode tolerance
tf = 365*24*3600                        # s, propagation end time in hours (change to one year)
t = np.linspace(0, tf, tf)          # propagate 5 orbits
rvVector1, rvVector2, rvVector3 = pos_vel_vecs(G, masses, rijo, ri_dot_mag, angles)

# integration parameters and ODE integration (mass 1)
y0_orb1 = rvVector1
sol_orb1 = odeint(twoBody, y0_orb1, t, args=(mu1, ), atol=tol)

# integration parameters and ODE integration (mass 2)
y0_orb2 = rvVector2
sol_orb2 = odeint(twoBody, y0_orb2, t, args=(mu2, ), atol=tol)

# integration parameters and ODE integration (mass 3)
y0_orb3 = rvVector3
sol_orb3 = odeint(twoBody, y0_orb3, t, args=(mu3, ), atol=tol)

# cartsian plot
fig = plt.figure(1)
ax = plt.axes(projection='3d')
ax.plot3D(sol_orb1[:, 0], sol_orb1[:,1], label='Orbit of m1')
ax.plot3D(sol_orb2[:, 0], sol_orb2[:,1], label='Orbit of m2')
ax.plot3D(sol_orb3[:, 0], sol_orb3[:,1], label='Orbit of m3')

ax.set_xlabel('rx (m)')
ax.set_ylabel('ry (m)')
ax.set_zlabel('rz (m)')
ax.set_title('ODE Propagation')
plt.legend()
plt.show()
'''


# PROBLEM 2: EXAMPLE 10.2 (pg. 587) ==============================

# givens
chi = 0.451027          # Lagrange's quintic eqn for scalar parameter
x12 = 1e9               # m, scalar distance
x23 = 4.510273e8        # m, relative dist.
x13 = 1.451027e9        # m, "

r1_dot_mag = 32.9901    # m/s, magnitude of inital velocity vector 1
r2_dot_mag = 150.903    # m/s, magnitude of inital velocity vector 2
r3_dot_mag = 233.844    # m/s, magnitude of inital velocity vector 3

# input parameters for all orbits (changed from problem 1)
rijo = np.array([x12, x23, x13])                            # m, initial relative distance
ri_dot_mag = np.array([r1_dot_mag, r2_dot_mag, r3_dot_mag]) # initial velocity vector magnitudes 
rvVector1, rvVector2, rvVector3 = pos_vel_vecs(G, masses, rijo, ri_dot_mag, angles)

# integration parameters and ODE integration (mass 1)
y0_orb1 = rvVector1
sol_orb1 = odeint(twoBody, y0_orb1, t, args=(mu1, ), atol=tol)

# integration parameters and ODE integration (mass 2)
y0_orb2 = rvVector2
sol_orb2 = odeint(twoBody, y0_orb2, t, args=(mu2, ), atol=tol)

# integration parameters and ODE integration (mass 3)
y0_orb3 = rvVector3
sol_orb3 = odeint(twoBody, y0_orb3, t, args=(mu3, ), atol=tol)

# cartsian plot
fig = plt.figure(2)
ax = plt.axes(projection='3d')
ax.plot3D(sol_orb1[:, 0], sol_orb1[:,1], label='Orbit of m1')
ax.plot3D(sol_orb2[:, 0], sol_orb2[:,1], label='Orbit of m2')
ax.plot3D(sol_orb3[:, 0], sol_orb3[:,1], label='Orbit of m3')

ax.set_xlabel('rx (m)')
ax.set_ylabel('ry (m)')
ax.set_zlabel('rz (m)')
ax.set_title('ODE Propagation')
plt.legend()
plt.show()'''


# PROBLEM 3: EXAMPLE 10.6 (pg. 605) ==============================















# GRAVEYARD =========================++++++++++++++++++++++++++

# from P1
# # cylindrical to cartesian coordinates
# def cyl2cart(r_vec, theta_vec, z_vec):

#     '''converts a vector from cylindrical to cartesian coordinates'''
#     rVector = np.zeros((len(r_vec),3))  # build empty output

#     for i in np.arange(len(r_vec)):

#         r = r_vec[i]
#         theta = theta_vec[i]
#         z = z_vec[i]

#         rVector[i,0] = r*np.cos(theta)
#         rVector[i,1] = r*np.sin(theta)
#         rVector[i,2] = z
#     return rVector

# def cyl2cart(r_vec, v_vec):

#     # init len + empty pos/vel vector
#     nstep = 2*len(r_vec)
#     rvVector = np.zeros((nstep,))
#     thetaDot = w

#     # cylindrical to cartesian coordinates
#     rvVector[0] = r*np.cos(theta)
#     rvVector[1] = r*np.sin(theta)
#     # rz will stay zero
#     rvVector[3] = rDot*np.cos(theta) - r*thetaDot*np.sin(theta)
#     rvVector[4] = rDot*np.sin(theta) + r*thetaDot*np.cos(theta)
#     # vz will stay zero

#     return rvVector

# rVector1 = cyl2cart(sol_orb1[:, 0], sol_orb1[:, 1], sol_orb1[:, 2])
# rVector2 = cyl2cart(sol_orb2[:, 0], sol_orb2[:, 1], sol_orb2[:, 2])
# rVector3 = cyl2cart(sol_orb3[:, 0], sol_orb3[:, 1], sol_orb3[:, 2])

# ax.plot3D(rVector1[:, 0], rVector1[:,1], rVector1[:,2], label='Orbit of m1')
# ax.plot3D(rVector2[:, 0], rVector2[:,1], rVector2[:,2], label='Orbit of m2')
# ax.plot3D(rVector3[:, 0], rVector3[:,1], rVector3[:,2], label='Orbit of m3')

# from P2
'''
# assuming that A stays constant for all three masses
A = (-G*(m1+m2) + G*m3*(1/chi**2 - 1/(1+chi)**2)) / x12**3

x1 =  G*(m2/x12**2 + m3/x13**2) / A
x2 =  G*(m3/x23**2 - m1/x12**2) / A
x3 = -G*(m1/x13**2 + m2/x23**2) / A

# angular velocity
w = r1_dot_mag*np.sin(gamma) / x1       # 1/s, angular acceleration (v/r)

 # build position and velocity vectors
r1_vec = np.array([ x1, 0, 0 ])                 # eqn. 10.43
r2_vec = np.array([ x2, 0, 0 ]) 
r3_vec = np.array([ x3, 0, 0 ]) 
v1_vec = np.array([ r1_dot_mag, x1*w, 0 ])      # eqn. 10.18
v2_vec = np.array([ r2_dot_mag, x2*w, 0 ])
v3_vec = np.array([ r3_dot_mag, x3*w, 0 ])

# transform vectors from cylindrical to cartesian coordinates
rVector1 = cyl2cart(sol_orb1[:, 0], sol_orb1[:, 1], sol_orb1[:, 2])
rVector2 = cyl2cart(sol_orb2[:, 0], sol_orb2[:, 1], sol_orb2[:, 2])
rVector3 = cyl2cart(sol_orb3[:, 0], sol_orb3[:, 1], sol_orb3[:, 2])

# integration parameters and ODE integration (mass 1)
y0_orb1 = np.append(r1_vec, v1_vec)
sol_orb1 = odeint(twoBody, y0_orb1, t, args=(mu1, ), atol=tol)

# integration parameters and ODE integration (mass 2)
y0_orb2 = np.append(r2_vec, v2_vec)
sol_orb2 = odeint(twoBody, y0_orb2, t, args=(mu2, ), atol=tol)

# integration parameters and ODE integration (mass 3)
y0_orb3 = np.append(r3_vec, v3_vec)
sol_orb3 = odeint(twoBody, y0_orb3, t, args=(mu3, ), atol=tol)

# cartsian plot
fig = plt.figure(1)
ax = plt.axes(projection='3d')
ax.plot3D(rVector1[:, 0], rVector1[:,1], rVector1[:,2], label='Orbit of m1')
ax.plot3D(rVector2[:, 0], rVector2[:,1], rVector2[:,2], label='Orbit of m2')
ax.plot3D(rVector3[:, 0], rVector3[:,1], rVector3[:,2], label='Orbit of m3')

ax.set_xlabel('rx (m)')
ax.set_ylabel('ry (m)')
# ax.set_zlabel('rz (m)')
ax.set_title('ODE Propagation')
plt.legend()
plt.show()'''




# Collinear solution (second invariant shape for 3BP)

'''
# initial position and velocity vectors function
def init_position_velocity_vec_P2(xijo):

    # updated relative position vectors
    x12o = xijo[0]; x13o = xijo[1]; x23o = xijo[2]

    # initial radial distances, eqn. 10.11
    x1o = np.sqrt(m2**2 * x12o**2 + m3**2 * x13o**2 + 2*m2*m3*x12o*x13o*np.cos(alpha)) / M
    x2o = np.sqrt(m1**2 * x12o**2 + m3**2 * x23o**2 - 2*m1*m3*x12o*x23o*np.cos(alpha)) / M
    x3o = np.sqrt(m1**2 * x13o**2 + m2**2 * x23o**2 + 2*m1*m2*x13o*x23o*np.cos(alpha)) / M

    # transfer to radial distances (m)
    x1 = x1o        # r1(t), eqn. 10.12
    x2 = x2o        # r2(t), eqn. 10.13
    x3 = x3o        # r3(t), eqn. 10.14

    w = r1_dot_mag*np.sin(gamma) / x1       # 1/s, angular acceleration (v/r)
    
    # eqn. 10.35, rearranged (m/s)
    r1_dot = np.sqrt(( r1_dot_mag**2 - x1**2 * w**2 )) 
    r2_dot = np.sqrt(( r2_dot_mag**2 - x2**2 * w**2 )) 
    r3_dot = np.sqrt(( r3_dot_mag**2 - x3**2 * w**2 )) 

    # build position and velocity vectors
    r1_vec = np.array([ x1, 0, 0 ])             # eqn. 10.43
    r2_vec = np.array([ x2, 0, 0 ]) 
    r3_vec = np.array([ x3, 0, 0 ]) 
    v1_vec = np.array([ r1_dot, x1*w, 0 ])      # eqn. 10.18
    v2_vec = np.array([ r2_dot, x2*w, 0 ])
    v3_vec = np.array([ r3_dot, x3*w, 0 ])
    
    #print('Check:', np.cross(r1_vec,r2_vec)) # should be zero

    return r1_vec, v1_vec, r2_vec, v2_vec, r3_vec, v3_vec


xijo = np.array([x12, x23, x13])
r1_vec, v1_vec, r2_vec, v2_vec, r3_vec, v3_vec = init_position_velocity_vec_P2(xijo)

# integration parameters and ODE integration (mass 1)
y0_orb1 = np.append(r1_vec, v1_vec)
sol_orb1 = odeint(twoBody, y0_orb1, t, args=(mu1, ), atol=tol)

# integration parameters and ODE integration (mass 2)
y0_orb2 = np.append(r2_vec, v2_vec)
sol_orb2 = odeint(twoBody, y0_orb2, t, args=(mu2, ), atol=tol)

# integration parameters and ODE integration (mass 3)
y0_orb3 = np.append(r3_vec, v3_vec)
sol_orb3 = odeint(twoBody, y0_orb3, t, args=(mu3, ), atol=tol)

# plot in the format of (theta, r)
rNorm1 = np.zeros((tf,)); theta1 = np.zeros((tf,))
rNorm2 = np.zeros((tf,)); theta2 = np.zeros((tf,))
rNorm3 = np.zeros((tf,)); theta3 = np.zeros((tf,))
for t in np.arange(tf):
    rNorm1[t] = np.linalg.norm(sol_orb1[t, 0:3])
    # theta1[t] = np.arctan2(sol_orb1[t, 4], sol_orb1[t, 3])
    rNorm2[t] = np.linalg.norm(sol_orb2[t, 0:3])
    # theta2[t] = np.arctan2(sol_orb2[t, 4], sol_orb2[t, 3])
    rNorm3[t] = np.linalg.norm(sol_orb3[t, 0:3])
    # theta3[t] = np.arctan2(sol_orb3[t, 4], sol_orb3[t, 3])

# polar orbit plot
fig = plt.figure(2)
ax = plt.subplot(111, projection = 'polar')
ax.plot(sol_orb1[:, 4], rNorm1, label='tbd')
ax.plot(sol_orb2[:, 4], rNorm2, label='tbd')
ax.plot(sol_orb3[:, 4], rNorm3, label='tbd')


ax.set_xlabel('rx (m)')
ax.set_ylabel('ry (m)')
# ax.set_zlabel('rz (m)')
ax.set_title('ODE Propagation')
plt.legend()
plt.show()'''
