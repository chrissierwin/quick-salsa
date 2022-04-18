# A624 HW 3, due April 21, 2022

import numpy as np
from scipy.integrate import odeint
from numpy.polynomial.polynomial import polyfit
from mpl_toolkits import mplot3d 
import matplotlib.pyplot as plt
from scipy.stats import uniform #norm

# integrate Gauss' Variational eqns to update orbit
def Gauss_Vari(y, t, t0):

    # grab inital/updated variables
    O = y[0]; i = y[1]; w = y[2]; a = y[3]; e = y[4]; M0 = y[5]

    # calculate true anomaly
    mu = 3.986004415e14                   # universal gravitational parameter, m^3/s^2
    n = np.sqrt(mu/(a**3))                # 1/s, mean motion
    M = M0 + n*(t-t0)
    f = M + (2*e - (e**3)/4)*np.sin(M) + (5*(e**2) * np.sin(2*M))/4 + (13*(e**3) * np.sin(3*M))/12

    # constants / givens
    p = a*(1 - e*e)                       # m, semi-latus rectum
    r = p/(1 + e*np.cos(f))               # m, radius
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

    yDot = np.array([dOdt, didt, dwdt, dadt, dedt, dM0dt])

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

    rVector = np.array([rx, ry, rz]) # km
    vVector = np.array([vx, vy, vz]) # km/s
    rvVector = np.array([rx, ry, rz, vx, vy, vz])

    return rvVector

# propogate orbit using Keplerian parameters
def keplerian_propagation(t0, tf, mu, a, e, i, O, w, sol_array):

    '''EDIT***takes t0: initial time (s), tf: final time (s), mu: gravitational constant (km^3/s^2), a: semi-major axis (km),
       e: eccentricity, i: inclination (rad), w: argument of perigee (rad), RAAN: (rad), sol_array: empty array. 
       Returns keplerian propagation in sol_array.'''

    for t in np.arange(t0, tf):

        O, i, w, a, e, M0 = updated_elems[t]   # grab orbital elements at current time
        # calculate true anomaly
        mu = 3.986004415e14                   # universal gravitational parameter, m^3/s^2
        n = np.sqrt(mu/(a**3))                # 1/s, mean motion
        M = M0 + n*(t-t0)
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

# integrate mean oscullation eqns to update orbit
def mean_osculation(y, t):

    # grab inital/updated variables
    O = y[0]; i = y[1]; w = y[2]; a = y[3]; e = y[4]; M0 = y[5]

    # constants / givens
    mu = 3.986004415e14                   # universal gravitational parameter, m^3/s^2
    n = np.sqrt(mu/(a**3))                # 1/s, mean motion
    p = a*(1 - e*e)                       # m, semi-latus rectum
    eta = np.sqrt(1 - e*e)                # unitless parameter
    r_eq = 6378.1370e3                    # m, e mean equatorial radius of the Earth
    J2 = 1082.63e-6                       # variation in Earth's Oblateness

    # Gauss' variational equations
    dOdt = -3*J2/2 *n*(r_eq/p)**2 * np.cos(i)
    didt = 0
    dwdt =  3*J2/4 *n*(r_eq/p)**2 * (5*np.cos(i)**2 - 1)
    dadt = 0
    dedt = 0
    dM0dt = 3*J2/4 *n*(r_eq/p)**2 * eta*(3*np.cos(i)**2 - 1)

    yDot = np.array([dOdt, didt, dwdt, dadt, dedt, dM0dt])

    return yDot


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
w0 = np.deg2rad(270) # rad, initial arg. of perigee
a0 = 7378e3          # m, initial semi-major axis
e0 = 0.01            # initial eccentricity
f0 = np.pi/10        # rad, initial true anomaly
M0 = np.arctan2((np.sqrt(1-e0*e0*np.sin(f0))/(1+e0*np.cos(f0))), (e0+np.cos(f0)/(1+e0*np.cos(f0)))) - e0*(np.sqrt(1-e0*e0*np.sin(f0))/(1+e0*np.cos(f0)))

# update elements over time
P = 2*np.pi*np.sqrt(a0**3/mu)
y0 = np.array([O0, i0, w0, a0, e0, M0])      # initial orital elements
t0 = 0                                       # s, initial time
tf = int(8*P)                                # s, final time
t = np.linspace(t0, tf, tf)                  # s, integrate across 8 orbital periods
updated_elems = odeint(Gauss_Vari, y0, t, args = (t0,), rtol=1e-6, atol=1e-6)

# propogate orbits
emp_sol = [] # empty array for solution
kep_sol = keplerian_propagation(t0, tf, mu, a0, e0, i0, O0, w0, emp_sol)
kep_sol = np.reshape(kep_sol, (len(kep_sol), 6))
'''
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
'''

# time varying element line fit
b0, m0 = polyfit(t, updated_elems[:,0], 1)
b1, m1 = polyfit(t, updated_elems[:,1], 1)
b2, m2 = polyfit(t, updated_elems[:,2], 1)
b3, m3 = polyfit(t, updated_elems[:,3], 1)
b4, m4 = polyfit(t, updated_elems[:,4], 1)
b5, m5 = polyfit(t, updated_elems[:,5], 1)

line0 = b0+m0*t
line1 = b1+m1*t
line2 = b2+m2*t
line3 = b3+m3*t
line4 = b4+m4*t
line5 = b5+m5*t

# calulate mean (mu) and standard deviation (std)
mu0 = np.mean(line0)
mu1 = np.mean(line1)
mu2 = np.mean(line2)
mu3 = np.mean(line3)
mu4 = np.mean(line4)
mu5 = np.mean(line5)

std0 = np.std(line0)   # std of RAAN
std1 = np.std(line1)   # std of inclination
std2 = np.std(line2)   # std of arg. of peri.
std3 = np.std(line3)   # std of semi-major axis
std4 = np.std(line4)   # std of eccentricity
std5 = np.std(line5)   # std of initial mean anomaly

# time varying elements
# s = 0.25
fig, ax1 = plt.subplots(2,3)
ax1[0,0].plot(t, updated_elems[:,0])#, s=s)
ax1[0,0].plot(t, line0, color='orange')
ax1[0,0].set_xlabel('Time (s)')
ax1[0,0].set_ylabel('\u03A9 (rad)')
ax1[0,0].set_title('Osculating RAAN')
ax1[0,0].text(0, 0.78538, 'mean: ' + str(mu0))
ax1[0,0].text(0, 0.783, 'std: ' + str(std0))

ax1[0,1].plot(t, updated_elems[:,1])#, s=s)
ax1[0,1].plot(t, line1, color='orange')
ax1[0,1].set_xlabel('Time (s)')
ax1[0,1].set_ylabel('i (rad)')
ax1[0,1].set_title('Osculating Inclination')
ax1[0,1].text(0, 0.523845, 'mean: ' + str(mu1))
ax1[0,1].text(0, 0.52382, 'std: ' + str(std1))

ax1[0,2].plot(t, updated_elems[:,2])#, s=s)
ax1[0,2].plot(t, line2, color='orange')
ax1[0,2].set_xlabel('Time (s)')
ax1[0,2].set_ylabel('w (rad)')
ax1[0,2].set_title('Osculating Argument of Perigee')
ax1[0,2].text(0, 4.844, 'mean: ' + str(mu2))
ax1[0,2].text(0, 4.83, 'std: ' + str(std2))

ax1[1,0].plot(t, updated_elems[:,3])#, s=s)
ax1[1,0].plot(t, line3, color='orange')
ax1[1,0].set_xlabel('Time (s)')
ax1[1,0].set_ylabel('a (m)')
ax1[1,0].set_title('Osculating Semi-Major Axis')
ax1[1,0].text(0, 7.38e6, 'mean: ' + str(mu3))
ax1[1,0].text(0, 7.3798e6, 'std: ' + str(std3))

ax1[1,1].plot(t, updated_elems[:,4])#, s=s)
ax1[1,1].plot(t, line4, color='orange')
ax1[1,1].set_xlabel('Time (s)')
ax1[1,1].set_ylabel('e')
ax1[1,1].set_title('Osculating Eccentricity')
ax1[1,1].text(0, 0.010057, 'mean: ' + str(mu4))
ax1[1,1].text(0, 0.0100, 'std: ' + str(std4))

ax1[1,2].plot(t, updated_elems[:,5])#, s=s)
ax1[1,2].plot(t, line5, color='orange')
ax1[1,2].set_xlabel('Time (s)')
ax1[1,2].set_ylabel('M0 (rad)')
ax1[1,2].set_title('Osculating Mean Anomaly')
ax1[1,2].text(0, 1.0144, 'mean: ' + str(mu5))
ax1[1,2].text(0, 1.001, 'std: ' + str(std5))

plt.tight_layout()
plt.show()

# find the x-intersection of perturbations and the mean
idx0 = np.argwhere(np.diff(np.sign(updated_elems[:,0] - line0))).flatten()
idx1 = np.argwhere(np.diff(np.sign(updated_elems[:,1] - line1))).flatten()
idx2 = np.argwhere(np.diff(np.sign(updated_elems[:,2] - line2))).flatten()
idx3 = np.argwhere(np.diff(np.sign(updated_elems[:,3] - line3))).flatten()
idx4 = np.argwhere(np.diff(np.sign(updated_elems[:,4] - line4))).flatten()
idx5 = np.argwhere(np.diff(np.sign(updated_elems[:,5] - line5))).flatten()
# now, locate the y-intersection value
y_int0 = updated_elems[idx0[0],0]
y_int1 = updated_elems[idx1[0],1]
y_int2 = updated_elems[idx2[0],2]
y_int3 = updated_elems[idx3[0],3]
y_int4 = updated_elems[idx4[0],4]
y_int5 = updated_elems[idx5[0],5]


# PROBLEM 2 =======================================
''' Using the standard deviation of the LS fit of problem 1 for each orbital element, 
generate a set of uniform sample points (of 1 standard deviation) about each of the 6 
orbital elements. Carryout the orbital propagation using the initial conditions thus 
obtained. Compute the numerical average of the elements obtained at each instant of 
time. Perform LS fits of each of the trajectories and calculate the mean and std 
(for every time step used in the integration process).'''

# calculate uniform distribution using mean and std
nd0 = uniform(mu0-std0, mu0+std0)
nd1 = uniform(mu1-std1, mu1+std1)
nd2 = uniform(mu2-std2, mu2+std2)
nd3 = uniform(mu3-std3, mu3+std3)
nd4 = uniform(mu4-std4, mu4+std4)
nd5 = uniform(mu5-std5, mu5+std5)

# x-values for dist axis
x0 = np.linspace(mu0-std0-5, mu0+std0+5, 100)
x1 = np.linspace(mu1-std1-5, mu1+std1+5, 100)
x2 = np.linspace(mu2-std2-10, mu2+std2+10, 100)
x3 = np.linspace(mu3-std3-100, mu3+std3+100, 100)
x4 = np.linspace(mu4-std4-5, mu4+std4+5, 100)
x5 = np.linspace(mu5-std5-5, mu5+std5+5, 100)

# generate a set of uniform sample points
sample0 = nd0.pdf(x0)   # 1 std RAAN
sample1 = nd1.pdf(x1)   # 1 std inclination
sample2 = nd2.pdf(x2)   # 1 std arg. of peri.
sample3 = nd3.pdf(x3)   # 1 std semi-major axis
sample4 = nd4.pdf(x4)   # 1 std eccentricity
sample5 = nd5.pdf(x5)   # 1 std mean anomaly (M0)

'''
# plot distribution
fig, ax2 = plt.subplots(2,3)
ax2[0,0].plot(x0, sample0)
ax2[0,0].set_xlabel('x values about mean')
ax2[0,0].set_ylabel('Probability')
ax2[0,0].set_title('RAAN Distribution')

ax2[0,1].plot(x1, sample1)
ax2[0,1].set_xlabel('x values about mean')
ax2[0,1].set_ylabel('Probability')
ax2[0,1].set_title('Inclination Distribution')

ax2[0,2].plot(x2, sample2) 
ax2[0,2].set_xlabel('x values about mean')
ax2[0,2].set_ylabel('Probability')
ax2[0,2].set_title('Arg. of Perigee Distribution')

ax2[1,0].plot(x3, sample3)
ax2[1,0].set_xlabel('x values about mean')
ax2[1,0].set_ylabel('Probability')
ax2[1,0].set_title('Semi-Major Axis Distribution')

ax2[1,1].plot(x4, sample4)
ax2[1,1].set_xlabel('x values about mean')
ax2[1,1].set_ylabel('Probability')
ax2[1,1].set_title('Eccentricity Distribution')

ax2[1,2].plot(x5, sample5)
ax2[1,2].set_xlabel('x values about mean')
ax2[1,2].set_ylabel('Probability')
ax2[1,2].set_title('Mean Anomaly Distribution')

plt.tight_layout()
plt.show()'''

# PLUS ONE STD
# new initial elements using old initials + 1std
O0_plu = O0 + sample0[int(len(sample0)/2)]     # rad, initial RAAN
i0_plu = i0 + sample1[int(len(sample1)/2)]     # rad, initial inclination
w0_plu = w0 + sample3[int(len(sample3)/2)]     # rad, initial arg. of perigee
a0_plu = a0 + sample2[int(len(sample2)/2)]     # m, initial semi-major axis
e0_plu = e0 + sample4[int(len(sample4)/2)]     # initial eccentricity
M0_plu = M0 + sample5[int(len(sample5)/2)]     # rad, initial true anomaly

# update elements over time
P_plu = 2*np.pi*np.sqrt(a0_plu**3/mu)                                 # orbital period
y0_plu = np.array([O0_plu, i0_plu, w0_plu, a0_plu, e0_plu, M0_plu])   # initial orital elements
t0 = 0                                                                # s, initial time
tf_plu = int(8*P_plu)                                                 # s, final time
t_plu = np.linspace(t0, tf_plu, tf_plu)                               # s, integrate across 8 orbital periods
updated_elems_plu = odeint(Gauss_Vari, y0_plu, t_plu, args = (t0,), rtol=1e-6, atol=1e-6)

# propogate orbits
emp_sol_plu = [] # empty array for solution
kep_sol_plu = keplerian_propagation(t0, tf, mu, a0_plu, e0_plu, i0_plu, O0_plu, w0_plu, emp_sol_plu)
kep_sol_plu = np.reshape(kep_sol_plu, (len(kep_sol_plu), 6))
'''
# plot orbital position
fig = plt.figure('Orbital Propagation + 1std')
ax = plt.axes(projection='3d')
ax.plot3D(kep_sol_plu[:, 0], kep_sol_plu[:,1], kep_sol_plu[:,2], label='Orbit')

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = rEarth*np.cos(u)*np.sin(v)
y = rEarth*np.sin(u)*np.sin(v)
z = rEarth*np.cos(v)
ax.plot_wireframe(x, y, z, color="g", label='Earth')

ax.set_xlabel('rx (m)')
ax.set_ylabel('ry (m)')
ax.set_zlabel('rz (m)')
ax.set_title('Orbital Propagation using + 1 \u03C3')
plt.legend()
plt.show()'''


# MINUS ONE STD
# new initial elements using old initials + 1std
O0_min = O0 - sample0[int(len(sample0)/2)]     # rad, initial RAAN
i0_min = i0 - sample1[int(len(sample1)/2)]     # rad, initial inclination
w0_min = w0 - sample3[int(len(sample3)/2)]     # rad, initial arg. of perigee
a0_min = a0 - sample2[int(len(sample2)/2)]     # m, initial semi-major axis
e0_min = e0 - sample4[int(len(sample4)/2)]     # initial eccentricity
M0_min = M0 - sample5[int(len(sample5)/2)]     # rad, initial true anomaly

# update elements over time
P_min = 2*np.pi*np.sqrt(a0_min**3/mu)                                 # orbital period
y0_min = np.array([O0_min, i0_min, w0_min, a0_min, e0_min, M0_min])   # initial orital elements
t0 = 0                                                                # s, initial time
tf_min = int(8*P_min)                                                 # s, final time
t_min = np.linspace(t0, tf_min, tf_min)                               # s, integrate across 8 orbital periods
updated_elems_min = odeint(Gauss_Vari, y0_min, t_min, args = (t0,), rtol=1e-6, atol=1e-6)

# propogate orbits
emp_sol_min = [] # empty array for solution
kep_sol_min = keplerian_propagation(t0, tf, mu, a0_min, e0_min, i0_min, O0_min, w0_min, emp_sol_min)
kep_sol_min = np.reshape(kep_sol_min, (len(kep_sol_min), 6))
'''
# plot orbital position
fig = plt.figure('Orbital Propagation - 1std')
ax = plt.axes(projection='3d')
ax.plot3D(kep_sol_min[:, 0], kep_sol_min[:,1], kep_sol_min[:,2], label='Orbit')

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = rEarth*np.cos(u)*np.sin(v)
y = rEarth*np.sin(u)*np.sin(v)
z = rEarth*np.cos(v)
ax.plot_wireframe(x, y, z, color="g", label='Earth')

ax.set_xlabel('rx (m)')
ax.set_ylabel('ry (m)')
ax.set_zlabel('rz (m)')
ax.set_title('Orbital Propagation using - 1 \u03C3')
plt.legend()
plt.show()'''
'''
# PLUS 1STD CASE
# time varying element line fit (numerical average) of +1std
b0_plu, m0_plu = polyfit(t_plu, updated_elems_plu[:,0], 1)
b1_plu, m1_plu = polyfit(t_plu, updated_elems_plu[:,1], 1)
b2_plu, m2_plu = polyfit(t_plu, updated_elems_plu[:,2], 1)
b3_plu, m3_plu = polyfit(t_plu, updated_elems_plu[:,3], 1)
b4_plu, m4_plu = polyfit(t_plu, updated_elems_plu[:,4], 1)
b5_plu, m5_plu = polyfit(t_plu, updated_elems_plu[:,5], 1)

# calulate mean (mu) and standard deviation (std) for +1std
mu0_plu = np.mean(b0_plu+m0_plu*t_plu)
mu1_plu = np.mean(b1_plu+m1_plu*t_plu)
mu2_plu = np.mean(b2_plu+m2_plu*t_plu)
mu3_plu = np.mean(b3_plu+m3_plu*t_plu)
mu4_plu = np.mean(b4_plu+m4_plu*t_plu)
mu5_plu = np.mean(b5_plu+m5_plu*t_plu)

std0_plu = np.std(b0_plu+m0_plu*t_plu)   # std of RAAN
std1_plu = np.std(b1_plu+m1_plu*t_plu)   # std of inclination
std2_plu = np.std(b2_plu+m2_plu*t_plu)   # std of arg. of peri.
std3_plu = np.std(b3_plu+m3_plu*t_plu)   # std of semi-major axis
std4_plu = np.std(b4_plu+m4_plu*t_plu)   # std of eccentricity
std5_plu = np.std(b5_plu+m5_plu*t_plu)   # std of initial mean anomaly
'''
'''
# time varying elements
s = 0.25
fig, ax3 = plt.subplots(2,3)
ax3[0,0].scatter(t_plu, updated_elems_plu[:,0], s=s)
ax3[0,0].plot(t_plu, b0_plu+m0_plu*t_plu, color='orange')
ax3[0,0].set_xlabel('Time (s)')
ax3[0,0].set_ylabel('\u03A9 (rad)')
ax3[0,0].set_title('Osculating RAAN')
ax3[0,0].text(0, 2.123, 'mean: ' + str(mu0_plu))
ax3[0,0].text(0, 2.12, 'std: ' + str(std0_plu))

ax3[0,1].scatter(t_plu, updated_elems_plu[:,1], s=s)
ax3[0,1].plot(t_plu, b1_plu+m1_plu*t_plu, color='orange')
ax3[0,1].set_xlabel('Time (s)')
ax3[0,1].set_ylabel('i (rad)')
ax3[0,1].set_title('Osculating Inclination')
ax3[0,1].text(0, 2.43399, 'mean: ' + str(mu1_plu))
ax3[0,1].text(0, 2.43396, 'std: ' + str(std1_plu))

ax3[0,2].scatter(t_plu, updated_elems_plu[:,2], s=s)
ax3[0,2].plot(t_plu, b2_plu+m2_plu*t_plu, color='orange')
ax3[0,2].set_xlabel('Time (s)')
ax3[0,2].set_ylabel('w (rad)')
ax3[0,2].set_title('Osculating Argument of Perigee')
ax3[0,2].text(0, 4.792, 'mean: ' + str(mu2_plu))
ax3[0,2].text(0, 4.782, 'std: ' + str(std2_plu))

ax3[1,0].scatter(t_plu, updated_elems_plu[:,3], s=s)
ax3[1,0].plot(t_plu, b3_plu+m3_plu*t_plu, color='orange')
ax3[1,0].set_xlabel('Time (s)')
ax3[1,0].set_ylabel('a (m)')
ax3[1,0].set_title('Osculating Semi-Major Axis')
ax3[1,0].text(0, 7.379e6, 'mean: ' + str(mu3_plu))
ax3[1,0].text(0, 7.3787e6, 'std: ' + str(std3_plu))

ax3[1,1].scatter(t_plu, updated_elems_plu[:,4], s=s)
ax3[1,1].plot(t_plu, b4_plu+m4_plu*t_plu, color='orange')
ax3[1,1].set_xlabel('Time (s)')
ax3[1,1].set_ylabel('e')
ax3[1,1].set_title('Osculating Eccentricity')
ax3[1,1].text(0, 0.0109, 'mean: ' + str(mu4_plu))
ax3[1,1].text(0, 0.01085, 'std: ' + str(std4_plu))

ax3[1,2].scatter(t_plu, updated_elems_plu[:,5], s=s)
ax3[1,2].plot(t_plu, b5_plu+m5_plu*t_plu, color='orange')
ax3[1,2].set_xlabel('Time (s)')
ax3[1,2].set_ylabel('M0 (rad)')
ax3[1,2].set_title('Osculating Mean Anomaly')
ax3[1,2].text(0, 2.128, 'mean: ' + str(mu5_plu))
ax3[1,2].text(0, 2.12, 'std: ' + str(std5_plu))

plt.tight_layout()
plt.show()'''


# MINUS 1STD CASE
# time varying element line fit (numerical average) of -1std
b0_min, m0_min = polyfit(t_min, updated_elems_min[:,0], 1)
b1_min, m1_min = polyfit(t_min, updated_elems_min[:,1], 1)
b2_min, m2_min = polyfit(t_min, updated_elems_min[:,2], 1)
b3_min, m3_min = polyfit(t_min, updated_elems_min[:,3], 1)
b4_min, m4_min = polyfit(t_min, updated_elems_min[:,4], 1)
b5_min, m5_min = polyfit(t_min, updated_elems_min[:,5], 1)

# calulate mean (mu) and standard deviation (std) for -1std
mu0_min = np.mean(b0_min+m0_min*t_min)
mu1_min = np.mean(b1_min+m1_min*t_min)
mu2_min = np.mean(b2_min+m2_min*t_min)
mu3_min = np.mean(b3_min+m3_min*t_min)
mu4_min = np.mean(b4_min+m4_min*t_min)
mu5_min = np.mean(b5_min+m5_min*t_min)

std0_min = np.std(b0_min+m0_min*t_min)   # std of RAAN
std1_min = np.std(b1_min+m1_min*t_min)   # std of inclination
std2_min = np.std(b2_min+m2_min*t_min)   # std of arg. of peri.
std3_min = np.std(b3_min+m3_min*t_min)   # std of semi-major axis
std4_min = np.std(b4_min+m4_min*t_min)   # std of eccentricity
std5_min = np.std(b5_min+m5_min*t_min)   # std of initial mean anomaly
'''
# time varying elements
s = 0.25
fig, ax4 = plt.subplots(2,3)
ax4[0,0].scatter(t_min, updated_elems_min[:,0], s=s)
ax4[0,0].plot(t_min, b0_min+m0_min*t_min, color='orange')
ax4[0,0].set_xlabel('Time (s)')
ax4[0,0].set_ylabel('\u03A9 (rad)')
ax4[0,0].set_title('Osculating RAAN')
ax4[0,0].text(0, -.5051, 'mean: ' + str(mu0_min))
ax4[0,0].text(0, -.5058, 'std: ' + str(std0_min))

ax4[0,1].scatter(t_min, updated_elems_min[:,1], s=s)
ax4[0,1].plot(t_min, b1_min+m1_min*t_min, color='orange')
ax4[0,1].set_xlabel('Time (s)')
ax4[0,1].set_ylabel('i (rad)')
ax4[0,1].set_title('Osculating Inclination')
ax4[0,1].text(0, -1.386247, 'mean: ' + str(mu1_min))
ax4[0,1].text(0, -1.386257, 'std: ' + str(std1_min))

ax4[0,2].scatter(t_min, updated_elems_min[:,2], s=s)
ax4[0,2].plot(t_min, b2_min+m2_min*t_min, color='orange')
ax4[0,2].set_xlabel('Time (s)')
ax4[0,2].set_ylabel('w (rad)')
ax4[0,2].set_title('Osculating Argument of Perigee')
ax4[0,2].text(0, 4.715, 'mean: ' + str(mu2_min))
ax4[0,2].text(0, 4.705, 'std: ' + str(std2_min))

ax4[1,0].scatter(t_min, updated_elems_min[:,3], s=s)
ax4[1,0].plot(t_min, b3_min+m3_min*t_min, color='orange')
ax4[1,0].set_xlabel('Time (s)')
ax4[1,0].set_ylabel('a (m)')
ax4[1,0].set_title('Osculating Semi-Major Axis')
ax4[1,0].text(0, 7.3933e6, 'mean: ' + str(mu3_min))
ax4[1,0].text(0, 7.3925e6, 'std: ' + str(std3_min))

ax4[1,1].scatter(t_min, updated_elems_min[:,4], s=s)
ax4[1,1].plot(t_min, b4_min+m4_min*t_min, color='orange')
ax4[1,1].set_xlabel('Time (s)')
ax4[1,1].set_ylabel('e')
ax4[1,1].set_title('Osculating Eccentricity')
ax4[1,1].text(0, 0.0126, 'mean: ' + str(mu4_min))
ax4[1,1].text(0, 0.0124, 'std: ' + str(std4_min))

ax4[1,2].scatter(t_min, updated_elems_min[:,5], s=s)
ax4[1,2].plot(t_min, b5_min+m5_min*t_min, color='orange')
ax4[1,2].set_xlabel('Time (s)')
ax4[1,2].set_ylabel('M0 (rad)')
ax4[1,2].set_title('Osculating Mean Anomaly')
ax4[1,2].text(0, -.18, 'mean: ' + str(mu5_min))
ax4[1,2].text(0, -.195, 'std: ' + str(std5_min))

plt.tight_layout()
plt.show()'''


# PROBLEM 3 =======================================
'''Pick a particular point in the time axis of each mean to osculating element plots 
where the mean and osculating elements coincide. Using this point, propagate the averaged 
elements provided by the equations.'''
# following eqn is available for 0std
# idx = np.argwhere(np.diff(np.sign(f - g))).flatten()


# givens
O0_mean = b0     # rad, initial RAAN
i0_mean = b1     # rad, initial inclination
w0_mean = b2     # rad, initial arg. of perigee
a0_mean = b3     # m, initial semi-major axis
e0_mean = b4     # initial eccentricity
M0_mean = b5     # rad, initial mean anomaly

# update elements over time
P_mean = 2*np.pi*np.sqrt(a0_mean**3/mu)
y0_mean = np.array([O0_mean, i0_mean, w0_mean, a0_mean, e0_mean, M0_mean])      # initial orital elements
t0 = 0                                                                          # s, initial time
tf_mean = int(8*P)                                                              # s, final time
t_mean = np.linspace(t0, tf_mean, tf_mean)                                      # s, integrate across 8 orbital periods
updated_elems_mean = odeint(mean_osculation, y0_mean, t_mean, rtol=1e-6, atol=1e-6)

# propogate orbits
emp_sol_mean = [] # empty array for solution
kep_sol_mean = keplerian_propagation(t0, tf_mean, mu, a0_mean, e0_mean, i0_mean, O0_mean, w0_mean, emp_sol_mean)
kep_sol_mean = np.reshape(kep_sol_mean, (len(kep_sol_mean), 6))

# plot orbital position
fig = plt.figure('Mean Oscullating Orbital Propagation')
ax = plt.axes(projection='3d')
ax.plot3D(kep_sol_mean[:, 0], kep_sol_mean[:,1], kep_sol_mean[:,2], label='Orbit')

# draw sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = rEarth*np.cos(u)*np.sin(v)
y = rEarth*np.sin(u)*np.sin(v)
z = rEarth*np.cos(v)
ax.plot_wireframe(x, y, z, color="g", label='Earth')

ax.set_xlabel('rx (m)')
ax.set_ylabel('ry (m)')
ax.set_zlabel('rz (m)')
ax.set_title('Mean Oscullating Orbital Propagation')
plt.legend()
plt.show()

# IC QUESTIONS:
# Q3: set mean_RAAN(0)= b (if y=mx+b from LS fit)
# ^ if so, how do we add in the use of some t1 to find the mean / account for missing time?
# Q2: Compute the numerical average of the elements obtained at each instant of time. Is this a LS fit or oscillatory avg of 0 and +-1 std?