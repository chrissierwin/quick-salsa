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
    O = y[0]; i = y[1]; a = y[2]; w = y[3]; e = y[4]; M0 = y[5]

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

        O, i, a, w, e, M = updated_elems[t]   # grab orbital elements at current time

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
    O = y[0]; i = y[1]; a = y[2]; w = y[3]; e = y[4]; M0 = y[5]

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
a0 = 7378e3          # m, initial semi-major axis
w0 = np.deg2rad(270) # rad, initial arg. of perigee
e0 = 0.01            # initial eccentricity
f0 = np.pi/10        # rad, initial true anomaly
M0 = np.arctan2((np.sqrt(1-e0*e0*np.sin(f0))/(1+e0*np.cos(f0))), (e0+np.cos(f0)/(1+e0*np.cos(f0)))) - e0*(np.sqrt(1-e0*e0*np.sin(f0))/(1+e0*np.cos(f0)))

# update elements over time
P = 2*np.pi*np.sqrt(a0**3/mu)
y0 = np.array([O0, i0, a0, w0, e0, M0])      # initial orital elements
t0 = 0                                       # s, initial time
tf = int(8*P)                                # s, final time
t = np.linspace(t0, tf, tf)                  # s, integrate across 8 orbital periods
updated_elems = odeint(Gauss_Vari, y0, t, args = (t0,), rtol=1e-6, atol=1e-6)
'''
# propogate orbits
emp_sol = [] # empty array for solution
kep_sol = keplerian_propagation(t0, tf, mu, a0, e0, i0, O0, w0, emp_sol)
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
'''

# time varying element line fit
b0, m0 = polyfit(t, updated_elems[:,0], 1)
b1, m1 = polyfit(t, updated_elems[:,1], 1)
b2, m2 = polyfit(t, updated_elems[:,2], 1)
b3, m3 = polyfit(t, updated_elems[:,3], 1)
b4, m4 = polyfit(t, updated_elems[:,4], 1)
b5, m5 = polyfit(t, updated_elems[:,5], 1)

# calulate mean (mu) and standard deviation (std)
mu0 = np.mean(b0+m0*t)
mu1 = np.mean(b1+m1*t)
mu2 = np.mean(b2+m2*t)
mu3 = np.mean(b3+m3*t)
mu4 = np.mean(b4+m4*t)
mu5 = np.mean(b5+m5*t)

std0 = np.std(b0+m0*t)   # std of RAAN
std1 = np.std(b0+m1*t)   # std of inclination
std2 = np.std(b0+m2*t)   # std of arg. of peri.
std3 = np.std(b0+m3*t)   # std of semi-major axis
std4 = np.std(b0+m4*t)   # std of eccentricity
std5 = np.std(b0+m5*t)   # std of initial mean anomaly
'''
# time varying elements
s = 0.25
fig, ax1 = plt.subplots(2,3)
ax1[0,0].scatter(t, updated_elems[:,0], s=s)
ax1[0,0].plot(t, b0+m0*t, color='orange')
ax1[0,0].set_xlabel('Time (s)')
ax1[0,0].set_ylabel('\u03A9 (rad)')
ax1[0,0].set_title('Osculating RAAN')
ax1[0,0].text(0, 0.78538, 'std: ' + str(std0))

ax1[0,1].scatter(t, updated_elems[:,1], s=s)
ax1[0,1].plot(t, b1+m1*t, color='orange')
ax1[0,1].set_xlabel('Time (s)')
ax1[0,1].set_ylabel('i (rad)')
ax1[0,1].set_title('Osculating Inclination')
ax1[0,1].text(0, 0.523588, 'std: ' + str(std1))

ax1[0,2].scatter(t, updated_elems[:,2], s=s)
ax1[0,2].plot(t, b2+m2*t, color='orange')
ax1[0,2].set_xlabel('Time (s)')
ax1[0,2].set_ylabel('w (rad)')
ax1[0,2].set_title('Osculating Argument of Perigee')
ax1[0,2].text(0, 7377999.786, 'std: ' + str(std2))

ax1[1,0].scatter(t, updated_elems[:,3], s=s)
ax1[1,0].plot(t, b3+m3*t, color='orange')
ax1[1,0].set_xlabel('Time (s)')
ax1[1,0].set_ylabel('a (m)')
ax1[1,0].set_title('Osculating Semi-Major Axis')
ax1[1,0].text(0, -43, 'std: ' + str(std3))

ax1[1,1].scatter(t, updated_elems[:,4], s=s)
ax1[1,1].plot(t, b4+m4*t, color='orange')
ax1[1,1].set_xlabel('Time (s)')
ax1[1,1].set_ylabel('e')
ax1[1,1].set_title('Osculating Eccentricity')
ax1[1,1].text(0, 0.008, 'std: ' + str(std4))

ax1[1,2].scatter(t, updated_elems[:,5], s=s)
ax1[1,2].plot(t, b5+m5*t, color='orange')
ax1[1,2].set_xlabel('Time (s)')
ax1[1,2].set_ylabel('M0 (rad)')
ax1[1,2].set_title('Osculating Mean Anomaly')
ax1[1,2].text(0, 0.76, 'std: ' + str(std5))

plt.tight_layout()
plt.show()'''



# PROBLEM 2 =======================================
''' Using the standard deviation of the LS fit of problem 1 for each orbital element, 
generate a set of uniform sample points (of 1 standard deviation) about each of the 6 
orbital elements. Carryout the orbital propagation using the initial conditions thus 
obtained. Compute the numerical average of the elements obtained at each instant of 
time. Perform LS fits of each of the trajectories and calculate the mean and std 
(for every time step used in the integration process).'''

# calulate mean and standard deviation
# mu0 = sum(updated_elems[:,0]) / len(updated_elems[:,0])
# mu1 = sum(updated_elems[:,1]) / len(updated_elems[:,1])
# mu2 = sum(updated_elems[:,2]) / len(updated_elems[:,2])
# mu3 = sum(updated_elems[:,3]) / len(updated_elems[:,3])
# mu4 = sum(updated_elems[:,4]) / len(updated_elems[:,4])
# mu5 = sum(updated_elems[:,5]) / len(updated_elems[:,5])

# # print(mu0)
# # print(mu1)
# # print(mu2)
# # print(mu3)
# # print(mu4)
# # print(mu5)

# std0 = np.std(updated_elems[:,0])   # std of RAAN
# std1 = np.std(updated_elems[:,1])   # std of inclination
# std2 = np.std(updated_elems[:,2])   # std of arg. of peri.
# std3 = np.std(updated_elems[:,3])   # std of semi-major axis
# std4 = np.std(updated_elems[:,4])   # std of eccentricity
# std5 = np.std(updated_elems[:,5])   # std of initial mean anomaly

# calculate uniform distribution using mean and std
nd0 = uniform(mu0, std0)
nd1 = uniform(mu1, std1)
nd2 = uniform(mu2, std2)
nd3 = uniform(mu3, std3)
nd4 = uniform(mu4, std4)
nd5 = uniform(mu5, std5)

# x-values for dist axis
x0 = np.linspace(mu0-.00005, mu0+.00005, 100)
x1 = np.linspace(mu1-.00005, mu1+.00005, 200)
x2 = np.linspace(mu2-.25, mu2+.25, 200)
x3 = np.linspace(mu3+80, mu3-80, 100)
x4 = np.linspace(mu4-.0005, mu4+.0005, 200)
x5 = np.linspace(mu5-.05, mu5+.05, 1000)

# generate a set of uniform sample points
sample0 = nd0.pdf(x0)   # 1 std RAAN
sample1 = nd1.pdf(x1)   # 1 std inclination
sample2 = nd2.pdf(x2)   # 1 std arg. of peri.
sample3 = nd3.pdf(x3)   # 1 std semi-major axis
sample4 = nd4.pdf(x4)   # 1 std eccentricity
sample5 = nd5.pdf(x5)   # 1 std mean anomaly (M0)



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
plt.show()

# print(int(len(sample0)/2))
# print(sample0[int(len(sample0)/2)])
# print(nd0[0])
# new initial elements based on sample
O0 = sample0[int(len(sample0)/2)]     # rad, initial RAAN
i0 = sample1[int(len(sample1)/2)]     # rad, initial inclination
a0 = sample2[int(len(sample2)/2)]     # m, initial semi-major axis
w0 = sample3[int(len(sample3)/2)]     # rad, initial arg. of perigee
e0 = sample4[int(len(sample4)/2)]     # initial eccentricity
f0 = sample5[int(len(sample5)/2)]     # rad, initial true anomaly
M0 = np.arctan2((np.sqrt(1-e0*e0*np.sin(f0))/(1+e0*np.cos(f0))), (e0+np.cos(f0)/(1+e0*np.cos(f0)))) - e0*(np.sqrt(1-e0*e0*np.sin(f0))/(1+e0*np.cos(f0)))

# NORMALIZE VALUES TO BE WITHIN NORMAL PARAMETERS FOR INITIAL VALUES!!!
print(O0); print(i0); print(a0); print(w0); print(e0); print(f0)

'''
# update elements over time
P = 2*np.pi*np.sqrt(a0**3/mu)
y0 = np.array([O0, i0, a0, w0, e0, M0])      # initial orital elements
t0 = 0                                       # s, initial time
tf = int(8*P)                                # s, final time
t = np.linspace(t0, tf, tf)                  # s, integrate across 8 orbital periods
updated_elems = odeint(Gauss_Vari, y0, t, args = (t0,), rtol=1e-6, atol=1e-6)

# propogate orbits
emp_sol = [] # empty array for solution
kep_sol = keplerian_propagation(t0, tf, mu, a0, e0, i0, O0, w0, emp_sol)
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
ax.set_title('Orbital Propagation using +\u03C3')
plt.legend()'''
# plt.show()