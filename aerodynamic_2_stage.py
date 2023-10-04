"""
second_stage_rocket.py

Created on Tues Oct 3 22:43:36 2023

Edited by Mark Graham

Original work done by Carlos Montalvo
"""

#### Import used modules
# %%
import numpy as np                  ### numeric python
import matplotlib.pyplot as mpl     ### matlab style plotting
import scipy.integrate as sci       ### integration toolbox

### Defining Constant Parameters

## Gravitational Constant
G = 6.6742*10**-11

## Planet Constants
# Earth
#Rplanet = 6357000.0                 ### meters
#mplanet = 5.972e24                  ### kilograms
#name = 'Earth'

#Kerbin
Rplanet = 600000                    ### meters
mplanet = 5.2915158*10**22          ### kilograms
name = 'Kerbin'

## Rocket 
weighttons = 5.3                    ### tons
mass0 = weighttons * 2000 / 2.2     ### kilograms
max_thrust = 167970.0               ### Newtons
Isp1 = 250.0                        ### seconds
Isp2 = 400.0                        ### seconds
tMECO = 38.0                        ### seconds
tSep1 = 2.0                         ### seconds
tonsMass1 = 0.2                     ### ton
tMass1 = tonsMass1 * 2000 / 2.2     ### kilograms
t2start = 261.0                     ### seconds
t2end = t2start + 17.5              ### seconds
x0 = Rplanet                        ### meters
z0 = 0.0                            ### meters
velz0 = 0.0                         ### meters / second
velx0 = 0.0                         ### meters / second
r0 = 200000.0 + Rplanet
period = 12000 #2*np.pi/np.sqrt(G*mplanet)*r0**(3.0/2.0) *1.5
D = 0.85
#S = np.pi*(D/2.0)**2
CD = 0.1

### Aerodynamics Class
class Aerodynamics():
    def __init__(self, name):
        ### Import Aero Model
        self.name = name
        if name == 'Kerbin':
            ### use kerbin aero model
            data = np.loadtxt('kerbin_aerodynamics.txt')
            self.altitude = data[:,0]
            self.density = data[:,3]
            self.rhos = self.density[0]
        elif name == 'Earth':
            ### use earth aero model
            self.beta = 0.1354/1000.0
            self.rhos = 1.225
    def getDensity(self, altitude):
        if self.name == 'Kerbin':
            rho = np.interp(altitude, self.altitude, self.density)
        elif self.name == 'Earth':
            rho = self.rhos * np.exp(-self.beta * altitude)
        return rho

### Create a Global Aerodynamics instance.
aeroModel = Aerodynamics(name)

## Gravitational Acceleration
def gravity(x, z):
    global Rplanet, mplanet

    r = np.sqrt(x**2 + z**2)

    if r < Rplanet:
            accelx = 0.0
            accelz = 0.0
    else:
        accelx = G*mplanet / (r**3)*x
        accelz = G*mplanet / (r**3)*z
    
    return np.array([accelx, accelz]), r

def propulsion(t):
    global max_thrust, Isp, tMECO, ve
    if t < tMECO:
        theta = 10.0 * np.pi / 180.0

        ### Firing the main thruster
        thrustF = max_thrust

        ### Compute Exit Velocity
        ve = Isp1 * 9.81            ### meters / second

        ### mdot
        mdot = -thrustF/ve
    if t > tMECO and t < (tMECO + tSep1):
        theta = 0.0
        thrustF = 0.0
        mdot = -tMass1 / tSep1
    if t > (tMECO + tSep1):
        theta = 0.0
        thrustF = 0.0
        mdot = 0.0
    if t > (t2start) and t < (t2end): 
        theta = 90 * np.pi / 180.0
        thrustF = max_thrust
        ve = Isp2 * 9.81            ### meters / second
        mdot = -thrustF / ve
    if t > t2end:
        theta = 0.0
        thrustF = 0.0
        mdot = 0.0

    thrustx = thrustF * np.cos(theta)
    thrustz = thrustF * np.sin(theta)

    

    return np.asarray([thrustx, thrustz]), mdot

### Equations of motion:
## Force = mass * acceleration = mass * zddot
# z = altitude from the center of the planet along the north pole
# x = alititude from center along equator through Africa
# zdot = velocity along z
# xdot = velocity along x
# zddot = acceleration along z
# xddot = acceleration along x

## Second Order Differential Equation 
def Derivatives(state, t):
    global aeroModel, Rplanet
    #state vector
    x = state[0]
    z = state[1]
    velx = state[2]
    velz = state[3]
    mass = state[4]

    # Compute zdot
    zdot = velz
    xdot = velx

    ### Total Forces:
    ## Gravity
    accel, r = gravity(x, z)
    gravityF = -accel * mass

    ## Aerodynamics
    altitude = r - Rplanet
    rho = aeroModel.getDensity(altitude)
    V = np.sqrt(velz**2 + velx**2)
    qinf = (np.pi / 8.0) * rho * (D**2)  * abs(V)
    aeroF = -qinf * CD * np.asarray([velx, velz])

    ## Thrust
    thrustF, mdot = propulsion(t)
         
    Forces = gravityF + aeroF + thrustF

    # Compute Acceleration
    if mass > 0:
        ddot = Forces/mass
    else:
        ddot = 0.0
        mdot = 0.0 

    # Compute the statedot
    statedot = np.asarray([xdot, zdot, ddot[0], ddot[1], mdot])

    return statedot


##### Main Script Below #####

### Test Surface Gravity
# print(f'Surface Gravity (m/s^2) = {gravity(0, Rplanet)}')

### Populate Intial Condition Vector
stateinitial = np.asarray([x0, z0, velx0, velz0, mass0])

## Time window 
""" period = 2*np.pi/np.sqrt(G*mplanet)*r0**(3.0/2.0) *1.5 """
# Over 30 secs, give 1000 data points
tout = np.linspace(0, period, 1000)

### Numerical Integration Call
stateout = sci.odeint(Derivatives, stateinitial, tout)

### Rename Variables
xout = stateout[:,0]
zout = stateout[:,1]
altitude = np.sqrt(xout**2 + zout**2) - Rplanet
velxout = stateout[:,2]
velzout = stateout[:,3]
velout = np.sqrt(velxout**2 + velzout**2)
massout = stateout[:,4]

### Plots:

## Air Density
test_altitude = np.linspace(0, 100000, 100)
test_rho = aeroModel.getDensity(test_altitude)
mpl.figure()
mpl.plot(test_altitude, test_rho, 'b-')
mpl.xlabel('Altitude (meters [m])')
mpl.ylabel('Air Density (kg / m ^3)')
mpl.grid()
# mpl.savefig('air_density.png')

## Altitude
mpl.figure() 
mpl.plot(tout, altitude)
mpl.xlabel('Time (seconds [sec])')
mpl.ylabel('Altitude (meters [m])')
mpl.grid()
# mpl.savefig('orbit.png')

## Velocity
mpl.figure()
mpl.plot(tout, velout)
mpl.xlabel('Time (seconds [sec])')
mpl.ylabel('Total Speed (meters per second [m/s])')
mpl.grid()
# mpl.savefig('velocity.png')

## Mass 
mpl.figure()
mpl.plot(tout, massout)
mpl.xlabel('Time (seconds [sec])')
mpl.ylabel('Mass (kilograms [kgs])')
mpl.grid()
# mpl.savefig('mass.png')

## 2D Orbit
mpl.figure()
mpl.plot(xout, zout, 'r-', label='Orbit')
mpl.plot(xout[0], zout[0], 'g*')
theta = np.linspace(0, 2 * np.pi, 1000)
xplanet = Rplanet*np.sin(theta)
yplanet = Rplanet*np.cos(theta)
mpl.plot(xplanet, yplanet, 'b-', label='Planet')
mpl.grid()
mpl.legend()
# mpl.savefig('orbit.png')
# %%
