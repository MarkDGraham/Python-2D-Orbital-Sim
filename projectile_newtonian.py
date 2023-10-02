"""
projectile_newtonian.py

Created on Mon Oct 2 11:42:21 2023

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
Rplanet = 6357000.0                 ### meters
mplanet = 5.972e24                  ### kilograms

#Kerbin
RKerbin = 600000                    ### meters
mkerbin = 5.2915158*10**22          ### kilograms

## Rocket (hobby rocket scale)
mass = 640.0/1000.0 #kgs

## Gravitational Acceleration
def gravity(z):
    global Rplanet, mplanet

    r = np.sqrt(z**2)

    if r < Rplanet:
            accel = 0.0
    else:
        accel = G*mplanet / (r**3)*r
    
    return accel

### Equations of motion:
## Force = mass * acceleration = mass * zddot
# z = altitude of the surface
# zdot = velocity
# zddot = acceleration

## Second Order Differential Equation 
def Derivatives(state, t):
    global mass
    #state vector
    z = state[0]
    velz = state[1]

    # Compute zdot
    zdot = velz

    ### Total Forces:
    ## Gravity
    gravityF = -gravity(z)*mass

    ## Aerodynamics
    aeroF = 0.0

    ## Thrust
    thrustF = 0.0

    Forces = gravityF + aeroF + thrustF

    # Compute zddot
    zddot = Forces/mass

    # Compute the statedot
    statedot = np.asarray([zdot, zddot])

    return statedot


##### Main Script Below #####

### Test Surface Gravity
print(f'Surface Gravity (m/s^2) = {gravity(Rplanet)}')

### Initial Conditions:
z0 = Rplanet                         ### meters
velz0 = 25*331.0                     ### meters / second
stateinitial = np.asarray([z0, velz0])


## Time window 
# Over 30 secs, give 1000 data points
tout = np.linspace(0, 340, 1000)

### Numerical Integration Call
stateout = sci.odeint(Derivatives, stateinitial, tout)

### Rename Variables
zout = stateout[:,0]
altitude = zout - Rplanet
velzout = stateout[:,1]

### Plot
### Altitude 
mpl.plot(tout, altitude)
mpl.xlabel('Time (seconds [sec])')
mpl.ylabel('Altitude (meters [m])')
mpl.grid()

### Velocity
mpl.figure()
mpl.plot(tout, velzout)
mpl.xlabel('Time (seconds [sec])')
mpl.ylabel('Normal Speed (meters per second [m/s])')
mpl.grid()