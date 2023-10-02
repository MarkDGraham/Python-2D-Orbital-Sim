"""
projectile.py

Created on Sun Oct 1 20:59:46 2023

Edited by Mark Graham

Original work done by Carlos Montalvo
"""

#### Import used modules
# %%
import numpy as np                  ### numeric python
import matplotlib.pyplot as mpl     ### matlab style plotting
import scipy.integrate as sci       ### integration toolbox

### Defining Constant Parameters

## Rocket (hobby rocket scale)
mass = 640.0/1000.0 #kgs

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
    gravity = -9.81*mass

    ## Aerodynamics
    aero = 0.0

    ## Thrust
    thrust = 0.0

    Forces = gravity + aero + thrust

    # Compute zddot
    zddot = Forces/mass

    # Compute the statedot
    statedot = np.asarray([zdot, zddot])

    return statedot


##### Main Script Below #####

### Initial Conditions:
z0 = 0.0                            ### meters
velz0 = 164.0                        ### meters / second
stateinitial = np.asarray([z0, velz0])


## Time window 
# Over 30 secs, give 1000 data points
tout = np.linspace(0, 35, 1000)

### Numerical Integration Call
stateout = sci.odeint(Derivatives, stateinitial, tout)

### Rename Variables
zout = stateout[:,0]
velzout = stateout[:,1]

### Plot

# %%
### Altitude 
mpl.plot(tout, zout)
mpl.xlabel('Time (seconds [sec])')
mpl.ylabel('Altitude (meters [m])')
mpl.grid()
# %%
### Velocity
mpl.figure()
mpl.plot(tout, velzout)
mpl.xlabel('Time (seconds [sec])')
mpl.ylabel('Normal Speed (meters per second [m/s])')
mpl.grid()
