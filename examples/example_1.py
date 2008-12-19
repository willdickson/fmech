#!/usr/bin/python
#
# example_1.py - This example demostrates how a simple controller can be
# can be incorporated into the flight simulation. 
#
# Will Dickson 03/17/2006
# ------------------------------------------------------------------------
import ode
from fmech.sim_util import read_config
from fmech import sim_fly, sim_log
import ctl_examp

# Read configuration file
config = read_config('config.txt')

# Get simulation parameters
dt  = config['time_step'] 
t_stop = config['stop_time']
grav_const = config['grav_const']

# Create world and simulated fly
t = 0.0
world = ode.World()
world.setGravity((0.0, 0.0, -grav_const))
fly = sim_fly(world, config)

# Set CM poition to origin
fly.setPosition((0,0,0))

# Initialize log
log = sim_log(t, config)
cnt = 0
print 'Running simulation ... '

# Get state variables for controller initialization
ang_vel = fly.getAngularVel()
x, y, z = fly.getPosition()
fvel = fly.getForwardVel()

# Initialize controller
fvel_ctl = ctl_examp.fvel_ctl(fvel, ang_vel, z)
                              
# Control signal low pass filter
ctl_sig_filter = ctl_examp.ctl_sig_filter((0.0,0.0,0.0,0.0))

while t < t_stop:

    # Change controller set point
    if t > 300.0:
        fvel_ctl.set_pt = 0.15
        
    # Get required state variables for controller
    ang_vel = fly.getAngularVel()
    x, y, z = fly.getPosition()
    fvel = fly.getForwardVel()

    # Set control singal 
    fvel_sig = fvel_ctl.get_sig(fvel, ang_vel, z, dt)
    fly.ctl_sig = ctl_sig_filter.update(fvel_sig,dt)
     
    # Update the fly
    fly.update(t,dt)

    # Update the world
    t += dt
    world.step(dt)

    # Update log
    log.update(fly,t)
       
    # Silly progress message
    if 100*t/t_stop > cnt+1:
        cnt+=1
        print 'cnt %d, t %1.2f, fvel %1.2f'%(cnt, t, fvel)
    
        
print 'Done'

# Write data to log
log.write()

# Clean up
del fly
del world
ode.CloseODE()
