#!/usr/bin/python
#
# exmaple_0.py - A very simple example demonstrating the use os fmech for
# insect flight simulation. No controller is used in this example.  
#
# Will Dickson 03/17/2006
# ------------------------------------------------------------------------
import ode
from fmech.sim_util import read_config
from fmech import sim_fly, sim_log

# Read configuration file
config = read_config('config.txt')

# Get simulation parameters
dt = config['time_step']
t_stop = config['stop_time'] 
grav_const = config['grav_const']

# Create world and simulated fly
t = 0.0
world = ode.World()
world.setGravity((0.0, 0.0, -grav_const))
fly = sim_fly(world, config)

# Initialize log
log = sim_log(t, config)

print 'Running simulation ... '
cnt = 0

while t < t_stop:
     
    # Update the fly
    fly.update(t,dt)

    # Update the world
    t += dt
    world.step(dt)

    # Update log
    log.update(fly,t)
       
    # Progress message
    if 100*t/t_stop > cnt+1:
        cnt+=1
        print 'cnt %d, t %1.1f'%(cnt,t)
        
print 'Done'

# Write data to log
log.write()

# Clean up
del fly
del world
ode.CloseODE()
