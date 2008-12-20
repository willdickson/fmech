#!/usr/bin/env python
"""
-----------------------------------------------------------------------
fmech
Copyright (C) William Dickson, 2008.
  
wbd@caltech.edu
www.willdickson.com

Released under the LGPL Licence, Version 3

This file is part of fmech.

fmech is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
    
fmech is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with fmech.  If not, see <http://www.gnu.org/licenses/>.

------------------------------------------------------------------------   
plot_utils.py

Purpose: exmaple utility demonstrating how to extract and plot data
from the log file generated during a simulation.

Author: William Dickson 
------------------------------------------------------------------------
"""
import sys
import math, cPickle
import matplotlib.pylab as pylab
from fmech import *
from fmech.sim_util import rad2deg

# Get log file
log_file = sys.argv[1]

# Read log file
print 'Reading %s ... '%(log_file,),
sys.stdout.flush()
log_fid = open(log_file, 'r')
log = cPickle.load(log_fid)
log_fid.close()
print 'done'

# Unpack time
t = log.t

# Unpack kine angles
l_ang = [x for x,y in log.wing_angles]
r_ang = [y for x,y in log.wing_angles]

# Left wing stroke, deviation and rotation angles  
l_phi = [rad2deg(x) for x,y,z in l_ang]
l_theta = [rad2deg(y) for x,y,z in l_ang]
l_alpha = [rad2deg(z) for x,y,z in l_ang]

# Right wing stroke, deviation and rotation angles  
r_phi = [rad2deg(x) for x,y,z in r_ang]
r_theta = [rad2deg(y) for x,y,z in r_ang]
r_alpha = [rad2deg(z) for x,y,z in r_ang]

# Unpack control signals
y_sig = [x for x,y,z,w in log.ctl_sig]
p_sig = [y for x,y,z,w in log.ctl_sig]
r_sig = [z for x,y,z,q in log.ctl_sig]
t_sig = [w for x,y,z,w in log.ctl_sig]

# Unpack period
period = [x for x in log.period]

# Unpack haltere forces
haltere_l = [x for x,y in log.haltere_forces]
haltere_r = [y for x,y in log.haltere_forces]

pylab.figure(1)
pylab.subplot(3,1,1)
pylab.title('wing kinematics')
pylab.plot(t,l_phi,'b')
pylab.plot(t,r_phi,'r')
pylab.ylabel('phi')

pylab.subplot(3,1,2)
pylab.plot(t,l_alpha, 'b')
pylab.plot(t,r_alpha, 'r')
pylab.ylabel('alpha')

pylab.subplot(3,1,3)
pylab.plot(t,l_theta, 'b')
pylab.plot(t,r_theta, 'r')
pylab.ylabel('theta')
pylab.xlabel('t (ms)')

pylab.figure(2)
pylab.plot(t,period)
pylab.ylabel('wing beat period (ms)')
pylab.xlabel('t (ms)')

pylab.figure(3)
pylab.subplot(2,1,1)
pylab.title('control signals')
pylab.plot(t,y_sig,'b')
pylab.plot(t,p_sig,'r')
pylab.plot(t,r_sig,'g')
pylab.ylabel('y,p,r signals')

pylab.subplot(2,1,2)
pylab.plot(t,t_sig,'b')
pylab.ylabel('throttle')
pylab.xlabel('t (ms)')

pylab.figure(4)
pylab.plot(t,haltere_l, 'b')
pylab.plot(t,haltere_r, 'r')
pylab.title('Haltere Forces')
pylab.xlabel('t (ms)')
pylab.ylabel('F  (mg*mm)/ms')
pylab.show()
