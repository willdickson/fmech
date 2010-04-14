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
fmech.py 

Purpose: base classes for ODE based articulated rigid-body insect
flight simulation.

classes:
  sim_fly       - articulated rigid-body insect model
  angle_driver  - wing kinematics driver    
  wing_aero     - wing aerodynamics model
  body_aero     - body aerodynamics model   
  sim_log       - data logger 

Author: William Dickson 

------------------------------------------------------------------------
"""
import ode, math, cPickle
from mirtich import get_mass_props
from sim_util import *

__all__ = ['sim_fly', 'sim_log']

class sim_fly(object):
    """ 
    Simulated fly class encapsulating the ODE articulated rigid-body
    ODE model

    The class represents a fly, or other effectively two winged
    insect, as an articulated rigid body. Included are wing kinematics
    drivers, an aerodynamics model, and simple controllers.
    
    """
    def update(self, t, dt):
        """
        Update the fly's state for the next time step -- sets the wing
        angular velocities and the wing/body forces.

        Inputs:
            t  = current simulation time (ms)
            dt = update time step (ms)
        """
        # Wing and body forces
        self.updateWingKine(t, dt)
        self.updateWingForce(dt)
        self.updateBodyForce()
        # Udpate haltere forces - for logging
        self.getHaltereForce()
        
    def updateBodyForce(self):
        """
        Update forces and torques acting on the fly's body.
        """
        # Get forces and moments
        forces, moments = self.aero_body_model.get_forces()
        norm_force, para_force = forces
        y_moment, p_moment, r_moment = moments

        # Store for data logging
        self.body_norm_force = norm_force
        self.body_para_force = para_force
        self.body_y_moment = y_moment
        self.body_p_moment = p_moment
        self.body_r_moment = r_moment
        
        # Apply forces and moments
        if self.force_flag:
            self.body.addForce(norm_force)
            self.body.addForce(para_force)
            self.body.addRelTorque((y_moment, p_moment, r_moment))

    def updateWingForce(self,dt):
        """
        Update the forces acting on the wings

        Inputs:
            dt = update time step (ms)
        """
        # Get left and right wing forces
        r_qs_forces = self.aero_r_wing.get_force(dt)
        l_qs_forces = self.aero_l_wing.get_force(dt)
        self.r_wing_forces = r_qs_forces
        self.l_wing_forces = l_qs_forces

        # Apply forces
        if self.force_flag:
            for pos, force in r_qs_forces:
                self.r_wing.addForceAtPos(force, pos)
            for pos, force in l_qs_forces:
                self.l_wing.addForceAtPos(force, pos)
            
    def updateWingKine(self, t, dt):
        """
        Update the angular velocity of the wings simulation based on
        the current time and the update step.

        Inputs:
            t  = current simulation time (ms)
            dt = update time step (ms)

        """
        # Get current kinematic wing angles and angle rates
        l_ang, r_ang = [], []
        for i in range(0,3):
            l_ang.append(self.l_amotor.getAngle(i))
            r_ang.append(self.r_amotor.getAngle(i))
        l_ang, r_ang = tuple(l_ang), tuple(r_ang)
        
        # Store current wing angles 
        self.wing_angles = (l_ang, r_ang)

        # Get next set of wing angles
        l_ang_new, r_ang_new = self.angle_driver.get_angles(dt,self.ctl_sig)

        # Determine wing angular velocities using wing controllers 
        l_ang_vel = self.l_wing_ctl.get_ang_vel(l_ang,l_ang_new,dt)
        r_ang_vel = self.r_wing_ctl.get_ang_vel(r_ang,r_ang_new,dt)

        # Store current wing angular velocities
        self.wing_ang_vel = (l_ang_vel, r_ang_vel)

        # Set the angular velocity of the wings
        self.setWingAngularVel(l_ang_vel,r_ang_vel)        
        
    def setWingAngularVel(self, l_wing_vel, r_wing_vel):
        """
        Set the angular velocities of the right and left wings

        Inputs:
            l_wing_vel = left wing angular velocities 
            r_wing_vel = right wing angular velocities

        Note, both l_wing_vel and r_wing_vel are 3-tuples specifying
        the angular velocity for the stroke position angle, deviation
        angle, and rotation angle of the wing. 
            
        """
        # Set right wing angular velocities
        self.r_amotor.setParam(ode.ParamVel,  r_wing_vel[0]) 
        self.r_amotor.setParam(ode.ParamVel2, r_wing_vel[1])
        self.r_amotor.setParam(ode.ParamVel3, r_wing_vel[2])
        # Set left wing angular velocities
        self.l_amotor.setParam(ode.ParamVel,  l_wing_vel[0])
        self.l_amotor.setParam(ode.ParamVel2, l_wing_vel[1])
        self.l_amotor.setParam(ode.ParamVel3, l_wing_vel[2])

    def setAngularVel(self, v):
        """
        Set angular velocity of simulated fly.

        Input:
            v = angular velcoity 3-tuple
        """
        self.body.setAngularVel(v)
        self.l_wing.setAngularVel(v)
        self.r_wing.setAngularVel(v)        

    def getAngularVel(self):
        """
        Gets the angular velocity of the fly. Currently returns the
        angular velocity of the fly's body which isn't really correct,
        but is probably what you want for control. I should really
        change this ... it is a simple calculation ... However, this
        function of for control purposes only. 
        """
        return self.body.getAngularVel()

    def setLinearVel(self, v):
        """
        Set linear velocity of simulated fly

        Input:
            v = linear velocity (vx,vy,vz)
        """
        self.body.setLinearVel(v)
        self.l_wing.setLinearVel(v)
        self.r_wing.setLinearVel(v)

    def getWingAngles(self):
        """
        Access function for wing angles

        Output:
           r_wing_ang, l_wing_ang = left and right wing angles
        
        """
        return self.wing_angles

    def getWingAngularVel(self):
        """
        Access function for wing angular velocities

        Output:
            r_wing_angvel, l_wing_angvel = left and right wing angular velocities
        """
        return self.wing_angvel

    def rotate(self, ax, ang):
        """
        Rotate the fly about the systems center of mass (body + wings)
        using the given axis and angle.

        Inputs:
            ax  = rotation axis
            ang = rotation angle
        
        """
        # Get position system (body+wings) center of mass
        cm_pos = self.getPosition()
        
        # Get body and wing center of mass positions
        body_pos = self.body.getPosition()
        r_wing_pos = self.r_wing.getPosition()
        l_wing_pos = self.l_wing.getPosition()

        # Shift position vectors so that system cm is at the origin
        body_pos = vec_sub(body_pos, cm_pos)
        r_wing_pos = vec_sub(r_wing_pos, cm_pos)
        l_wing_pos = vec_sub(l_wing_pos, cm_pos)

        # Rotate position vectors
        body_pos = rotate_vec(body_pos,ax,ang)
        r_wing_pos = rotate_vec(r_wing_pos,ax,ang)
        l_wing_pos = rotate_vec(l_wing_pos,ax,ang)

        # Shift back to cm position
        body_pos = vec_add(body_pos, cm_pos)
        r_wing_pos = vec_add(r_wing_pos, cm_pos)
        l_wing_pos = vec_add(l_wing_pos, cm_pos)

        # Set new positions
        self.body.setPosition(body_pos)
        self.r_wing.setPosition(r_wing_pos)
        self.l_wing.setPosition(l_wing_pos)
        
        # Get orientation quaternions (tuples) for body and wings
        body_q = self.body.getQuaternion()
        r_wing_q = self.r_wing.getQuaternion()
        l_wing_q = self.l_wing.getQuaternion()

        # Rotate orientation quaternions  w.r.t the given axis and angle
        body_q = rotate_quat(body_q, ax, ang)
        r_wing_q = rotate_quat(r_wing_q, ax, ang)
        l_wing_q = rotate_quat(l_wing_q, ax, ang)
        
        # Set body and wings to new orientation quaterions
        self.body.setQuaternion(body_q)
        self.r_wing.setQuaternion(r_wing_q)
        self.l_wing.setQuaternion(l_wing_q)
        
    def translate(self, v):
        """
        Translates fly's center of mass (body+wings) using the vector v.
        For example if w is the original position vector for the fly's
        center of mass them w+v will be the new position of the cener of
        mass.

        Input:
            v = translation vector (x,y,z) 
        
        """
        # Translate body
        body_pos = self.body.getPosition()
        body_pos_new = vec_add(body_pos,v) 
        self.body.setPosition(body_pos_new)
        # Translate right wing
        r_wing_pos = self.r_wing.getPosition()
        r_wing_pos_new = vec_add(r_wing_pos,v) 
        self.r_wing.setPosition(r_wing_pos_new)
        # Translate left wing
        l_wing_pos = self.l_wing.getPosition()
        l_wing_pos_new = vec_add(l_wing_pos,v) 
        self.l_wing.setPosition(l_wing_pos_new)

    def setPosition(self, pos):
        """
        Set the position of the center of mass of the fly.

        Input:
            pos = new position (x,y,z)
        """
        cm_pos = self.getPosition()
        v = vec_sub(pos,cm_pos) 
        self.translate(v)
        
    def getPosition(self):
        """
        Returns the position of the center of mass of the fly
        """
        # Get Body position and mass
        body_pos = self.body.getPosition()
        body_mass = self.body.getMass().mass
        # Get right wing position and mass
        r_wing_pos = self.r_wing.getPosition()
        r_wing_mass = self.r_wing.getMass().mass
        # Get left wing position and mass
        l_wing_pos = self.l_wing.getPosition()
        l_wing_mass = self.l_wing.getMass().mass
        # Compute center of mass
        mass_total = body_mass + r_wing_mass + l_wing_mass
        c1 = body_mass/mass_total
        c2 = r_wing_mass/mass_total
        c3 = l_wing_mass/mass_total
        cm_pos = (
            c1*body_pos[0] + c2*r_wing_pos[0] + c3*l_wing_pos[0],
            c1*body_pos[1] + c2*r_wing_pos[1] + c3*l_wing_pos[1],
            c1*body_pos[2] + c2*r_wing_pos[2] + c3*l_wing_pos[2]
            )
        return cm_pos

    def getLinearVel(self):
        """
        Returns the linear velocity of the center of mass of the fly
        """
        # Get Body position and mass
        body_vel = self.body.getLinearVel()
        body_mass = self.body.getMass().mass
        # Get right wing position and mass
        r_wing_vel = self.r_wing.getLinearVel()
        r_wing_mass = self.r_wing.getMass().mass
        # Get left wing position and mass
        l_wing_vel = self.l_wing.getLinearVel()
        l_wing_mass = self.l_wing.getMass().mass
        # Compute center of mass
        mass_total = body_mass + r_wing_mass + l_wing_mass
        c1 = body_mass/mass_total
        c2 = r_wing_mass/mass_total
        c3 = l_wing_mass/mass_total
        cm_vel = (
            c1*body_vel[0] + c2*r_wing_vel[0] + c3*l_wing_vel[0],
            c1*body_vel[1] + c2*r_wing_vel[1] + c3*l_wing_vel[1],
            c1*body_vel[2] + c2*r_wing_vel[2] + c3*l_wing_vel[2]
            )
        return cm_vel

    def getForwardVel(self):
        """
        Get forwad velocity of the fly - for use with forward velocity
        controller. Doesn't work for all orientations, but its good
        enough for a simple control demonstrations. 
        """
        vel = self.getLinearVel()
        # Get forward xy plane reference vector
        p0 = self.body.getRelPointPos((0.0,0.0,0.0))
        p1 = self.body.getRelPointPos((1.0,0.0,1.0))
        body_fw = vec_sub(p1,p0)
        fref = vec_sub(body_fw, (0.0,0.0,body_fw[2]))
        fref = vec2unit(fref)
        # Get component of velocity in 
        fvel = vec_dot(vel, fref)
        return fvel


    def getSideVel(self):
        vel = self.getLinearVel()
        # Get side xy plane reference vector
        p0 = self.body.getRelPointPos((0.0,0.0,0.0))
        p1 = self.body.getRelPointPos((0.0,1.0,0.0))
        body_side = vec_sub(p1,p0)
        sref = vec_sub(body_side, (0.0,0.0,body_side[2]))
        sref = vec2unit(sref)
        # Get component of velocity in 
        fvel = vec_dot(vel, sref)
        return fvel

    def getHaltereForce(self):
        """
        Return the Lateral component of the forces acting on the left and right
        halteres. 
        """
        body_ang_vel = self.body.getAngularVel()
        l_ang, r_ang = self.wing_angles
        l_ang_vel, r_ang_vel = self.wing_ang_vel
        # Compute left haltere force
        temp0 = -2.0*self.haltere_mass*l_ang_vel[0]*self.haltere_length
        temp1 = body_ang_vel[1]*math.cos(self.haltere_angle)*math.cos(l_ang[0])
        temp2 = body_ang_vel[0]*math.sin(l_ang[0])
        temp3 = body_ang_vel[2]*math.sin(self.haltere_angle)*math.cos(l_ang[0])
        l_haltere_force = temp0*(temp1 - temp2 + temp3)
        # Compute right haltere force
        temp0 = -2.0*self.haltere_mass*r_ang_vel[0]*self.haltere_length
        temp1 = body_ang_vel[1]*math.cos(self.haltere_angle)*math.cos(r_ang[0])
        temp2 = body_ang_vel[0]*math.sin(r_ang[0])
        temp3 = body_ang_vel[2]*math.sin(self.haltere_angle)*math.cos(r_ang[0])
        r_haltere_force = temp0*(temp1 + temp2 - temp3)
        # Store values for logging
        self.haltere_forces = l_haltere_force, r_haltere_force
        return l_haltere_force, r_haltere_force
        
    def set2Config(self, config):
        """
        Initializes fly according the parameters in the configuration
        dictionary 'config'.

        Notes:
        
        1.) Initial position is with the body center of gravity at the world
            frame origin.

        2.) The fly is facing the down the positive x axis.

        3.) The left and right wings are parallel to the y axis with the left
            wing pointing down the positive y axis and the right wing pointing
            down the negative y axis.

        4.) The vertical axis is the z axis

        5.) Wing's are initialized leading edge up in the (phi,alpha,theta) =
            (0,0,0) position. 
        
        """
        # Set body mass properties
        body_stl_file = config['body_stl_file']
        body_density = config['body_density']
        body_mass, body_cg, body_it = get_mass_props(body_stl_file, body_density)        
        mass = ode.Mass()
        p = (body_mass,) + body_cg + body_it
        mass.setParameters(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9])
        mass.translate((-p[1],-p[2],-p[3]))
        self.body.setMass(mass)
        
        # Rotate body for correct body angle
        self.body_angle = deg2rad( config['body_angle'] ) 
        q = axis_angle2quat((0,1,0),self.body_angle)  
        self.body.setQuaternion(q)

        # Set mass properties of wings
        wing_stl_file = config['wing_stl_file']
        wing_density = config['wing_density']
        wing_mass, wing_cg, wing_it = get_mass_props(wing_stl_file, wing_density)
        p = (wing_mass,) + wing_cg + wing_it
        mass.setParameters(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9])       
        mass.translate((-p[1],-p[2],-p[3]))
        self.r_wing.setMass(mass)
        self.l_wing.setMass(mass)
        
        # Set positions of left and right wings
        joint2body_cg = config['joint2body_cg']
        joint2wing_cg = config['joint2wing_cg']
        l_joint_pos = joint2body_cg
        r_joint_pos = (l_joint_pos[0], -l_joint_pos[1], l_joint_pos[2])
        l_joint_pos = self.body.getRelPointPos(l_joint_pos)
        r_joint_pos = self.body.getRelPointPos(r_joint_pos)
        l_wing_cg = [x-y for x,y in zip(l_joint_pos, joint2wing_cg)]
        r_wing_cg = [l_wing_cg[0], -l_wing_cg[1], l_wing_cg[2]]
        self.r_wing.setPosition(r_wing_cg)
        self.l_wing.setPosition(l_wing_cg)
        
        # Rotate right wing by 180 deg
        q = axis_angle2quat((0,0,1),math.pi)
        self.r_wing.setQuaternion(q)   
        
        # Set linear and angular velocities of body and wings to zero
        v0 = (0.0,0.0,0.0)
        self.setLinearVel(v0)
        self.setAngularVel(v0)

        # Attach and anchor wing joints 
        self.r_joint.attach(self.body, self.r_wing)
        self.l_joint.attach(self.body, self.l_wing)
        self.r_joint.setAnchor(r_joint_pos)
        self.l_joint.setAnchor(l_joint_pos)

        # Initialize joint angular motors
        self.r_amotor.attach(self.body, self.r_wing)
        self.l_amotor.attach(self.body, self.l_wing)
        self.r_amotor.setMode(ode.AMotorEuler)
        self.l_amotor.setMode(ode.AMotorEuler)
              
        # Set angular motor axes - AMotorEuler mode
        self.stroke_angle =  deg2rad(config['stroke_angle'])
        # ---------------------------------------------------
        #ang = (0.5*math.pi - self.stroke_angle) # old method
        # ---------------------------------------------------
        ang = self.stroke_angle - 0.5*math.pi + self.body_angle
        r_axis_0 = rotate_vec((0, 0, 1),(0, 1, 0),ang)
        r_axis_2 = rotate_vec((0,1, 0),(0, 1, 0),ang)
        l_axis_0 = rotate_vec((0, 0,-1),(0, 1, 0),ang)
        l_axis_2 = rotate_vec((0,1, 0),(0, 1, 0),ang)
        self.r_amotor.setAxis(0,1, r_axis_0)
        self.r_amotor.setAxis(2,2, r_axis_2)
        self.l_amotor.setAxis(0,1, l_axis_0)
        self.l_amotor.setAxis(2,2, l_axis_2)

        # Set force limits 
        self.r_amotor.setParam(ode.ParamFMax, ode.Infinity)
        self.l_amotor.setParam(ode.ParamFMax, ode.Infinity)
        self.r_amotor.setParam(ode.ParamFMax2, ode.Infinity)
        self.l_amotor.setParam(ode.ParamFMax2, ode.Infinity)
        self.r_amotor.setParam(ode.ParamFMax3, ode.Infinity)
        self.l_amotor.setParam(ode.ParamFMax3, ode.Infinity)

        # Default wing beat period and period limits
        start_tc = config['start_tc']
        period_limits = config['period_limits']
        ypr_sig_limits = config['ypr_sig_limits']

        # Setup wing kinematics driver - adjust rotation angles to account for stroke angle 
        angles = read_wing_kine(config['wing_kine_file'])
        angles = [ (x,y,z-ang) for x,y,z in angles]
        self.angle_driver = angle_driver(angles,period_limits,ypr_sig_limits,t_ninety=start_tc)
        self.kine_angles =((0.0,0.0,0.0),(0.0,0.0,0.0)) 

        # Default control signal
        self.ctl_sig = (0.0,0.0,0.0,0.0)
        
        # Setup wing controllers
        wing_ctl_gain = config['wing_ctl_gain']
        self.l_wing_ctl = wing_controller(gain=wing_ctl_gain)
        self.r_wing_ctl = wing_controller(gain=wing_ctl_gain)
        
        # Setup wing aerodynamics
        self.force_flag = config['force_flag']
        self.aero_r_wing = wing_aero(self.r_wing, config)
        self.aero_l_wing = wing_aero(self.l_wing, config)
        self.r_wing_forces = []
        self.l_wing_forces = []

        # Setup body_aerodynamics
        self.aero_body_model = body_aero(self.body, config) 
        self.body_norm_force = 0.0
        self.body_para_force = 0.0

        # Haltere information
        self.haltere_length = config['haltere_length']
        self.haltere_mass = config['haltere_mass']
        self.haltere_angle = deg2rad(config['haltere_angle'])
        self.haltere_forces = 0.0, 0.0
        
    def __init__(self, world, config):
        """
        Create articulated rigid body fly object for dynamic simulation.

        The fly's body is initialized so that it's center of mass is
        (0,0,0) in world coordinates.

        Inputs:
            world = ODE simulation world object
            config = simulation configuration dictionary
        
        """
        # Create fly body
        self.body = ode.Body(world) 

        # Create wings
        self.r_wing = ode.Body(world)
        self.l_wing = ode.Body(world)

        # Create wing joints
        self.r_joint = ode.BallJoint(world)
        self.l_joint = ode.BallJoint(world)

        # Create angular motors
        self.r_amotor = ode.AMotor(world)
        self.l_amotor = ode.AMotor(world)        
        # Set body according to configuration
        self.set2Config(config)

        
class angle_driver:
    """
    A simple class for setting the wing positions based on a set of
    baseline wing kinematics
    """

    def __init__(self, angles, period_limits=(4.0,6.0), ypr_sig_limits=(20.0,10.0,20.0), t_ninety = 3.5):
        """
        Initializes wing driver.
         angles   = the data set which defines the baseline wing kinematics.
         period_limits =  mimimum and maximum allowed wing beat periods (clamp)
         ypr_sig_limits = maximum absolute values of the y,p,r signals (clamp)
         t_ninety = startup response - time to reach 90% of full stroke amplitude.
                    This stops the startup process from being too violent. 
        """
        self.angles = angles # baseline kinematics data
        self.t_ninety = t_ninety # startup time constant
        self.x = 0.0 # Current position on kinematics
        self.t = 0.0 # time (ms)
        self.period_limits = period_limits
        self.ypr_sig_limits = map(deg2rad, ypr_sig_limits)
        self.period = 0.5*(period_limits[0] + period_limits[1])
      
    def get_angles(self, dt, ctl_sig):
        """
        Get angular position of the wing given the time step dt and
        the ctl_signal.

        Inputs:
            dt      = time step (ms)
            ctl_sig = wing deformation control signal 4-tuple
                      ctl_sig[0] = yaw mode signal
                      ctl_sig[1] = pitch mode signal
                      ctl_sig[2] = roll mode signal
                      ctl_sig[3] = throttle control signal

        Outputs:
            l_ang_new, r_ang_new = new wing angles

        Note: currently the throttle control signal only sets the wing
        beat period. This could be altered so that the throttle
        adjusts both the wing beat period and the stroke amplitude
        which would be more realistic. 
        """
        # Extract the four basic control signals:
        y_sig, p_sig, r_sig, t_sig = ctl_sig

        # Set throttle - currentl just sets wing beat period
        period_max = self.period_limits[1]
        period_min = self.period_limits[0]
        temp0 = 0.5*(period_max + period_min)
        temp1 = 0.5*(period_min - period_max)
        period = temp0 + t_sig*temp1
        
        # Clamp wing beat period within limits
        period = min([period, period_max])
        period = max([period, period_min])
        self.period = period # For logging

        # Update position with in kinematics and time
        self.x += len(self.angles)*dt/period
        self.t += dt
        
        # Get neighboring indices
        s1 = math.floor(self.x) 
        s2 = math.ceil(self.x)
        if s1 == s2:
            s2 = s1 + 1
        n1 = int(s1 % len(self.angles))
        n2 = int(s2 % len(self.angles))
          
        # Get angles at current t - by interpolation
        kine_ang = []
        for i in range(0,3):
            a = (self.angles[n2][i] - self.angles[n1][i])/(s2-s1)
            b = self.angles[n1][i] - a*s1
            ang = a*self.x  + b
            kine_ang.append(ang)

        # Include startup function - prevents kicks during wing start up
        ang_base = []
        tau = self.t_ninety/math.sqrt(-math.log(1.0 - 0.9))
        for i in range(0,3):
            ang = (1 - math.exp(-(self.t/tau)**2))*kine_ang[i] 
            ang_base.append(ang)

        # Clamp the p,r,y signals
        y_sig = min([y_sig,  self.ypr_sig_limits[2]])
        y_sig = max([y_sig, -self.ypr_sig_limits[2]])
        p_sig = min([p_sig,  self.ypr_sig_limits[0]])
        p_sig = max([p_sig, -self.ypr_sig_limits[0]])
        r_sig = min([r_sig,  self.ypr_sig_limits[1]])
        r_sig = max([r_sig, -self.ypr_sig_limits[1]])
            
        # 1st apply pitch and yaw deformation
        l_ang_new = (1 + y_sig)*ang_base[0] + p_sig, ang_base[1], ang_base[2]
        r_ang_new = (1 - y_sig)*ang_base[0] + p_sig, ang_base[1], ang_base[2]

        # 2nd rotate strokes planes for roll deformation
        l_ang_new = rotate_vec(l_ang_new, (0.0,0.0,1.0), -r_sig)
        r_ang_new = rotate_vec(r_ang_new, (0.0,0.0,1.0),  r_sig)

        # Adjust rotation angle
        l_ang_new = l_ang_new[0], l_ang_new[1], l_ang_new[2] + r_sig
        r_ang_new = r_ang_new[0], r_ang_new[1], r_ang_new[2] - r_sig
        
        return l_ang_new, r_ang_new
            
 
class wing_controller:
    """
    A simple proportional controller which to keeps the wings tracking
    the desired kinematics
    """

    def __init__(self, gain = 0.1, max_vel = (100.0, 100.0, 300.0)):
        """
        Initialize wing controller.

        keyword args:
            gain = proportional gain for wing controller
            max_vel = maximum allowable wing velocities (deg/s)
        """
        self.gain = gain
        self.max_vel = [ deg2rad(x) for x in max_vel]
        self.min_vel = [-deg2rad(x) for x in max_vel]

    def get_ang_vel(self, ang_meas, ang_setp, dt):
        """
        Get wing angular velocity of the wing given the current
        (measured) wing position and the next desired (set point) wing
        position. 

        Inputs:
            ang_meas = measured or actual wing position 
            ang_setp = set point or desired wing position

        Outputs:
            ang_vel  = tuple of wing angular velocities - stoke position,
                       deviation, rotation. 
        """
        # Compute wing velocity
        vel_list = []
        for i in range(0,3):
            err = (ang_setp[i] - ang_meas[i])/dt
            vel = self.gain*err 
            # Apply max velocity clamp
            vel = min([vel, self.max_vel[i]])
            vel = max([vel, self.min_vel[i]])
            vel_list.append(vel)
        return tuple(vel_list)

        
class body_aero:
    """
    This class repsents the body aerodynamics model.  
    """

    def __init__(self, body, config):
        self.body = body
        self.force_params = config['body_force_params']
        self.air_density = config['air_density']
        self.body_ref_area = config['body_ref_area']
        self.body_ref_length = config['body_ref_length']
        n = len(config['body_moment_params'])
        assert n%2==0, 'body moment coefficient must be even'
        moment_param_sin = config['body_moment_params'][:n/2]
        moment_param_cos = config['body_moment_params'][n/2:]
        self.moment_param = zip(moment_param_sin,moment_param_cos)

    def get_forces(self):
        """
        Compute the aerodynamic forces on the body
        """
        # Get body and free stream velocity
        body_vel = self.body.getLinearVel()
        free_stream = vec_scalar_mul(-1.0, body_vel)

        # Get reference vectors (body vector is longitudal body axis)
        body_vector = self.body.getRelPointPos((0.0,0.0,1.0))
        body_pos = self.body.getPosition()
        body_vector = vec_sub(body_vector, body_pos) # Unit vector
        
        # Get parallel and normal components of free stream velocity
        norm_vel = vec_comp(free_stream, body_vector)
        para_vel = vec_remove_comp(free_stream, body_vector)

        # Compute normal and parallel drag forces
        vel_mag = vec_mag(free_stream)
        norm_vv = vec_scalar_mul(vel_mag, norm_vel)
        para_vv = vec_scalar_mul(vel_mag, para_vel)
        force_const = 0.5*self.air_density*self.body_ref_area
        scal_norm = force_const*self.force_params[0]
        scal_para = force_const*self.force_params[1]
        force_norm = vec_scalar_mul(scal_norm, norm_vv)
        force_para = vec_scalar_mul(scal_para, para_vv)

        # Compute dorsoventral body vector
        dorso_vector  = self.body.getRelPointPos((1.0,0.0,0.0))
        dorso_vector = vec_sub(dorso_vector, body_pos) # Unit vector

        # Get angle between longitudinal body axis and velocity (delta)
        delta = vec_angle(body_vector,body_vel)
        
        # Get side slip angle (beta)
        body_vel_proj = vec_dot(body_vector,body_vel)
        vel_longi = vec_scalar_mul(body_vel_proj, body_vector)
        vel_dorso_trans = vec_sub(body_vel,vel_longi) 
        beta = vec_angle(vel_dorso_trans,dorso_vector)
        
        # Compute moment fit functions 
        f_pos = 0.0
        f_neg = 0.0
        for i,param in enumerate(self.moment_param):
            n = i+1
            a,b = param
            f_pos += a*math.sin(n*delta) + b*math.cos(n*delta)
            f_neg += a*math.sin(-n*delta) + b*math.cos(-n*delta)

        # Get pitch and yaw moments coefficients (roll moment is small)
        temp0 = 0.5*(1+math.cos(2.0*beta))
        temp1 = 0.5*(1-math.cos(2.0*beta)) 
        coeff_p = temp0*(f_pos*temp0 + f_neg*temp1)
        coeff_y = temp1*(f_pos*temp0 + f_neg*temp1)
    
        # Get pitch, roll and yaw moments
        coeff2moment = 0.5*self.air_density*self.body_ref_area*self.body_ref_length*(vel_mag**2) 
        moment_p = coeff_p*coeff2moment
        moment_y = coeff_y*coeff2moment
        moment_r = 0.0

        # Pack results into tuples
        forces = force_norm, force_para
        moments = moment_y, moment_p, moment_r        
        return forces, moments 

class wing_aero:
    """
    A simple class for representing the wing aerodynamics model. 
    """

    def __init__(self, wing, config):
        self.wing = wing
        self.qs_lift_param = config['wing_qs_lift_param']
        self.qs_drag_param = config['wing_qs_drag_param']
        self.rotation_coeff = config['wing_rotation_coeff']
        self.num_element = config['wing_num_elem']
        self.mean_chord, self.length, self.chord_data = get_wing_info(config['wing_chord_file'])
        self.joint2wing_cg = config['joint2wing_cg']
        self.chord_ref_pt = (0,0,1)
        self.span_ref_pt = (0,1,0)
        self.norm_ref_pt = (-1,0,0)
        self.elem_data = self.get_elem() # (elem_pos, elem_len)
        self.air_density = config['air_density']
        self.kine_derivs = _kine_derivs(self.elem_data)
        
    def get_qs_force(self, elem, elem_len):
        """
        Computes the quasi-steady (steady-state) force for a single wing element
        """
        # Get element velocity 
        elem_vel_vec = self.wing.getRelPointVel(elem)
        # Get refernce vectors
        chord_vec, span_vec, norm_vec = self.get_ref_vec()
        # Get wing element width
        elem_width = self.length/float(self.num_element)
        # Calculate chordwise component of velocity + magnitude and unit vector
        chord_vel_vec = vec_remove_comp(elem_vel_vec,span_vec)
        chord_vel_mag = vec_mag(chord_vel_vec)
        chord_vel_unit = vec2unit(chord_vel_vec)
        # Compute angle of attack and force coefficients
        aoa = vec_angle(chord_vel_vec, chord_vec)
        # Get lift and drag coefficients
        lift_coeff = get_qslift_coeff(aoa, self.qs_lift_param)
        drag_coeff = get_qsdrag_coeff(aoa, self.qs_drag_param)
        # Compute steady state lift and drag magnitudes for elements
        lift_mag = 0.5*self.air_density*lift_coeff*elem_len*elem_width*(chord_vel_mag)**2
        drag_mag = 0.5*self.air_density*drag_coeff*elem_len*elem_width*(chord_vel_mag)**2
        # Compute steady state lift, drag, and total force vectors
        drag_vec = vec_scalar_mul(-drag_mag,chord_vel_unit) 
        lift_dir = vec_remove_comp(chord_vec, chord_vel_unit)
        lift_unit = vec2unit(lift_dir)
        lift_vec = vec_scalar_mul(lift_mag,lift_unit)
        force_vec = vec_add(lift_vec, drag_vec)
        return force_vec

    def get_rot_force(self, elem, elem_len, dt):
        """
        Computes the rotational component of the forces on the wing.

        This sucks, something is wrong with the rotational theory. It
        doesn't really make sense for for a general case such as
        this. So I've fudged it so that it works reasonably well, but
        it is not great.
        """
        # Get element velocity 
        elem_vel_vec = self.wing.getRelPointVel(elem)
        # Get refernce vectors
        chord_vec, span_vec, norm_vec = self.get_ref_vec()
        # Get wing element width
        elem_width = self.length/float(self.num_element)
        # Calculate chordwise component of velocity + magnitude 
        chord_vel_vec = vec_remove_comp(elem_vel_vec,span_vec)
        # Get magnitude of normal velocity  
        norm_vel_mag = math.fabs(vec_dot(chord_vel_vec, norm_vec))
        # Get direction vector for rotational force
        if vec_dot(norm_vec, chord_vel_vec) >= 0:
            rot_force_unit = norm_vec
        else:
            rot_force_unit = vec_scalar_mul(-1.0,norm_vec)
        # Compute angle of attack and force coefficients
        aoa = vec_angle(chord_vel_vec, chord_vec)
        # Compute vel_mag*d(aoa)/dt
        vel_daoa = self.kine_derivs.get_vel_daoa(elem, aoa, norm_vel_mag, dt)
        # Compute rotatinal force magnitude and vector
        rot_force_mag = self.rotation_coeff*self.air_density*vel_daoa*elem_width*elem_len**2
        rot_force_vec = vec_scalar_mul(rot_force_mag, rot_force_unit)
        return rot_force_vec
        
#     def get_am_force(self):
#         """
#         Computes the added mass component of the forces on the
#         wing. I'm leaving this out for now as this component is
#         usually more trouble than it is worth. I'll try to re-work the
#         code and make this component more robust and then re-insert
#         the code.
#         """
#         pass
    
    def get_force(self, dt):
        """
        Iterates over wing/blade elements and computes the aerodynamic forces
        on the wing.
        """
        # Array for force positions and vectors
        elem_force = []
        # Compute force vector for each element
        for elem, elem_len in self.elem_data:
            # Get element position vector and velocity in world frame
            elem_pos = self.wing.getRelPointPos(elem)
            # Get quasi-steady (steady-state) force component
            qs_force_vec = self.get_qs_force(elem, elem_len)
            # Get rotational force component
            rot_force_vec = self.get_rot_force(elem, elem_len, dt)
            # Get total force
            force_vec = vec_add(qs_force_vec, rot_force_vec)
            # Add force to list 
            elem_force.append((elem_pos, force_vec))
        return elem_force
        
    def get_elem(self):
        """
        Get location of wing elements in wing frame and local chord length
        """
        elem = []
        elem_len = []
        for n in range(0,self.num_element):
            # Get element location in wing frame
            s = float(n+0.5)/float(self.num_element)
            x = self.joint2wing_cg[0]
            y = self.joint2wing_cg[1] + s*self.length
            z = self.joint2wing_cg[2]
            elem.append((x,y,z))
            # Get local chord length
            chord = self.mean_chord*interp_list(self.chord_data, s)
            elem_len.append(chord)
        return zip(elem, elem_len)

    def get_ref_vec(self):
        """
        Get refernce vectors for force calculations. 
        """
        # Get reference points in world coordinates
        wing_cg = self.wing.getRelPointPos((0,0,0))
        chord_ref = self.wing.getRelPointPos(self.chord_ref_pt)
        span_ref = self.wing.getRelPointPos(self.span_ref_pt)
        norm_ref = self.wing.getRelPointPos(self.norm_ref_pt)
        
        # Get reference vectors in world coordinates
        chord_vec = vec_sub(chord_ref, wing_cg)
        span_vec = vec_sub(span_ref, wing_cg)
        norm_vec = vec_sub(norm_ref, wing_cg)
        return chord_vec, span_vec, norm_vec

class _kine_derivs:
    """
    Kinemaic derivative data used in the force calculations. This is
    kind of kludgey.
    """
    def __init__(self, elem_data, filter_tc = 0.1):
        # Initialize data dictionaries
        self.aoa_last = {}
        self.vel_last = {}
        self.vel_daoa = {}
        self.vel_daoa_filter = {}
        for elem, elem_lem in elem_data:
            self.aoa_last[elem] = 0.0
            self.vel_last[elem] = 0.0
            self.vel_daoa[elem] = 0.0
            self.vel_daoa_filter[elem] = Filt1(0.0,filter_tc)  
    def get_vel_daoa(self, elem, aoa, vel, dt):
        # Compute vel times d(aoa)/dt - using low pass filter for smoothing
        temp0 = (aoa*vel - self.aoa_last[elem]*self.vel_last[elem])/dt
        temp1 = aoa*(vel - self.vel_last[elem])/dt
        temp3 = temp0 - temp1
        self.vel_daoa[elem] = self.vel_daoa_filter[elem].update(temp3,dt)
        # Reset last states
        self.aoa_last[elem] = aoa
        self.vel_last[elem] = vel
        return self.vel_daoa[elem]
        
class sim_log:
    """
    Simple class for logging data during the simulation. This is
    definitely not the way to go for long simulations as everything is
    stored in memory - not logged to disk. However, it is useful
    shorter simulations and debugging.
    """
    def __init__(self, t, log_config):

        # Set log file, update interval, and update time
        self.log_file = log_config['log_file']
        self.log_dt = log_config['log_dt']
        self.t_last_log = t - self.log_dt
        self.t = []
        
        # Initialize body/wings position and quaternion lists
        self.body_p = []
        self.body_q = []
        self.l_wing_p = []
        self.l_wing_q = []
        self.r_wing_p = []
        self.r_wing_q = []
        
        # Initialize system cm and velocity list
        self.fly_cm = []
        self.fly_vel = []
        self.fly_fvel = []
        
        # Initialize wing angle, control signal and period lists
        self.wing_angles = []
        self.ctl_sig = []
        self.period = []
        
        # Initialize wing for lists
        self.r_wing_forces = []
        self.l_wing_forces = []

        # Initialize body force and moment lists
        self.body_norm_forces =[]
        self.body_para_forces =[]
        self.body_forces = []
        self.body_p_moment = []
        self.body_r_moment = []
        self.body_y_moment = []

        # Haltere forces
        self.haltere_forces = []
        
    def update(self, fly, t):
        """
        Update the log for the given fly at the given time point. 
        """
        if t - self.t_last_log >= self.log_dt:
            
            # Time point
            self.t.append(t)
            # Body position and orientation
            self.body_p.append(fly.body.getPosition())
            self.body_q.append(fly.body.getQuaternion())
        
            # Left wing position and orientation
            self.l_wing_p.append(fly.l_wing.getPosition())
            self.l_wing_q.append(fly.l_wing.getQuaternion())
        
            # Right wing position and orientation
            self.r_wing_p.append(fly.r_wing.getPosition())
            self.r_wing_q.append(fly.r_wing.getQuaternion())
            self.t_last_log = t

            # Fly center of mass
            self.fly_cm.append(fly.getPosition())
            self.fly_vel.append(fly.getLinearVel())
            self.fly_fvel.append(fly.getForwardVel())

            # Wing angles, control values, wing beat period
            self.wing_angles.append(fly.wing_angles)
            self.ctl_sig.append(fly.ctl_sig)            
            self.period.append(fly.angle_driver.period)

            # Wing forces
            self.r_wing_forces.append(fly.r_wing_forces)
            self.l_wing_forces.append(fly.l_wing_forces)

            # Body forces and moments
            self.body_norm_forces.append(fly.body_norm_force)
            self.body_para_forces.append(fly.body_para_force)
            total_force = vec_add(fly.body_norm_force,fly.body_para_force)
            self.body_forces.append(total_force)
            self.body_p_moment.append(fly.body_p_moment)
            self.body_r_moment.append(fly.body_r_moment)
            self.body_y_moment.append(fly.body_y_moment)

            # Haltere forces
            self.haltere_forces.append(fly.haltere_forces)
            
    def write(self):
        log_fid = open(self.log_file,'w')
        cPickle.dump(self, log_fid)
        log_fid.close()
