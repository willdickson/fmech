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
sim_util.py 

Purpose: utility functions for fmech.

Author: William Dickson 
------------------------------------------------------------------------
"""
import math, stl_tools, string
import os.path
try:
    import cgkit.cgtypes as cgtypes # cgkit 2
except ImportError, err:
    import cgtypes # cgkit 1

class Filt1:
    """
    Encapsulates simple 1st order filter.
    """
    def __init__(self, x0, tc):
        """
        Initialize 1st order filter.
            x0 = initial value,
            tc  = time constant
        """
        self.val = x0
        self.tc = tc

    def update(self, x, dt):
        """
        Update filter. x = new value, dt = time step
        """
        c0 = self.tc/(self.tc+dt)
        c1 = dt/(self.tc+dt)
        if type(self.val)==list or type(self.val)==tuple:
            self.val = [c0*a + c1*b for a,b in zip(self.val, x)]
        else:
            self.val = c0*self.val + c1*x
        return self.val

class Integrator:
    """
    Simple integrator with clamp for use in control systems
    """
    def __init__(self, x0, clamp = None):
        self.val = x0
        self.clamp = clamp

    def update(self, x):
        self.val += x
        if self.clamp:
            self.val = min([self.val,  self.clamp])
            self.val = max([self.val, -self.clamp])
        return self.val

def deg2rad(angle_deg):
    """
    Converts and angle from degrees to radians
    """
    return angle_deg*math.pi/180.0


def rad2deg(angle_rad):
    """
    Converts an angle from radians to degrees
    """
    return angle_rad*180.0/math.pi


def quat2mat(q):
    """
    Convert quaterion to rotation matrix
    """
    qq1 = 2*q[1]*q[1]
    qq2 = 2*q[2]*q[2]
    qq3 = 2*q[3]*q[3]

    R00 = 1 - qq2 - qq3;
    R01 = 2*(q[1]*q[2] - q[0]*q[3]);
    R02 = 2*(q[1]*q[3] + q[0]*q[2]);
    R10 = 2*(q[1]*q[2] + q[0]*q[3]);
    R11 = 1 - qq1 - qq3;
    R12 = 2*(q[2]*q[3] - q[0]*q[1]);
    R20 = 2*(q[1]*q[3] - q[0]*q[2]);
    R21 = 2*(q[2]*q[3] + q[0]*q[1]);
    R22 = 1 - qq1 - qq2;
    return R00, R01, R02, R10, R11, R12, R20, R21, R22

def axis_angle2quat(ax, ang):
    """
    Get quaterion from axis and angle (rad)

    """
    ax_len = math.sqrt(float(ax[0]**2 + ax[1]**2 + ax[2]**2))
    if ax_len > 0.0:
        q0 = math.cos(0.5*ang)
        q1 = ax[0]*math.sin(0.5*ang)/ax_len
        q2 = ax[1]*math.sin(0.5*ang)/ax_len
        q3 = ax[2]*math.sin(0.5*ang)/ax_len
    else:
        q0 = 1.0
        q1 = 0.0
        q2 = 0.0
        q3 = 0.0
    return q0, q1, q2, q3


def euler2quat():
    """
    Convet zyx euler angles to quaternions
    """
    pass

def quat2euler(q):
    """
    Convert quaternions to zyx euler angles
    """
    pass

def quat2tuple(q):
    """
    Convert cgtypes quaternion to tuple
    """
    return q.w, q.x, q.y, q.z


def vec_add(v,w):
    """
    Add tuples as  vectors
    """
    return tuple([x+y for x,y in zip(v,w)])

def vec_sub(v,w):
    """
    Subract tuples as  vectors
    """
    return tuple([x-y for x,y in zip(v,w)])

def vec_mag(v):
    """
    Compute vector magnitude - simple minded and slow
    """
    val = 0
    for x in v:
        val+=x**2
    return math.sqrt(val)

def vec_dot(v,w):
    """
    Compute vector dot product
    """
    val = 0
    for x,y in zip(v,w):
        val+= x*y
    return val

def vec_cross(v,w):
    """
    Computes the cross product of two vectors
    """
    return v[1]*w[2]-v[2]*w[1], -v[0]*w[2]+v[2]*w[0], v[0]*w[1]-v[1]*w[0]

def vec_comp(v,w):
    """
    Computes the component of the vector v in vector w direction.
    """
    v_dot_w = vec_dot(v,w)
    w_unit = vec2unit(w)
    return vec_scalar_mul(v_dot_w, w_unit)

def vec_remove_comp(v,w):
    """
    Remove the component of the vector v in the vector w direction
    """
    w_comp = vec_comp(v,w)
    return vec_sub(v,w_comp)

def vec2unit(v):
    """
    Computes the unit vector in v direction
    """
    mag = vec_mag(v)
    if mag != 0:
        return tuple([x/mag for x in v])
    else:
        return tuple([0 for x in v])

def vec_scalar_mul(scal, vec):
    """
    multiply elements of vector v by scalar scal
    """
    return tuple([scal*x for x in vec])

def vec_angle(v,w):
    """
    get the angle between two vectors
    """
    v_mag = vec_mag(v)
    w_mag = vec_mag(w)
    if v_mag==0 or w_mag==0:
        return 0.0
    else:
        denom = vec_dot(v,w)
        num = vec_mag(vec_cross(v,w))
        return math.atan2( num, denom )

def qrotate_vec(v,q):
    """
    Rotate vector using a quaternion.

    v = input vector (tuple)
    q = quaternion (tuple)

    """
    # Get cgtypes rotation quaternion and its inverse
    rot_q = cgtypes.quat(q)
    rot_q_inv = rot_q.inverse()
    # Rotate vector
    vq = cgtypes.quat(0.0, v[0], v[1], v[2])
    vq_new = rot_q*vq*rot_q_inv
    return vq_new.x, vq_new.y, vq_new.z

def rotate_vec(v, ax, ang):
    """
    Rotate a vector by a given angle about a given axis.

    v = input vector
    ax = rotation axis
    ang = rotation angle (radians)

    """
    # Get rotation quaterion and its inverse
    rot_q = cgtypes.quat(ang, ax)
    rot_q_inv = rot_q.inverse()
    # Rotate vector
    vq = cgtypes.quat(0.0,v[0],v[1],v[2])
    vq_new = rot_q*vq*rot_q_inv
    return vq_new.x, vq_new.y, vq_new.z

def rotate_quat(q, ax, ang):
    """
    Rotate quaternion orientation by given axis and angle

    q = input quaternion
    ax rotation axis
    ang = rotation angle

    """

    # Convert quaternion to cgtypes quaternion
    q1 = cgtypes.quat(q)
    q_rot = cgtypes.quat(ang,ax)
    q2 = q_rot*q1
    return quat2tuple(q2)

def get_wing_length(stl_file):
    """
    Get the wing length from the polyhedral model in the stl file
    """
    wing_facets = stl_tools.read_stl('wing.stl')
    wing_length = stl_tools.get_extent(wing_facets,1)
    return wing_length

def get_wing_info(chord_file):
    """
    Get the mean chord and wing chord data from the wing chord file
    """
    fid = open(chord_file, 'r')
    # Read mean chord
    line = string.split(fid.readline())
    mean_chord = float(line[0])
    line = string.split(fid.readline())
    length = float(line[0])
    # Read chord data
    chord_data = []
    for line in fid.readlines():
        line = string.split(line)
        chord_data.append(float(line[1]))
    fid.close()
    return mean_chord, length, chord_data

def read_wing_kine(filename):
    """
    Read wing kinematics data file. Note, angles are assumed to be in degrees
    and are convered to radians.
    """
    angles = []
    fid = open(filename, 'r')
    for line in fid.readlines():
        line = string.split(line)
        ang0 = deg2rad(float(line[0]))
        ang1 = deg2rad(float(line[1]))
        ang2 = deg2rad(float(line[2]))
        angles.append((ang0, ang2, -ang1))
    fid.close()
    return angles

def interp_list(the_list, x):
    """
    Interpolate values in the list. Note, elements of the list
    [l[0], l[1], ..., l[n-1]] are assumed to form (x,y) pairs
    (0,l[0]) ... (n-1,l[n-1]). The list is then interpolated using
    the given x which must be in the range [0,1].
    """
    if x < 0 or x > 1:
        raise RuntimeError, 'interpolation location must be within [0,1]'
    s = float(x*(len(the_list)-1))
    s0, s1 = math.floor(s), math.ceil(s)
    n0, n1 = int(s0), int(s1)
    if s==s0:
        y = the_list[n0]
    else:
        a = (the_list[n1] - the_list[n0])/(s1-s0)
        b = the_list[n1] - a*s1
        y = a*s + b
    return y

def ang_vel_test(fly, t, dt, amp, period, t_set):
    """
    Sample angular velocity function
    """
    # Get startup function
    tau = t_set/math.sqrt(-math.log(1 - 0.9))
    startup = 1.0 - math.exp(-(t/tau)**2)

    # Get the current wing angles
    r_ang0 = fly.r_amotor.getAngle(0)
    r_ang1 = fly.r_amotor.getAngle(1)
    r_ang2 = fly.r_amotor.getAngle(2)
    l_ang0 = fly.l_amotor.getAngle(0)
    l_ang1 = fly.l_amotor.getAngle(1)
    l_ang1 = fly.l_amotor.getAngle(2)

    # Get next angle - from kinematics
    r_ang0_next = startup*amp*math.cos(2*math.pi*(t+dt)/period)
    r_ang1_next = startup*amp*math.sin(2*math.pi*(t+dt)/period)
    l_ang0_next = startup*amp*math.cos(2*math.pi*(t+dt)/period)
    l_ang1_next = startup*amp*math.sin(2*math.pi*(t+dt)/period)

    # Compute angular velocities
    r_ang_vel = 0.9*(r_ang0_next - r_ang0)/dt, 0.9*(r_ang1_next - r1)/dt, 0
    l_ang_vel = 0.9*(l_ang0_next - l_ang0)/dt, 0.9*(l_ang1_next - l_ang1)/dt, 0
    return l_ang_vel, r_ang_vel

def get_qslift_coeff(aoa, qs_params):
    """
    Compute the quasi-steady lift coefficient
    """
    return qs_params[0]*math.sin(aoa)*math.cos(aoa) + qs_params[1]

def get_qsdrag_coeff(aoa, qs_params):
    """
    Compute the quasi-steady lift coefficient
    """
    return qs_params[0]*math.sin(aoa)*math.sin(aoa) + qs_params[1]

def read_config(filename=None):
    """
    Reads simulation configuration file - replace missing or items set
    equal to default entrees with the default values.
    """
    # Defualt files which require package data path
    def_files_special = [
        'body_stl_file',
        'wing_stl_file',
        'wing_kine_file',
        'wing_chord_file'
        ]

    # Read configuration file
    if filename == None:
        sym_config = {}
    else:
        sym_config = eval_read(filename)

    # Read default configuration file
    data_dir, junk = os.path.split(__file__)
    data_dir = os.path.join(data_dir, 'data')
    def_config_file = os.path.join(data_dir, 'default_config.txt')
    def_config = eval_read(def_config_file)

    # Replace missing or values set to default with the default values
    for key, val in def_config.iteritems():
        if not sym_config.has_key(key) or sym_config[key] == 'default':
            if key in def_files_special:
                # Add package data path to special file names
                val = os.path.join(data_dir, val)
            sym_config[key] = val
    return sym_config

def eval_read(filename):
    """
    Create dictionary by reading a file line by line and setting the
    key of the dictionary to value on the left of equals sign and the
    value equal to eval of the stuff on the right of the equals sign.
    """
    the_dict = {'default':'default'}
    execfile(filename, {}, the_dict)
    del the_dict['default']
    return the_dict

# Testing ----------------------------------------------------------------------------
if __name__ == '__main__':
    filename = 'sim_config.txt'
    sym_config = read_config(filename)
    print sym_config


