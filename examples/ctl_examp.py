#
# ctl_examp.py - simple examples of controllers which can be used
# with the GUFM simulation. All exmaples use direct state feedback.
#
# classes:
#
#    ctl_sig_filter - low pass filter for control signals
#    pry_ctl        - proportional pitch/roll/yaw rate controller
#    alt_ctl        - proportional derivative altitude controller
#    vz_ctl         - proportional vertical velocity controller
#    fvel_ctl       - foward velocity+altitide  controller (uses pry_ctl
#                     and alt_ctl as nested controllers)
#
# Author Will Dickson 03/14/2006
# --------------------------------------------------------------------------
import math
from fmech.sim_util import Filt1, Integrator

class ctl_sig_filter:
    """
    Simple set of 1st order filters for the control signals. These are
    useful as changing the control signals too quickly can lead to
    some strange and unrealistic things - such as a kink in the wing
    kinematics. Ideally I will deal with this in the sim_fly class,
    but this kludge will have to do for now. 
    """
    def __init__(self, ctl_sig, tc = (5.0,5.0,5.0,5.0)):
        self.tc = tc
        self.filters = [Filt1(x,tc[i]) for i,x in enumerate(ctl_sig)]

    def update(self,ctl_sig,dt):
        """
        Update filter set and return filtered values
        """
        filt_ctl_sig = [self.filters[i].update(x,dt) for i,x in enumerate(ctl_sig)]
        return tuple(filt_ctl_sig)
    
class ypr_ctl:
    """
    This class implements a simple proportional yaw/pitch/roll rate
    controller.
    """

    def __init__(self, ang_vel, p_gain = 3.0, r_gain = 1.0, y_gain = 1.0, av_filt_tc = 0.0):
        """
        Initialize controller.

        Input:
            ang_vel = current angular velocity of the fly
            
        keyword args:
            p_gain = pitch rate controller gain
            r_gain = roll rate controller gain
            y_gain = yaw rate controller gain
            av_filt_tc = angular velocity filter time constant.            
        """
        self.set_pt = ang_vel  
        self.p_gain = p_gain
        self.r_gain = r_gain
        self.y_gain = y_gain
        self.filter = Filt1(ang_vel, av_filt_tc)

    def get_sig(self, ang_vel, dt):
        """
        Get wing deformation control signals.

        Inputs:
            ang_vel = angular velocity of the fly
            dt = update time step

        Output:
            y_sig, p_sig, r_sig = yaw, pitch, roll wing
            deformation control signals
        """
        ang_vel_filt = self.filter.update(ang_vel, dt)
        # Get control values
        y_sig = self.y_gain*(self.set_pt[0] - ang_vel_filt[0]) # Yaw
        p_sig = self.p_gain*(self.set_pt[1] - ang_vel_filt[1]) # Pitch
        r_sig = self.r_gain*(self.set_pt[2] - ang_vel_filt[2]) # Roll
        return y_sig, p_sig, r_sig    

class alt_ctl:
    """
    This class implements a simple proportional derivative altitude
    controller.
    """

    def __init__(self, alt, prop_gain = 1.5, deriv_gain = 15.0, filt_tc = 0.0):
        """
        Initialize altitude controller.

        Input:
            alt = current altitude (mm)

        Keyword Args:
            prop_gain  = proportional gain
            deriv_gain = derivative gain
            filt_tc = altitude filter time constant
        """
        self.prop_gain = prop_gain 
        self.deriv_gain = deriv_gain 
        self.filter = Filt1(alt,filt_tc) 
        self.set_pt = alt
        self.err = 0.0 

    def get_sig(self, alt, dt):
        """
        Get altitude control signal given the current altitude and the
        update time step. 

        Inputs:
            alt = current altitud of the fly
            dt  = update time step

        Ouput:
            sig = throttle control signal
        """
        alt_filt = self.filter.update(alt,dt)
        # Compute altitude error and rate of change
        new_err =  self.set_pt - alt_filt
        derr_dt = (new_err - self.err)/dt
        # Update stored error 
        self.err = new_err
        # Compute throttle control signal
        t_sig = self.prop_gain*self.err + self.deriv_gain*derr_dt
        return t_sig

class vz_ctl:
    """
    The class encapsulates a simple vertical velocity (vz) controller.
    This is currently implemented as a proportional controller - however
    it might work better as proportional integral. 
    """

    def __init__(self, vz, prop_gain=20.0, filt_tc = 0.0):
        """
        Initialize vertical velocity controller.

        Input:
           vz = current (estimated) vertical velocity

        keyword args:
           prop_gain = proportional gain
           filt_tc = filter time constant
           
        """
        self.prop_gain = prop_gain
        self.set_pt = vz
        self.filter = Filt1(vz,filt_tc)

    def get_sig(self, vz, dt):
        """
        Compute throttle control signal.

        Input:
            vz = vertical velocity
            dt = update interval

        Output:
            sig = throttle control signal
        """
        vz_filt = self.filter.update(vz,dt)
        err = self.set_pt - vz_filt
        sig = self.prop_gain*err
        return sig

class fvel_ctl:
    """
    This class implements a simple pitch rate based forward velocity
    and altitude controller. It uses alt_ctl and pry_ctl as nested
    controllers. 
    """

    def __init__(self, fvel, ang_vel, alt,
                 prop_gain = 0.1,
                 deriv_gain = 10.0,
                 integ_gain = 0.0000,
                 filt_tc = 1.0,
                 deriv_clamp = 0.1,
                 integ_clamp = 2e3
                 ):
        """
        Initialize controller.

        Inputs:
           fvel = current forward velocity
           ang_vel = angular velocity of the fly
           alt = altitude of the fly

        keyword args:
           prop_gain = proportional gain
           deriv_gain = derivative gain
           integ_gain = integral gain 
           filt_tc = filter time constant
           deriv_clamp = derivative clamp
           integ_clamp = integrator clamp
           
        """
        self.err = 0.0
        self.integ_err = 0.0
        self.set_pt = fvel 
        self.err_filter = Filt1(fvel,filt_tc) 
        self.prop_gain = prop_gain
        self.deriv_gain = deriv_gain
        self.deriv_clamp = deriv_clamp
        self.integ_gain = integ_gain
        self.integ = Integrator(0.0, integ_clamp)

        # Sub controllers - altitude yaw/pitch/roll rate
        self.alt_ctl = alt_ctl(alt)
        self.ypr_ctl = ypr_ctl(ang_vel) 
  
    def get_sig(self, fvel, ang_vel, alt, dt):
        """
        Compute yaw, pitch, roll and thottle control signals.

        Inputs:
           fvel = forward velocity
           ang_vel = angular velocity
           alt = altitude
           dt = time step

        Output:
            y_sig, p_sig, r_sig, t_sig = yaw, pitch, roll, throttle
            control signals
        """
        # Compute error and low pass filter
        err = self.set_pt - fvel
        ferr_last = self.err_filter.val
        ferr = self.err_filter.update(err,dt)
        self.err = ferr

        # Get derivative error
        dferr_dt = (ferr - ferr_last)/dt
        self.derr_dt = dferr_dt
        
        # Clamp derivative of error term
        if self.deriv_clamp:
            dferr_dt = min([dferr_dt,  self.deriv_clamp])
            dferr_dt = max([dferr_dt, -self.deriv_clamp])
        
        # Get integral error 
        integ_err = self.integ.update(err)
        self.integ_err = integ_err
        
        # Compute pitch rate set point
        p_set_pt = self.prop_gain*ferr + self.deriv_gain*dferr_dt + self.integ_gain*integ_err
        self.ypr_ctl.set_pt = (0.0,p_set_pt, 0.0)
        
        # Get pitch/roll/yaw mode control signals using pry controller
        y_sig, p_sig, r_sig = self.ypr_ctl.get_sig(ang_vel,dt)
        
        # Get throttle control signal using atlitude controller 
        t_sig = self.alt_ctl.get_sig(alt,dt)
        return y_sig, p_sig, r_sig, t_sig 


class svel_ctl:
    """
    Simple proportional derivative side-slip controller - used to
    prevent side slip during long flights
    """

    def __init__(self, svel, prop_gain = 0.01, deriv_gain = 0.5, filt_tc = 1.0):
        self.prop_gain = prop_gain
        self.deriv_gain = deriv_gain
        self.set_pt = svel
        self.filter = Filt1(0.0,filt_tc)
        self.err = 0.0

    def get_sig(self, svel, dt):
        self.filter.update(svel,dt)
        err = self.set_pt - svel
        err_filt = self.filter.update(err,dt)
        derr_dt = (err_filt - self.err)/dt
        self.err = err_filt
        sig = self.prop_gain*err + self.deriv_gain*derr_dt
        return -sig
