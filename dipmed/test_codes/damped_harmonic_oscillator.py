#! /usr/bin/env python                                                             
# -*- coding: utf-8 -*-

#########################################################################################
from numpy import *
import pylab as p
from scipy.integrate import odeint

#########################################################################################
# Definition of the equations and parameters 
Omega = 1.0              # Natural frequency
b = 0.1                  # Damping factor
theta_0 = pi*1.0/2.0     # Initial angle
v_0     = 0.0            # Initial angular velocity

def dX_dt(X, t=0):
    """ Returns the variations of both the angle and velocity. """
    return array([ X[1], - b*X[1] -Omega**2*sin(X[0])])
    
#########################################################################################
# The pendulums equations of motion are integrated and the dynamics is plotted
t = linspace(0, 100*Omega, 1000)           # time
X0 = array([theta_0,v_0])                  # initials conditions

################### Solution of the ode set #############################################
X, infodict = odeint(dX_dt, X0, t, full_output=True)
infodict['message']                     # >>> 'Integration successful.'

################### The result is plotted ###############################################
theta, v = X.T
time = t/Omega
fs=16
f1 = p.figure()
# Angle
p.subplot(2,1,1)
p.title(r"Time evolution of a damped pendulum",fontsize=fs)
p.plot(time, theta, 'r-', label=r"$\theta$")
p.plot(time, theta_0*exp(-b**1.3*t), 'g-', label=r"damping fit")
p.ylabel(r"$\theta$ (rad)",fontsize=fs)
p.grid()
p.legend(loc='best')

# Velocity
p.subplot(2,1,2)
p.plot(time, v    , 'b-', label=r"$v$")
p.grid()
p.xlabel(r"Time (in units of $\Omega$)",fontsize=fs)
p.ylabel(r"Angular velocity (rad/s)",fontsize=fs)

#########################################################################################
#p.savefig('fig.pdf')
p.show()
################### End of file #########################################################
