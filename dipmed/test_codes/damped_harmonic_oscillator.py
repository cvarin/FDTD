#! /usr/bin/env python                                                             
# -*- coding: utf-8 -*-

#########################################################################################
from numpy import *
import pylab as p
from scipy.integrate import odeint

#########################################################################################
# The equations of motion are integrated and the dynamics is plotted
tmax = 100
numpoints = 300

omega = 1.0          # Natural frequency
b = 0.00001          # Damping factor
x_0 = 0.0            # initials conditions
v_0 = 0.0
X0 = array([x_0,v_0]) 

t = linspace(0,tmax,numpoints)
dt = float(tmax)/float(numpoints)

#########################################################################################
# Source
def F_sine_gauss(t):
     return exp(-(t-40.0)**2/200.0)*cos(t*omega/2.0)

#########################################################################################
F = F_sine_gauss

################### Solution of the ode set #############################################
def dX_dt(X,t=0):
    return array([X[1], - b*X[1] -omega**2*F(t)])
X, infodict = odeint(dX_dt,X0,t,full_output=True)
infodict['message']

################### The result is plotted ###############################################
x,v = X.T
time = t/omega
fs=16
p.figure(figsize=(16.0,8.5))
p.title(r"Time evolution of a damped pendulum",fontsize=fs)
p.plot(time,x/x.max(),'r-',label=r"Adaptive integrator")
p.plot(time,F(time)/F(time).max(),'k',label="source")
p.xlabel(r"Time (in units of $\omega$)",fontsize=fs)
p.ylabel(r"$x$",fontsize=1.5*fs)
p.grid(True)
p.legend(loc='best')

#########################################################################################
#p.savefig('fig.pdf')
p.show()
################### End of file #########################################################
