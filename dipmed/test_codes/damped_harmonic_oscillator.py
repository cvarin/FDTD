#! /usr/bin/env python                                                             
# -*- coding: utf-8 -*-

#########################################################################################
# The equation for the damped and driven harmonic oscillator:
# d^2 x/dT^2 + (gamma/omega_0) d x/dt + x = f(T)
# where T = omega_0 t. It can be rewritten with two equations
# d x/dt = v
# d v/dt = f(T) - x - (gamma/omega_0) v
# Ref. : http://en.wikipedia.org/wiki/Harmonic_oscillator

#########################################################################################
from numpy import *
import pylab as p
from scipy.integrate import odeint

#########################################################################################
# The equations of motion are integrated and the dynamics is plotted
tmax = 100
numpoints = 1000

b = 0.15         # Damping factor, effectively gamma/omega_0
x_0 = 0.0        # initials conditions
v_0 = 0.0
X0 = array([x_0,v_0]) 

t = linspace(0,tmax,numpoints)
dt = float(tmax)/float(numpoints)

#########################################################################################
# Source
def F_sine_gauss(t):
     return exp(-(t-40.0)**2/200.0)*cos(t/2.0)
     
def F_constant(t):
     return 10.0;

#########################################################################################
F = F_sine_gauss
#F = F_constant

################### Solution of the ode set #############################################
def dX_dt(X,t=0):
    return array([X[1], F(t) - X[0] - b*X[1]])
X, infodict = odeint(dX_dt,X0,t,full_output=True)
infodict['message']

#########################################################################################
x_fd = zeros(numpoints)
common = 1.0/(1 + 0.5*b*dt)
a1 = (2.0 - dt**2)*common
a2 = -(1.0 - 0.5*b*dt)*common
a3 = dt**2*common
for i in range(2,numpoints):
     x_fd[i] = a1*x_fd[i-1] + a2*x_fd[i-2] + a3*F(t[i])

################### The result is plotted ###############################################
x,v = X.T
fs=16
p.figure(figsize=(16.0,8.5))
p.title(r"Time evolution of a damped oscillator",fontsize=fs)
p.plot(t,x,'r-',label=r"Adaptive integrator")
p.plot(t,x_fd,'b.',label=r"Finite differences")
p.plot(t,F(t),'k',label="source")
p.xlabel(r"Time (in units of $\omega$)",fontsize=fs)
p.ylabel(r"$x$",fontsize=1.5*fs)
p.grid(True)
p.legend(loc='best')

#########################################################################################
#p.savefig('fig.pdf')
p.show()
################### End of file #########################################################
