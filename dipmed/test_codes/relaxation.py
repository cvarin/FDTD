#! /usr/bin/env python                                                             
# -*- coding: utf-8 -*-

#########################################################################################
from numpy import linspace,squeeze,zeros
import pylab as p
from scipy.integrate import odeint
from scipy.special import erf,exp,sin

#########################################################################################
# The equations of motion are integrated and the dynamics is plotted
tmax = 200
numpoints = 100
t = linspace(0, tmax, numpoints)
dt = float(tmax)/float(numpoints)

X0 = 0.0
b = 1.0/10.0
A = 5.0

#########################################################################################
# Different sources
def F_erf(t):
     t1 = (t-20.0)/1.0e-1
     t2 = (t-50.0)/1.0e-1
     return A*(erf(t1) - erf(t2))
     
#########################################################################################
def F_gauss(t):
     dt = 5.0
     return A*exp(-(t-30)**2/dt**2)
     
#########################################################################################
def F_sine(t):
     omega = 7.0e-1
     return A*sin(t*omega)

#########################################################################################
def F_gauss_sine(t):
     dt = 10.0
     omega = 1.5
     return A*exp(-(t-40)**2/dt**2)*sin(t*omega)

#########################################################################################
F = F_erf
    
################### Solution of the ode set #############################################
def dX_dt(X, t=0):
    """ Returns the variations of both the angle and velocity. """
    return - b*X + b*F(t)
X, infodict = odeint(dX_dt, X0, t, full_output=True)
infodict['message']                     # >>> 'Integration successful.'

################### Solution with finite differences ####################################
x_fd = zeros(numpoints)
for i in range(1,numpoints):
     gam = 2.0/(dt*b)
     x_fd[i] = F(t[i])/(1+gam) - (1-gam)/(1+gam)*x_fd[i-1]
     
x_fd2 = zeros(numpoints)
for i in range(1,numpoints):
     gam = 2.0/(dt*b)
     x_fd2[i] = (F(t[i]) + F(t[i-1]))/2.0/(1+gam) - (1-gam)/(1+gam)*x_fd2[i-1]
     
x_fd3 = zeros(numpoints)
for i in range(1,numpoints):
     gam = 2.0/(dt*b)
     x_fd3[i] = F( 0.5*(t[i] + t[i-1]) )/(1+gam) - (1-gam)/(1+gam)*x_fd3[i-1]

################### The result is plotted ###############################################
x = X.T
fs=16
p.figure(figsize=(16.0,8.5))
p.plot(t,x_fd/x_fd.max(),'.', label=r"Source at $E^{n+1}$")
p.plot(t,x_fd2/x_fd2.max(),'.',label=r"Source at $(E^{n+1} + E^n)/2$")
p.plot(t,x_fd3/x_fd3.max(),'.',label=r"Source at $E^{n+1/2}$")
p.plot(t,squeeze(x)/squeeze(x).max(), '-', label=r"adaptive")
p.plot(t,F(t)/F(t).max(),'k',label="source")
#p.axhline(y=1.0/exp(1.0),ls='--',label=r"$1/e$")
#p.axvline(x=t[squeeze(x).argmax()],ls='--')
p.ylabel(r"amplitude (a.u.)",fontsize=fs)
p.grid(True)
p.legend(loc='best')
p.xlabel(r"Time",fontsize=fs)
p.ylim(-0.2,1.2)

#########################################################################################
#p.savefig('fig.pdf')
p.show()
################### End of file #########################################################
