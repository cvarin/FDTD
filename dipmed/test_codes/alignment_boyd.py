#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pylab import *
from scipy.special import erf
from scipy.integrate import quadrature

###############################################################################
# Evaluate the solution of Eq. (4.4.13) of R. W. Boyd, Nonlinear Optics, 3rd Ed.
tol = 1.0e-3
nint = 50
N = 50
J = linspace(0,25,N)
f_series = 1.0/3.0 + 4.0*J/45.0 + 8.0*J**2/945

###############################################################################
# Analytical solution (see alignement_boyd.wxmx)
f = (1j*sqrt(J)*(exp(J)/J-(sqrt(pi)*erf(1j*sqrt(J)))/(2.0*1j*sqrt(J)*J)))/(sqrt(pi)*erf(1j*sqrt(J)))

###############################################################################
# Numerical solution
def func_norm(theta,J,):
     return sin(theta)*exp(J*cos(theta)**2)
def func_avg(theta,J,):
     return cos(theta)**2*func_norm(theta,J)
def cos_square_avg(J):
     global tol
     global nint
     norm,error1 = quadrature(func_norm,0,pi,args=(J,),tol=tol,maxiter=nint)
     avg,error2 = quadrature(func_avg,0,pi,args=(J,),tol=tol,maxiter=nint)
     avg *= 1.0/norm
     return avg
f_num = zeros(N)
for i in range(0,N):
     f_num[i] = cos_square_avg(J[i])

###############################################################################
# Plots
figure(figsize=(16,8.5))
subplot(2,1,1)
plot(J,f-1.0/3.0,label="Analytical")
plot(J,f_series-1.0/3.0,label="Second order")
plot(J,f_num-1.0/3.0,'.',label="Numerical")
axhline(y=2.0/3.0,ls='--')
ylabel(r"$\left<\cos^2\theta\right>-\frac{1}{3}$",size=18)
grid(True)
ylim(0,0.8)
legend(loc="lower right")

subplot(2,1,2)
plot(J,(f_num-f)/f)
grid(True)
xlabel(r"$J=\frac{1}{2}(\alpha_3-\alpha_1)E^2/kT$",size=18)
ylabel(r"(num - ana)/ana",size=18)

###############################################################################
show()

###############################################################################