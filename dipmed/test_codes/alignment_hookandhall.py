#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pylab import *
from scipy.integrate import quadrature

###############################################################################
# Evaluate Eqs. (9.23) of Hook and Hall, Solide state physics, 2nd Ed.
tol = 1.5e-8
nint = 100
N = 40
x = linspace(0.0,25,N)
f = (1.0/tanh(x)-1.0/x)
f_series = 1.0/3.0*x

###############################################################################
# Integrate numerically Eq. (9.22) of Hook and Hall, Solide state physics, 2nd Ed.
def func_norm(u,x,):
     return exp(u*x)
def func_avg(u,x,):
     return u*func_norm(u,x)
def p_avg(x):
     global tol
     global nint
     norm,error1 = quadrature(func_norm,-1,1,args=(x,),tol=tol,maxiter=nint)
     avg,error2 = quadrature(func_avg,-1,1,args=(x,),tol=tol,maxiter=nint)
     avg *= 1.0/norm
     return avg
f_num = zeros(N)
for i in range(0,N):
     f_num[i] = p_avg(x[i])

###############################################################################
# Plot
figure(figsize=(16,8.5))
subplot(2,1,1)
plot(x,f,label="Analytical")
plot(x,f_series,label="First order")
plot(x,f_num,'.',label="Numerical")
ylabel(r"$\left<\cos\theta\right>$",size=18)
axhline(y=1,ls='--')
grid(True)
ylim(0,1.2)
legend(loc="lower right")

subplot(2,1,2)
plot(x,(f_num-f)/f)
grid(True)
xlabel(r"$x=p_0E/kT$",size=18)
ylabel(r"(num - ana)/ana",size=18)

###############################################################################
# Check the denominator
#norm_ana = exp(x)/x-exp(-x)/x
#norm_num = zeros(N)
#for i in range(0,N):
     #norm_num[i],error = quadrature(func_norm,-1,1,args=(x[i],),tol=tol,maxiter=nint)
#figure()
#plot(x,1.0/norm_ana,label="analytical")
#plot(x,1.0/norm_num,".",label="numerical")
#xlabel(r"$x$",size=18)
#ylabel(r"$\left[\int_{-1}^1\exp(ux)du\right]^{-1}$",size=16)
#grid(True)
#legend()

###############################################################################
show()
###############################################################################
