#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pylab import *
from scipy.integrate import quadrature
from scipy.special import erf

#########################################################################################
# Integrals to find <p>
def U(theta,p0,a1,a3,Ex):
     return -p0*Ex*cos(theta)-0.5*(a1 + (a3 - a1)*cos(theta)**2)*Ex**2
     
def Boltzmann(theta,p0,a1,a3,Ex,T):
     return exp(-U(theta,p0,a1,a3,Ex)/(kB*T))
     
def kernel(theta,p0,a1,a3,Ex,T):
     return sin(theta)*Boltzmann(theta,p0,a1,a3,Ex,T)

def func_cos_avg(theta,p0,a1,a3,Ex,T):
     return cos(theta)*kernel(theta,p0,a1,a3,Ex,T)
     
def func_cos_sqr_avg(theta,p0,a1,a3,Ex,T):
     return cos(theta)**2*kernel(theta,p0,a1,a3,Ex,T)

def cos_square_avg(p0,a1,a3,Ex,T,tol=1.0e-3,maxiter=50):
     args=(p0,a1,a3,Ex,T)
     norm,error_norm = quadrature(kernel,0,pi,args=args,tol=tol,maxiter=maxiter)
     cos_square_avg,error_cos_sqr = quadrature(func_cos_sqr_avg,0,pi,args=args,tol=tol,maxiter=maxiter)
     cos_square_avg *= 1.0/norm
     return cos_square_avg

def p_avg(p0,a1,a3,Ex,T,tol=1.0e-3,maxiter=50):
     args=(p0,a1,a3,Ex,T)
     norm,error_norm = quadrature(kernel,0,pi,args=args,tol=tol,maxiter=maxiter)
     cos_avg,error_cos = quadrature(func_cos_avg,0,pi,args=args,tol=tol,maxiter=maxiter)
     cos_square_avg,error_cos_sqr = quadrature(func_cos_sqr_avg,0,pi,args=args,tol=tol,maxiter=maxiter)
     cos_avg *= 1.0/norm
     cos_square_avg *= 1.0/norm
     return p0*cos_avg + (a1 + (a3 - a1)*cos_square_avg)*Ex 

#########################################################################################
# With and without pemanent dipole moment...
#kB = 1.3806503e-23
kB = 1.0
T = 1.0
p0 = 1.0
a1 = 0.05
a3 = 0.1

N = 50
MAX = 15
x = linspace(0.0,MAX,N)
E = linspace(0.0,MAX,N)

#########################################################################################
# Numerical solutions
f_num1 = zeros(N)
f_num2 = zeros(N)
f_num3 = zeros(N)
f_num4 = zeros(N)
f_num5 = zeros(N)
f_num6 = zeros(N)
for i in range(N):
     f_num1[i] = p_avg(p0,a1,a3,E[i],T,1.0e-3,100)
     f_num2[i] = p_avg(p0,0.0,0.0,E[i],T,1.0e-3,100)
     f_num3[i] = p_avg(0.0,a1,a3,E[i],T,1.0e-3,100)
     f_num4[i] = p_avg(0.0,0.0,a3,E[i],T,1.0e-3,100)
     f_num5[i] = p_avg(0.0,a3,a3,E[i],T,1.0e-3,100)
     f_num6[i] = p_avg(p0,a3,a3,E[i],T,1.0e-3,100)

#########################################################################################
# Analytical solutions
# Hook and Hall, Eq. (9.23)
f_ana1 = p0*(1.0/tanh(p0*x)-1.0/(p0*x)) 
# Alignement_boyd.wxmx
J = 0.5*(a3 - a1)*E**2/(kB*T)
cos_sqr_avg = (1j*sqrt(J)*(exp(J)/J-(sqrt(pi)*erf(1j*sqrt(J)))/(2.0*1j*sqrt(J)*J)))/(sqrt(pi)*erf(1j*sqrt(J)))
f_ana2 = (a1 + (a3 - a1)*cos_sqr_avg)*E

#########################################################################################
figure(figsize=(16.0,8.5))
plot(E,f_num5,label=r"Non-polar - isotropic response (p$_0$= 0, a$_\bot$= a$_\parallel$)")
plot(E,f_num3,label=r"Non-polar - anisotropic response (p$_0$= 0, a$_\bot\neq$ a$_\parallel$)")
#plot(E,f_num4,label="p = 0.0, a1 = 0.0, a3 = 0.1")
plot(E,f_num2,label=r"Polar - no electronic response, (p$_0\neq$0, a$_\bot$= a$_\parallel$= 0)")
plot(E,f_num1,"k",label=r"Polar - anisotropic response, (p$_0\neq$0, a$_\bot\neq$ a$_\parallel$)")
#plot(E,f_num6,label=r"Polar and isotropic, (p=1.0, a$_\bot$= a$_\parallel$= 0.1)")
plot(x,f_ana1,'.',label="Analytical [Jackson, 3rd ed., Eq. (4.80)]")
plot(x,f_ana2,'.',label="Analytical [Boyd, 3rd ed., Eqs. (4.4.10) and (4.4.11)]")
#axhline(y=1.0,ls="--")
#ylim(0,1.2)
xlim(0,MAX)
xlabel(r"$E$",size=18)
ylabel(r"$<p>$",size=18)
legend(loc="upper left")
grid(True)

#########################################################################################
show()

#########################################################################################
