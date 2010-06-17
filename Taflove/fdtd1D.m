%***********************************************************************
%     1-D FDTD code
%***********************************************************************
%
%     Program author: Susan C. Hagness
%                     Department of Electrical and Computer Engineering
%                     University of Wisconsin-Madison
%                     1415 Engineering Drive
%                     Madison, WI 53706-1691
%                     hagness@engr.wisc.edu
%
%     Copyright 2005
%
%     This MATLAB M-file implements a finite-difference time-domain
%     solution of Maxwell's curl equations over a one-dimensional space
%     lattice comprised of uniform grid cells.
%
%     To illustrate the algorithm, a 1-GHz sinusoidal wave propagating 
%     in a nonpermeable lossy medium (epsr=1.0, sigma=5.0e-3 S/m) is 
%     modeled.  The simplified finite difference system for nonpermeable
%     media (discussed in Section 3.6.6 of the text) is implemented.
%
%     The grid resolution (dx = 1.5 cm) is chosen to provide 20
%     samples per wavelength.  The Courant factor S=c*dt/dx is set to
%     the stability limit (S=1).  In 1-D, this is the "magic time step."
%     The total number of time steps (nmax=240) corresponds to a physical 
%     time of 12 ns.
%
%     The grid is terminated with electric-field components at the far-left
%     (i=1) and far-right (i=ie) boundaries.  The sinusoidal wave is
%     launched by an electric-field hard-source condition at i=1 (see
%     Eq. 5.1 in the text).  The simplest radiation boundary condition
%     for plane wave propagation is used to update the electric field
%     at i=ie:  
%
%                      Ez(ie,n+1) = Ez(ie-1,n)
%
%     To execute this M-file, type "fdtd1D" at the MATLAB prompt.
%
%     This code has been tested in the following Matlab environments:
%     Matlab version 6.1.0.450 Release 12.1 (May 18, 2001)
%     Matlab version 6.5.1.199709 Release 13 Service Pack 1 (August 4, 2003)
%     Matlab version 7.0.0.19920 R14 (May 6, 2004)
%     Matlab version 7.0.1.24704 R14 Service Pack 1 (September 13, 2004)
%     Matlab version 7.0.4.365 R14 Service Pack 2 (January 29, 2005)  
%
%***********************************************************************

clear

%***********************************************************************
%     Fundamental constants
%***********************************************************************

cc=2.99792458e8;            %speed of light in free space
muz=4.0*pi*1.0e-7;          %permeability of free space
epsz=1.0/(cc*cc*muz);       %permittivity of free space

freq=1.0e+9;                %frequency of source excitation
lambda=cc/freq;             %wavelength of source excitation
omega=2.0*pi*freq;

%***********************************************************************
%     Grid parameters
%***********************************************************************

ie=200;                     %number of Ez samples in grid
ih=ie-1;                    %number of Hy samples in grid

dx=lambda/20.0;             %space increment of 1-D lattice
dt=dx/cc;                   %time step (S=1.0)
omegadt=omega*dt;

nmax=round(12.0e-9/dt);     %total number of time steps

x=linspace(dx,ih*dx,ih);

%***********************************************************************
%     Material parameters
%***********************************************************************

eps=1.0;
sig=5.0e-3;

%***********************************************************************
%     Updating coefficients for space region with nonpermeable media
%***********************************************************************

scfact=dt/muz/dx;

ca=(1.0-(dt*sig)/(2.0*epsz*eps))/(1.0+(dt*sig)/(2.0*epsz*eps));
cb=scfact*(dt/epsz/eps/dx)/(1.0+(dt*sig)/(2.0*epsz*eps));

%***********************************************************************
%     Field arrays
%***********************************************************************

ez(1:ie)=0.0;
hy(1:ih)=0.0;

%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************

for n=1:nmax

%***********************************************************************
%     Update electric fields
%***********************************************************************

ez(1)=scfact*sin(omegadt*n);

rbc=ez(ih);
ez(2:ih)=ca*ez(2:ih)+cb*(hy(2:ih)-hy(1:ih-1));
ez(ie)=rbc;

%***********************************************************************
%     Update magnetic fields
%***********************************************************************

hy(1:ih)=hy(1:ih)+ez(2:ie)-ez(1:ih);

%***********************************************************************
%     Visualize fields
%***********************************************************************

rtime=num2str(round(n*dt/1.0e-9));

subplot(2,1,1),plot(x,ez(1:ih)/scfact,'r'),axis([0 3 -1 1]);
title(['time = ',rtime,' ns']); ylabel('E_z');
subplot(2,1,2),plot(x,hy,'b'),axis([0 3 -3.0e-3 3.0e-3]);
xlabel('x (meters)');ylabel('H_y');

pause(0.05)

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************

end
