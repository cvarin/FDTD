%***********************************************************************
%     3-D FDTD code with UPML absorbing boundary conditions
%***********************************************************************
%
%     Program author: Keely J. Willis, Graduate Student
%                     UW Computational Electromagnetics Laboratory
%                           Director: Susan C. Hagness
%                     Department of Electrical and Computer Engineering
%                     University of Wisconsin-Madison
%                     1415 Engineering Drive
%                     Madison, WI 53706-1691
%                     kjwillis@wisc.edu
%
%     Copyright 2005
%
%     This MATLAB M-file implements the finite-difference time-domain
%     solution of Maxwell's curl equations over a three-dimensional
%     Cartesian space lattice comprised of uniform cubic grid cells.
%     
%     The dimensions of the computational domain are 8.2 cm
%     (x-direction), 3.4 cm (y-direction), and 3.2 cm (z-direction).  
%     The grid is terminated with UPML absorbing boundary conditions.
%
%     An electric current source comprised of two collinear Jz components
%     (realizing a Hertzian dipole) excites a radially propagating wave.  
%     The current source is located in the center of the grid.  The 
%     source waveform is a differentiated Gaussian pulse given by 
%          J(t)=J0*(t-t0)*exp(-(t-t0)^2/tau^2), 
%     where tau=50 ps.  The FWHM spectral bandwidth of this zero-dc-
%     content pulse is approximately 7 GHz. The grid resolution 
%     (dx = 2 mm) was chosen to provide at least 10 samples per 
%     wavelength up through 15 GHz.
%
%     To execute this M-file, type "fdtd3D_UPML" at the MATLAB prompt.  
%
%     This code has been tested in the following Matlab environments:
%     Matlab version 6.1.0.450 Release 12.1 (May 18, 2001)
%     Matlab version 6.5.1.199709 Release 13 Service Pack 1 (August 4, 2003)
%     Matlab version 7.0.0.19920 R14 (May 6, 2004)
%     Matlab version 7.0.1.24704 R14 Service Pack 1 (September 13, 2004)
%     Matlab version 7.0.4.365 R14 Service Pack 2 (January 29, 2005)
%
%     Note: if you are using Matlab version 6.x, you may wish to make
%     one or more of the following modifications to this code: 
%       --uncomment line numbers 485 and 486
%       --comment out line numbers 552 and 561
%
%***********************************************************************

clear

%***********************************************************************
%     Fundamental constants
%***********************************************************************

cc=2.99792458e8;
muz=4.0*pi*1.0e-7;
epsz=1.0/(cc*cc*muz);
etaz=sqrt(muz/epsz);

%***********************************************************************
%     Material parameters 
%***********************************************************************

mur=1.0;
epsr=1.0;
eta=etaz*sqrt(mur/epsr);

%***********************************************************************
%     Grid parameters
%
%     Each grid size variable name describes the number of sampled points 
%     for a particular field component in the direction of that component.
%     Additionally, the variable names indicate the region of the grid 
%     for which the dimension is relevant.  For example, ie_tot is the 
%     number of sample points of Ex along the x-axis in the total 
%     computational grid, and jh_bc is the number of sample points of Hy 
%     along the y-axis in the y-normal UPML regions.
%
%***********************************************************************

ie=41;          % Size of main grid
je=17;
ke=16;
ih=ie+1;
jh=je+1;   
kh=ke+1;   

upml=10;        % Thickness of PML boundaries
ih_bc=upml+1;
jh_bc=upml+1;
kh_bc=upml+1;

ie_tot=ie+2*upml;          % Size of total computational domain
je_tot=je+2*upml;        
ke_tot=ke+2*upml;        
ih_tot=ie_tot+1;
jh_tot=je_tot+1;          
kh_tot=ke_tot+1;          

is=round(ih_tot/2);         % Location of z-directed current source
js=round(jh_tot/2);
ks=round(ke_tot/2);

%***********************************************************************
%     Fundamental grid parameters
%***********************************************************************

delta=0.002;
dt=delta*sqrt(epsr*mur)/(2.0*cc);
nmax=100;

%***********************************************************************
%     Differentiated Gaussian pulse excitation
%***********************************************************************

rtau=50.0e-12;
tau=rtau/dt;
ndelay=3*tau;
J0=-1.0*epsz;

%***********************************************************************
%     Initialize field arrays
%***********************************************************************

ex=zeros(ie_tot,jh_tot,kh_tot);
ey=zeros(ih_tot,je_tot,kh_tot);
ez=zeros(ih_tot,jh_tot,ke_tot);
dx=zeros(ie_tot,jh_tot,kh_tot);
dy=zeros(ih_tot,je_tot,kh_tot);
dz=zeros(ih_tot,jh_tot,ke_tot);

hx=zeros(ih_tot,je_tot,ke_tot);
hy=zeros(ie_tot,jh_tot,ke_tot);
hz=zeros(ie_tot,je_tot,kh_tot);
bx=zeros(ih_tot,je_tot,ke_tot);
by=zeros(ie_tot,jh_tot,ke_tot);
bz=zeros(ie_tot,je_tot,kh_tot);

%***********************************************************************
%     Initialize update coefficient arrays
%***********************************************************************

C1ex=zeros(size(ex));
C2ex=zeros(size(ex));
C3ex=zeros(size(ex));
C4ex=zeros(size(ex));
C5ex=zeros(size(ex));
C6ex=zeros(size(ex));

C1ey=zeros(size(ey));
C2ey=zeros(size(ey));
C3ey=zeros(size(ey));
C4ey=zeros(size(ey));
C5ey=zeros(size(ey));
C6ey=zeros(size(ey));

C1ez=zeros(size(ez));
C2ez=zeros(size(ez));
C3ez=zeros(size(ez));
C4ez=zeros(size(ez));
C5ez=zeros(size(ez));
C6ez=zeros(size(ez));

D1hx=zeros(size(hx));
D2hx=zeros(size(hx));
D3hx=zeros(size(hx));
D4hx=zeros(size(hx));
D5hx=zeros(size(hx));
D6hx=zeros(size(hx));

D1hy=zeros(size(hy));
D2hy=zeros(size(hy));
D3hy=zeros(size(hy));
D4hy=zeros(size(hy));
D5hy=zeros(size(hy));
D6hy=zeros(size(hy));

D1hz=zeros(size(hz));
D2hz=zeros(size(hz));
D3hz=zeros(size(hz));
D4hz=zeros(size(hz));
D5hz=zeros(size(hz));
D6hz=zeros(size(hz));

%***********************************************************************
%     Update coefficients, as described in Section 7.8.2.
%
%     In order to simplify the update equations used in the time-stepping
%     loop, we implement UPML update equations throughout the entire
%     grid.  In the main grid, the electric-field update coefficients of 
%     Equations 7.91a-f and the correponding magnetic field update
%     coefficients extracted from Equations 7.89 and 7.90 are simplified
%     for the main grid (free space) and calculated below.
%
%***********************************************************************

C1=1.0;
C2=dt;
C3=1.0;
C4=1.0/2.0/epsr/epsr/epsz/epsz;
C5=2.0*epsr*epsz;
C6=2.0*epsr*epsz;

D1=1.0;
D2=dt;
D3=1.0;
D4=1.0/2.0/epsr/epsz/mur/muz;
D5=2.0*epsr*epsz;
D6=2.0*epsr*epsz;

%***********************************************************************
%     Initialize main grid update coefficients
%***********************************************************************

C1ex(:,jh_bc:jh_tot-upml,:)=C1;     
C2ex(:,jh_bc:jh_tot-upml,:)=C2;
C3ex(:,:,kh_bc:kh_tot-upml)=C3;
C4ex(:,:,kh_bc:kh_tot-upml)=C4;
C5ex(ih_bc:ie_tot-upml,:,:)=C5;
C6ex(ih_bc:ie_tot-upml,:,:)=C6;

C1ey(:,:,kh_bc:kh_tot-upml)=C1;
C2ey(:,:,kh_bc:kh_tot-upml)=C2;
C3ey(ih_bc:ih_tot-upml,:,:)=C3;
C4ey(ih_bc:ih_tot-upml,:,:)=C4;
C5ey(:,jh_bc:je_tot-upml,:)=C5;
C6ey(:,jh_bc:je_tot-upml,:)=C6;

C1ez(ih_bc:ih_tot-upml,:,:)=C1;
C2ez(ih_bc:ih_tot-upml,:,:)=C2;
C3ez(:,jh_bc:jh_tot-upml,:)=C3;
C4ez(:,jh_bc:jh_tot-upml,:)=C4;
C5ez(:,:,kh_bc:ke_tot-upml)=C5;
C6ez(:,:,kh_bc:ke_tot-upml)=C6;

D1hx(:,jh_bc:je_tot-upml,:)=D1;
D2hx(:,jh_bc:je_tot-upml,:)=D2;
D3hx(:,:,kh_bc:ke_tot-upml)=D3;
D4hx(:,:,kh_bc:ke_tot-upml)=D4;
D5hx(ih_bc:ih_tot-upml,:,:)=D5;
D6hx(ih_bc:ih_tot-upml,:,:)=D6;

D1hy(:,:,kh_bc:ke_tot-upml)=D1;
D2hy(:,:,kh_bc:ke_tot-upml)=D2;
D3hy(ih_bc:ie_tot-upml,:,:)=D3;
D4hy(ih_bc:ie_tot-upml,:,:)=D4;
D5hy(:,jh_bc:jh_tot-upml,:)=D5;
D6hy(:,jh_bc:jh_tot-upml,:)=D6;

D1hz(ih_bc:ie_tot-upml,:,:)=D1;
D2hz(ih_bc:ie_tot-upml,:,:)=D2;
D3hz(:,jh_bc:je_tot-upml,:)=D3;
D4hz(:,jh_bc:je_tot-upml,:)=D4;
D5hz(:,:,kh_bc:kh_tot-upml)=D5;
D6hz(:,:,kh_bc:kh_tot-upml)=D6;

%***********************************************************************
%     Fill in PML regions
% 
%     PML theory describes a continuous grading of the material properties
%     over the PML region.  In the FDTD grid it is necessary to discretize
%     the grading by averaging the material properties over a grid cell 
%     width centered on each field component.  As an example of the 
%     implementation of this averaging, we take the integral of the 
%     continuous sigma(x) in the PML region
%   
%         sigma_i = integral(sigma(x))/delta
%   
%     where the integral is over a single grid cell width in x, and is 
%     bounded by x1 and x2.  Applying this to the polynomial grading of 
%     Equation 7.60a produces
%
%         sigma_i = (x2^(m+1)-x1^(m+1))*sigmam/(delta*(m+1)*d^m)
%
%     where sigmam is the maximum value of sigma as described by Equation 
%     7.62. 
%         
%     The definitions of x1 and x2 depend on the position of the field 
%     component within the grid cell.  We have either
%
%         x1 = (i-0.5)*delta,  x2 = (i+0.5)*delta
%  
%     or
%  
%         x1 = (i)*delta,      x2 = (i+1)*delta
%
%     where i varies over the PML region.
%  
%***********************************************************************

rmax=exp(-16);  %desired reflection error, designated as R(0) in Equation 7.62 

orderbc=4;      %order of the polynomial grading, designated as m in Equation 7.60a,b

%   x-varying material properties
delbc=upml*delta;
sigmam=-log(rmax)*(orderbc+1.0)/(2.0*eta*delbc); 
sigfactor=sigmam/(delta*(delbc^orderbc)*(orderbc+1.0));
kmax=1;
kfactor=(kmax-1.0)/delta/(orderbc+1.0)/delbc^orderbc;

for i=1:upml
    
    % Coefficients for field components in the center of the grid cell
    x1=(upml-i+1)*delta;
    x2=(upml-i)*delta;
    sigma=sigfactor*(x1^(orderbc+1)-x2^(orderbc+1));
    ki=1+kfactor*(x1^(orderbc+1)-x2^(orderbc+1));
    facm=(2*epsr*epsz*ki-sigma*dt);
    facp=(2*epsr*epsz*ki+sigma*dt);

    C5ex(i,:,:)=facp;
    C5ex(ie_tot-i+1,:,:)=facp;
    C6ex(i,:,:)=facm;
    C6ex(ie_tot-i+1,:,:)=facm;
    D1hz(i,:,:)=facm/facp;
    D1hz(ie_tot-i+1,:,:)=facm/facp;
    D2hz(i,:,:)=2.0*epsr*epsz*dt/facp;
    D2hz(ie_tot-i+1,:,:)=2.0*epsr*epsz*dt/facp;
    D3hy(i,:,:)=facm/facp;
    D3hy(ie_tot-i+1,:,:)=facm/facp;
    D4hy(i,:,:)=1.0/facp/mur/muz;
    D4hy(ie_tot-i+1,:,:)=1.0/facp/mur/muz;

    % Coefficients for field components on the grid cell boundary
    x1=(upml-i+1.5)*delta;
    x2=(upml-i+0.5)*delta;
    sigma=sigfactor*(x1^(orderbc+1)-x2^(orderbc+1));
    ki=1.0+kfactor*(x1^(orderbc+1)-x2^(orderbc+1));
    facm=(2.0*epsr*epsz*ki-sigma*dt);
    facp=(2.0*epsr*epsz*ki+sigma*dt);

    C1ez(i,:,:)=facm/facp;
    C1ez(ih_tot-i+1,:,:)=facm/facp;
    C2ez(i,:,:)=2.0*epsr*epsz*dt/facp;
    C2ez(ih_tot-i+1,:,:)=2.0*epsr*epsz*dt/facp;
    C3ey(i,:,:)=facm/facp;
    C3ey(ih_tot-i+1,:,:)=facm/facp;
    C4ey(i,:,:)=1.0/facp/epsr/epsz;
    C4ey(ih_tot-i+1,:,:)=1.0/facp/epsr/epsz;
    D5hx(i,:,:)=facp;
    D5hx(ih_tot-i+1,:,:)=facp;
    D6hx(i,:,:)=facm;
    D6hx(ih_tot-i+1,:,:)=facm;
    
end

%   PEC walls
C1ez(1,:,:)=-1.0;
C1ez(ih_tot,:,:)=-1.0;
C2ez(1,:,:)=0.0;
C2ez(ih_tot,:,:)=0.0;
C3ey(1,:,:)=-1.0;
C3ey(ih_tot,:,:)=-1.0;
C4ey(1,:,:)=0.0;
C4ey(ih_tot,:,:)=0.0;

%   y-varying material properties
delbc=upml*delta;
sigmam=-log(rmax)*epsr*epsz*cc*(orderbc+1.0)/(2.0*delbc); 
sigfactor=sigmam/(delta*(delbc^orderbc)*(orderbc+1.0));
kmax=1.0;
kfactor=(kmax-1.0)/delta/(orderbc+1.0)/delbc^orderbc;

for j=1:upml
    
    % Coefficients for field components in the center of the grid cell
    y1=(upml-j+1)*delta;
    y2=(upml-j)*delta;
    sigma=sigfactor*(y1^(orderbc+1)-y2^(orderbc+1));
    ki=1+kfactor*(y1^(orderbc+1)-y2^(orderbc+1));
    facm=(2*epsr*epsz*ki-sigma*dt);
    facp=(2*epsr*epsz*ki+sigma*dt);
    
    C5ey(:,j,:)=facp;
    C5ey(:,je_tot-j+1,:)=facp;
    C6ey(:,j,:)=facm;
    C6ey(:,je_tot-j+1,:)=facm;
    D1hx(:,j,:)=facm/facp;
    D1hx(:,je_tot-j+1,:)=facm/facp;
    D2hx(:,j,:)=2*epsr*epsz*dt/facp;
    D2hx(:,je_tot-j+1,:)=2*epsr*epsz*dt/facp;
    D3hz(:,j,:)=facm/facp;
    D3hz(:,je_tot-j+1,:)=facm/facp;
    D4hz(:,j,:)=1/facp/mur/muz;
    D4hz(:,je_tot-j+1,:)=1/facp/mur/muz;
    
    % Coefficients for field components on the grid cell boundary
    y1=(upml-j+1.5)*delta;
    y2=(upml-j+0.5)*delta;
    sigma=sigfactor*(y1^(orderbc+1)-y2^(orderbc+1));
    ki=1+kfactor*(y1^(orderbc+1)-y2^(orderbc+1));
    facm=(2*epsr*epsz*ki-sigma*dt);
    facp=(2*epsr*epsz*ki+sigma*dt);    
     
    C1ex(:,j,:)=facm/facp;
    C1ex(:,jh_tot-j+1,:)=facm/facp;
    C2ex(:,j,:)=2*epsr*epsz*dt/facp;
    C2ex(:,jh_tot-j+1,:)=2*epsr*epsz*dt/facp;
    C3ez(:,j,:)=facm/facp;
    C3ez(:,jh_tot-j+1,:)=facm/facp;
    C4ez(:,j,:)=1/facp/epsr/epsz;
    C4ez(:,jh_tot-j+1,:)=1/facp/epsr/epsz;   
    D5hy(:,j,:)=facp;
    D5hy(:,jh_tot-j+1,:)=facp;
    D6hy(:,j,:)=facm;
    D6hy(:,jh_tot-j+1,:)=facm;

end

%   PEC walls
C1ex(:,1,:)=-1;
C1ex(:,jh_tot,:)=-1;
C2ex(:,1,:)=0;
C2ex(:,jh_tot,:)=0;
C3ez(:,1,:)=-1;
C3ez(:,jh_tot,:)=-1;
C4ez(:,1,:)=0;
C4ez(:,jh_tot,:)=0;   

%   z-varying material properties
delbc=upml*delta;
sigmam=-log(rmax)*epsr*epsz*cc*(orderbc+1)/(2*delbc); 
sigfactor=sigmam/(delta*(delbc^orderbc)*(orderbc+1));
kmax=1;
kfactor=(kmax-1)/delta/(orderbc+1)/delbc^orderbc;

for k=1:upml

    % Coefficients for field components in the center of the grid cell
    z1=(upml-k+1)*delta;
    z2=(upml-k)*delta;
    sigma=sigfactor*(z1^(orderbc+1)-z2^(orderbc+1));
    ki=1+kfactor*(z1^(orderbc+1)-z2^(orderbc+1));
    facm=(2*epsr*epsz*ki-sigma*dt);
    facp=(2*epsr*epsz*ki+sigma*dt);
    
    C5ez(:,:,k)=facp;
    C5ez(:,:,ke_tot-k+1)=facp;
    C6ez(:,:,k)=facm;
    C6ez(:,:,ke_tot-k+1)=facm;
    D1hy(:,:,k)=facm/facp;
    D1hy(:,:,ke_tot-k+1)=facm/facp;
    D2hy(:,:,k)=2*epsr*epsz*dt/facp;
    D2hy(:,:,ke_tot-k+1)=2*epsr*epsz*dt/facp;
    D3hx(:,:,k)=facm/facp;
    D3hx(:,:,ke_tot-k+1)=facm/facp;
    D4hx(:,:,k)=1/facp/mur/muz;
    D4hx(:,:,ke_tot-k+1)=1/facp/mur/muz;
    
    % Coefficients for field components on the grid cell boundary
    z1=(upml-k+1.5)*delta;
    z2=(upml-k+0.5)*delta;
    sigma=sigfactor*(z1^(orderbc+1)-z2^(orderbc+1));
    ki=1+kfactor*(z1^(orderbc+1)-z2^(orderbc+1));
    facm=(2*epsr*epsz*ki-sigma*dt);
    facp=(2*epsr*epsz*ki+sigma*dt);
    
    C1ey(:,:,k)=facm/facp;
    C1ey(:,:,kh_tot-k+1)=facm/facp;
    C2ey(:,:,k)=2*epsr*epsz*dt/facp;
    C2ey(:,:,kh_tot-k+1)=2*epsr*epsz*dt/facp;
    C3ex(:,:,k)=facm/facp;
    C3ex(:,:,kh_tot-k+1)=facm/facp;
    C4ex(:,:,k)=1/facp/epsr/epsz;
    C4ex(:,:,kh_tot-k+1)=1/facp/epsr/epsz;
    D5hz(:,:,k)=facp;
    D5hz(:,:,kh_tot-k+1)=facp;
    D6hz(:,:,k)=facm;
    D6hz(:,:,kh_tot-k+1)=facm;

end

%   PEC walls
C1ey(:,:,1)=-1;
C1ey(:,:,kh_tot)=-1;
C2ey(:,:,1)=0;
C2ey(:,:,kh_tot)=0;
C3ex(:,:,1)=-1;
C3ex(:,:,kh_tot)=-1;
C4ex(:,:,1)=0;
C4ex(:,:,kh_tot)=0;

%figure
%set(gcf,'DoubleBuffer','on')

%***********************************************************************
%     Begin time stepping loop
%***********************************************************************

for n=1:nmax
    
    % Update magnetic field
    bstore=bx;
    bx(2:ie_tot,:,:)=D1hx(2:ie_tot,:,:).*  bx(2:ie_tot,:,:)-...
                     D2hx(2:ie_tot,:,:).*((ez(2:ie_tot,2:jh_tot,:)-ez(2:ie_tot,1:je_tot,:))-...
                                          (ey(2:ie_tot,:,2:kh_tot)-ey(2:ie_tot,:,1:ke_tot)))./delta;
    hx(2:ie_tot,:,:)= D3hx(2:ie_tot,:,:).*hx(2:ie_tot,:,:)+...
                      D4hx(2:ie_tot,:,:).*(D5hx(2:ie_tot,:,:).*bx(2:ie_tot,:,:)-...
                                           D6hx(2:ie_tot,:,:).*bstore(2:ie_tot,:,:));
    bstore=by;
    by(:,2:je_tot,:)=D1hy(:,2:je_tot,:).*  by(:,2:je_tot,:)-...
                     D2hy(:,2:je_tot,:).*((ex(:,2:je_tot,2:kh_tot)-ex(:,2:je_tot,1:ke_tot))-...
                                          (ez(2:ih_tot,2:je_tot,:)-ez(1:ie_tot,2:je_tot,:)))./delta;
    hy(:,2:je_tot,:)= D3hy(:,2:je_tot,:).*hy(:,2:je_tot,:)+...
                      D4hy(:,2:je_tot,:).*(D5hy(:,2:je_tot,:).*by(:,2:je_tot,:)-...
                                           D6hy(:,2:je_tot,:).*bstore(:,2:je_tot,:));
    bstore=bz;
    bz(:,:,2:ke_tot)=D1hz(:,:,2:ke_tot).*  bz(:,:,2:ke_tot)-...
                     D2hz(:,:,2:ke_tot).*((ey(2:ih_tot,:,2:ke_tot)-ey(1:ie_tot,:,2:ke_tot))-...
                                          (ex(:,2:jh_tot,2:ke_tot)-ex(:,1:je_tot,2:ke_tot)))./delta;
    hz(:,:,2:ke_tot)= D3hz(:,:,2:ke_tot).*hz(:,:,2:ke_tot)+...
                      D4hz(:,:,2:ke_tot).*(D5hz(:,:,2:ke_tot).*bz(:,:,2:ke_tot)-...
                                           D6hz(:,:,2:ke_tot).*bstore(:,:,2:ke_tot));
    
    % Update electric field
    dstore=dx;
    dx(:,2:je_tot,2:ke_tot)=C1ex(:,2:je_tot,2:ke_tot).*  dx(:,2:je_tot,2:ke_tot)+...
                            C2ex(:,2:je_tot,2:ke_tot).*((hz(:,2:je_tot,2:ke_tot)-hz(:,1:je_tot-1,2:ke_tot))-...
                                                        (hy(:,2:je_tot,2:ke_tot)-hy(:,2:je_tot,1:ke_tot-1)))./delta;
    ex(:,2:je_tot,2:ke_tot)=C3ex(:,2:je_tot,2:ke_tot).*ex(:,2:je_tot,2:ke_tot)+...
                            C4ex(:,2:je_tot,2:ke_tot).*(C5ex(:,2:je_tot,2:ke_tot).*dx(:,2:je_tot,2:ke_tot)-...
                                                        C6ex(:,2:je_tot,2:ke_tot).*dstore(:,2:je_tot,2:ke_tot));
    dstore=dy;
    dy(2:ie_tot,:,2:ke_tot)=C1ey(2:ie_tot,:,2:ke_tot).*  dy(2:ie_tot,:,2:ke_tot)+...
                            C2ey(2:ie_tot,:,2:ke_tot).*((hx(2:ie_tot,:,2:ke_tot)-hx(2:ie_tot,:,1:ke_tot-1))-...
                                                        (hz(2:ie_tot,:,2:ke_tot)-hz(1:ie_tot-1,:,2:ke_tot)))./delta;
    ey(2:ie_tot,:,2:ke_tot)=C3ey(2:ie_tot,:,2:ke_tot).*ey(2:ie_tot,:,2:ke_tot)+...
                            C4ey(2:ie_tot,:,2:ke_tot).*(C5ey(2:ie_tot,:,2:ke_tot).*dy(2:ie_tot,:,2:ke_tot)-...
                                                        C6ey(2:ie_tot,:,2:ke_tot).*dstore(2:ie_tot,:,2:ke_tot));
    dstore=dz;
    dz(2:ie_tot,2:je_tot,:)=C1ez(2:ie_tot,2:je_tot,:).*  dz(2:ie_tot,2:je_tot,:)+...
                            C2ez(2:ie_tot,2:je_tot,:).*((hy(2:ie_tot,2:je_tot,:)-hy(1:ie_tot-1,2:je_tot,:))-...
                                                        (hx(2:ie_tot,2:je_tot,:)-hx(2:ie_tot,1:je_tot-1,:)))./delta;
    dz(is,js,ks:ks+1)=dz(is,js,ks:ks+1)+J0*(n-ndelay)*exp(-((n-ndelay)^2/tau^2));
    ez(2:ie_tot,2:je_tot,:)=C3ez(2:ie_tot,2:je_tot,:).*ez(2:ie_tot,2:je_tot,:)+...
                            C4ez(2:ie_tot,2:je_tot,:).*(C5ez(2:ie_tot,2:je_tot,:).*dz(2:ie_tot,2:je_tot,:)-...
                                                        C6ez(2:ie_tot,2:je_tot,:).*dstore(2:ie_tot,2:je_tot,:));

    %***********************************************************************
    %     Visualize fields
    %***********************************************************************

    timestep=int2str(n);
    tview(:,:)=squeeze(ez(ih_bc:upml+ie,jh_bc:upml+je,ks));
    sview(:,:)=squeeze(ez(ih_bc:upml+ie,js,kh_bc:upml+ke));
    
    subplot('position',[0.15 0.57 0.7 0.35])
    imagesc(tview');
    caxis([-0.2 0.2]);
    colorbar;
    axis image; axis xy;
    title(['E_z(i,j,k=k_s_o_u_r_c_e), time step = ',timestep]);
    xlabel('i coordinate');
    ylabel('j coordinate');
    
    subplot('position',[0.15 0.08 0.7 0.35])
    imagesc(sview');
    caxis([-0.2 0.2]);
    colorbar;
    axis image; axis xy;
    title(['E_z(i,j=j_s_o_u_r_c_e,k), time step = ',timestep]);
    xlabel('i coordinate');
    ylabel('k coordinate');
    
    pause(0.05)
    
end
