%***********************************************************************
%     2-D FDTD TE code with PML absorbing boundary conditions
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
%     This MATLAB M-file implements the finite-difference time-domain
%     solution of Maxwell's curl equations over a two-dimensional
%     Cartesian space lattice comprised of uniform square grid cells.
%
%     To illustrate the algorithm, a 6-cm-diameter metal cylindrical 
%     scatterer in free space is modeled. The source excitation is 
%     a Gaussian pulse with a carrier frequency of 5 GHz.
%
%     The grid resolution (dx = 3 mm) was chosen to provide 20 samples
%     per wavelength at the center frequency of the pulse (which in turn
%     provides approximately 10 samples per wavelength at the high end
%     of the excitation spectrum, around 10 GHz).
%
%     The computational domain is truncated using the perfectly matched
%     layer (PML) absorbing boundary conditions.  The formulation used 
%     in this code is based on the original split-field Berenger PML. 
%     Exponential time stepping is implemented in the PML regions. 
%     The PML regions are labeled as shown in the following diagram: 
%
%            ----------------------------------------------
%           |  |                BACK PML                |  |
%            ----------------------------------------------
%           |L |                                       /| R|
%           |E |                                (ib,jb) | I|
%           |F |                                        | G|
%           |T |                                        | H|
%           |  |                MAIN GRID               | T|
%           |P |                                        |  |
%           |M |                                        | P|
%           |L | (1,1)                                  | M|
%           |  |/                                       | L|
%            ----------------------------------------------
%           |  |                FRONT PML               |  |
%            ----------------------------------------------
%
%     To execute this M-file, type "fdtd2D" at the MATLAB prompt.  
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
%       --uncomment line numbers 418 and 419
%       --comment out line numbers 610, 619, and 628
%
%***********************************************************************

clear

%***********************************************************************
%     Fundamental constants
%***********************************************************************

cc=2.99792458e8;            %speed of light in free space
muz=4.0*pi*1.0e-7;          %permeability of free space
epsz=1.0/(cc*cc*muz);       %permittivity of free space
etaz=sqrt(muz/epsz);

freq=5.0e+9;                %center frequency of source excitation
lambda=cc/freq;             %center wavelength of source excitation
omega=2.0*pi*freq;          

%***********************************************************************
%     Grid parameters
%***********************************************************************

ie=100;           %number of grid cells in x-direction
je=50;            %number of grid cells in y-direction

ib=ie+1;
jb=je+1;

is=15;            %location of z-directed hard source
js=je/2;          %location of z-directed hard source

dx=3.0e-3;        %space increment of square lattice
dt=dx/(2.0*cc);   %time step

nmax=300;         %total number of time steps

iebc=8;           %thickness of left and right PML region
jebc=8;           %thickness of front and back PML region
rmax=1.0e-7;
orderbc=2;
ibbc=iebc+1;
jbbc=jebc+1;
iefbc=ie+2*iebc;
jefbc=je+2*jebc;
ibfbc=iefbc+1;
jbfbc=jefbc+1;

%***********************************************************************
%     Material parameters
%***********************************************************************

media=2;

eps=[1.0 1.0];
sig=[0.0 1.0e+7];
mur=[1.0 1.0];
sim=[0.0 0.0];

%***********************************************************************
%     Wave excitation
%***********************************************************************

rtau=160.0e-12;
tau=rtau/dt;
delay=3*tau;

source=zeros(1,nmax);
for n=1:7.0*tau
  source(n)=sin(omega*(n-delay)*dt)*exp(-((n-delay)^2/tau^2));
end

%***********************************************************************
%     Field arrays
%***********************************************************************

ex=zeros(ie,jb);           %fields in main grid 
ey=zeros(ib,je);
hz=zeros(ie,je);

exbcf=zeros(iefbc,jebc);   %fields in front PML region
eybcf=zeros(ibfbc,jebc);
hzxbcf=zeros(iefbc,jebc);
hzybcf=zeros(iefbc,jebc);

exbcb=zeros(iefbc,jbbc);   %fields in back PML region
eybcb=zeros(ibfbc,jebc);
hzxbcb=zeros(iefbc,jebc);
hzybcb=zeros(iefbc,jebc);

exbcl=zeros(iebc,jb);      %fields in left PML region
eybcl=zeros(iebc,je);
hzxbcl=zeros(iebc,je);
hzybcl=zeros(iebc,je);

exbcr=zeros(iebc,jb);      %fields in right PML region
eybcr=zeros(ibbc,je);
hzxbcr=zeros(iebc,je);
hzybcr=zeros(iebc,je);

%***********************************************************************
%     Updating coefficients
%***********************************************************************

for i=1:media
  eaf  =dt*sig(i)/(2.0*epsz*eps(i));
  ca(i)=(1.0-eaf)/(1.0+eaf);
  cb(i)=dt/epsz/eps(i)/dx/(1.0+eaf);
  haf  =dt*sim(i)/(2.0*muz*mur(i));
  da(i)=(1.0-haf)/(1.0+haf);
  db(i)=dt/muz/mur(i)/dx/(1.0+haf);
end

%***********************************************************************
%     Geometry specification (main grid)
%***********************************************************************

%     Initialize entire main grid to free space

caex(1:ie,1:jb)=ca(1);     
cbex(1:ie,1:jb)=cb(1);

caey(1:ib,1:je)=ca(1);
cbey(1:ib,1:je)=cb(1);

dahz(1:ie,1:je)=da(1);
dbhz(1:ie,1:je)=db(1);

%     Add metal cylinder

diam=20;          % diameter of cylinder: 6 cm
rad=diam/2.0;     % radius of cylinder: 3 cm
icenter=4*ie/5;   % i-coordinate of cylinder's center
jcenter=je/2;     % j-coordinate of cylinder's center

for i=1:ie
for j=1:je
  dist2=(i+0.5-icenter)^2 + (j-jcenter)^2;
  if dist2 <= rad^2 
     caex(i,j)=ca(2);
     cbex(i,j)=cb(2);
  end
  dist2=(i-icenter)^2 + (j+0.5-jcenter)^2;
  if dist2 <= rad^2 
     caey(i,j)=ca(2);
     cbey(i,j)=cb(2);
  end
end
end

%***********************************************************************
%     Fill the PML regions
%
%     PML theory describes a continuous grading of the material properties
%     over the PML region.  In the FDTD grid it is necessary to discretize
%     the grading by averaging the material properties over a grid cell 
%     width centered on each field component.  As an example of the 
%     implementation of this averaging, we take the integral of the 
%     continuous sigma(x) in the PML region
%   
%         sigma_i = integral(sigma(x))/dx
%   
%     where the integral is over a single grid cell width in x, and is 
%     bounded by x1 and x2.  Applying this to the polynomial grading of 
%     Equation 7.60a produces
%
%         sigma_i = sigmam/(dx*(m+1)*d^m)*(x2^(m+1)-x1^(m+1))
%
%     where sigmam is the maximum value of sigma as described by Equation 
%     7.62.  In the code below, the coefficient in the expression for
%     sigma_i is denoted "bcfactor".
%         
%     The definitions of x1 and x2 depend on the position of the field 
%     component within the grid cell.  We have either
%
%         x1 = (i-0.5)*dx,  x2 = (i+0.5)*dx
%  
%     or
%  
%         x1 = (i)*dx,      x2 = (i+1)*dx
%
%     where i varies over the PML region.
%  
%***********************************************************************

delbc=iebc*dx;
sigmam=-log(rmax)*(orderbc+1)/(2*etaz*delbc);   %Eq. 7.62 in text
bcfactor=sigmam/(dx*(orderbc+1)*(delbc^orderbc));  

%     FRONT region 

caexbcf(1:iefbc,1)=1.0;
cbexbcf(1:iefbc,1)=0.0;
for j=2:jebc
  y1=(jebc-j+1.5)*dx;
  y2=(jebc-j+0.5)*dx;
  sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
  ca1=exp(-sigmay*dt/epsz);
  cb1=(1.0-ca1)/(sigmay*dx);
  caexbcf(1:iefbc,j)=ca1;
  cbexbcf(1:iefbc,j)=cb1;
end
sigmay = bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmay*dt/epsz);
cb1=(1-ca1)/(sigmay*dx);
caex(1:ie,1)=ca1;
cbex(1:ie,1)=cb1;
caexbcl(1:iebc,1)=ca1;
cbexbcl(1:iebc,1)=cb1;
caexbcr(1:iebc,1)=ca1;
cbexbcr(1:iebc,1)=cb1;

for j=1:jebc
  y1=(jebc-j+1)*dx;
  y2=(jebc-j)*dx;
  sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
  sigmays=sigmay*(muz/epsz);
  da1=exp(-sigmays*dt/muz);
  db1=(1-da1)/(sigmays*dx);
  dahzybcf(1:iefbc,j)=da1;
  dbhzybcf(1:iefbc,j)=db1;
  caeybcf(1:ibfbc,j)=ca(1);
  cbeybcf(1:ibfbc,j)=cb(1);
  dahzxbcf(1:iefbc,j)=da(1);
  dbhzxbcf(1:iefbc,j)=db(1);
end

%     BACK region 

caexbcb(1:iefbc,jbbc)=1.0;
cbexbcb(1:iefbc,jbbc)=0.0;
for j=2:jebc
  y1=(j-0.5)*dx;
  y2=(j-1.5)*dx;
  sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
  ca1=exp(-sigmay*dt/epsz);
  cb1=(1-ca1)/(sigmay*dx);
  caexbcb(1:iefbc,j)=ca1;
  cbexbcb(1:iefbc,j)=cb1;
end
sigmay = bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmay*dt/epsz);
cb1=(1-ca1)/(sigmay*dx);
caex(1:ie,jb)=ca1;
cbex(1:ie,jb)=cb1;
caexbcl(1:iebc,jb)=ca1;
cbexbcl(1:iebc,jb)=cb1;
caexbcr(1:iebc,jb)=ca1;
cbexbcr(1:iebc,jb)=cb1;

for j=1:jebc
  y1=j*dx;
  y2=(j-1)*dx;
  sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
  sigmays=sigmay*(muz/epsz);
  da1=exp(-sigmays*dt/muz);
  db1=(1-da1)/(sigmays*dx);
  dahzybcb(1:iefbc,j)=da1;
  dbhzybcb(1:iefbc,j)=db1;
  caeybcb(1:ibfbc,j)=ca(1);
  cbeybcb(1:ibfbc,j)=cb(1);
  dahzxbcb(1:iefbc,j)=da(1);
  dbhzxbcb(1:iefbc,j)=db(1);
end

%     LEFT region 

caeybcl(1,1:je)=1.0;
cbeybcl(1,1:je)=0.0;
for i=2:iebc
  x1=(iebc-i+1.5)*dx;
  x2=(iebc-i+0.5)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  ca1=exp(-sigmax*dt/epsz);
  cb1=(1-ca1)/(sigmax*dx);
  caeybcl(i,1:je)=ca1;
  cbeybcl(i,1:je)=cb1;
  caeybcf(i,1:jebc)=ca1;
  cbeybcf(i,1:jebc)=cb1;
  caeybcb(i,1:jebc)=ca1;
  cbeybcb(i,1:jebc)=cb1;
end
sigmax=bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmax*dt/epsz);
cb1=(1-ca1)/(sigmax*dx);
caey(1,1:je)=ca1;
cbey(1,1:je)=cb1;
caeybcf(iebc+1,1:jebc)=ca1;
cbeybcf(iebc+1,1:jebc)=cb1;
caeybcb(iebc+1,1:jebc)=ca1;
cbeybcb(iebc+1,1:jebc)=cb1;

for i=1:iebc
  x1=(iebc-i+1)*dx;
  x2=(iebc-i)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  sigmaxs=sigmax*(muz/epsz);
  da1=exp(-sigmaxs*dt/muz);
  db1=(1-da1)/(sigmaxs*dx);
  dahzxbcl(i,1:je)=da1;
  dbhzxbcl(i,1:je)=db1;
  dahzxbcf(i,1:jebc)=da1;
  dbhzxbcf(i,1:jebc)=db1;
  dahzxbcb(i,1:jebc)=da1;
  dbhzxbcb(i,1:jebc)=db1;
  caexbcl(i,2:je)=ca(1);
  cbexbcl(i,2:je)=cb(1);
  dahzybcl(i,1:je)=da(1);
  dbhzybcl(i,1:je)=db(1);
end

%     RIGHT region 

caeybcr(ibbc,1:je)=1.0;
cbeybcr(ibbc,1:je)=0.0;
for i=2:iebc
  x1=(i-0.5)*dx;
  x2=(i-1.5)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  ca1=exp(-sigmax*dt/epsz);
  cb1=(1-ca1)/(sigmax*dx);
  caeybcr(i,1:je)=ca1;
  cbeybcr(i,1:je)=cb1;
  caeybcf(i+iebc+ie,1:jebc)=ca1;
  cbeybcf(i+iebc+ie,1:jebc)=cb1;
  caeybcb(i+iebc+ie,1:jebc)=ca1;
  cbeybcb(i+iebc+ie,1:jebc)=cb1;
end
sigmax=bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmax*dt/epsz);
cb1=(1-ca1)/(sigmax*dx);
caey(ib,1:je)=ca1;
cbey(ib,1:je)=cb1;
caeybcf(iebc+ib,1:jebc)=ca1;
cbeybcf(iebc+ib,1:jebc)=cb1;
caeybcb(iebc+ib,1:jebc)=ca1;
cbeybcb(iebc+ib,1:jebc)=cb1;

for i=1:iebc
  x1=i*dx;
  x2=(i-1)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  sigmaxs=sigmax*(muz/epsz);
  da1=exp(-sigmaxs*dt/muz);
  db1=(1-da1)/(sigmaxs*dx);
  dahzxbcr(i,1:je) = da1;
  dbhzxbcr(i,1:je) = db1;
  dahzxbcf(i+ie+iebc,1:jebc)=da1;
  dbhzxbcf(i+ie+iebc,1:jebc)=db1;
  dahzxbcb(i+ie+iebc,1:jebc)=da1;
  dbhzxbcb(i+ie+iebc,1:jebc)=db1;
  caexbcr(i,2:je)=ca(1);
  cbexbcr(i,2:je)=cb(1);
  dahzybcr(i,1:je)=da(1);
  dbhzybcr(i,1:je)=db(1);
end

%figure
%set(gcf,'DoubleBuffer','on')

%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************

for n=1:nmax

%***********************************************************************
%     Update electric fields (EX and EY) in main grid
%***********************************************************************

ex(:,2:je)=caex(:,2:je).*ex(:,2:je)+...
           cbex(:,2:je).*(hz(:,2:je)-hz(:,1:je-1));

ey(2:ie,:)=caey(2:ie,:).*ey(2:ie,:)+...
           cbey(2:ie,:).*(hz(1:ie-1,:)-hz(2:ie,:));

%***********************************************************************
%     Update EX in PML regions
%***********************************************************************

%     FRONT

exbcf(:,2:jebc)=caexbcf(:,2:jebc).*exbcf(:,2:jebc)-...  
  cbexbcf(:,2:jebc).*(hzxbcf(:,1:jebc-1)+hzybcf(:,1:jebc-1)-...
                      hzxbcf(:,2:jebc)-hzybcf(:,2:jebc));
ex(1:ie,1)=caex(1:ie,1).*ex(1:ie,1)-...
  cbex(1:ie,1).*(hzxbcf(ibbc:iebc+ie,jebc)+...
                hzybcf(ibbc:iebc+ie,jebc)-hz(1:ie,1));
 
%     BACK

exbcb(:,2:jebc-1)=caexbcb(:,2:jebc-1).*exbcb(:,2:jebc-1)-...
  cbexbcb(:,2:jebc-1).*(hzxbcb(:,1:jebc-2)+hzybcb(:,1:jebc-2)-...
                        hzxbcb(:,2:jebc-1)-hzybcb(:,2:jebc-1));
ex(1:ie,jb)=caex(1:ie,jb).*ex(1:ie,jb)-...
  cbex(1:ie,jb).*(hz(1:ie,jb-1)-hzxbcb(ibbc:iebc+ie,1)-...
                 hzybcb(ibbc:iebc+ie,1));
 
%     LEFT

exbcl(:,2:je)=caexbcl(:,2:je).*exbcl(:,2:je)-...
  cbexbcl(:,2:je).*(hzxbcl(:,1:je-1)+hzybcl(:,1:je-1)-...
                    hzxbcl(:,2:je)-hzybcl(:,2:je));
exbcl(:,1)=caexbcl(:,1).*exbcl(:,1)-...
  cbexbcl(:,1).*(hzxbcf(1:iebc,jebc)+hzybcf(1:iebc,jebc)-...
                 hzxbcl(:,1)-hzybcl(:,1));
exbcl(:,jb)=caexbcl(:,jb).*exbcl(:,jb)-...
  cbexbcl(:,jb).*(hzxbcl(:,je)+hzybcl(:,je)-...
                  hzxbcb(1:iebc,1)-hzybcb(1:iebc,1));
 
%     RIGHT

exbcr(:,2:je)=caexbcr(:,2:je).*exbcr(:,2:je)-...
  cbexbcr(:,2:je).*(hzxbcr(:,1:je-1)+hzybcr(:,1:je-1)-...
                    hzxbcr(:,2:je)-hzybcr(:,2:je));
exbcr(:,1)=caexbcr(:,1).*exbcr(:,1)-...
  cbexbcr(:,1).*(hzxbcf(1+iebc+ie:iefbc,jebc)+...
                 hzybcf(1+iebc+ie:iefbc,jebc)-...
                 hzxbcr(:,1)-hzybcr(:,1));
exbcr(:,jb)=caexbcr(:,jb).*exbcr(:,jb)-...
  cbexbcr(:,jb).*(hzxbcr(:,je)+hzybcr(:,je)-...
                  hzxbcb(1+iebc+ie:iefbc,1)-...
                  hzybcb(1+iebc+ie:iefbc,1));
 
%***********************************************************************
%     Update EY in PML regions
%***********************************************************************

%     FRONT

eybcf(2:iefbc,:)=caeybcf(2:iefbc,:).*eybcf(2:iefbc,:)-...
  cbeybcf(2:iefbc,:).*(hzxbcf(2:iefbc,:)+hzybcf(2:iefbc,:)-...
                       hzxbcf(1:iefbc-1,:)-hzybcf(1:iefbc-1,:));
 
%     BACK

eybcb(2:iefbc,:)=caeybcb(2:iefbc,:).*eybcb(2:iefbc,:)-...
  cbeybcb(2:iefbc,:).*(hzxbcb(2:iefbc,:)+hzybcb(2:iefbc,:)-...
                       hzxbcb(1:iefbc-1,:)-hzybcb(1:iefbc-1,:));
 
%     LEFT

eybcl(2:iebc,:)=caeybcl(2:iebc,:).*eybcl(2:iebc,:)-...
  cbeybcl(2:iebc,:).*(hzxbcl(2:iebc,:)+hzybcl(2:iebc,:)-...
                      hzxbcl(1:iebc-1,:)-hzybcl(1:iebc-1,:));
ey(1,:)=caey(1,:).*ey(1,:)-...
  cbey(1,:).*(hz(1,:)-hzxbcl(iebc,:)-hzybcl(iebc,:));
 
%     RIGHT

eybcr(2:iebc,:)=caeybcr(2:iebc,:).*eybcr(2:iebc,:)-...
  cbeybcr(2:iebc,:).*(hzxbcr(2:iebc,:)+hzybcr(2:iebc,:)-...
                      hzxbcr(1:iebc-1,:)-hzybcr(1:iebc-1,:));
ey(ib,:)=caey(ib,:).*ey(ib,:)-...
  cbey(ib,:).*(hzxbcr(1,:)+hzybcr(1,:)- hz(ie,:));


%***********************************************************************
%     Update magnetic fields (HZ) in main grid
%***********************************************************************

hz(1:ie,1:je)=dahz(1:ie,1:je).*hz(1:ie,1:je)+... 
              dbhz(1:ie,1:je).*(ex(1:ie,2:jb)-ex(1:ie,1:je)+...
                                ey(1:ie,1:je)-ey(2:ib,1:je));

hz(is,js)=source(n);


%***********************************************************************
%     Update HZX in PML regions
%***********************************************************************

%     FRONT

hzxbcf(1:iefbc,:)=dahzxbcf(1:iefbc,:).*hzxbcf(1:iefbc,:)-...
  dbhzxbcf(1:iefbc,:).*(eybcf(2:ibfbc,:)-eybcf(1:iefbc,:));
 
%     BACK
 
hzxbcb(1:iefbc,:)=dahzxbcb(1:iefbc,:).*hzxbcb(1:iefbc,:)-...
  dbhzxbcb(1:iefbc,:).*(eybcb(2:ibfbc,:)-eybcb(1:iefbc,:));
 
%     LEFT
 
hzxbcl(1:iebc-1,:)=dahzxbcl(1:iebc-1,:).*hzxbcl(1:iebc-1,:)-...
  dbhzxbcl(1:iebc-1,:).*(eybcl(2:iebc,:)-eybcl(1:iebc-1,:));
hzxbcl(iebc,:)=dahzxbcl(iebc,:).*hzxbcl(iebc,:)-...
  dbhzxbcl(iebc,:).*(ey(1,:)-eybcl(iebc,:));
 
%     RIGHT
 
hzxbcr(2:iebc,:)=dahzxbcr(2:iebc,:).*hzxbcr(2:iebc,:)-...
  dbhzxbcr(2:iebc,:).*(eybcr(3:ibbc,:)-eybcr(2:iebc,:));
hzxbcr(1,:)=dahzxbcr(1,:).*hzxbcr(1,:)-...
  dbhzxbcr(1,:).*(eybcr(2,:)-ey(ib,:));
 
%***********************************************************************
%     Update HZY in PML regions
%***********************************************************************

%     FRONT
 
hzybcf(:,1:jebc-1)=dahzybcf(:,1:jebc-1).*hzybcf(:,1:jebc-1)-...
  dbhzybcf(:,1:jebc-1).*(exbcf(:,1:jebc-1)-exbcf(:,2:jebc));
hzybcf(1:iebc,jebc)=dahzybcf(1:iebc,jebc).*hzybcf(1:iebc,jebc)-...
  dbhzybcf(1:iebc,jebc).*(exbcf(1:iebc,jebc)-exbcl(1:iebc,1));
hzybcf(iebc+1:iebc+ie,jebc)=...
  dahzybcf(iebc+1:iebc+ie,jebc).*hzybcf(iebc+1:iebc+ie,jebc)-...
  dbhzybcf(iebc+1:iebc+ie,jebc).*(exbcf(iebc+1:iebc+ie,jebc)-...
                                  ex(1:ie,1));
hzybcf(iebc+ie+1:iefbc,jebc)=...
  dahzybcf(iebc+ie+1:iefbc,jebc).*hzybcf(iebc+ie+1:iefbc,jebc)-...
  dbhzybcf(iebc+ie+1:iefbc,jebc).*(exbcf(iebc+ie+1:iefbc,jebc)-...
                                   exbcr(1:iebc,1));

%     BACK
 
hzybcb(1:iefbc,2:jebc)=dahzybcb(1:iefbc,2:jebc).*hzybcb(1:iefbc,2:jebc)-...
  dbhzybcb(1:iefbc,2:jebc).*(exbcb(1:iefbc,2:jebc)-exbcb(1:iefbc,3:jbbc));
hzybcb(1:iebc,1)=dahzybcb(1:iebc,1).*hzybcb(1:iebc,1)-...
  dbhzybcb(1:iebc,1).*(exbcl(1:iebc,jb)-exbcb(1:iebc,2));
hzybcb(iebc+1:iebc+ie,1)=...
  dahzybcb(iebc+1:iebc+ie,1).*hzybcb(iebc+1:iebc+ie,1)-...
  dbhzybcb(iebc+1:iebc+ie,1).*(ex(1:ie,jb)-exbcb(iebc+1:iebc+ie,2));
hzybcb(iebc+ie+1:iefbc,1)=...
  dahzybcb(iebc+ie+1:iefbc,1).*hzybcb(iebc+ie+1:iefbc,1)-...
  dbhzybcb(iebc+ie+1:iefbc,1).*(exbcr(1:iebc,jb)-...
                                exbcb(iebc+ie+1:iefbc,2));
 
%     LEFT
 
hzybcl(:,1:je)=dahzybcl(:,1:je).*hzybcl(:,1:je)-...
  dbhzybcl(:,1:je).*(exbcl(:,1:je)-exbcl(:,2:jb));
 
%     RIGHT
 
hzybcr(:,1:je)=dahzybcr(:,1:je).*hzybcr(:,1:je)-...
  dbhzybcr(:,1:je).*(exbcr(:,1:je)-exbcr(:,2:jb));

%***********************************************************************
%     Visualize fields
%***********************************************************************

timestep=int2str(n);

subplot(3,1,1),imagesc(ex');
shading flat;
caxis([-80.0 80.0]);
axis([1 ie 1 jb]);
colorbar;
axis image; axis xy
axis off;
title(['Ex at time step = ',timestep]);

subplot(3,1,2),imagesc(ey');
shading flat;
caxis([-80.0 80.0]);
axis([1 ib 1 je]);
colorbar;
axis image; axis xy
axis off;
title(['Ey at time step = ',timestep]);

subplot(3,1,3),imagesc(hz');
shading flat;
caxis([-0.2 0.2]);
axis([1 ie 1 je]);
colorbar;
axis image; axis xy
axis off;
title(['Hz at time step = ',timestep]);
pause(0.01)

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************

end

