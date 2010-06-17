!***********************************************************************
!     3-D FDTD code with CPML absorbing boundary conditions
!***********************************************************************
!
!     Program author: Jamesina J. Simpson, Assistant Professor
!                     National Science Foundation Fellow
!                     Department of Electrical and Computer Engineering 
!                     University of New Mexico
!
!                 Address:
!                     MSC01 1100
!                     1 Universty of New Mexico
!                     Albuquerque, NM 87131
!                 Email:
!                     simpson@ece.unm.edu
!
!     Copyright 2005
!
!     This FORTRAN 90 code implements the finite-difference time-domain
!     solution of Maxwell's curl equations over a three-dimensional
!     Cartesian space lattice.  The grid is terminated by CPML absorbing
!     boundary conditions.  However, in a straightforward manner it can
!     be altered to have PEC planes at any of the outer grid boundaries.  
!     Also, the number of grid cells as well as the thickness of 
!     the PML in each Cartesian direction can be varied independently.
!
!     For illustrative purposes, the code supplied here models the PML 
!     numerical experiment described in Section 7.11.2.  Three output 
!     files are generated having the following data:
!
!       (1) Ez at the source for each time step;
!       (2) Ey at the probe point for each time step;  and
!       (3) The plane of Ez values 1 mm above the source in the k-
!           direction recorded at the time step, "record_grid".  This 
!           data can be viewed in Matlab using the following commands:
!
!             >>  load view_grid.dat
!             >>  image = reshape(view_grid,Imax,Jmax);
!             >>  pcolor(image'); shading flat
!
!     The relative error (equation 7.135) graphed in Fig. 7.6 on 
!     page 318 can be reproduced by comparing the output file 
!     "probe.dat" with that of a much larger benchmark grid having 
!     Imax, Jmax, and Kmax increased to the values mentioned in the text
!     (201 x 276 x 176).
!
!     This code has been tested using the Intel Fortran Compiler 8.0 for
!     Linux.  The executable can be created by typing
!     "ifort fdtd3D_CPML.f90" at the command prompt.
!
!     This program is not guaranteed to be free of defects or bugs.
!     Please report any bugs that may exist to Dr. Simpson at:
!                                                    simpson@ece.unm.edu
!
!
!***********************************************************************


   PROGRAM fdtd3D_CPML

   IMPLICIT NONE

!  ..................................
!  Input Fundamental Constants (MKS units)
   REAL, PARAMETER ::                            &
      pi = 3.14159265358979, C = 2.99792458E8, &
      muO = 4.0 * pi * 1.0E-7, epsO = 1.0/(C*C*muO)

!  ..................................
!  Specify Material Relative Permittivity and Conductivity
   REAL, PARAMETER::                      &
      epsR = 1.0, sigM1 = 0.0   ! free space

!  ..................................
!  Specify Number of Time Steps and Grid Size Parameters
   INTEGER, PARAMETER ::                                     &
      nmax = 2100, &  ! total number of time steps
      Imax = 51, Jmax = 126, Kmax = 26  

!  ..................................
!  Specify Grid Cell Size in Each Direction and Calculate the 
!  Resulting Courant-Stable Time Step
   REAL, PARAMETER ::                                        &
      dx = 1.0E-3, dy = 1.0E-3, dz = 1.0E-3,  &  ! cell size in each direction
      dt = 0.99 / (C*(1.0/dx**2+1.0/dy**2+1.0/dz**2)**0.5)
                                                 ! time step increment

!  ..................................
!  Specify the Impulsive Source (See Equation 7.134)
   REAL, PARAMETER ::                                        &
      tw = 53.0E-12, tO = 4.0*tw  

!  ..................................
!  Specify the Time Step at which the Grid is Recorded for Visualization
   INTEGER, PARAMETER ::                                        &
      record_grid = 300 

!  ..................................
!  Specify the PEC Plate Boundaries and the Source/Recording Points
   INTEGER, PARAMETER ::                                    &
      istart = (Imax-1)/2-11, iend = istart+24, jstart = Jmax/2-49,   &
      jend = jstart + 99, kstart = Kmax/2, kend = kstart,      &
      isource = istart, jsource = jstart, ksource = kstart,     &
      irecv1 = iend, jrecv1 = jend+1, krecv1 = kend  ! Ey at probe point

!  ..................................
!  Specify the CPML Thickness in Each Direction (Value of Zero 
!  Corresponds to No PML, and the Grid is Terminated with a PEC)
   INTEGER, PARAMETER ::                        &
      ! PML thickness in each direction 
      nxPML_1 = 11, nxPML_2 = nxPML_1, nyPML_1 = nxPML_1,      &
      nyPML_2 = nxPML_1, nzPML_1 = nxPML_1, nzPML_2 = nxPML_2
!  ..................................
!  Specify the CPML Order and Other Parameters
   INTEGER, PARAMETER ::                        & 
      m = 3, ma = 1 
   REAL, PARAMETER  ::                     &
      sig_x_max = 0.75 * (0.8*(m+1)/(dx*(muO/epsO*epsR)**0.5)),   &
      sig_y_max = 0.75 * (0.8*(m+1)/(dy*(muO/epsO*epsR)**0.5)),   &
      sig_z_max = 0.75 * (0.8*(m+1)/(dz*(muO/epsO*epsR)**0.5)),   &
      alpha_x_max = 0.24,   &
      alpha_y_max = alpha_x_max, alpha_z_max = alpha_x_max, &
      kappa_x_max = 15.0, &
      kappa_y_max = kappa_x_max, kappa_z_max = kappa_x_max

   INTEGER ::                                                &
	i,j,ii,jj,k,kk,n

   REAL  ::                                                   &
      source, P1, P2

!     TM components
   REAL,DIMENSION(Imax, Jmax, Kmax-1)  ::                      &
      Ez, CA, CB, sig, eps

   REAL,DIMENSION(Imax-1, Jmax, Kmax-1)  ::                      &
      Hy

   REAL,DIMENSION(Imax,Jmax-1, Kmax-1)  ::                      &
      Hx

   REAL  ::                        &
      DA, DB

!     TE components
   REAL,DIMENSION(Imax-1, Jmax-1, Kmax)  ::                      &
      Hz

   REAL,DIMENSION(Imax-1, Jmax, Kmax)  ::                      &
      Ex

   REAL,DIMENSION(Imax,Jmax-1, Kmax)  ::                      &
      Ey

!  PML
   REAL ,DIMENSION(nxPML_1,Jmax,Kmax) ::                       &
      psi_Ezx_1

   REAL ,DIMENSION(nxPML_2,Jmax,Kmax) ::                       &
      psi_Ezx_2

   REAL ,DIMENSION(nxPML_1-1,Jmax,Kmax) ::                       &
      psi_Hyx_1

   REAL ,DIMENSION(nxPML_2-1,Jmax,Kmax) ::                       &
      psi_Hyx_2

   REAL ,DIMENSION(Imax,nyPML_1,Kmax) ::                       &
      psi_Ezy_1                               

   REAL ,DIMENSION(Imax,nyPML_2,Kmax) ::                       &
      psi_Ezy_2

   REAL ,DIMENSION(Imax,nyPML_1-1,Kmax) ::                       &
      psi_Hxy_1                               

   REAL ,DIMENSION(Imax,nyPML_2-1,Kmax) ::                       &
      psi_Hxy_2

   REAL ,DIMENSION(Imax,Jmax-1,nzPML_1-1) ::                       &
      psi_Hxz_1

   REAL ,DIMENSION(Imax,Jmax-1,nzPML_2-1) ::                       &
      psi_Hxz_2

   REAL ,DIMENSION(Imax-1,Jmax,nzPML_1-1) ::                       &
      psi_Hyz_1

   REAL ,DIMENSION(Imax-1,Jmax,nzPML_2-1) ::                       &
      psi_Hyz_2

   REAL ,DIMENSION(Imax-1,Jmax,nzPML_1) ::                       &
      psi_Exz_1

   REAL ,DIMENSION(Imax-1,Jmax,nzPML_2) ::                       &
      psi_Exz_2

   REAL ,DIMENSION(Imax,Jmax-1,nzPML_1) ::                       &
      psi_Eyz_1

   REAL ,DIMENSION(Imax,Jmax-1,nzPML_2) ::                       &
      psi_Eyz_2

   REAL ,DIMENSION(nxPML_1-1,Jmax-1,Kmax) ::                       &
      psi_Hzx_1
   REAL ,DIMENSION(nxPML_1,Jmax-1,Kmax) ::                       &
      psi_Eyx_1

   REAL ,DIMENSION(nxPML_2-1,Jmax-1,Kmax) ::                       &
      psi_Hzx_2
   REAL ,DIMENSION(nxPML_2,Jmax-1,Kmax) ::                       &
      psi_Eyx_2

   REAL ,DIMENSION(Imax-1,nyPML_1-1,Kmax) ::                       &
      psi_Hzy_1                               
   REAL ,DIMENSION(Imax-1,nyPML_1,Kmax) ::                       &
      psi_Exy_1                               

   REAL ,DIMENSION(Imax-1,nyPML_2-1,Kmax) ::                       &
      psi_Hzy_2
   REAL ,DIMENSION(Imax-1,nyPML_2,Kmax) ::                       &
      psi_Exy_2

   REAL ,DIMENSION(nxPML_1) ::                       &
      be_x_1, ce_x_1, alphae_x_PML_1, sige_x_PML_1, kappae_x_PML_1
   REAL ,DIMENSION(nxPML_1-1) ::                       &
      bh_x_1, ch_x_1, alphah_x_PML_1, sigh_x_PML_1, kappah_x_PML_1

   REAL ,DIMENSION(nxPML_2) ::                       &
      be_x_2, ce_x_2, alphae_x_PML_2, sige_x_PML_2, kappae_x_PML_2
   REAL ,DIMENSION(nxPML_2-1) ::                       &
      bh_x_2, ch_x_2, alphah_x_PML_2, sigh_x_PML_2, kappah_x_PML_2

   REAL ,DIMENSION(nyPML_1) ::                       &
      be_y_1, ce_y_1, alphae_y_PML_1, sige_y_PML_1, kappae_y_PML_1
   REAL ,DIMENSION(nyPML_1-1) ::                       &
      bh_y_1, ch_y_1, alphah_y_PML_1, sigh_y_PML_1, kappah_y_PML_1

   REAL ,DIMENSION(nyPML_2) ::                       &
      be_y_2, ce_y_2, alphae_y_PML_2, sige_y_PML_2, kappae_y_PML_2
   REAL ,DIMENSION(nyPML_2-1) ::                       &
      bh_y_2, ch_y_2, alphah_y_PML_2, sigh_y_PML_2, kappah_y_PML_2

   REAL ,DIMENSION(nzPML_1) ::                       &
      be_z_1, ce_z_1, alphae_z_PML_1, sige_z_PML_1, kappae_z_PML_1
   REAL ,DIMENSION(nzPML_1-1) ::                       &
      bh_z_1, ch_z_1, alphah_z_PML_1, sigh_z_PML_1, kappah_z_PML_1

   REAL ,DIMENSION(nzPML_2) ::                       &
      be_z_2, ce_z_2, alphae_z_PML_2, sige_z_PML_2, kappae_z_PML_2
   REAL ,DIMENSION(nzPML_2-1) ::                       &
      bh_z_2, ch_z_2, alphah_z_PML_2, sigh_z_PML_2, kappah_z_PML_2

!     denominators for the update equations
   REAL,DIMENSION(Imax-1)  ::                      &
      den_ex, den_hx

   REAL,DIMENSION(Jmax-1)  ::                      &
      den_ey, den_hy

   REAL,DIMENSION(Kmax-1)  ::                      &
      den_ez, den_hz


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  OPEN OUTPUT FILES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   OPEN (UNIT = 30, FILE = "source.dat")
   OPEN (UNIT = 31, FILE = "probe.dat")
   OPEN (UNIT = 33, FILE = "view_grid.dat")

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  INITIALIZE VARIABLES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   P1 = 0.0
   P2 = 0.0

   Ez(:,:,:) = 0.0
   Hz(:,:,:) = 0.0
   Ex(:,:,:) = 0.0
   Hx(:,:,:) = 0.0
   Ey(:,:,:) = 0.0
   Hy(:,:,:) = 0.0
   sig(:,:,:) = sigM1
   eps(:,:,:) = epsR*epsO

   psi_Exy_1(:,:,:) = 0.0
   psi_Exy_2(:,:,:) = 0.0
   psi_Exz_1(:,:,:) = 0.0
   psi_Exz_2(:,:,:) = 0.0
   psi_Eyx_1(:,:,:) = 0.0
   psi_Eyx_2(:,:,:) = 0.0
   psi_Eyz_1(:,:,:) = 0.0
   psi_Eyz_2(:,:,:) = 0.0
   psi_Ezy_1(:,:,:) = 0.0
   psi_Ezy_2(:,:,:) = 0.0
   psi_Ezx_1(:,:,:) = 0.0
   psi_Ezx_2(:,:,:) = 0.0
   psi_Hxy_1(:,:,:) = 0.0
   psi_Hxy_2(:,:,:) = 0.0
   psi_Hxz_1(:,:,:) = 0.0
   psi_Hxz_2(:,:,:) = 0.0
   psi_Hyx_1(:,:,:) = 0.0
   psi_Hyx_2(:,:,:) = 0.0
   psi_Hyz_1(:,:,:) = 0.0
   psi_Hyz_2(:,:,:) = 0.0
   psi_Hzy_1(:,:,:) = 0.0
   psi_Hzy_2(:,:,:) = 0.0
   psi_Hzx_1(:,:,:) = 0.0
   psi_Hzx_2(:,:,:) = 0.0

   write(*,*)"Imax: ", Imax
   write(*,*)"Jmax: ", Jmax
   write(*,*)"Kmax: ", Kmax
   write(*,*)"dt: ", dt
   write(*,*)"nmax: ", nmax
   write(*,*)"max time: ", nmax*dt
   write(*,*)"record grid after ", record_grid, "dt"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  SET CPML PARAMETERS IN EACH DIRECTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO i = 1,nxPML_1
      sige_x_PML_1(i) = sig_x_max * ( (nxPML_1 - i) / (nxPML_1 - 1.0) )**m
      alphae_x_PML_1(i) = alpha_x_max*((i-1.0)/(nxPML_1-1.0))**ma
      kappae_x_PML_1(i) = 1.0+(kappa_x_max-1.0)*   &
                                 ((nxPML_1 - i) / (nxPML_1 - 1.0))**m
      be_x_1(i) = EXP(-(sige_x_PML_1(i) / kappae_x_PML_1(i) +   &
                                 alphae_x_PML_1(i))*dt/epsO)
      if ((sige_x_PML_1(i) == 0.0) .and.        &
         (alphae_x_PML_1(i) == 0.0) .and. (i == nxPML_1)) then
         ce_x_1(i) = 0.0
      else
         ce_x_1(i) = sige_x_PML_1(i)*(be_x_1(i)-1.0)/       &
            (sige_x_PML_1(i)+kappae_x_PML_1(i)*alphae_x_PML_1(i)) &
            / kappae_x_PML_1(i)
      endif
   ENDDO
   DO i = 1,nxPML_1-1
      sigh_x_PML_1(i) = sig_x_max * ( (nxPML_1 - i - 0.5)/(nxPML_1-1.0))**m
      alphah_x_PML_1(i) = alpha_x_max*((i-0.5)/(nxPML_1-1.0))**ma
      kappah_x_PML_1(i) = 1.0+(kappa_x_max-1.0)*   &
                            ((nxPML_1 - i - 0.5) / (nxPML_1 - 1.0))**m
      bh_x_1(i) = EXP(-(sigh_x_PML_1(i) / kappah_x_PML_1(i) +   &
                                 alphah_x_PML_1(i))*dt/epsO)
      ch_x_1(i) = sigh_x_PML_1(i)*(bh_x_1(i)-1.0)/      &
                  (sigh_x_PML_1(i)+kappah_x_PML_1(i)*alphah_x_PML_1(i)) &
                  / kappah_x_PML_1(i)
   ENDDO

   DO i = 1,nxPML_2
      sige_x_PML_2(i) = sig_x_max * ( (nxPML_2 - i) / (nxPML_2 - 1.0) )**m
      alphae_x_PML_2(i) = alpha_x_max*((i-1.0)/(nxPML_2-1.0))**ma
      kappae_x_PML_2(i) = 1.0+(kappa_x_max-1.0)*   &
                                 ((nxPML_2 - i) / (nxPML_2 - 1.0))**m
      be_x_2(i) = EXP(-(sige_x_PML_2(i) / kappae_x_PML_2(i) +   &
                                 alphae_x_PML_2(i))*dt/epsO)
      if ((sige_x_PML_2(i) == 0.0) .and.        &
         (alphae_x_PML_2(i) == 0.0) .and. (i == nxPML_2)) then
         ce_x_2(i) = 0.0
      else
         ce_x_2(i) = sige_x_PML_2(i)*(be_x_2(i)-1.0)/       &
            (sige_x_PML_2(i)+kappae_x_PML_2(i)*alphae_x_PML_2(i)) &
            / kappae_x_PML_2(i)
      endif
   ENDDO
   DO i = 1,nxPML_2-1
      sigh_x_PML_2(i) = sig_x_max * ( (nxPML_2 - i - 0.5)/(nxPML_2-1.0))**m
      alphah_x_PML_2(i) = alpha_x_max*((i-0.5)/(nxPML_2-1.0))**ma
      kappah_x_PML_2(i) = 1.0+(kappa_x_max-1.0)*   &
                            ((nxPML_2 - i - 0.5) / (nxPML_2 - 1.0))**m
      bh_x_2(i) = EXP(-(sigh_x_PML_2(i) / kappah_x_PML_2(i) +   &
                                 alphah_x_PML_2(i))*dt/epsO)
      ch_x_2(i) = sigh_x_PML_2(i)*(bh_x_2(i)-1.0)/      &
                  (sigh_x_PML_2(i)+kappah_x_PML_2(i)*alphah_x_PML_2(i)) &
                  / kappah_x_PML_2(i)
   ENDDO

   DO j = 1,nyPML_1
      sige_y_PML_1(j) = sig_y_max * ( (nyPML_1 - j ) / (nyPML_1 - 1.0) )**m
      alphae_y_PML_1(j) = alpha_y_max*((j-1)/(nyPML_1-1.0))**ma
      kappae_y_PML_1(j) = 1.0+(kappa_y_max-1.0)*   &
                                 ((nyPML_1 - j) / (nyPML_1 - 1.0))**m
      be_y_1(j) = EXP(-(sige_y_PML_1(j) / kappae_y_PML_1(j) +   &
                                 alphae_y_PML_1(j))*dt/epsO)
      if ((sige_y_PML_1(j) == 0.0) .and.        &
         (alphae_y_PML_1(j) == 0.0) .and. (j == nyPML_1)) then
         ce_y_1(j) = 0.0
      else
         ce_y_1(j) = sige_y_PML_1(j)*(be_y_1(j)-1.0)/       &
            (sige_y_PML_1(j)+kappae_y_PML_1(j)*alphae_y_PML_1(j)) &
            / kappae_y_PML_1(j)
      endif
   ENDDO
   DO j = 1,nyPML_1-1
      sigh_y_PML_1(j) = sig_y_max * ( (nyPML_1 - j - 0.5)/(nyPML_1-1.0))**m
      alphah_y_PML_1(j) = alpha_y_max*((j-0.5)/(nyPML_1-1.0))**ma
      kappah_y_PML_1(j) = 1.0+(kappa_y_max-1.0)*   &
                            ((nyPML_1 - j - 0.5) / (nyPML_1 - 1.0))**m
      bh_y_1(j) = EXP(-(sigh_y_PML_1(j) / kappah_y_PML_1(j) +   &
                                 alphah_y_PML_1(j))*dt/epsO)
      ch_y_1(j) = sigh_y_PML_1(j)*(bh_y_1(j)-1.0)/      &
                  (sigh_y_PML_1(j)+kappah_y_PML_1(j)*alphah_y_PML_1(j)) &
                  / kappah_y_PML_1(j)
   ENDDO
   DO j = 1,nyPML_2
      sige_y_PML_2(j) = sig_y_max * ( (nyPML_2 - j ) / (nyPML_2 - 1.0) )**m
      alphae_y_PML_2(j) = alpha_y_max*((j-1)/(nyPML_2-1.0))**ma
      kappae_y_PML_2(j) = 1.0+(kappa_y_max-1.0)*   &
                                 ((nyPML_2 - j) / (nyPML_2 - 1.0))**m
      be_y_2(j) = EXP(-(sige_y_PML_2(j) / kappae_y_PML_2(j) +   &
                                 alphae_y_PML_2(j))*dt/epsO)
      if ((sige_y_PML_2(j) == 0.0) .and.        &
         (alphae_y_PML_2(j) == 0.0) .and. (j == nyPML_2)) then
         ce_y_2(j) = 0.0
      else
         ce_y_2(j) = sige_y_PML_2(j)*(be_y_2(j)-1.0)/       &
            (sige_y_PML_2(j)+kappae_y_PML_2(j)*alphae_y_PML_2(j)) &
            / kappae_y_PML_2(j)
      endif
   ENDDO
   DO j = 1,nyPML_2-1
      sigh_y_PML_2(j) = sig_y_max * ( (nyPML_2 - j - 0.5)/(nyPML_2-1.0))**m
      alphah_y_PML_2(j) = alpha_y_max*((j-0.5)/(nyPML_2-1.0))**ma
      kappah_y_PML_2(j) = 1.0+(kappa_y_max-1.0)*   &
                            ((nyPML_2 - j - 0.5) / (nyPML_2 - 1.0))**m
      bh_y_2(j) = EXP(-(sigh_y_PML_2(j) / kappah_y_PML_2(j) +   &
                                 alphah_y_PML_2(j))*dt/epsO)
      ch_y_2(j) = sigh_y_PML_2(j)*(bh_y_2(j)-1.0)/      &
                  (sigh_y_PML_2(j)+kappah_y_PML_2(j)*alphah_y_PML_2(j)) &
                  / kappah_y_PML_2(j)
   ENDDO

   DO k = 1,nzPML_1
      sige_z_PML_1(k) = sig_z_max * ( (nzPML_1 - k ) / (nzPML_1 - 1.0) )**m
      alphae_z_PML_1(k) = alpha_z_max*((k-1)/(nzPML_1-1.0))**ma
      kappae_z_PML_1(k) = 1.0+(kappa_z_max-1.0)*   &
                                 ((nzPML_1 - k) / (nzPML_1 - 1.0))**m
      be_z_1(k) = EXP(-(sige_z_PML_1(k) / kappae_z_PML_1(k) +   &
                                 alphae_z_PML_1(k))*dt/epsO)
      if ((sige_z_PML_1(k) == 0.0) .and.        &
         (alphae_z_PML_1(k) == 0.0) .and. (k == nzPML_1)) then
         ce_z_1(k) = 0.0
      else
         ce_z_1(k) = sige_z_PML_1(k)*(be_z_1(k)-1.0)/       &
            (sige_z_PML_1(k)+kappae_z_PML_1(k)*alphae_z_PML_1(k)) &
            / kappae_z_PML_1(k)
      endif
   ENDDO
   DO k = 1,nzPML_1-1
      sigh_z_PML_1(k) = sig_z_max * ( (nzPML_1 - k - 0.5)/(nzPML_1-1.0))**m
      alphah_z_PML_1(k) = alpha_z_max*((k-0.5)/(nzPML_1-1.0))**ma
      kappah_z_PML_1(k) = 1.0+(kappa_z_max-1.0)*   &
                            ((nzPML_1 - k - 0.5) / (nzPML_1 - 1.0))**m
      bh_z_1(k) = EXP(-(sigh_z_PML_1(k) / kappah_z_PML_1(k) +   &
                                 alphah_z_PML_1(k))*dt/epsO)
      ch_z_1(k) = sigh_z_PML_1(k)*(bh_z_1(k)-1.0)/      &
                  (sigh_z_PML_1(k)+kappah_z_PML_1(k)*alphah_z_PML_1(k)) &
                  / kappah_z_PML_1(k)
   ENDDO

   DO k = 1,nzPML_2
      sige_z_PML_2(k) = sig_z_max * ( (nzPML_2 - k ) / (nzPML_2 - 1.0) )**m
      alphae_z_PML_2(k) = alpha_z_max*((k-1)/(nzPML_2-1.0))**ma
      kappae_z_PML_2(k) = 1.0+(kappa_z_max-1.0)*   &
                                 ((nzPML_2 - k) / (nzPML_2 - 1.0))**m
      be_z_2(k) = EXP(-(sige_z_PML_2(k) / kappae_z_PML_2(k) +   &
                                 alphae_z_PML_2(k))*dt/epsO)
      if ((sige_z_PML_2(k) == 0.0) .and.        &
         (alphae_z_PML_2(k) == 0.0) .and. (k == nzPML_2)) then
         ce_z_2(k) = 0.0
      else
         ce_z_2(k) = sige_z_PML_2(k)*(be_z_2(k)-1.0)/       &
            (sige_z_PML_2(k)+kappae_z_PML_2(k)*alphae_z_PML_2(k)) &
            / kappae_z_PML_2(k)
      endif
   ENDDO
   DO k = 1,nzPML_2-1
      sigh_z_PML_2(k) = sig_z_max * ( (nzPML_2 - k - 0.5)/(nzPML_2-1.0))**m
      alphah_z_PML_2(k) = alpha_z_max*((k-0.5)/(nzPML_2-1.0))**ma
      kappah_z_PML_2(k) = 1.0+(kappa_z_max-1.0)*   &
                            ((nzPML_2 - k - 0.5) / (nzPML_2 - 1.0))**m
      bh_z_2(k) = EXP(-(sigh_z_PML_2(k) / kappah_z_PML_2(k) +   &
                                 alphah_z_PML_2(k))*dt/epsO)
      ch_z_2(k) = sigh_z_PML_2(k)*(bh_z_2(k)-1.0)/      &
                  (sigh_z_PML_2(k)+kappah_z_PML_2(k)*alphah_z_PML_2(k)) &
                  / kappah_z_PML_2(k)
   ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  FILL IN UPDATING COEFFICIENTS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DA = 1.0
   DB = (dt/(muO)) 

   DO i = 1,Imax
      DO j = 1,Jmax
         DO k = 1,Kmax-1
            CA(i,j,k) = (1.0 - sig(i,j,k)*dt / (2.0*eps(i,j,k))) /        &
               (1.0 + sig(i,j,k) * dt / (2.0*eps(i,j,k)))
            CB(i,j,k) = (dt/(eps(i,j,k))) /                            &
               (1.0 + sig(i,j,k)*dt / (2.0*eps(i,j,k)))
            ENDDO
	ENDDO
   ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  FILL IN DENOMINATORS FOR FIELD UPDATES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ii = nxPML_2-1
   DO i = 1,Imax-1
      if (i <= nxPML_1-1) then
         den_hx(i) = 1.0/(kappah_x_PML_1(i)*dx)
      elseif (i >= Imax+1-nxPML_2) then
         den_hx(i) = 1.0/(kappah_x_PML_2(ii)*dx)
         ii = ii-1
      else
         den_hx(i) = 1.0/dx
      endif
   ENDDO
   jj = nyPML_2-1
   DO j = 1,Jmax-1
      if (j <= nyPML_1-1) then
         den_hy(j) = 1.0/(kappah_y_PML_1(j)*dy)
      elseif (j >= Jmax+1-nyPML_2) then
         den_hy(j) = 1.0/(kappah_y_PML_2(jj)*dy)
         jj = jj-1
      else
         den_hy(j) = 1.0/dy
      endif
   ENDDO
   kk = nzPML_2-1
   DO k = 1,Kmax-1
      if (k <= nzPML_1-1) then
         den_hz(k) = 1.0/(kappah_z_PML_1(k)*dz)
      elseif (k >= Kmax+1-nzPML_2) then
         den_hz(k) = 1.0/(kappah_z_PML_2(kk)*dz)
         kk = kk - 1
      else
         den_hz(k) = 1.0/dz
      endif
   ENDDO
   ii = nxPML_2
   DO i = 1,Imax-1
      if (i <= nxPML_1) then
         den_ex(i) = 1.0/(kappae_x_PML_1(i)*dx)
      elseif (i >= Imax+1-nxPML_2) then
         den_ex(i) = 1.0/(kappae_x_PML_2(ii)*dx)
         ii = ii-1
      else
         den_ex(i) = 1.0/dx
      endif
   ENDDO
   jj = nyPML_2
   DO j = 1,Jmax-1
      if (j <= nyPML_1) then
         den_ey(j) = 1.0/(kappae_y_PML_1(j)*dy)
      elseif (j >= Jmax+1-nyPML_2) then
         den_ey(j) = 1.0/(kappae_y_PML_2(jj)*dy)
         jj = jj-1
      else
         den_ey(j) = 1.0/dy
      endif
   ENDDO
   kk = nzPML_2
   DO k = 1,Kmax-1
      if (k <= nzPML_1) then
         den_ez(k) = 1.0/(kappae_z_PML_1(k)*dz)
      elseif (k >= Kmax+1-nzPML_2) then
         den_ez(k) = 1.0/(kappae_z_PML_2(kk)*dz)
         kk = kk - 1
      else
         den_ez(k) = 1.0/dz
      endif
   ENDDO

!.:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  BEGIN TIME STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  write(*,*)"begin time-stepping"
  DO n = 1,nmax
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Hx
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 1,Kmax-1
      DO i = 1,Imax-1
	   DO j = 1,Jmax-1
	      Hx(i,j,k) = DA * Hx(i,j,k) + DB *       &
			( (Ez(i,j,k) - Ez(i,j+1,k))*den_hy(j)  +    &
			  (Ey(i,j,k+1) - Ey(i,j,k))*den_hz(k) )
	   ENDDO
	ENDDO 
      DO i = 1,Imax-1
!.....................................................................
!  PML for bottom Hx, j-direction
!.....................................................................
         DO j = 1,nyPML_1-1
  	      psi_Hxy_1(i,j,k) = bh_y_1(j)*psi_Hxy_1(i,j,k)                 &
	 			+ ch_y_1(j) *(Ez(i,j,k) - Ez(i,j+1,k))/dy
            Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxy_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hx, j-direction
!.....................................................................
         jj = nyPML_2-1
         DO j = Jmax+1-nyPML_2,Jmax-1
  	      psi_Hxy_2(i,jj,k) = bh_y_2(jj)*psi_Hxy_2(i,jj,k)       &
	 			+ ch_y_2(jj) *(Ez(i,j,k) -       &
                                Ez(i,(j+1),k))/dy
            Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxy_2(i,jj,k)
            jj = jj-1
         ENDDO
	ENDDO
   ENDDO
   DO i = 1,Imax-1
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Hx, k-direction
!.....................................................................
         DO k = 1,nzPML_1-1
  	      psi_Hxz_1(i,j,k) = bh_z_1(k)*psi_Hxz_1(i,j,k)            &
	 			+ ch_z_1(k) *(Ey(i,j,k+1) - Ey(i,j,k))/dz
            Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxz_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hx, k-direction
!.....................................................................
         kk = nzPML_2-1
         DO k = Kmax+1-nzPML_2,Kmax-1
  	      psi_Hxz_2(i,j,kk) = bh_z_2(kk)*psi_Hxz_2(i,j,kk)             &
	 			+ ch_z_2(kk) *(Ey(i,j,k+1) -       &
                                Ey(i,j,k))/dz
            Hx(i,j,k) = Hx(i,j,k) + DB*psi_Hxz_2(i,j,kk)
           kk = kk-1
         ENDDO
	ENDDO
   ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Hy
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 1,Kmax-1
      DO i = 1,Imax-1
	   DO j = 1,Jmax-1
            Hy(i,j,k) = DA * Hy(i,j,k) + DB *             &
			( (Ez(i+1,j,k) - Ez(i,j,k))*den_hx(i) +      &
			  (Ex(i,j,k) - Ex(i,j,k+1))*den_hz(k) )
	   ENDDO 
      ENDDO
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Hy, i-direction
!.....................................................................
         DO i = 1,nxPML_1-1
	      psi_Hyx_1(i,j,k) = bh_x_1(i)*psi_Hyx_1(i,j,k)                   &
				+ ch_x_1(i)*(Ez(i+1,j,k) - Ez(i,j,k))/dx
	      Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyx_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hy, i-direction
!.....................................................................
         ii = nxPML_2-1
         DO i = Imax+1-nxPML_2,Imax-1
	      psi_Hyx_2(ii,j,k) = bh_x_2(ii)*psi_Hyx_2(ii,j,k)     &
				+ ch_x_2(ii)*(Ez(i+1,j,k) -    &
                                Ez(i,j,k))/dx
	      Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyx_2(ii,j,k)
            ii = ii-1
	   ENDDO
     ENDDO
   ENDDO
   DO i = 1,Imax-1
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Hy, k-direction
!.....................................................................
         DO k = 1,nzPML_1-1
	      psi_Hyz_1(i,j,k) = bh_z_1(k)*psi_Hyz_1(i,j,k)          &
				+ ch_z_1(k)*(Ex(i,j,k) - Ex(i,j,k+1))/dz
	      Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyz_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hy, k-direction
!.....................................................................
         kk = nzPML_2-1
         DO k = Kmax+1-nzPML_2,Kmax-1
	    psi_Hyz_2(i,j,kk) = bh_z_2(kk)*psi_Hyz_2(i,j,kk)               &
				+ ch_z_2(kk)*(Ex(i,j,k) -    &
                                Ex(i,j,k+1))/dz
	    Hy(i,j,k) = Hy(i,j,k) + DB*psi_Hyz_2(i,j,kk)
            kk = kk-1
         ENDDO
     ENDDO
   ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Hz
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 2,Kmax-1
      DO i = 1,Imax-1
        DO j = 1,Jmax-1
            Hz(i,j,k) = DA * Hz(i,j,k) + DB       &
                  * ((Ey(i,j,k) - Ey(i+1,j,k))*den_hx(i) +        &
			    (Ex(i,j+1,k) - Ex(i,j,k))*den_hy(j))
	   ENDDO
      ENDDO
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Hz, x-direction
!.....................................................................
         DO i = 1,nxPML_1-1
   	      psi_Hzx_1(i,j,k) = bh_x_1(i)*psi_Hzx_1(i,j,k)                 &
	 			+ ch_x_1(i) *(Ey(i,j,k) - Ey(i+1,j,k))/dx
	      Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzx_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hz, x-direction
!.....................................................................
         ii = nxPML_2-1
         DO i = Imax+1-nxPML_2,Imax-1
   	      psi_Hzx_2(ii,j,k) = bh_x_2(ii)*psi_Hzx_2(ii,j,k)            &
	 			+ ch_x_2(ii) *(Ey(i,j,k) -       &
                                Ey(i+1,j,k))/dx
	      Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzx_2(ii,j,k)
            ii = ii-1
	   ENDDO
      ENDDO
      DO i = 1,Imax-1
!.....................................................................
!  PML for bottom Hz, y-direction
!.....................................................................
         DO j = 1,nyPML_1-1
            psi_Hzy_1(i,j,k) = bh_y_1(j)*psi_Hzy_1(i,j,k)                   &
				+ ch_y_1(j)*(Ex(i,j+1,k) - Ex(i,j,k))/dy
	      Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzy_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Hz, y-direction
!.....................................................................
         jj = nyPML_2-1
         DO j = Jmax+1-nyPML_2,Jmax-1
            psi_Hzy_2(i,jj,k) = bh_y_2(jj)*psi_Hzy_2(i,jj,k)               &
				+ ch_y_2(jj)*(Ex(i,j+1,k) -    &
                                Ex(i,j,k))/dy
	      Hz(i,j,k) = Hz(i,j,k) + DB*psi_Hzy_2(i,jj,k)
            jj = jj-1
	   ENDDO
      ENDDO
   ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Ex
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 2,Kmax-1
      DO i = 1,Imax-1
	   DO j = 2,Jmax-1
              IF (i >= istart .and. i <= iend .and. j >= jstart .and.  &
                   j <= jend .and. k >= kstart .and. k <= kend) THEN
                 Ex(i,j,k) = 0.0
              ELSE
	         Ex(i,j,k) = CA(i,j,k) * Ex(i,j,k) + CB(i,j,k) *       &
  			( (Hz(i,j,k) - Hz(i,j-1,k))*den_ey(j)  +    &
			  (Hy(i,j,k-1) - Hy(i,j,k))*den_ez(k) )
              ENDIF
	   ENDDO
	ENDDO 
      DO i = 1,Imax-1
!.....................................................................
!  PML for bottom Ex, j-direction
!.....................................................................
         DO j = 2,nyPML_1
  	      psi_Exy_1(i,j,k) = be_y_1(j)*psi_Exy_1(i,j,k)                 &
	 			+ ce_y_1(j) *(Hz(i,j,k) - Hz(i,j-1,k))/dy
            Ex(i,j,k) = Ex(i,j,k) + CB(i,j,k)*psi_Exy_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ex, j-direction
!.....................................................................
         jj = nyPML_2
         DO j = Jmax+1-nyPML_2,Jmax-1
  	      psi_Exy_2(i,jj,k) = be_y_2(jj)*psi_Exy_2(i,jj,k)       &
	 			+ ce_y_2(jj) *(Hz(i,j,k) -       &
                                Hz(i,(j-1),k))/dy
            Ex(i,j,k) = Ex(i,j,k) + CB(i,j,k)*psi_Exy_2(i,jj,k)
            jj = jj-1
         ENDDO
	ENDDO
   ENDDO
   DO i = 1,Imax-1
      DO j = 2,Jmax-1
!.....................................................................
!  PML for bottom Ex, k-direction
!.....................................................................
         DO k = 2,nzPML_1
  	      psi_Exz_1(i,j,k) = be_z_1(k)*psi_Exz_1(i,j,k)                 &
	 			+ ce_z_1(k) *(Hy(i,j,k-1) - Hy(i,j,k))/dz
            Ex(i,j,k) = Ex(i,j,k) + CB(i,j,k)*psi_Exz_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ex, k-direction
!.....................................................................
         kk = nzPML_2
         DO k = Kmax+1-nzPML_2,Kmax-1
  	      psi_Exz_2(i,j,kk) = be_z_2(kk)*psi_Exz_2(i,j,kk)             &
	 			+ ce_z_2(kk) *(Hy(i,j,k-1) -       &
                                Hy(i,j,k))/dz
            Ex(i,j,k) = Ex(i,j,k) + CB(i,j,k)*psi_Exz_2(i,j,kk)
            kk = kk-1
         ENDDO
	ENDDO
   ENDDO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Ey
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 2,Kmax-1
      DO i = 2,Imax-1
	   DO j = 1,Jmax-1
              IF (i >= istart .and. i <= iend .and. j >= jstart .and. &
                   j <= jend .and. k >= kstart .and. k <= kend) THEN
                 Ey(i,j,k) = 0.0
              ELSE
                 Ey(i,j,k) = CA(i,j,k) * Ey(i,j,k) + CB(i,j,k) *    &
			( (Hz(i-1,j,k) - Hz(i,j,k))*den_ex(i) +         &
			  (Hx(i,j,k) - Hx(i,j,k-1))*den_ez(k) )
              ENDIF
	   ENDDO 
      ENDDO
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Ey, i-direction
!.....................................................................
         DO i = 2,nxPML_1
	      psi_Eyx_1(i,j,k) = be_x_1(i)*psi_Eyx_1(i,j,k)                   &
				+ ce_x_1(i)*(Hz(i-1,j,k) - Hz(i,j,k))/dx
	      Ey(i,j,k) = Ey(i,j,k) + CB(i,j,k)*psi_Eyx_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ey, i-direction
!.....................................................................
         ii = nxPML_2
         DO i = Imax+1-nxPML_2,Imax-1
	      psi_Eyx_2(ii,j,k) = be_x_2(ii)*psi_Eyx_2(ii,j,k)              &
				+ ce_x_2(ii)*(Hz(i-1,j,k) -    &
                               Hz(i,j,k))/dx
	      Ey(i,j,k) = Ey(i,j,k) + CB(i,j,k)*psi_Eyx_2(ii,j,k)
            ii = ii-1
	   ENDDO
     ENDDO
   ENDDO
   DO i = 2,Imax-1
      DO j = 1,Jmax-1
!.....................................................................
!  PML for bottom Ey, k-direction
!.....................................................................
         DO k = 2,nzPML_1
	      psi_Eyz_1(i,j,k) = be_z_1(k)*psi_Eyz_1(i,j,k)                   &
				+ ce_z_1(k)*(Hx(i,j,k) - Hx(i,j,k-1))/dz
	      Ey(i,j,k) = Ey(i,j,k) + CB(i,j,k)*psi_Eyz_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ey, k-direction
!.....................................................................
         kk = nzPML_2
         DO k = Kmax+1-nzPML_2,Kmax-1
	    psi_Eyz_2(i,j,kk) = be_z_2(kk)*psi_Eyz_2(i,j,kk)               &
				+ ce_z_2(kk)*(Hx(i,j,k) -    &
                                Hx(i,j,k-1))/dz
	    Ey(i,j,k) = Ey(i,j,k) + CB(i,j,k)*psi_Eyz_2(i,j,kk)
            kk = kk-1
         ENDDO
     ENDDO
   ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  UPDATE Ez
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DO k = 1,Kmax-1
      DO i = 2,Imax-1
         DO j = 2,Jmax-1
            Ez(i,j,k) = CA(i,j,k) * Ez(i,j,k) + CB(i,j,k)       &
                  * ((Hy(i,j,k) - Hy(i-1,j,k))*den_ex(i) +        &
			    (Hx(i,j-1,k) - Hx(i,j,k))*den_ey(j))
	   ENDDO
      ENDDO
      DO j = 2,Jmax-1
!.....................................................................
!  PML for bottom Ez, x-direction
!.....................................................................
         DO i = 2,nxPML_1
   	      psi_Ezx_1(i,j,k) = be_x_1(i)*psi_Ezx_1(i,j,k)             &
	 			+ ce_x_1(i) *(Hy(i,j,k) - Hy(i-1,j,k))/dx
	      Ez(i,j,k) = Ez(i,j,k) + CB(i,j,k)*psi_Ezx_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ez, x-direction
!.....................................................................
         ii = nxPML_2
         DO i = Imax+1-nxPML_2,Imax-1
   	      psi_Ezx_2(ii,j,k) = be_x_2(ii)*psi_Ezx_2(ii,j,k)       &
	 			+ ce_x_2(ii) *(Hy(i,j,k) -       &
                                Hy(i-1,j,k))/dx
	      Ez(i,j,k) = Ez(i,j,k) + CB(i,j,k)*psi_Ezx_2(ii,j,k)
            ii = ii-1
	   ENDDO
      ENDDO
      DO i = 2,Imax-1
!.....................................................................
!  PML for bottom Ez, y-direction
!.....................................................................
         DO j = 2,nyPML_1
            psi_Ezy_1(i,j,k) = be_y_1(j)*psi_Ezy_1(i,j,k)            &
				+ ce_y_1(j)*(Hx(i,j-1,k) - Hx(i,j,k))/dy
	      Ez(i,j,k) = Ez(i,j,k) + CB(i,j,k)*psi_Ezy_1(i,j,k)
         ENDDO
!.....................................................................
!  PML for top Ez, y-direction
!.....................................................................
         jj = nyPML_2
         DO j = Jmax+1-nyPML_2,Jmax-1
            psi_Ezy_2(i,jj,k) = be_y_2(jj)*psi_Ezy_2(i,jj,k)         &
				+ ce_y_2(jj)*(Hx(i,j-1,k) -    &
                                Hx(i,j,k))/dy
	      Ez(i,j,k) = Ez(i,j,k) + CB(i,j,k)*psi_Ezy_2(i,jj,k)
            jj = jj-1
	   ENDDO
      ENDDO
   ENDDO

!-----------------------------------------------------------------------
!   SOURCE
!-----------------------------------------------------------------------
   i = isource
   j = jsource
   k = ksource
   source = -2.0*((n*dt-tO)/tw) * exp(-((n*dt-tO)/tw)**2.0)
   Ez(i,j,k) = Ez(i,j,k) - CB(i,j,k)*source

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  RECORD GRID FOR VISUALIZATION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IF (n == record_grid) then
      DO j = 1,Jmax
         DO i = 1,Imax      
            write(33,*)Ez(i,j,ksource)
	   ENDDO
	ENDDO
   CLOSE(UNIT = 33)
   endif   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  WRITE TO OUTPUT FILES
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
   P1 = Ez(isource,jsource,ksource)
   write(30,*)P1      
   P2 = Ey(irecv1,jrecv1,krecv1)
   write(31,*)P2      

   IF (mod(n,10) == 0) then
    WRITE(*,*)n, " of ", nmax
   ENDIF

   ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  END TIME STEP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!.:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:.
    WRITE(*,*)"done time-stepping"

!-----------------------------------------------------------------------
! CLOSE OUTPUT FILES
!-----------------------------------------------------------------------   
   CLOSE(UNIT = 30)
   CLOSE(UNIT = 31)

   END PROGRAM fdtd3D_CPML
!cccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5cccgtgccc6ccccccccc7cc

