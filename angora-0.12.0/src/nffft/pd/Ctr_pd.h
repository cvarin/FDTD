/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CTR_PD_H
#define CTR_PD_H

//Declaration of the abstract base class "Ctr_pd" for a PHASOR-DOMAIN near-field-to-far-field transformer
//This class takes a 2D range of observation directions, defined parametrically by a 2D matrix of size (num_of_dirs_1,num_of_dirs_2).
//The values in this matrix are converted to spherical longitude and azimuth angles THETA and PHI, using the direction specification variable "directionspec". This variable provides all the necessary information for converting the direction-parameter matrix (dir1,dir2) to angle values (THETA,PHI).
//Currently, only two parametrizations are supported:
// 1) Direct THETA and PHI parametrization : (dir1,dir2) = (THETA,PHI)  [in degrees]
// 2) Direction-cosine [wrt the x and y axes] parametrization : (dir1,dir2) = (sin(theta)cos(phi),sin(theta)sin(phi))  [the direction specification variable "directionspec" also specifies whether the observation direction is in the upper (+z) or lower (-z) hemisphere]

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

#include <fstream>

//only the declaration of Cpointsources needed: use forward declaration
class Cpointsources;
//only the declaration of Ctfsf needed: use forward declaration
class Ctfsf;

// when this is defined, the MPI 2.0 construct "MPI_IN_PLACE" is used in MPI_Reduce commands for eliminating the buffer arrays
#define USE_MPI_IN_PLACE

extern int OriginX,OriginY,OriginZ;


class TrDataType_pd
{
 public:
	 TrDataType_pd(const Array<double,1>& mylambda, const Array<double,1>& mydir1, const Array<double,1>& mydir2, const double& mymax_s, const string& mydirectionspec):
	   	lambda(mylambda),
	   	dir1(mydir1),
	   	dir2(mydir2),
		directionspec(mydirectionspec),
		max_s(mymax_s),
		scale_with_wavelength(false),
		NFFFTMarginBackX(3),
		NFFFTMarginFrontX(3),
		NFFFTMarginLeftY(3),
		NFFFTMarginRightY(3),
		NFFFTMarginLowerZ(3),
		NFFFTMarginUpperZ(3),
		NFFFTOriginX(OriginX),
		NFFFTOriginY(OriginY),
		NFFFTOriginZ(OriginZ),
	   	PointSourcesPtr(NULL),	//if left NULL, theoretical far-field is left 0
	   	TFSFPtr(NULL)	//if left NULL, no PW data is written
	   {};

	 const Array<double,1> lambda;	//wavelengths at which freq-domain NFFFT is evaluated (in meters)
	 const Array<double,1> dir1;	//first parametrization variable of the 2D observation direction range
	 const Array<double,1> dir2;	//second parametrization variable of the 2D observation direction range
	 const string directionspec;	//specifies the manner in which (dir1,dir2) is mapped to spherical angles (THETA,PHI)
	 const double max_s; //maximum absolute direction cosine (sx^2+sy^2)^1/2 [only used if "directionspec" is "dircosx-dircosy-upper" or "dircosx-dircosy-lower"]
	 bool scale_with_wavelength; //are the direction cosines scaled linearly with the wavelength? (useful for optical imaging)
	 int NFFFTMarginBackX,NFFFTMarginFrontX,NFFFTMarginLeftY,NFFFTMarginRightY,NFFFTMarginLowerZ,NFFFTMarginUpperZ;	//distances (in cells) of the NFFFT box from the PML boundary
	 double NFFFTOriginX,NFFFTOriginY,NFFFTOriginZ;	 //the index of the cell assigned as the origin of the NFFFT box
	 //(its rear-left-lower corner is the origin)
	 //(can be fractional, therefore defined as "double")
	 const Cpointsources* PointSourcesPtr;	//points to a non-changeable Cpointsources object
	 const Ctfsf* TFSFPtr;	//points to a non-changeable CTFSF object
};

class Ctr_pd
{
 public:
	 Ctr_pd(const TrDataType_pd& MyData, const string& FileName, const int& Index);	//constructor
 	 virtual ~Ctr_pd(){};	//virtual destructor for cleaning up some objects (e.g. ofstream object)
						//this is needed when deleting derived objects using a base class pointer

	 void UpdateFarField(const int& n);		//Do some auxiliary updates using the current E,H on the virtual surface
	 virtual void ConstructFarField()=0;	//Do the necessary post-processing steps to obtain the final far field

 protected:
	 void UpdateElectric(const int& n);	//Update far-field arrays using the H values on the virtual surface
	 void UpdateMagnetic(const int& n);	//Update far-field arrays using the E values on the virtual surface

	 void UpdatePWPhasors(const int& n); //Update the phasors for the scattered PW components

 	 virtual void ConstructPotential_A(const int& l, const int& d1, const int& d2)=0;	//computes the theta and phi components of the magnetic potential A
 	 virtual void ConstructPotential_F(const int& l, const int& d1, const int& d2)=0;	//computes the theta and phi components of the electric potential F

	 virtual complex<double> TheoreticalFarFieldTheta(const int& l, const int& d1, const int& d2)=0;	//the theoretical theta-component of the far field (l: wavenumber, d1: first direction parameter, d2: second direction parameter)
	 virtual complex<double> TheoreticalFarFieldPhi(const int& l, const int& d1, const int& d2)=0;	//the theoretical phi-component of the far field (l: wavenumber, d1: first direction parameter, d2: second direction parameter)

	 void GatherAndWriteFarField();		//gathers the far-field arrays from the nodes and writes them into file

	 double GridVelocity(const double& w, const double& dir1, const double& dir2); //calculates the grid velocity v(dir1,dir2) at angular frequency w, normalized by c

     //functions that compute for some direction-related parameters (angles, direction cosines, etc)
     // Although these are direction-related parameters, their definition may depend on the wavelength; hence the l-dependence.
	 //spherical observation angles
	 double THETA(const double& sx, const double& sy, const double& sz);//theta angle corresponding to a given set of x,y,z direction cosines
	 double THETA(const int& l, const int& d1, const int& d2);
	 double PHI(const double& sx, const double& sy);//phi angle corresponding to a given set of x,y direction cosines (it is independent of the z-direction cosine)
	 double PHI(const int& l, const int& d1, const int& d2);
	 //boolean flag specifying whether there is a valid observation angle
	 //if not, the far fields are left at 0, and no post processing is done for those angles
	 bool observation_angle_exists(const int& l, const int& d1, const int& d2);
	 //Direction sines and cosines
	 double SinT(const int& l, const int& d1, const int& d2);
	 double CosT(const int& l, const int& d1, const int& d2);
	 double SinP(const int& l, const int& d1, const int& d2);
	 double CosP(const int& l, const int& d1, const int& d2);
	 double SinTCosP(const int& l, const int& d1, const int& d2);
	 double SinTSinP(const int& l, const int& d1, const int& d2);
	 double CosTSinP(const int& l, const int& d1, const int& d2);
	 double CosTCosP(const int& l, const int& d1, const int& d2);

	 int L,D1,D2;	//sizes of the lambda, 1st direction, and 2nd direction arrays

	 const int TransformerIndex;	//index of the current plane-wave in the NFFFT object Cnffft (1st, 2nd,.. etc.)

	 const TrDataType_pd Data;	//data type that holds the transformer data

	 //wavelength array (defined for FREE SPACE)
	 Array<double,1> lambda;
	 //observation direction parameters
	 Array<double,1> dir1,dir2;
	 //maximum (absolute) direction cosine
	 double max_s;

	 //temporary variables that hold the x and y direction cosines for an observation angle
	 double sx,sy;

     //smallest wavelength in the wavelength array
     double min_lambda;

	 //temporal frequency array
	 Array<double,1> ww;	//in sec^-1

	 //the index of the cell assigned as the origin of the NFFFT box
	 //(its rear-left-lower corner is the origin)
	 double FarFieldOriginX,FarFieldOriginY,FarFieldOriginZ;

	 //Indices of the cells immediately inside the NFFFT box
	 int SurfaceBackX,SurfaceFrontX,SurfaceLeftY,SurfaceRightY,SurfaceLowerZ,SurfaceUpperZ;

	 int l,d1,d2;	//wavenumber and angle counters

	//complex exponential Fourier kernels for the Fourier transform integrals (1st dim: frequency, 2nd: time step)
	//phasor_Mt is for the magnetic currents (coincident in time with the E-field, no extra delay necessary)
	//phasor_Jt is for the electric currents (coincident in time with the H-field, which is 0.5t ahead of E)
	Array<complex<double>,2> phasors_Jt,phasors_Mt;

	//Phasor storage arrays (first dim: frequency, second&third dims: spatial index along the 2D section of the NFFFT surface)
	Array<complex<double>,3> J_x_lower,J_x_upper,J_x_left,J_x_right;
	Array<complex<double>,3> J_y_lower,J_y_upper,J_y_back,J_y_front;
	Array<complex<double>,3> J_z_left,J_z_right,J_z_back,J_z_front;
	Array<complex<double>,3> M_x_lower,M_x_upper,M_x_left,M_x_right;
	Array<complex<double>,3> M_y_lower,M_y_upper,M_y_back,M_y_front;
	Array<complex<double>,3> M_z_left,M_z_right,M_z_back,M_z_front;

	 // Do the faces of the NFFFT surface "pass through" this node?
	 // A face of the NFFFT surface is defined to be "passing through the node" when at least one of the field components on the face "belongs" to one of the cells in the node. Note that the E-field components on the E-shell are right on the NFFFT box, and the H-field components on the H-shell are half-cell away from the NFFFT box. Thus, the definition of "passing through the node" is different for each field component (Ex, Hy, etc.) on the NFFFT surface. For example, if an NFFFT surface coincides with the lower boundary of node A, the H-field components half-cell below the NFFFT surface belong to the node below node A, while the E-field components on the NFFFT surface belong to node A itself.
	 // "Belonging to a cell" is defined as usual as being at the back/left/lower side of a cell (either face-wise or side-wise).
	 // The "downside" (!) of this definition is that the frontmost, righmost and uppermost E or H-field components in the grid don't belong to any cell, which is not a big deal.
	 // The NFFFT box is not supposed to touch the rightmost and uppermost boundary of the FDTD grid anyway.
	 bool Back_Ey_PassesThroughNode,Back_Ez_PassesThroughNode,Back_Hy_PassesThroughNode,Back_Hz_PassesThroughNode;
	 bool Front_Ey_PassesThroughNode,Front_Ez_PassesThroughNode,Front_Hy_PassesThroughNode,Front_Hz_PassesThroughNode;
	 bool Left_Ex_PassesThroughNode,Left_Ez_PassesThroughNode,Left_Hx_PassesThroughNode,Left_Hz_PassesThroughNode;
	 bool Right_Ex_PassesThroughNode,Right_Ez_PassesThroughNode,Right_Hx_PassesThroughNode,Right_Hz_PassesThroughNode;
	 bool Lower_Ex_PassesThroughNode,Lower_Ey_PassesThroughNode,Lower_Hx_PassesThroughNode,Lower_Hy_PassesThroughNode;
	 bool Upper_Ex_PassesThroughNode,Upper_Ey_PassesThroughNode,Upper_Hx_PassesThroughNode,Upper_Hy_PassesThroughNode;

	 // the back-front, left-right and lower-upper limits of the indices of the field components on the NFFFT surface
	 int Ex_backlimit,Ey_backlimit,Ez_backlimit,Hx_backlimit,Hy_backlimit,Hz_backlimit;
	 int Ex_frontlimit,Ey_frontlimit,Ez_frontlimit,Hx_frontlimit,Hy_frontlimit,Hz_frontlimit;
	 int Ex_leftlimit,Ey_leftlimit,Ez_leftlimit,Hx_leftlimit,Hy_leftlimit,Hz_leftlimit;
	 int Ex_rightlimit,Ey_rightlimit,Ez_rightlimit,Hx_rightlimit,Hy_rightlimit,Hz_rightlimit;
	 int Ex_lowerlimit,Ey_lowerlimit,Ez_lowerlimit,Hx_lowerlimit,Hy_lowerlimit,Hz_lowerlimit;
	 int Ex_upperlimit,Ey_upperlimit,Ez_upperlimit,Hx_upperlimit,Hy_upperlimit,Hz_upperlimit;
	 // the limits of the indices of these components that "belong" to the current node (see definition above)
	 int Ex_backlimit_in_node,Ey_backlimit_in_node,Ez_backlimit_in_node,Hx_backlimit_in_node,Hy_backlimit_in_node,Hz_backlimit_in_node;
	 int Ex_frontlimit_in_node,Ey_frontlimit_in_node,Ez_frontlimit_in_node,Hx_frontlimit_in_node,Hy_frontlimit_in_node,Hz_frontlimit_in_node;
	 int Ex_leftlimit_in_node,Ey_leftlimit_in_node,Ez_leftlimit_in_node,Hx_leftlimit_in_node,Hy_leftlimit_in_node,Hz_leftlimit_in_node;
	 int Ex_rightlimit_in_node,Ey_rightlimit_in_node,Ez_rightlimit_in_node,Hx_rightlimit_in_node,Hy_rightlimit_in_node,Hz_rightlimit_in_node;
	 int Ex_lowerlimit_in_node,Ey_lowerlimit_in_node,Ez_lowerlimit_in_node,Hx_lowerlimit_in_node,Hy_lowerlimit_in_node,Hz_lowerlimit_in_node;
	 int Ex_upperlimit_in_node,Ey_upperlimit_in_node,Ez_upperlimit_in_node,Hx_upperlimit_in_node,Hy_upperlimit_in_node,Hz_upperlimit_in_node;

	 //indices of the field components at the back,front,left,right,lower,upper parts of the NFFFT surface
	 int Back_E_shell_index, Front_E_shell_index, Left_E_shell_index, Right_E_shell_index, Lower_E_shell_index, Upper_E_shell_index;
	 int Back_H_shell_index, Front_H_shell_index, Left_H_shell_index, Right_H_shell_index, Lower_H_shell_index, Upper_H_shell_index;

	 //In the following arrays, first dimension is wavelength, second and third are direction parameters (dir1 and dir2, respectively)
	 //Final far-field arrays
	 Array<complex<double>,3> E_theta;	//Theta-component of the radiated electric field (normalized by 1/r, advanced by r+rOffset)
	 Array<complex<double>,3> E_phi;	//Phi component of the radiated electric field (normalized by 1/r, advanced by r+rOffset)
	 //Auxiliary far-field potential arrays
	 Array<complex<double>,3> A_theta;	//Theta-component of the magnetic potential A
	 Array<complex<double>,3> A_phi;	//Phi component of the magnetic potential A
	 Array<complex<double>,3> F_theta;	//Theta-component of the electric potential F
	 Array<complex<double>,3> F_phi;	//Phi component of the electric potential F
#ifndef USE_MPI_IN_PLACE
	 Array<complex<double>,3> E_theta_buf;	//buffer array for MPI_REDUCE
	 Array<complex<double>,3> E_phi_buf;	//buffer array for MPI_REDUCE

	 Array<complex<double>,3> A_theta_buf;	//buffer array for MPI_REDUCE
	 Array<complex<double>,3> A_phi_buf;	//buffer array for MPI_REDUCE
	 Array<complex<double>,3> F_theta_buf;	//buffer array for MPI_REDUCE
	 Array<complex<double>,3> F_phi_buf;	//buffer array for MPI_REDUCE
#endif

#ifndef MPI_DISABLE
	 MPI_Status Status;
#endif
	 ofstream FarFieldFile;		//far-field file
	 const string FarFieldFileName; //name of the image output file

	 double a;	//interpolation parameter

	 double Jx,Jy,Jz,Mx,My,Mz;

	 //Scattered plane-wave information (from the Ctfsf class)
	 //These are not initialized in Ctr_pd, since they are resized as necessary within Ctfsf
	 Array<double,1> PW_THETA,PW_PHI;	//scattering angles
	 Array<double,1> E_x_array,E_y_array;	//the x and y components of the E-field
	 Array<double,1> origindelay_array;		//amount of time by which "field_array" is delayed compared to the field at the transformer origin
	 Array<complex<double>,2> E_x_phasor_array,E_y_phasor_array;	//phasor of E_theta or E_phi of the scattered PW at the origin (first dim: frequency, second dim: plane-wave index)
};

//inline function definitions (inlined into the code whenever they are called from a .cpp that includes this .h file)
inline double Ctr_pd::THETA(const double& sx, const double& sy, const double& sz)
{//theta angle corresponding to the given x,y,z direction cosines
	return atan2(sqrt(pow2(sx)+pow2(sy)),sz);//guaranteed to be between [0,pi] because sqrt(sx^2+sy^2) is always positive
}

inline double Ctr_pd::PHI(const double& sx, const double& sy)
{//phi angle corresponding to the given x,y,z direction cosines
	//is sy is positive, atan2 is between [0,pi], so nothing is needed
	//is sy is negative, atan2 is between [-pi,0], so we make it positive by adding 2*pi
	return ((sy<0)?(atan2(sy,sx)+2*M_PI):atan2(sy,sx));
}
#endif
