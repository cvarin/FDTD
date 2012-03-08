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

//Definition of the abstract base class "Ctr_pd" for a PHASOR-DOMAIN near-field-to-far-field transformer

#include "headers.h"

#include "Ctr_pd.h"

//definition of Ctfsf needed
#include "tfsf/Ctfsf.h"

#include "float_comp.h"

//routines for getting/setting the initial time value in the simulation
#include "time_axis.h"


extern double courant,dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;
extern int NSTEPS;

extern int OriginX,OriginY,OriginZ;

extern double c_upper,c_lower,epsilon_r_upper,epsilon_r_lower,mu_r_upper,mu_r_lower;

extern int GridIndex;
extern int rank;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern void MPI_exit(const int& exitcode);


Ctr_pd::Ctr_pd(const TrDataType_pd& MyData, const string& FileName, const int& Index)
		: Data(MyData), FarFieldFileName(FileName), TransformerIndex(Index)
{

    if ((Data.directionspec!="theta-phi")&&(Data.directionspec!="dircosx-dircosy-upper")&&(Data.directionspec!="dircosx-dircosy-lower"))
    {
        if (rank==0)
        {
            cout << "Error: Invalid direction specification (" << Data.directionspec << ") for phasor-domain near-field-to-far-field transformer (NFFFT) " << TransformerIndex << " in grid " << GridIndex << endl;
            exit(-1);
        }
    }

	//NFFFT wavelengths
	lambda.resize(Data.lambda.size());
	lambda = Data.lambda;

    //smallest wavelength in the wavelength array
    min_lambda = min(abs(lambda));

	//Far-field observation direction parameters
	dir1.resize(Data.dir1.size());
	dir2.resize(Data.dir2.size());
	dir1 = Data.dir1;
	dir2 = Data.dir2;
	max_s = Data.max_s;

	// sizes of the lambda, 1st direction, and 2nd direction arrays
	L = lambda.size();
	D1 = dir1.size();
	D2 = dir2.size();

	//construct temporal frequency array
	ww.resize(L);
	for (int l=0; l<L; l++)
	{
		ww(l) = 2*M_PI*c/lambda(l);	//note that lambda is defined for free space
	}

	//Origin of the far-field graphs (r is measured from this point)
	FarFieldOriginX = Data.NFFFTOriginX;
	FarFieldOriginY = Data.NFFFTOriginY;
	FarFieldOriginZ = Data.NFFFTOriginZ;

	//Indices of the cells immediately inside the NFFFT box
	SurfaceBackX=NPML+Data.NFFFTMarginBackX+1;
	SurfaceFrontX=NCELLS_X+NPML-Data.NFFFTMarginFrontX;
	SurfaceLeftY=NPML+Data.NFFFTMarginLeftY+1;
	SurfaceRightY=NCELLS_Y+NPML-Data.NFFFTMarginRightY;
	SurfaceLowerZ=NPML+Data.NFFFTMarginLowerZ+1;
	SurfaceUpperZ=NCELLS_Z+NPML-Data.NFFFTMarginUpperZ;

	//allocate and initialize the far-field arrays
	E_theta.resize(L,D1,D2);
	E_theta=0;
	E_phi.resize(L,D1,D2);
	E_phi=0;

	//allocate and initialize the potential arrays
	A_theta.resize(L,D1,D2);
	A_theta=0;
	A_phi.resize(L,D1,D2);
	A_phi=0;
	F_theta.resize(L,D1,D2);
	F_theta=0;
	F_phi.resize(L,D1,D2);
	F_phi=0;

#ifndef USE_MPI_IN_PLACE
	A_theta_buf.resize(L,D1,D2);
	A_theta_buf=0;
	A_phi_buf.resize(L,D1,D2);
	A_phi_buf=0;
	F_theta_buf.resize(L,D1,D2);
	F_theta_buf=0;
	F_phi_buf.resize(L,D1,D2);
	F_phi_buf=0;
	E_theta_buf.resize(L,D1,D2);
	E_theta_buf=0;
	E_phi_buf.resize(L,D1,D2);
	E_phi_buf=0;
#endif

	//allocate and fill the phasor arrays
	//temporal phasors
	phasors_Mt.resize(L,NSTEPS);
	phasors_Jt.resize(L,NSTEPS);
	phasors_Mt=0;
	phasors_Jt=0;
	for (int l=0;l<L;l++)
	{
		for (int n=0;n<NSTEPS;n++)
		{
			//The 1/(2pi) factor is to get rid of the same factor in the inverse transform, such that the result represents the true weight of the desired frequency component.
			phasors_Mt(l,n)=1.0/(2*M_PI)*exp(-ii*ww(l)*(n*dt+get_initial_time_value()));
			phasors_Jt(l,n)=1.0/(2*M_PI)*exp(-ii*ww(l)*((n+0.5)*dt+get_initial_time_value()));
		}
	}

	// the minimum and maximum indices of the field components on parts of the NFFFT surface, in transverse directions (e.g. z index on the left surface)
	// Ex
	Ex_backlimit = SurfaceBackX;
	Ex_frontlimit = SurfaceFrontX;
	Ex_leftlimit = SurfaceLeftY;
	Ex_rightlimit = SurfaceRightY+1;
	Ex_lowerlimit = SurfaceLowerZ;
	Ex_upperlimit = SurfaceUpperZ+1;
	// Ey
	Ey_backlimit = SurfaceBackX;
	Ey_frontlimit = SurfaceFrontX+1;
	Ey_leftlimit = SurfaceLeftY;
	Ey_rightlimit = SurfaceRightY;
	Ey_lowerlimit = SurfaceLowerZ;
	Ey_upperlimit = SurfaceUpperZ+1;
	// Ez
	Ez_backlimit = SurfaceBackX;
	Ez_frontlimit = SurfaceFrontX+1;
	Ez_leftlimit = SurfaceLeftY;
	Ez_rightlimit = SurfaceRightY+1;
	Ez_lowerlimit = SurfaceLowerZ;
	Ez_upperlimit = SurfaceUpperZ;
	// Hx
	Hx_backlimit = SurfaceBackX;
	Hx_frontlimit = SurfaceFrontX+1;
	Hx_leftlimit = SurfaceLeftY;
	Hx_rightlimit = SurfaceRightY;
	Hx_lowerlimit = SurfaceLowerZ;
	Hx_upperlimit = SurfaceUpperZ;
	// Hy
	Hy_backlimit = SurfaceBackX;
	Hy_frontlimit = SurfaceFrontX;
	Hy_leftlimit = SurfaceLeftY;
	Hy_rightlimit = SurfaceRightY+1;
	Hy_lowerlimit = SurfaceLowerZ;
	Hy_upperlimit = SurfaceUpperZ;
	// Hz
	Hz_backlimit = SurfaceBackX;
	Hz_frontlimit = SurfaceFrontX;
	Hz_leftlimit = SurfaceLeftY;
	Hz_rightlimit = SurfaceRightY;
	Hz_lowerlimit = SurfaceLowerZ;
	Hz_upperlimit = SurfaceUpperZ+1;

	// the limits of the indices of the field components in the current node
	// Ex
	Ex_backlimit_in_node = max(iback,Ex_backlimit);
	Ex_frontlimit_in_node = min(ifront,Ex_frontlimit);
	Ex_leftlimit_in_node = max(jleft,Ex_leftlimit);
	Ex_rightlimit_in_node = min(jright,Ex_rightlimit);
	Ex_lowerlimit_in_node = max(klower,Ex_lowerlimit);
	Ex_upperlimit_in_node = min(kupper,Ex_upperlimit);
	// Ey
	Ey_backlimit_in_node = max(iback,Ey_backlimit);
	Ey_frontlimit_in_node = min(ifront,Ey_frontlimit);
	Ey_leftlimit_in_node = max(jleft,Ey_leftlimit);
	Ey_rightlimit_in_node = min(jright,Ey_rightlimit);
	Ey_lowerlimit_in_node = max(klower,Ey_lowerlimit);
	Ey_upperlimit_in_node = min(kupper,Ey_upperlimit);
	// Ez
	Ez_backlimit_in_node = max(iback,Ez_backlimit);
	Ez_frontlimit_in_node = min(ifront,Ez_frontlimit);
	Ez_leftlimit_in_node = max(jleft,Ez_leftlimit);
	Ez_rightlimit_in_node = min(jright,Ez_rightlimit);
	Ez_lowerlimit_in_node = max(klower,Ez_lowerlimit);
	Ez_upperlimit_in_node = min(kupper,Ez_upperlimit);
	// Hx
	Hx_backlimit_in_node = max(iback,Hx_backlimit);
	Hx_frontlimit_in_node = min(ifront,Hx_frontlimit);
	Hx_leftlimit_in_node = max(jleft,Hx_leftlimit);
	Hx_rightlimit_in_node = min(jright,Hx_rightlimit);
	Hx_lowerlimit_in_node = max(klower,Hx_lowerlimit);
	Hx_upperlimit_in_node = min(kupper,Hx_upperlimit);
	// Hy
	Hy_backlimit_in_node = max(iback,Hy_backlimit);
	Hy_frontlimit_in_node = min(ifront,Hy_frontlimit);
	Hy_leftlimit_in_node = max(jleft,Hy_leftlimit);
	Hy_rightlimit_in_node = min(jright,Hy_rightlimit);
	Hy_lowerlimit_in_node = max(klower,Hy_lowerlimit);
	Hy_upperlimit_in_node = min(kupper,Hy_upperlimit);
	// Hz
	Hz_backlimit_in_node = max(iback,Hz_backlimit);
	Hz_frontlimit_in_node = min(ifront,Hz_frontlimit);
	Hz_leftlimit_in_node = max(jleft,Hz_leftlimit);
	Hz_rightlimit_in_node = min(jright,Hz_rightlimit);
	Hz_lowerlimit_in_node = max(klower,Hz_lowerlimit);
	Hz_upperlimit_in_node = min(kupper,Hz_upperlimit);

// 	// determine if the faces of the NFFFT box pass through the node
// 	// In the following, "belonging to a cell" is defined as usual as being at the back/left/lower side of a cell (either face-wise or side-wise).
// 	// E-shell is defined to be "passing through the node" when at least one of the E-field components on the box "belong" to one of the cells in the node.
// 	// Deciding whether an E-field component "belongs" to a cell inside a node is done by the following check:
// 	//		(nodestart<=E-field-index)&&(nodeend>=E-field-index)
// 	// where (nodestart,nodeend) are the lower and upper limits of the cell indices in the node (like jleft,jright)
// 	// and E-field-index is the index of the E-field component in the component array (e.g.from 1 to NCELLS+2*NPML+1)
// 	// Similar to above, deciding whether either one of a RANGE OF E-field components "belongs" to a cell inside a node is done by the following check:
// 	//		(nodestart<=upper-limit-of-E-field-indices)&&(nodeend>=lower-limit-of-E-field-indices)
// 	// where (upper,lower)-limit-of-E-field-indices are the upper and lower limits of the E-field-indices in the range.
// 	Lower_E_Shell_InNode = ((klower<=SurfaceLowerZ)&&(kupper>=SurfaceLowerZ)&&(jleft<=Shell_rightlimit)&&(jright>=Shell_leftlimit));
// 	Upper_E_Shell_InNode = ((klower<=SurfaceUpperZ+1)&&(kupper>=SurfaceUpperZ+1)&&(jleft<=Shell_rightlimit)&&(jright>=Shell_leftlimit));
// 	Left_E_Shell_InNode = ((jleft<=SurfaceLeftY)&&(jright>=SurfaceLeftY)&&(klower<=Shell_upperlimit)&&(kupper>=Shell_lowerlimit));
// 	Right_E_Shell_InNode = ((jleft<=SurfaceRightY+1)&&(jright>=SurfaceRightY+1)&&(klower<=Shell_upperlimit)&&(kupper>=Shell_lowerlimit));
// 	// H-shell is defined to be "passing through the node" when at least one of the H-field components outside the box "belong" to one of the cells in the node.
// 	Lower_H_Shell_InNode = ((klower<=SurfaceLowerZ-1)&&(kupper>=SurfaceLowerZ-1)&&(jleft<=Shell_rightlimit)&&(jright>=Shell_leftlimit));
// 	Upper_H_Shell_InNode = ((klower<=SurfaceUpperZ+1)&&(kupper>=SurfaceUpperZ+1)&&(jleft<=Shell_rightlimit)&&(jright>=Shell_leftlimit));
// 	Left_H_Shell_InNode = ((jleft<=SurfaceLeftY-1)&&(jright>=SurfaceLeftY-1)&&(klower<=Shell_upperlimit)&&(kupper>=Shell_lowerlimit));
// 	Right_H_Shell_InNode = ((jleft<=SurfaceRightY+1)&&(jright>=SurfaceRightY+1)&&(klower<=Shell_upperlimit)&&(kupper>=Shell_lowerlimit));

	//indices of the field components on parts of the NFFFT surface in the direction perpendicular to that part (e.g. the y index on the left surface)
	Back_E_shell_index = SurfaceBackX;
	Back_H_shell_index = SurfaceBackX-1;
	Front_E_shell_index = SurfaceFrontX+1;
	Front_H_shell_index = SurfaceFrontX+1;
	Left_E_shell_index = SurfaceLeftY;
	Left_H_shell_index = SurfaceLeftY-1;
	Right_E_shell_index = SurfaceRightY+1;
	Right_H_shell_index = SurfaceRightY+1;
	Lower_E_shell_index = SurfaceLowerZ;
	Lower_H_shell_index = SurfaceLowerZ-1;
	Upper_E_shell_index = SurfaceUpperZ+1;
	Upper_H_shell_index = SurfaceUpperZ+1;

	Back_Ey_PassesThroughNode = ((iback<=Back_E_shell_index)&&(ifront>=Back_E_shell_index)&&(Ey_leftlimit_in_node<=Ey_rightlimit_in_node)&&(Ey_lowerlimit_in_node<=Ey_upperlimit_in_node));
	Back_Ez_PassesThroughNode = ((iback<=Back_E_shell_index)&&(ifront>=Back_E_shell_index)&&(Ez_leftlimit_in_node<=Ez_rightlimit_in_node)&&(Ez_lowerlimit_in_node<=Ez_upperlimit_in_node));
	Back_Hy_PassesThroughNode = ((iback<=Back_H_shell_index)&&(ifront>=Back_H_shell_index)&&(Hy_leftlimit_in_node<=Hy_rightlimit_in_node)&&(Hy_lowerlimit_in_node<=Hy_upperlimit_in_node));
	Back_Hz_PassesThroughNode = ((iback<=Back_H_shell_index)&&(ifront>=Back_H_shell_index)&&(Hz_leftlimit_in_node<=Hz_rightlimit_in_node)&&(Hz_lowerlimit_in_node<=Hz_upperlimit_in_node));
	Front_Ey_PassesThroughNode = ((iback<=Front_E_shell_index)&&(ifront>=Front_E_shell_index)&&(Ey_leftlimit_in_node<=Ey_rightlimit_in_node)&&(Ey_lowerlimit_in_node<=Ey_upperlimit_in_node));
	Front_Ez_PassesThroughNode = ((iback<=Front_E_shell_index)&&(ifront>=Front_E_shell_index)&&(Ez_leftlimit_in_node<=Ez_rightlimit_in_node)&&(Ez_lowerlimit_in_node<=Ez_upperlimit_in_node));
	Front_Hy_PassesThroughNode = ((iback<=Front_H_shell_index)&&(ifront>=Front_H_shell_index)&&(Hy_leftlimit_in_node<=Hy_rightlimit_in_node)&&(Hy_lowerlimit_in_node<=Hy_upperlimit_in_node));
	Front_Hz_PassesThroughNode = ((iback<=Front_H_shell_index)&&(ifront>=Front_H_shell_index)&&(Hz_leftlimit_in_node<=Hz_rightlimit_in_node)&&(Hz_lowerlimit_in_node<=Hz_upperlimit_in_node));
	Left_Ex_PassesThroughNode = ((jleft<=Left_E_shell_index)&&(jright>=Left_E_shell_index)&&(Ex_backlimit_in_node<=Ex_frontlimit_in_node)&&(Ex_lowerlimit_in_node<=Ex_upperlimit_in_node));
	Left_Ez_PassesThroughNode = ((jleft<=Left_E_shell_index)&&(jright>=Left_E_shell_index)&&(Ez_backlimit_in_node<=Ez_frontlimit_in_node)&&(Ez_lowerlimit_in_node<=Ez_upperlimit_in_node));
	Left_Hx_PassesThroughNode = ((jleft<=Left_H_shell_index)&&(jright>=Left_H_shell_index)&&(Hx_backlimit_in_node<=Hx_frontlimit_in_node)&&(Hx_lowerlimit_in_node<=Hx_upperlimit_in_node));
	Left_Hz_PassesThroughNode = ((jleft<=Left_H_shell_index)&&(jright>=Left_H_shell_index)&&(Hz_backlimit_in_node<=Hz_frontlimit_in_node)&&(Hz_lowerlimit_in_node<=Hz_upperlimit_in_node));
	Right_Ex_PassesThroughNode = ((jleft<=Right_E_shell_index)&&(jright>=Right_E_shell_index)&&(Ex_backlimit_in_node<=Ex_frontlimit_in_node)&&(Ex_lowerlimit_in_node<=Ex_upperlimit_in_node));
	Right_Ez_PassesThroughNode = ((jleft<=Right_E_shell_index)&&(jright>=Right_E_shell_index)&&(Ez_backlimit_in_node<=Ez_frontlimit_in_node)&&(Ez_lowerlimit_in_node<=Ez_upperlimit_in_node));
	Right_Hx_PassesThroughNode = ((jleft<=Right_H_shell_index)&&(jright>=Right_H_shell_index)&&(Hx_backlimit_in_node<=Hx_frontlimit_in_node)&&(Hx_lowerlimit_in_node<=Hx_upperlimit_in_node));
	Right_Hz_PassesThroughNode = ((jleft<=Right_H_shell_index)&&(jright>=Right_H_shell_index)&&(Hz_backlimit_in_node<=Hz_frontlimit_in_node)&&(Hz_lowerlimit_in_node<=Hz_upperlimit_in_node));
	Lower_Ex_PassesThroughNode = ((klower<=Lower_E_shell_index)&&(kupper>=Lower_E_shell_index)&&(Ex_backlimit_in_node<=Ex_frontlimit_in_node)&&(Ex_leftlimit_in_node<=Ex_rightlimit_in_node));
	Lower_Ey_PassesThroughNode = ((klower<=Lower_E_shell_index)&&(kupper>=Lower_E_shell_index)&&(Ey_backlimit_in_node<=Ey_frontlimit_in_node)&&(Ey_leftlimit_in_node<=Ey_rightlimit_in_node));
	Lower_Hx_PassesThroughNode = ((klower<=Lower_H_shell_index)&&(kupper>=Lower_H_shell_index)&&(Hx_backlimit_in_node<=Hx_frontlimit_in_node)&&(Hx_leftlimit_in_node<=Hx_rightlimit_in_node));
	Lower_Hy_PassesThroughNode = ((klower<=Lower_H_shell_index)&&(kupper>=Lower_H_shell_index)&&(Hy_backlimit_in_node<=Hy_frontlimit_in_node)&&(Hy_leftlimit_in_node<=Hy_rightlimit_in_node));
	Upper_Ex_PassesThroughNode = ((klower<=Upper_E_shell_index)&&(kupper>=Upper_E_shell_index)&&(Ex_backlimit_in_node<=Ex_frontlimit_in_node)&&(Ex_leftlimit_in_node<=Ex_rightlimit_in_node));
	Upper_Ey_PassesThroughNode = ((klower<=Upper_E_shell_index)&&(kupper>=Upper_E_shell_index)&&(Ey_backlimit_in_node<=Ey_frontlimit_in_node)&&(Ey_leftlimit_in_node<=Ey_rightlimit_in_node));
	Upper_Hx_PassesThroughNode = ((klower<=Upper_H_shell_index)&&(kupper>=Upper_H_shell_index)&&(Hx_backlimit_in_node<=Hx_frontlimit_in_node)&&(Hx_leftlimit_in_node<=Hx_rightlimit_in_node));
	Upper_Hy_PassesThroughNode = ((klower<=Upper_H_shell_index)&&(kupper>=Upper_H_shell_index)&&(Hy_backlimit_in_node<=Hy_frontlimit_in_node)&&(Hy_leftlimit_in_node<=Hy_rightlimit_in_node));


	//allocate and fill the phasor storage arrays
	//Back and front surfaces
	//Jz is created by Hy
	if (Back_Hy_PassesThroughNode)
	{
		J_z_back.resize(Range(0,L-1),
						Range(Hy_leftlimit_in_node,Hy_rightlimit_in_node),
						Range(Hy_lowerlimit_in_node,Hy_upperlimit_in_node));
		J_z_back=0;
	}
	if (Front_Hy_PassesThroughNode)
	{
		J_z_front.resize(Range(0,L-1),
						 Range(Hy_leftlimit_in_node,Hy_rightlimit_in_node),
						 Range(Hy_lowerlimit_in_node,Hy_upperlimit_in_node));
		J_z_front=0;
	}
	//My is created by Ez
	if (Back_Ez_PassesThroughNode)
	{
		M_y_back.resize(Range(0,L-1),
						Range(Ez_leftlimit_in_node,Ez_rightlimit_in_node),
						Range(Ez_lowerlimit_in_node,Ez_upperlimit_in_node));
		M_y_back=0;
	}
	if (Front_Ez_PassesThroughNode)
	{
		M_y_front.resize(Range(0,L-1),
						Range(Ez_leftlimit_in_node,Ez_rightlimit_in_node),
						Range(Ez_lowerlimit_in_node,Ez_upperlimit_in_node));
		M_y_front=0;
	}
	//Jy is created by Hz
	if (Back_Hz_PassesThroughNode)
	{
		J_y_back.resize(Range(0,L-1),
						Range(Hz_leftlimit_in_node,Hz_rightlimit_in_node),
						Range(Hz_lowerlimit_in_node,Hz_upperlimit_in_node));
		J_y_back=0;
	}
	if (Front_Hz_PassesThroughNode)
	{
		J_y_front.resize(Range(0,L-1),
						Range(Hz_leftlimit_in_node,Hz_rightlimit_in_node),
						Range(Hz_lowerlimit_in_node,Hz_upperlimit_in_node));
		J_y_front=0;
	}
	//Mz is created by Ey
	if (Back_Ey_PassesThroughNode)
	{
		M_z_back.resize(Range(0,L-1),
						Range(Ey_leftlimit_in_node,Ey_rightlimit_in_node),
						Range(Ey_lowerlimit_in_node,Ey_upperlimit_in_node));
		M_z_back=0;
	}
	if (Front_Ey_PassesThroughNode)
	{
		M_z_front.resize(Range(0,L-1),
						Range(Ey_leftlimit_in_node,Ey_rightlimit_in_node),
						Range(Ey_lowerlimit_in_node,Ey_upperlimit_in_node));
		M_z_front=0;
	}

	//Left and right surfaces
	//Jx is created by Hz
	if (Left_Hz_PassesThroughNode)
	{
		J_x_left.resize(Range(0,L-1),
						Range(Hz_backlimit_in_node,Hz_frontlimit_in_node),
						Range(Hz_lowerlimit_in_node,Hz_upperlimit_in_node));
		J_x_left=0;
	}
	if (Right_Hz_PassesThroughNode)
	{
		J_x_right.resize(Range(0,L-1),
						Range(Hz_backlimit_in_node,Hz_frontlimit_in_node),
						Range(Hz_lowerlimit_in_node,Hz_upperlimit_in_node));
		J_x_right=0;
	}
	//Mz is created by Ex
	if (Left_Ex_PassesThroughNode)
	{
		M_z_left.resize(Range(0,L-1),
						Range(Ex_backlimit_in_node,Ex_frontlimit_in_node),
						Range(Ex_lowerlimit_in_node,Ex_upperlimit_in_node));
		M_z_left=0;
	}
	if (Right_Ex_PassesThroughNode)
	{
		M_z_right.resize(Range(0,L-1),
						Range(Ex_backlimit_in_node,Ex_frontlimit_in_node),
						Range(Ex_lowerlimit_in_node,Ex_upperlimit_in_node));
		M_z_right=0;
	}
	//Jz is created by Hx
	if (Left_Hx_PassesThroughNode)
	{
		J_z_left.resize(Range(0,L-1),
						Range(Hx_backlimit_in_node,Hx_frontlimit_in_node),
						Range(Hx_lowerlimit_in_node,Hx_upperlimit_in_node));
		J_z_left=0;
	}
	if (Right_Hx_PassesThroughNode)
	{
		J_z_right.resize(Range(0,L-1),
						Range(Hx_backlimit_in_node,Hx_frontlimit_in_node),
						Range(Hx_lowerlimit_in_node,Hx_upperlimit_in_node));
		J_z_right=0;
	}
	//Mx is created by Ez
	if (Left_Ez_PassesThroughNode)
	{
		M_x_left.resize(Range(0,L-1),
						Range(Ez_backlimit_in_node,Ez_frontlimit_in_node),
						Range(Ez_lowerlimit_in_node,Ez_upperlimit_in_node));
		M_x_left=0;
	}
	if (Right_Ez_PassesThroughNode)
	{
		M_x_right.resize(Range(0,L-1),
						Range(Ez_backlimit_in_node,Ez_frontlimit_in_node),
						Range(Ez_lowerlimit_in_node,Ez_upperlimit_in_node));
		M_x_right=0;
	}

	//Lower and upper surfaces
	//Jx is created by Hy
	if (Lower_Hy_PassesThroughNode)
	{
		J_x_lower.resize(Range(0,L-1),
						 Range(Hy_backlimit_in_node,Hy_frontlimit_in_node),
						 Range(Hy_leftlimit_in_node,Hy_rightlimit_in_node));
		J_x_lower=0;
	}
	if (Upper_Hy_PassesThroughNode)
	{
		J_x_upper.resize(Range(0,L-1),
						 Range(Hy_backlimit_in_node,Hy_frontlimit_in_node),
						 Range(Hy_leftlimit_in_node,Hy_rightlimit_in_node));
		J_x_upper=0;
	}
	//My is created by Ex
	if (Lower_Ex_PassesThroughNode)
	{
		M_y_lower.resize(Range(0,L-1),
						 Range(Ex_backlimit_in_node,Ex_frontlimit_in_node),
						 Range(Ex_leftlimit_in_node,Ex_rightlimit_in_node));
		M_y_lower=0;
	}
	if (Upper_Ex_PassesThroughNode)
	{
		M_y_upper.resize(Range(0,L-1),
						 Range(Ex_backlimit_in_node,Ex_frontlimit_in_node),
						 Range(Ex_leftlimit_in_node,Ex_rightlimit_in_node));
		M_y_upper=0;
	}
	//Jy is created by Hx
	if (Lower_Hx_PassesThroughNode)
	{
		J_y_lower.resize(Range(0,L-1),
						 Range(Hx_backlimit_in_node,Hx_frontlimit_in_node),
						 Range(Hx_leftlimit_in_node,Hx_rightlimit_in_node));
		J_y_lower=0;
	}
	if (Upper_Hx_PassesThroughNode)
	{
		J_y_upper.resize(Range(0,L-1),
						 Range(Hx_backlimit_in_node,Hx_frontlimit_in_node),
						 Range(Hx_leftlimit_in_node,Hx_rightlimit_in_node));
		J_y_upper=0;
	}
	//Mx is created by Ey
	if (Lower_Ey_PassesThroughNode)
	{
		M_x_lower.resize(Range(0,L-1),
						 Range(Ey_backlimit_in_node,Ey_frontlimit_in_node),
						 Range(Ey_leftlimit_in_node,Ey_rightlimit_in_node));
		M_x_lower=0;
	}
	if (Upper_Ey_PassesThroughNode)
	{
		M_x_upper.resize(Range(0,L-1),
						 Range(Ey_backlimit_in_node,Ey_frontlimit_in_node),
						 Range(Ey_leftlimit_in_node,Ey_rightlimit_in_node));
		M_x_upper=0;
	}


	//Get the propagation directions and associated delays (from the origin) of the scattered PWs from the Ctfsf object
	if (Data.TFSFPtr!=NULL)
	{//don't dereference TFSFPtr if it is NULL
	 // PW_THETA and PW_PHI are left uninitialized, with size 0
		Data.TFSFPtr->WriteScatteredPWDirections(PW_THETA,PW_PHI);
		Data.TFSFPtr->WriteScatteredPWDelaysFromOrigin(origindelay_array,FarFieldOriginX,FarFieldOriginY,FarFieldOriginZ);	//get the delay from the transformer origin
		//allocate the field phasor arrays
		E_x_phasor_array.resize(L,PW_THETA.size());		//could have also used "PW_PHI.size()", since they are supposed to be of the same size (num. of PWs)
		E_y_phasor_array.resize(L,PW_THETA.size());		//could have also used "PW_PHI.size()", since they are supposed to be of the same size (num. of PWs)
	}
}

bool Ctr_pd::observation_angle_exists(const int& l, const int& d1, const int& d2)
{//boolean value that indicates whether the observation angle exists and falls within the collection NA
    if (Data.directionspec=="theta-phi")
    {//dir1 represents THETA
//cout << l << "," << abs(lambda(l)/min_lambda) << endl;
//cout << l << endl;
    	if (Data.scale_with_wavelength)
    	{//direction cosines are scaled with lambda
            return lessThanOrEqual(abs(abs(lambda(l)/min_lambda)*sin(dir1(d1)*M_PI/180)),min(1.0,max_s),LIBSTD_DBL_EPSILON);
        }
        else
        {
        	//sin(theta) is always less than or equal to 1
        	return true;
//            return (abs(
//						sin(dir1(d1)*M_PI/180)
//					   )<min(1.0,max_s));
        }
    }
    else
    {//dir1 represents the x-direction cosine = sin(THETA)*cos(PHI)
     //dir2 represents the y-direction cosine = sin(THETA)*sin(PHI)
	//sqrt(dir1^2+dir2^2) = |sin(theta)| should be less than or equal to 1 for a valid sin(theta)
    // it should also be less than max_s (the maximum allowable direction cosine)
    	if (Data.scale_with_wavelength)
    	{//direction cosines are scaled with lambda
    	    return lessThanOrEqual(abs(lambda(l)/min_lambda)*sqrt(pow2(dir1(d1))+pow2(dir2(d2))),min(1.0,max_s),LIBSTD_DBL_EPSILON);
    	}
    	else
    	{
    	    return lessThanOrEqual(sqrt(pow2(dir1(d1))+pow2(dir2(d2))),min(1.0,max_s),LIBSTD_DBL_EPSILON);
    	}
	}
}

double Ctr_pd::THETA(const int& l, const int& d1, const int& d2)
{//spherical theta angle
//For efficiency, this function does NOT check if the observation angle actually exists
//If it does not, then the result can be anything, including NaN or inf
//It is the caller's responsibility to take this possibility into account
    if (Data.directionspec=="theta-phi")
    {//dir1 represents THETA
    	if (Data.scale_with_wavelength)
        {//return the new theta with scaled direction cosines
        	sx = abs(lambda(l)/min_lambda)*sin(dir1(d1)*M_PI/180)*cos(dir2(d2)*M_PI/180);
        	sy = abs(lambda(l)/min_lambda)*sin(dir1(d1)*M_PI/180)*sin(dir2(d2)*M_PI/180);
        	return THETA(sx,sy,((cos(dir1(d1)*M_PI/180)>0)?sqrt(1-pow2(sx)-pow2(sy)):-sqrt(1-pow2(sx)-pow2(sy))));
		}
		else
        {
            return dir1(d1)*M_PI/180;	//dir1 represents THETA  [in degrees, so convert to radian]
        }
    }
    else
    {//dir1 represents the x-direction cosine = sin(THETA)*cos(PHI)
     //dir2 represents the y-direction cosine = sin(THETA)*sin(PHI)
        if (Data.directionspec=="dircosx-dircosy-upper")
        {
            if (Data.scale_with_wavelength)
            {//scale the direction cosines with the wavelength, return the new theta
            	return asin(abs(lambda(l)/min_lambda)*sqrt(pow2(dir1(d1))+pow2(dir2(d2))));//between 0 and pi/2
			}
			else
            {
                return asin(sqrt(pow2(dir1(d1))+pow2(dir2(d2))));//between 0 and pi/2
            }
        }
        else if (Data.directionspec=="dircosx-dircosy-lower")
        {
            if (Data.scale_with_wavelength)
            {//scale the direction cosines with the wavelength, return the new theta
            	return (M_PI - asin(abs(lambda(l)/min_lambda)*sqrt(pow2(dir1(d1))+pow2(dir2(d2)))));//between pi/2 and pi
			}
			else
            {
                return (M_PI - asin(sqrt(pow2(dir1(d1))+pow2(dir2(d2)))));//between pi/2 and pi
            }
        }
    }
}

double Ctr_pd::PHI(const int& l, const int& d1, const int& d2)
{//spherical phi angle
//For efficiency, this function does NOT check if the observation angle actually exists
//If it does not, then the result can be anything, including NaN or inf
//It is the caller's responsibility to take this possibility into account
    double phi;

    if (Data.directionspec=="theta-phi")
    {
    	if (Data.scale_with_wavelength)
        {//return the new theta with scaled direction cosines
		 //(although this is unnecessary for z-centered spherical coordinates, it will be relevant for other orientations in the future)
        	sx = abs(lambda(l)/min_lambda)*sin(dir1(d1))*cos(dir2(d2));
        	sy = abs(lambda(l)/min_lambda)*sin(dir1(d1))*sin(dir2(d2));
        	return PHI(sx,sy);
        }
        else
        {
        	return dir2(d2)*M_PI/180;		//dir2 represents PHI  [in degrees, so convert to radian]
        }
    }
    else
    {//dir1 represents the x-direction cosine = sin(THETA)*cos(PHI)
     //dir2 represents the y-direction cosine = sin(THETA)*sin(PHI)
        if ((abs(dir1(d1))<LIBSTD_DBL_EPSILON*100)&&(abs(dir2(d2))<LIBSTD_DBL_EPSILON*100)) //is THETA=0?
        {
          return 0; //PHI is arbitrary, so assign PHI=0 (TODO:document this!!!) [theta^ unit angle = x^; phi^ unit angle = y^]
        }
        else
        {
            phi = atan2(dir2(d2),dir1(d1));
            //atan2 is between [-pi,pi]
            //if less than 0, make PHI positive by adding 2pi
            return ((phi<0)?(phi+2*M_PI):phi);
        }
    }
}

double Ctr_pd::SinT(const int& l, const int& d1, const int& d2)
{
    return sin(THETA(l,d1,d2));
}

double Ctr_pd::CosT(const int& l, const int& d1, const int& d2)
{
    return cos(THETA(l,d1,d2));
}

double Ctr_pd::SinP(const int& l, const int& d1, const int& d2)
{
    return sin(PHI(l,d1,d2));
}

double Ctr_pd::CosP(const int& l, const int& d1, const int& d2)
{
    return cos(PHI(l,d1,d2));
}

double Ctr_pd::SinTCosP(const int& l, const int& d1, const int& d2)
{
    return SinT(l,d1,d2)*CosP(l,d1,d2);
}

double Ctr_pd::SinTSinP(const int& l, const int& d1, const int& d2)
{
    return SinT(l,d1,d2)*SinP(l,d1,d2);
}

double Ctr_pd::CosTSinP(const int& l, const int& d1, const int& d2)
{
    return CosT(l,d1,d2)*SinP(l,d1,d2);
}

double Ctr_pd::CosTCosP(const int& l, const int& d1, const int& d2)
{
    return CosT(l,d1,d2)*CosP(l,d1,d2);
}

double Ctr_pd::GridVelocity(const double& w, const double& theta, const double& phi)
{//calculates the grid velocity v(theta,phi) at angular frequency w, normalized by c
	if (w!=0)
	{
		double lambda_0 = c_upper/(w/2/M_PI);	//free-space wavelength
		double dx_norm = dx/lambda_0;	//normalized spatial step
		double N_lambda = 1/dx_norm;	//steps per wavelength
		//propagation direction cosines
		double k_x = -cos(phi)*sin(theta);
		double k_y = -sin(phi)*sin(theta);
		double k_z = -cos(theta);
		double A = dx_norm*k_x/2;
		double B = dx_norm*k_y/2;
		double C = dx_norm*k_z/2;
		double S = courant/sqrt(3.0);
		double D = 1/pow(S,2)*pow(sin(M_PI*S/N_lambda),2);

		double k_i = 2*M_PI;	//initial guess for k
		double correction = 1e6;
		double threshold = 2*M_PI/1e10;	//error threshold for convergence

		int iter=1;
		while (abs(correction)>threshold)
		{//see Taflove, eqn. (4.16a)
			correction = -(pow(sin(A*k_i),2)+pow(sin(B*k_i),2)+pow(sin(C*k_i),2)-D)/(A*sin(2*A*k_i)+B*sin(2*B*k_i)+C*sin(2*C*k_i));
			k_i = k_i + correction;
		}
		return 2*M_PI/k_i;	//the grid velocity normalized by c
	}
	else
	{
		return 1;	//if w=0, then grid velocity is approx. c
	}
}
