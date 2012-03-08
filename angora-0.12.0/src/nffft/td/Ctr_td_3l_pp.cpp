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

//Defines the post-processing steps in the TIME_DOMAIN near-field-to-far-field transformer object "Ctr_td_3l" for 3-layered lossless media

#include "headers.h"

#include "Ctr_td_3l.h"

//for the definition of MaterialId
#include "material_id.h"

//definition of Cwf needed
#include "waveforms/Cwf.h"

//definition of Cpointsources needed
#include "pointsources/Cpointsources.h"

//routines for getting/setting the initial time value in the simulation
#include "time_axis.h"

extern double dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;

extern Array<ElectricMaterialIndexType_Z,1> Layering_e_z;

extern Array<double,1> eps_z;

extern double c_o;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif
extern int rank;

extern double DiffGaussian(double t, double tau);
extern double DoubleDiffGaussian(double t, double tau);
extern double TripleDiffGaussian(double t, double tau);


void Ctr_td_3l::ConstructFarField()
{
	Theta_J = dx*dx/(2*M_PI*c_o)*(-Electric_X_e*CosP - Electric_Y_e*SinP + Z_o*Electric_Z_e*SinT);
	Theta_M = dx*dx/(2*M_PI*c_o)*(Magnetic_X_e*SinP - Magnetic_Y_e*CosP);
	Phi_J = dx*dx/(2*M_PI*c_o)*(Electric_X_h*SinP*CosT - Electric_Y_h*CosP*CosT);
	Phi_M = dx*dx/(2*M_PI*c_o)*(Magnetic_X_h*CosP*CosT + Magnetic_Y_h*SinP*CosT - (1/Z_o)*Magnetic_Z_h*SinT*CosT);

	//It can be shown that for THETA>90, the preceding formula for the far field (Smith, Table 3.4) changes sign!
	if (CosT<0)
	{
		Theta_J *= -1;
		Theta_M *= -1;
		Phi_J *= -1;
		Phi_M *= -1;
	}

	//Now differentiate, bring both J and M forward half-step
	Range T(0,TotalSteps-2);
	Array<double,1> temparray(T);	//needed for safety, since array traversal order not known
	temparray = (Theta_J(T+1)-Theta_J(T))/dt;	//Two-step assignment for differentiation - step 1
	Theta_J(T) = temparray;						// - step 2
	Theta_J(TotalSteps-1)=0;	//zero the last element - will not be used
	temparray = (Phi_J(T+1)-Phi_J(T))/dt;		//Two-step assignment for differentiation - step 1
	Phi_J(T) = temparray;						// - step 2
	Phi_J(TotalSteps-1)=0;		//zero the last element - will not be used
	temparray = (Phi_M(T+1)-Phi_M(T))/dt;	//Two-step assignment for differentiation - step 1
	Phi_M(T) = temparray;						// - step 2
	Phi_M(TotalSteps-1)=0;	//zero the last element - will not be used
	temparray = (Theta_M(T+1)-Theta_M(T))/dt;		//Two-step assignment for differentiation - step 1
	Theta_M(T) = temparray;						// - step 2
	Theta_M(TotalSteps-1)=0;		//zero the last element - will not be used
	//Take into account the "0.5dt" lag of M with respect to J, bring M forward half-step
	T = Range(0,TotalSteps-3);
	temparray.resize(T);	//again needed for safety, since array traversal order not known
	temparray = (Phi_M(T+1)+Phi_M(T))/2;	//Two-step assignment for interpolation - step 1
	Phi_M(T+1) = temparray;					// - step 2
	Phi_M(0)=0;	//zero the first element - will not be used
	temparray = (Theta_M(T+1)+Theta_M(T))/2;		//Two-step assignment for interpolation - step 1
	Theta_M(T+1) = temparray;						// - step 2
	Theta_M(0)=0;		//zero the first element - will not be used
	//Align J with M:
	temparray = Theta_J(T);
	Theta_J(T+1) = temparray;
	Theta_J(0)=0;	//zero the first element - will not be used
	temparray = Phi_J(T);
	Phi_J(T+1) = temparray;
	Phi_J(0)=0;		//zero the first element - will not be used

	//Now, M has been brought forward ONE step, J is brought forward HALF step, and they coincide in time.
	//Both arrays are aligned with the initial M-waveform. By padding necessary zeros, no extra delaying
	//is introduced; therefore the total delay is still (r+r_offset)/c_o.
	// Before differentiation:
	//	J->		____x_______x_______x_______x___________
	//	M->		*_______*_______*_______*_______________>time
	// After differentiation:
	//	J->		________x_______x_______x_______0_______
	//	M->		____*_______*_______*_______0___________>time
	// After interpolation:
	//	J->		0_______x_______x_______0_______________
	//	M->		0_______*_______*_______0_______________>time

	//Now, they may be combined to yield the far-field waveforms:
	waveformTheta = Theta_J+Theta_M;
	waveformPhi = Phi_J+Phi_M;


	//Take the sum (MPI_SUM) of all partial far-field waveforms at each node and store in node 0
#ifndef MPI_DISABLE
	if (rank==0)
	{
		MPI_Reduce(MPI_IN_PLACE,waveformTheta.data(),waveformTheta.size(),MPI_DOUBLE,MPI_SUM,0,MPI_CartSubComm);
		MPI_Reduce(MPI_IN_PLACE,waveformPhi.data(),waveformPhi.size(),MPI_DOUBLE,MPI_SUM,0,MPI_CartSubComm);
	}
	else
	{
		MPI_Reduce(waveformTheta.data(),NULL,waveformTheta.size(),MPI_DOUBLE,MPI_SUM,0,MPI_CartSubComm);
		MPI_Reduce(waveformPhi.data(),NULL,waveformPhi.size(),MPI_DOUBLE,MPI_SUM,0,MPI_CartSubComm);
	}
#endif
	//if there is no MPI, there is only rank 0, which holds waveformTheta and waveformPhi anyway


	//write out the far field arrays
	if (rank==0)
	{
		const int NumOfWaveforms = 5;
		double Out[NumOfWaveforms];
		FarFieldFile.write((char*)&NumOfWaveforms,sizeof(NumOfWaveforms));
		for (int n=0; n<TotalSteps; n++)
		{
			Out[0] = n*dt+get_initial_time_value()-r_offset/c_o;
			Out[1] = waveformTheta(n);
			Out[2] = waveformPhi(n);
			Out[3] = TheoreticalFarFieldTheta(n);
			Out[4] = TheoreticalFarFieldPhi(n);
			FarFieldFile.write((char*)Out,NumOfWaveforms*sizeof(Out[0]));
		}
	}
}

double Ctr_td_3l::TheoreticalFarFieldTheta(const int& n)
//Theoretical theta-component of the radiated E-field due to a collection of electric dipoles (for which data is given in
// *PointSourcesPtr)
{
	theta_component = 0;
	if (Data.PointSourcesPtr!=NULL)	//calculate far-field if the pointer to the Cpointsource object is defined
	{
		for (int i=0; i<Data.PointSourcesPtr->NumberOfElectricDipoles(); i++)
		{
			Orientation = Data.PointSourcesPtr->ElectricDipole(i)->Orientation;
// 			WaveShape = Data.PointSourcesPtr->ElectricDipole(i)->WaveShape;
			x = Data.PointSourcesPtr->ElectricDipole(i)->x;
			y = Data.PointSourcesPtr->ElectricDipole(i)->y;
			z = Data.PointSourcesPtr->ElectricDipole(i)->z;
			SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX;
			SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY;
// 			delay = Data.PointSourcesPtr->ElectricDipole(i)->delay;
// 			j0 = Data.PointSourcesPtr->ElectricDipole(i)->j0;		//Maximum current moment (Smith 7.3)
// 			tau = Data.PointSourcesPtr->ElectricDipole(i)->tau;

/*			//the effect of differentiating the current moment waveform is absorbed into the maximum current moment j0
			if (WaveShape=="gaussian")
			{
				j0 *= exp(-0.5)/tau;
			}
			else if (WaveShape=="diffgaussian")
			{
				j0 *= exp(0.5)/tau;
			}
			else if (WaveShape=="doublediffgaussian")
			{
				j0 *= 1.3801190461/tau;
			}*/

			//begin calculation of theta_component
			if (Orientation=="x_directed")
			{
				temp=0;
				delay_lateral = (r_offset-dx*((SourceX+0.5)*SinTCosP+(SourceY)*SinTSinP))/c_o/dt;
				for (int impulse=0;impulse<V_e_i(z).size();impulse++)
				{
					delay_total = delay_lateral + V_e_i(z)(impulse).Delay;
// 					temp += V_e_i(z)(impulse).Amp*TheoreticalFarFieldWaveform((n-delay_total)*dt-delay*tau,tau);
					temp += V_e_i(z)(impulse).Amp*Data.PointSourcesPtr->ElectricDipole(i)->current_moment_waveform()->Derivative((n-delay_total)*dt+get_initial_time_value());
				}
// 				theta_component += (j0/(2*M_PI*c_o))*(-CosP)*temp;	//for x-directed dipole
				theta_component += (1/(2*M_PI*c_o))*(-CosP)*temp;	//for x-directed dipole
			}
			else if (Orientation=="y_directed")
			{
				temp=0;
				delay_lateral = (r_offset-dx*((SourceX)*SinTCosP+(SourceY+0.5)*SinTSinP))/c_o/dt;	//Lateral delay
				for (int impulse=0;impulse<V_e_i(z).size();impulse++)
				{
					delay_total = delay_lateral + V_e_i(z)(impulse).Delay;
// 					temp += V_e_i(z)(impulse).Amp*TheoreticalFarFieldWaveform((n-delay_total)*dt-delay*tau,tau);
					temp += V_e_i(z)(impulse).Amp*Data.PointSourcesPtr->ElectricDipole(i)->current_moment_waveform()->Derivative((n-delay_total)*dt+get_initial_time_value());
				}
// 				theta_component += (j0/(2*M_PI*c_o))*(-SinP)*temp;	//for y-directed dipole
				theta_component += (1/(2*M_PI*c_o))*(-SinP)*temp;	//for y-directed dipole
			}
			else if (Orientation=="z_directed")
			{
				temp=0;
				delay_lateral = (r_offset-dx*((SourceX)*SinTCosP+(SourceY)*SinTSinP))/c_o/dt;	//Lateral delay
				for (int impulse=0;impulse<V_e_v(z).size();impulse++)
				{
					delay_total = delay_lateral + V_e_v(z)(impulse).Delay;
// 					temp += V_e_v(z)(impulse).Amp*TheoreticalFarFieldWaveform((n-delay_total)*dt-delay*tau,tau);
					temp += V_e_v(z)(impulse).Amp*Data.PointSourcesPtr->ElectricDipole(i)->current_moment_waveform()->Derivative((n-delay_total)*dt+get_initial_time_value());
				}
				theta_component += (1/(2*M_PI*c_o))*(Z_o*SinT/(eps_z(Layering_e_z(z))/epsilon_o))*temp; //for z-directed dipole
			}
		}
	}

	if (CosT>=0)
		return theta_component;
	else
		return -theta_component;	//see note at the beginning about the sign reversal
}

double Ctr_td_3l::TheoreticalFarFieldPhi(const int& n)
//Theoretical phi-component of the radiated E-field due to a collection of electric dipoles (for which data is given in
// *PointSourcesPtr)
{
	phi_component = 0;

	if (Data.PointSourcesPtr!=NULL)	//calculate far-field if the pointer to the Cpointsource object is defined
	{
		for (int i=0; i<Data.PointSourcesPtr->NumberOfElectricDipoles(); i++)
		{
			Orientation = Data.PointSourcesPtr->ElectricDipole(i)->Orientation;
// 			WaveShape = Data.PointSourcesPtr->ElectricDipole(i)->WaveShape;
			x = Data.PointSourcesPtr->ElectricDipole(i)->x;
			y = Data.PointSourcesPtr->ElectricDipole(i)->y;
			z = Data.PointSourcesPtr->ElectricDipole(i)->z;
			SourceX = Data.PointSourcesPtr->ElectricDipole(i)->x-FarFieldOriginX;
			SourceY = Data.PointSourcesPtr->ElectricDipole(i)->y-FarFieldOriginY;
// 			delay = Data.PointSourcesPtr->ElectricDipole(i)->delay;
// 			j0 = Data.PointSourcesPtr->ElectricDipole(i)->j0;		//Maximum current moment (Smith 7.3)
// 			tau = Data.PointSourcesPtr->ElectricDipole(i)->tau;

// 			//the effect of differentiating the current moment waveform is absorbed into the maximum current moment j0
// 			if (WaveShape=="gaussian")
// 			{
// 				j0 *= exp(-0.5)/tau;
// 			}
// 			else if (WaveShape=="diffgaussian")
// 			{
// 				j0 *= exp(0.5)/tau;
// 			}
// 			else if (WaveShape=="doublediffgaussian")
// 			{
// 				j0 *= 1.3801190461/tau;
// 			}

			//begin calculation of phi_component
			if (Orientation=="x_directed")
			{
				temp=0;
				delay_lateral = (r_offset-dx*((SourceX+0.5)*SinTCosP+(SourceY)*SinTSinP))/c_o/dt;
				for (int impulse=0;impulse<V_h_i(z).size();impulse++)
				{
					delay_total = delay_lateral + V_h_i(z)(impulse).Delay;
// 					temp += V_h_i(z)(impulse).Amp*TheoreticalFarFieldWaveform((n-delay_total)*dt-delay*tau,tau);
					temp += V_h_i(z)(impulse).Amp*Data.PointSourcesPtr->ElectricDipole(i)->current_moment_waveform()->Derivative((n-delay_total)*dt+get_initial_time_value());
				}
// 				phi_component += (j0/(2*M_PI*c_o))*(SinP*CosT)*temp;	//for x-directed dipole
				phi_component += (1/(2*M_PI*c_o))*(SinP*CosT)*temp;	//for x-directed dipole
			}
			else if (Orientation=="y_directed")
			{
				temp=0;
				delay_lateral = (r_offset-dx*((SourceX)*SinTCosP+(SourceY+0.5)*SinTSinP))/c_o/dt;	//Lateral delay
				for (int impulse=0;impulse<V_h_i(z).size();impulse++)
				{
					delay_total = delay_lateral + V_h_i(z)(impulse).Delay;
// 					temp += V_h_i(z)(impulse).Amp*TheoreticalFarFieldWaveform((n-delay_total)*dt-delay*tau,tau);
					temp += V_h_i(z)(impulse).Amp*Data.PointSourcesPtr->ElectricDipole(i)->current_moment_waveform()->Derivative((n-delay_total)*dt+get_initial_time_value());
				}
// 				phi_component += (j0/(2*M_PI*c_o))*(-CosP*CosT)*temp;	//for y-directed dipole
				phi_component += (1/(2*M_PI*c_o))*(-CosP*CosT)*temp;	//for y-directed dipole
			}
			else if (Orientation=="z_directed")
			{
				//z-component does not contribute to the phi-component
			}
		}
	}

	if (CosT>=0)
		return phi_component;
	else
		return -phi_component;	//see note at the beginning about the sign reversal
}

// double Ctr_td_3l::TheoreticalFarFieldWaveform(double t, double tau)
// {//Return the appropriate far-field waveform shape (normalized), depending on the source current waveform
// 	if (WaveShape=="gaussian")
// 	{
// 		return DiffGaussian(t,tau);
// 	}
// 	else if (WaveShape=="diffgaussian")
// 	{
// 		return DoubleDiffGaussian(t,tau);
// 	}
// 	else if (WaveShape=="doublediffgaussian")
// 	{
// 		return TripleDiffGaussian(t,tau);
// 	}
// 	else return 0;
// }
