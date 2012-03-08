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
//This file includes the definitions of post-processing routines that are common to all phasor-domain transformers.

#include "headers.h"

#include "Ctr_pd.h"

// #define WRITE_VECTOR_POTENTIALS
#define WRITE_THEORETICAL_FIELD

//include the HDF5 header, if not disabled
#ifndef HDF5_DISABLE
#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif
#endif

#ifndef MPI_DISABLE
extern MPI_Comm MPI_CartSubComm;
#endif

extern int angora_version_major,angora_version_minor,angora_version_revision;

extern int rank;

extern void MPI_exit(const int& exitcode);


void Ctr_pd::GatherAndWriteFarField()
{
	// Finally, gather the far-field arrays from the nodes and write into file
#ifndef MPI_DISABLE
#ifdef USE_MPI_IN_PLACE
	//Take the sum (MPI_SUM) of all partial far-field waveforms at each node and store in node 0
	if (rank==0)
	{
#ifdef WRITE_VECTOR_POTENTIALS
		MPI_Reduce(MPI_IN_PLACE,A_theta.data(),A_theta.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
		MPI_Reduce(MPI_IN_PLACE,F_phi.data(),F_phi.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
		MPI_Reduce(MPI_IN_PLACE,A_phi.data(),A_phi.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
		MPI_Reduce(MPI_IN_PLACE,F_theta.data(),F_theta.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
#endif //WRITE_VECTOR_POTENTIALS
		MPI_Reduce(MPI_IN_PLACE,E_theta.data(),E_theta.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
		MPI_Reduce(MPI_IN_PLACE,E_phi.data(),E_phi.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
	}
	else
	{
#ifdef WRITE_VECTOR_POTENTIALS
		MPI_Reduce(A_theta.data(),NULL,A_theta.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
		MPI_Reduce(F_phi.data(),NULL,F_phi.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
		MPI_Reduce(A_phi.data(),NULL,A_phi.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
		MPI_Reduce(F_theta.data(),NULL,F_theta.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
#endif //WRITE_VECTOR_POTENTIALS
		MPI_Reduce(E_theta.data(),NULL,E_theta.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
		MPI_Reduce(E_phi.data(),NULL,E_phi.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
	}
#else //USE_MPI_IN_PLACE
	//Take the sum (MPI_SUM) of all partial far-field waveforms at each node and store in node 0
#ifdef WRITE_VECTOR_POTENTIALS
	MPI_Reduce(A_theta.data(),A_theta_buf.data(),A_theta_buf.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
	MPI_Reduce(F_phi.data(),F_phi_buf.data(),F_phi_buf.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
	MPI_Reduce(A_phi.data(),A_phi_buf.data(),A_phi_buf.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
	MPI_Reduce(F_theta.data(),F_theta_buf.data(),F_theta_buf.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
#endif //WRITE_VECTOR_POTENTIALS
	MPI_Reduce(E_theta.data(),E_theta_buf.data(),E_theta_buf.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);
	MPI_Reduce(E_phi.data(),E_phi_buf.data(),E_phi_buf.size(),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_CartSubComm);

	if (rank==0)
	{
#ifdef WRITE_VECTOR_POTENTIALS
		A_theta = A_theta_buf;
		F_phi = F_phi_buf;
		A_phi = A_phi_buf;
		F_theta = F_theta_buf;
#endif //WRITE_VECTOR_POTENTIALS
		E_theta = E_theta_buf;
		E_phi = E_phi_buf;
	}
#endif //USE_MPI_IN_PLACE
#endif //MPI_DISABLE
	//if there is no MPI, there is only rank 0, which holds all the far-field arrays anyway

	//write out the far field arrays
	if (rank==0)
	{
#ifndef HDF5_DISABLE
		//create the far-field file
		H5File far_field_file(FarFieldFileName.c_str(), H5F_ACC_TRUNC);
		// Default property list
		DSetCreatPropList plist;
		/** write the fdtd code version **/
		//dataspace for three values
		hsize_t dims=3;
		DataSpace dspace(1,&dims);
		DataSet dataset=far_field_file.createDataSet("angora_version", PredType::NATIVE_INT, dspace, plist);
		int version_array[3] = {angora_version_major,angora_version_minor,angora_version_revision};
		dataset.write(&version_array,PredType::NATIVE_INT);
		/** write the size of the wavenumber array **/
		dims=1;
		dspace = DataSpace(1,&dims);
		dataset=far_field_file.createDataSet("num_lambda", PredType::NATIVE_INT, dspace, plist);
		dataset.write(&L,PredType::NATIVE_INT);
		/** write the sizes of the direction parameter arrays **/
		string dir1_name,dir2_name;
		if (Data.directionspec=="theta-phi")
		{
			dir1_name = "theta";
			dir2_name = "phi";
		}
		else if ((Data.directionspec=="dircosx-dircosy-upper")||(Data.directionspec=="dircosx-dircosy-lower"))
		{
			dir1_name = "dircos_x";
			dir2_name = "dircos_y";
		}
		//dataspace for a one dimensional array
		dims=1;
		dspace = DataSpace(1,&dims);
		dataset=far_field_file.createDataSet("num_"+dir1_name, PredType::NATIVE_INT, dspace, plist);
		dataset.write(&D1,PredType::NATIVE_INT);
		dims=1;
		dspace = DataSpace(1,&dims);
		dataset=far_field_file.createDataSet("num_"+dir2_name, PredType::NATIVE_INT, dspace, plist);
		dataset.write(&D2,PredType::NATIVE_INT);
		/** write the wavenumber and direction arrays themselves **/
		//dataspace for the wavelength array
		hsize_t wv_range_length = lambda.size();
		dspace = DataSpace(1,&wv_range_length);
		dataset=far_field_file.createDataSet("lambda", PredType::NATIVE_DOUBLE, dspace, plist);
		dataset.write(lambda.data(),PredType::NATIVE_DOUBLE);
		//dataspace for the dir1 array
		hsize_t dir1_length = dir1.size();
		dspace = DataSpace(1,&dir1_length);
		dataset=far_field_file.createDataSet(dir1_name, PredType::NATIVE_DOUBLE, dspace, plist);
		dataset.write(dir1.data(),PredType::NATIVE_DOUBLE);
		//dataspace for the dir2 array
		hsize_t dir2_length = dir2.size();
		dspace = DataSpace(1,&dir2_length);
		dataset=far_field_file.createDataSet(dir2_name, PredType::NATIVE_DOUBLE, dspace, plist);
		dataset.write(dir2.data(),PredType::NATIVE_DOUBLE);
		/** write whether the direction cosines are scaled with wavelength **/
		dims=1;
		dspace = DataSpace(1,&dims);
		dataset=far_field_file.createDataSet("scale_with_wavelength", PredType::NATIVE_CHAR, dspace, plist);
		unsigned char temp;
		temp = Data.scale_with_wavelength;
		dataset.write(&temp,PredType::NATIVE_CHAR);
		/** write the number of scattered plane waves **/
		int num_of_PWs = PW_THETA.size();	//could have also used PW_PHI.size(), since it is supposed to be the same as PW_THETA.size() (=num. of PWs)
		dims = 1;
		dspace = DataSpace(1,&dims);
		dataset=far_field_file.createDataSet("num_pw", PredType::NATIVE_INT, dspace, plist);
		dataset.write(&num_of_PWs,PredType::NATIVE_INT);
		/** write the propagation directions of the scattered plane waves  **/
		if (num_of_PWs>0)
		{//don't write if there are no plane waves (num_of_PWs=0), otherwise PW_THETA(0) below will give runtime error
			dims = num_of_PWs;
			dspace = DataSpace(1,&dims);
			//theta angles
			dataset=far_field_file.createDataSet("pw_theta", PredType::NATIVE_DOUBLE, dspace, plist);
			dataset.write(PW_THETA.data(),PredType::NATIVE_DOUBLE);
			//phi angles
			dataset=far_field_file.createDataSet("pw_phi", PredType::NATIVE_DOUBLE, dspace, plist);
			dataset.write(PW_PHI.data(),PredType::NATIVE_DOUBLE);
		}

		/** write the far-field information **/
		//dataspace for the 3D arrays
		hsize_t farfielddims[] = {L, D1, D2};
		DataSpace dspace3D( 3, farfielddims );
		//create temporary arrays to hold the 3D data [less elegant, but more portable than creating strided dataspaces]
		Array<double,3> temp3D_1(L,D1,D2),temp3D_2(L,D1,D2);
		Array<double,2> temp2D_1(PW_THETA.size(),L),temp2D_2(PW_THETA.size(),L);
#ifdef WRITE_VECTOR_POTENTIALS
		dataset=far_field_file.createDataSet("A_theta_r", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = real(A_theta);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("A_theta_i", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = imag(A_theta);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("A_phi_r", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = real(A_phi);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("A_phi_i", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = imag(A_phi);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("F_theta_r", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = real(F_theta);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("F_theta_i", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = imag(F_theta);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("F_phi_r", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = real(F_phi);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("F_phi_i", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = imag(F_phi);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
#else
		dataset=far_field_file.createDataSet("E_theta_r", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = real(E_theta);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("E_theta_i", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = imag(E_theta);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("E_phi_r", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = real(E_phi);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("E_phi_i", PredType::NATIVE_DOUBLE, dspace3D, plist);
		temp3D_1 = imag(E_phi);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
#ifdef WRITE_THEORETICAL_FIELD
		for (int l=0; l<L; l++)
		{
			for (int d1=0; d1<D1; d1++)
			{
				for (int d2=0; d2<D2; d2++)
				{
					temp3D_1(l,d1,d2) = real(TheoreticalFarFieldTheta(l,d1,d2));
					temp3D_2(l,d1,d2) = imag(TheoreticalFarFieldTheta(l,d1,d2));
				}
			}
		}
		dataset=far_field_file.createDataSet("E_theta_th_r", PredType::NATIVE_DOUBLE, dspace3D, plist);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("E_theta_th_i", PredType::NATIVE_DOUBLE, dspace3D, plist);
		dataset.write(temp3D_2.data(),PredType::NATIVE_DOUBLE);
		for (int l=0; l<L; l++)
		{
			for (int d1=0; d1<D1; d1++)
			{
				for (int d2=0; d2<D2; d2++)
				{
					temp3D_1(l,d1,d2) = real(TheoreticalFarFieldPhi(l,d1,d2));
					temp3D_2(l,d1,d2) = imag(TheoreticalFarFieldPhi(l,d1,d2));
				}
			}
		}
		dataset=far_field_file.createDataSet("E_phi_th_r", PredType::NATIVE_DOUBLE, dspace3D, plist);
		dataset.write(temp3D_1.data(),PredType::NATIVE_DOUBLE);
		dataset=far_field_file.createDataSet("E_phi_th_i", PredType::NATIVE_DOUBLE, dspace3D, plist);
		dataset.write(temp3D_2.data(),PredType::NATIVE_DOUBLE);
#endif
#endif
		//deallocate the 3D temp arrays
		temp3D_1.resize(0,0,0);
		temp3D_2.resize(0,0,0);
		//write out the scattered PW phasors
		if (PW_THETA.size()>0)
		{
			//dataspace for the 2D arrays
			hsize_t pwphasordims[] = {PW_THETA.size(),L};
			DataSpace dspace2D( 2, pwphasordims);
			for (int l=0; l<L; l++)
			{
				for (int PW=0; PW<PW_THETA.size(); PW++)	//could have also used "PW_PHI.size()", since they are supposed to be of the same size (num. of PWs)
				{
					//write the theta-component of the E-field
					temp2D_1(PW,l) = real(E_x_phasor_array(l,PW));
					temp2D_2(PW,l) = imag(E_x_phasor_array(l,PW));
				}
			}
			dataset=far_field_file.createDataSet("E_x_pw_r", PredType::NATIVE_DOUBLE, dspace2D, plist);
			dataset.write(temp2D_1.data(),PredType::NATIVE_DOUBLE);
			dataset=far_field_file.createDataSet("E_x_pw_i", PredType::NATIVE_DOUBLE, dspace2D, plist);
			dataset.write(temp2D_2.data(),PredType::NATIVE_DOUBLE);
			for (int l=0; l<L; l++)
			{
				for (int PW=0; PW<PW_THETA.size(); PW++)	//could have also used "PW_PHI.size()", since they are supposed to be of the same size (num. of PWs)
				{
					//write the phi-component of the E-field
					temp2D_1(PW,l) = real(E_y_phasor_array(l,PW));
					temp2D_2(PW,l) = imag(E_y_phasor_array(l,PW));
				}
			}
			dataset=far_field_file.createDataSet("E_y_pw_r", PredType::NATIVE_DOUBLE, dspace2D, plist);
			dataset.write(temp2D_1.data(),PredType::NATIVE_DOUBLE);
			dataset=far_field_file.createDataSet("E_y_pw_i", PredType::NATIVE_DOUBLE, dspace2D, plist);
			dataset.write(temp2D_2.data(),PredType::NATIVE_DOUBLE);
		}

		far_field_file.close();
#else //#ifndef HDF5_DISABLE
		if (rank==0)
		{
			cout << "Error: Custom file format is obsolete for far-field output files." << endl;
		}
		MPI_exit(-1);
#endif //#ifndef HDF5_DISABLE
	}
}
