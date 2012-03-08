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

//Defines the field-value recorder class "Crecorder_fieldvalue"

#include "headers.h"

#include "Crecorder_fieldvalue.h"

//include the HDF5 header, if not disabled
#ifndef HDF5_DISABLE
#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif
#endif

//routines for getting/setting the initial time value in the simulation
#include "time_axis.h"

extern int angora_version_major,angora_version_minor,angora_version_revision;

extern double dt;
extern int OriginX,OriginY,OriginZ;
extern int NSTEPS;

extern Array<double,3> Ex,Ey,Ez;
extern Array<double,3> Hx,Hy,Hz;

extern int GridIndex;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern double max_field_value;
extern double dB_accuracy;


Crecorder_fieldvalue::Crecorder_fieldvalue(const int& xPos, const int& yPos, const int& zPos,
		 const string& mycomponent, const string& myscale, const string& myfilename,
		 const int& Index)
		 : relative_FieldPos_x(xPos), relative_FieldPos_y(yPos), relative_FieldPos_z(zPos), component(mycomponent), scale(myscale),
		 FieldValueFileName(myfilename), FieldValueRecorderIndex(Index)
{//constructor for the field value recorder
#ifdef HDF5_DISABLE
	throw AngoraDeveloperException("Custom file format is obsolete for field-value recorder output files.");
#endif
	//Extract the maximum and minimum values in the discretization of the fields
	if (scale=="dB")
	{
		FieldMax = 20.0*log10(abs(max_field_value)+1e-20);
		FieldMin = FieldMax+dB_accuracy;
	}
	else if (scale=="linear")
	{
		FieldMax = max_field_value;
		FieldMin = -max_field_value;
	}
	else if (scale=="absolute")
	{
		FieldMax = max_field_value;
		FieldMin = 0;
	}
	else
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//			InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,scale,
			"(valid arguments are \"dB\", \"linear\", or \"absolute\")");
	}

	if ((component!="E")&&(component!="Ex")&&(component!="Ey")&&(component!="Ez"))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//			InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,component,
			"(valid arguments are \"E\", \"Ex\", \"Ey\", or \"Ez\")");
	}

	//allocate the array
	field_value_array.resize(NSTEPS);
	field_value_array = 0;

	//initialize the counter
	counter = 0;

	FieldPos_x = relative_FieldPos_x + OriginX;
	FieldPos_y = relative_FieldPos_y + OriginY;
	FieldPos_z = relative_FieldPos_z + OriginZ;

	FieldComponentInNode = ((iback<=FieldPos_x)&&(ifront>=FieldPos_x)
		&&(jleft<=FieldPos_y)&&(jright>=FieldPos_y)
		&&(klower<=FieldPos_z)&&(kupper>=FieldPos_z));
}

void Crecorder_fieldvalue::RecordFieldValue()
{//records a single field value
	//act only if recorded position is within the node
	if (FieldComponentInNode)
	{
		if (component=="E")
		{
			field_value = sqrt(pow((Ex(FieldPos_x,FieldPos_y,FieldPos_z)+Ex(FieldPos_x,FieldPos_y+1,FieldPos_z)
									+Ex(FieldPos_x,FieldPos_y,FieldPos_z+1)+Ex(FieldPos_x,FieldPos_y+1,FieldPos_z+1))/4.0,2)
						+pow((Ey(FieldPos_x,FieldPos_y,FieldPos_z)+Ey(FieldPos_x+1,FieldPos_y,FieldPos_z)
								+Ey(FieldPos_x,FieldPos_y,FieldPos_z+1)+Ey(FieldPos_x+1,FieldPos_y,FieldPos_z+1))/4.0,2)
						+pow((Ez(FieldPos_x,FieldPos_y,FieldPos_z)+Ez(FieldPos_x,FieldPos_y+1,FieldPos_z)
								+Ez(FieldPos_x+1,FieldPos_y,FieldPos_z)+Ez(FieldPos_x+1,FieldPos_y+1,FieldPos_z))/4.0,2));
		}
		else if (component=="Ex")
		{
			field_value = Ex(FieldPos_x,FieldPos_y,FieldPos_z);
		}
		else if (component=="Ey")
		{
			field_value = Ey(FieldPos_x,FieldPos_y,FieldPos_z);
		}
		else if (component=="Ez")
		{
			field_value = Ez(FieldPos_x,FieldPos_y,FieldPos_z);
		}

		if (scale=="dB")
		{
			field_value = 20.0*log10(abs(field_value)+1e-20);
		}
		else if (scale=="absolute")
		{
			field_value = abs(field_value);
		}

		//store the field value (no need for discretization, since only one value is written)
		field_value_array(counter) = field_value;

		counter++;
	}
}

Crecorder_fieldvalue::~Crecorder_fieldvalue()
{//dtor writes the array into the file before quitting
	if (FieldComponentInNode)
	{
#ifndef HDF5_DISABLE
		//create the far-field file
		H5File field_value_file(FieldValueFileName.c_str(), H5F_ACC_TRUNC);
		// Default property list
		DSetCreatPropList plist;
		//number of values written at one time
		hsize_t chunksize;
		//dataspace for three values
		chunksize=3;
		DataSpace dspace(1,&chunksize);
		DataSet dataset=field_value_file.createDataSet("angora_version", PredType::NATIVE_INT, dspace, plist);
		int version_array[3] = {angora_version_major,angora_version_minor,angora_version_revision};
		dataset.write(&version_array,PredType::NATIVE_INT);
		//write number of steps
		chunksize = 1;
		dspace = DataSpace(1,&chunksize);
		dataset=field_value_file.createDataSet("num_time_steps", PredType::NATIVE_INT, dspace, plist);
		dataset.write(&NSTEPS,PredType::NATIVE_INT);
		//write time step
		dataset=field_value_file.createDataSet("time_step", PredType::NATIVE_DOUBLE, dspace, plist);
		dataset.write(&dt,PredType::NATIVE_DOUBLE);
		//write initial time value
		double initial_time_value = get_initial_time_value();
		dataset=field_value_file.createDataSet("initial_time_value", PredType::NATIVE_DOUBLE, dspace, plist);
		dataset.write(&initial_time_value,PredType::NATIVE_DOUBLE);
		//write the field-value array
		chunksize=field_value_array.size();
		dspace = DataSpace(1,&chunksize);
		dataset=field_value_file.createDataSet("field_values", PredType::NATIVE_DOUBLE, dspace, plist);
		dataset.write(field_value_array.data(),PredType::NATIVE_DOUBLE);
#endif //#ifndef HDF5_DISABLE
	}
}
