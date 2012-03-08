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

//Function that reads some global variables from the config file

#include "headers.h"

#include "read_global.h"

extern string config_filename;

extern double courant,dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML,NSTEPS;
extern int OriginX,OriginY,OriginZ;

extern double max_field_value;	//maximum electric field value encountered in grid
extern bool max_field_value_set_in_configfile;		//is the max. field value fixed in the config file?
extern double dB_accuracy;
extern bool dB_accuracy_set_in_configfile;


void read_global(const Config& fdtdconfig, const Config& validsettings)
//Reads some global variables
{
	//Read and override maximum field value, if specified in the config file
	if (!read_optional_value_from_group<double>(fdtdconfig.getRoot(),"max_field_value",max_field_value))
	{
		max_field_value_set_in_configfile = false;	//max. field value set in config file
	}
	//Read and override dB accuracy, if specified in the config file
	if (!read_optional_value_from_group<double>(fdtdconfig.getRoot(),"dB_accuracy",dB_accuracy))
	{
		dB_accuracy_set_in_configfile = false;	//dB accuracy set in config file
	}


	//Read courant stability factor
	read_value_from_group<double>(fdtdconfig.getRoot(),"courant",courant);

	//Read grid spacing
	read_value_from_group<double>(fdtdconfig.getRoot(),"dx",dx);

	dt = courant*dx/c/1.73205;		//Time step

	//Read number of grid cells in the x,y,z directions
	read_value_from_group<int>(fdtdconfig.getRoot(),"NCELLS_X",NCELLS_X);
	read_value_from_group<int>(fdtdconfig.getRoot(),"NCELLS_Y",NCELLS_Y);
	read_value_from_group<int>(fdtdconfig.getRoot(),"NCELLS_Z",NCELLS_Z);

	//Read PML thickness
	read_value_from_group<int>(fdtdconfig.getRoot(),"NPML",NPML);

	//Read number of time steps
	read_value_from_group<int>(fdtdconfig.getRoot(),"NSTEPS",NSTEPS);

	//Read x index of the origin cell
	try{read_value_from_group<int>(fdtdconfig.getRoot(),"OriginX",OriginX);}
	catch (AngoraSettingNotFoundException&)
	{
		OriginX = (NCELLS_X+2*NPML)/2+1;
	}
	try{read_value_from_group<int>(fdtdconfig.getRoot(),"OriginY",OriginY);}
	catch (AngoraSettingNotFoundException&)
	{
		OriginY = (NCELLS_Y+2*NPML)/2+1;
	}
	try{read_value_from_group<int>(fdtdconfig.getRoot(),"OriginZ",OriginZ);}
	catch (AngoraSettingNotFoundException&)
	{
		OriginZ = (NCELLS_Z+2*NPML)/2+1;
	}
}
