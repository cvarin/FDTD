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

//Includes the routine that reads the base and output directory info from the config file

#include "headers.h"

#include "read_basedirs.h"

//For file-directory  manipulations
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern bool check_mode;

extern string currentworkdir;
extern string angora_basepath;
extern string OutputDir,InputDir;

extern const string default_OutputDir,default_InputDir;

extern void add_slash_to_path(string& path);

extern bool is_absolute_path(string& path);

extern int create_path(const string& path);


void read_basedirs(const Config& fdtdconfig, const Config& validsettings)
{
	//Read the base path for file operations
	try{read_value_from_group<string>(fdtdconfig.getRoot(),"angora_basepath",angora_basepath);}
	catch (AngoraSettingNotFoundException&)
	{
		angora_basepath = currentworkdir;
	}
	//add slash to path if necessary
	add_slash_to_path(angora_basepath);
	//create base directory if it does not exist
	if (!check_mode)
	{
		if (create_path(angora_basepath)<0)
		{
			/** throw exception **/
		}
	}

	//Read the output base path
	try{read_value_from_group<string>(fdtdconfig.getRoot(),"output_dir",OutputDir);}
	catch (AngoraSettingNotFoundException&)
	{
		OutputDir = default_OutputDir;
	}
	//add slash to path if necessary
	add_slash_to_path(OutputDir);
	//if path is not absolute, prepend the base path to get full path
	if (!is_absolute_path(OutputDir))
	{
		//prepend the base path to get full path
		OutputDir = angora_basepath + OutputDir;
	}
	//create output directory if it does not exist
	if (!check_mode)
	{
		if (create_path(OutputDir)<0)
		{
			/** throw exception **/
		}
	}

	//Read the input base path
	try{read_value_from_group<string>(fdtdconfig.getRoot(),"input_dir",InputDir);}
	catch (AngoraSettingNotFoundException&)
	{
		InputDir = default_InputDir;
	}
	//add slash to path if necessary
	add_slash_to_path(InputDir);
	//if path is not absolute, prepend the base path to get full path
	if (!is_absolute_path(InputDir))
	{
		InputDir = angora_basepath + InputDir;
	}
}
