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

#include "headers.h"

#include "read_logging_data.h"

//definition of Cestimator needed
#include "Cestimator.h"

extern bool check_mode;

extern string OutputDir;
extern string LogOutputDir;
extern const string default_LogOutputDir;
extern const string default_LogFileName;

extern int NSTEPS;

extern void add_slash_to_path(string& path);
extern bool is_absolute_path(string& path);
extern int create_path(const string& path);


void read_logging_data(Cestimator &Estimator, const Config& fdtdconfig)
{
	//Read the log output path
	try{read_value_from_group<string>(fdtdconfig.getRoot(),"log_output_dir",LogOutputDir);}
	catch (AngoraSettingNotFoundException&)
	{
		LogOutputDir = default_LogOutputDir;
	}
	//add slash to path if necessary
	add_slash_to_path(LogOutputDir);
	//if path is not absolute, prepend the base path to get full path
	if (!is_absolute_path(LogOutputDir))
	{
		//prepend the output base path to get full path
		LogOutputDir = OutputDir + LogOutputDir;
	}
	//create output directory if it does not exist
	if (!check_mode)
	{
		if (create_path(LogOutputDir)<0)
		{
			/** throw exception **/
		}
	}
	//the name of the log file
	string logfilename,fulllogfilename;
	if (!read_optional_value_from_group(fdtdconfig.getRoot(),"log_file_name",logfilename))
	{
		logfilename = default_LogFileName;
	}
	fulllogfilename = LogOutputDir+logfilename;

	//initialize the estimator
	Estimator.Initialize(fulllogfilename,NSTEPS);
}
