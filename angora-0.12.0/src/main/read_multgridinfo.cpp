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

//Includes the routine that reads the multiple grid information

#include "headers.h"

#include "read_multgridinfo.h"

extern const int default_number_of_runs;
extern int number_of_runs;
extern Array<bool,1> grid_is_enabled;

extern int rank;


void read_multgridinfo(const Config& fdtdconfig, const Config& validsettings)
{
	//Read number of independent simulation runs
	try{read_value_from_group<int>(fdtdconfig.getRoot(),"number_of_runs",number_of_runs);}
	catch (AngoraSettingNotFoundException&)
	{
		number_of_runs = default_number_of_runs;
	}

	//Read the disabled grid indices
	// individually-disabled grids
	grid_is_enabled.resize(number_of_runs);
	grid_is_enabled = true;
	if (fdtdconfig.exists("disabled_runs"))
	{
		Setting& disabled_runs = fdtdconfig.lookup("disabled_runs");
		if (!disabled_runs.isArray())
		{
			throw AngoraInvalidSettingTypeException(disabled_runs,"should be an array");
		}
		if (disabled_runs.getLength()!=0)
		{
			if (disabled_runs[0].getType()!=Setting::TypeInt)
			{
				throw AngoraInvalidSettingTypeException(disabled_runs,"should be an array of integers");
			}
			//mark the disabled grids
			for (int i=0;i<disabled_runs.getLength();i++)
			{
				int disabled_run_index = disabled_runs[i];
				if ((disabled_run_index>=0)&&(disabled_run_index<number_of_runs))
				{
					grid_is_enabled(disabled_run_index) = false;
				}
			}
		}
	}
	// disabled grid ranges
	if (fdtdconfig.exists("disabled_run_range"))
	{
		Setting& disabled_run_range = fdtdconfig.lookup("disabled_run_range");
		if (!disabled_run_range.isArray())
		{
			throw AngoraInvalidSettingTypeException(disabled_run_range,"should be an array");
		}
		if (disabled_run_range.getLength()!=0)
		{
			if (disabled_run_range[0].getType()!=Setting::TypeInt)
			{
				throw AngoraInvalidSettingTypeException(disabled_run_range,"should be an array of integers");
			}
			if (disabled_run_range.getLength()!=2)
			{
				throw AngoraInvalidSettingValueException(disabled_run_range, "should be an array with 2 components");
			}
			//mark the disabled grids
			for (int disabled_run_index=disabled_run_range[0];disabled_run_index<=(int)disabled_run_range[1];disabled_run_index++)
			{
				if ((disabled_run_index>=0)&&(disabled_run_index<number_of_runs))
				{
					grid_is_enabled(disabled_run_index) = false;
				}
			}
		}
	}
}
