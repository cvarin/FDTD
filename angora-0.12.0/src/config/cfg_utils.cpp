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

//some utilities for config-file usage

#include "headers.h"

#include "cfg_utils.h"

extern int GridIndex;
extern int rank;

extern void MPI_exit(const int& exitcode);


void CheckAngoraGroupSetting(const Setting& mygroup, const Config& validsettings)
{//check the group "mygroup" for invalid settings
	//first, check if the setting is actually a group
	if (!mygroup.isGroup())
	{
		throw AngoraInvalidSettingTypeException(mygroup,"should be a group");
	}
	//get the corresponding group with all the valid settings
	try{
	string path_to_group = mygroup.getPath();
	//now, change all the integer indices in "path_to_group" to 0, since there is only one setting in the template file
	size_t open_paranthesis_pos,close_paranthesis_pos;
	open_paranthesis_pos=path_to_group.find("[");
	while (open_paranthesis_pos!=string::npos)
	{
	close_paranthesis_pos=path_to_group.find("]",open_paranthesis_pos+1);
	path_to_group.replace(open_paranthesis_pos+1,close_paranthesis_pos-open_paranthesis_pos-1,"0");
	close_paranthesis_pos -= close_paranthesis_pos-open_paranthesis_pos-2; //the string is shortened, so the position of "]" changes
	open_paranthesis_pos=path_to_group.find("[",close_paranthesis_pos+1);
	}
	//all indices changed to [0].
	//look up the corresponding group in the template config file
	Setting& group_with_valid_settings = validsettings.lookup(path_to_group);
	for (int setting_index = 0; setting_index<mygroup.getLength(); setting_index++)
	{
		string this_setting_name = mygroup[setting_index].getName();
		if (((this_setting_name!="enabled_for_runs")||mygroup.isRoot())&&(!group_with_valid_settings.exists(this_setting_name)))
		{
			throw AngoraInvalidSettingNameException(mygroup[setting_index]);
		}
	}
	}
	catch (SettingNotFoundException& exc)
	{
		throw AngoraDeveloperException("Angora is seeking an invalid configuration option named "+string(exc.getPath()));
	}
}

bool SettingEnabledForGrid(const Setting& mySetting)
{//determines if the setting "mySetting" is enabled for the current grid
	bool enabled_for_this_grid = false;
	if (mySetting.exists("enabled_for_runs"))
	{
		Setting& EnabledGrids = mySetting["enabled_for_runs"];
		if (EnabledGrids.getType()!=Setting::TypeArray)
		{
			throw AngoraInvalidSettingTypeException(EnabledGrids,"should be array");
		}
		if (EnabledGrids.getLength()==0)
		{
			throw AngoraSettingEmptyException(EnabledGrids);
		}
		if ((EnabledGrids[0].getType()!=Setting::TypeInt)&&(EnabledGrids[0].getType()!=Setting::TypeInt64))
		{
			throw AngoraInvalidSettingTypeException(EnabledGrids,"should be an array of unsigned integers");
		}
		//if the current grid index is in list EnabledGrids, "enabled_for_this_grid" is true
		for (int grid=0;grid<EnabledGrids.getLength();grid++)
		{
			if ((int)EnabledGrids[grid]==GridIndex)
			{
				enabled_for_this_grid = true;
			}
		}
	}
	else
	{
		enabled_for_this_grid = true;
	}

	return enabled_for_this_grid;
}

template<>
void read_value_from_group(const Setting& parent_group, const string& setting_name, string& setting_variable)
{
	try{setting_variable = (const char*)parent_group[setting_name];}
	catch (SettingNotFoundException& exc)
	{
		throw AngoraSettingNotFoundException(setting_name,parent_group);
	}
	catch (SettingTypeException& exc)
	{
		throw AngoraInvalidSettingTypeException(parent_group[setting_name],"should be of type string");
	}
}

Setting& read_list_from_group(const Setting& parent_group, const string& list_name)
{
	try{
	Setting& list_setting = parent_group[list_name];
	if (!list_setting.isList())
	{
		throw AngoraInvalidSettingTypeException(parent_group[list_name],"should be a list");
	}
	} //try block
	catch (SettingNotFoundException& exc)
	{
		throw AngoraSettingNotFoundException(list_name,parent_group);
	}
	return parent_group[list_name];
}

Setting& read_array_from_group(const Setting& parent_group, const string& array_name)
{
	try{
	Setting& array_setting = parent_group[array_name];
	if (!array_setting.isArray())
	{
		throw AngoraInvalidSettingTypeException(parent_group[array_name],"should be an array");
	}
	} //try block
	catch (SettingNotFoundException& exc)
	{
		throw AngoraSettingNotFoundException(array_name,parent_group);
	}
	return parent_group[array_name];
}

Setting& read_group_from_group(const Setting& parent_group, const string& group_name)
{
	try{
	Setting& group_setting = parent_group[group_name];
	if (!group_setting.isGroup())
	{
		throw AngoraInvalidSettingTypeException(parent_group[group_name],"should be a group");
	}
	} //try block
	catch (SettingNotFoundException& exc)
	{
		throw AngoraSettingNotFoundException(group_name,parent_group);
	}
	return parent_group[group_name];
}
