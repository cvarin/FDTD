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

// reading the config options that specify the material definitions

#include "headers.h"

#include "read_mat.h"

#include "Cmats.h"

#include "material.h"


void read_mat(Cmats &Materials, const Config& fdtdconfig, const Config& validsettings)
{
	string material_setting_path = "Materials";
	if (fdtdconfig.exists(material_setting_path))
	{
		Setting& Materiallistsettings = read_list_from_group(fdtdconfig.getRoot(),material_setting_path);

		int num_of_newmaterials = Materiallistsettings.getLength();
		for (int newmaterialindex=0; newmaterialindex<num_of_newmaterials; newmaterialindex++)
		{
			Setting& Newmaterialsettings = Materiallistsettings[newmaterialindex];	//go to the newmaterialindex'th new material setting
			//check group for invalid settings
			CheckAngoraGroupSetting(Newmaterialsettings,validsettings);

			if (SettingEnabledForGrid(Materiallistsettings))		//apply only if enabled for this grid
			{
				double rel_permeability;
				try{read_value_from_group<double>(Newmaterialsettings,"rel_permeability",rel_permeability);}
				catch (AngoraSettingNotFoundException&)
				{
					rel_permeability = 1; //default: mu_r = 1
				}

				double electric_conductivity;
				try{read_value_from_group<double>(Newmaterialsettings,"electric_conductivity",electric_conductivity);}
				catch (AngoraSettingNotFoundException&)
				{
					//May implement better warning system later
					electric_conductivity = 0; //default: sigma = 0
				}

				double magnetic_conductivity;
				try{read_value_from_group<double>(Newmaterialsettings,"magnetic_conductivity",magnetic_conductivity);}
				catch (AngoraSettingNotFoundException&)
				{
					//May implement better warning system later
					magnetic_conductivity = 0; //default: sigma = 0
				}

				//this should be specified
				double rel_permittivity;
				read_value_from_group<double>(Newmaterialsettings,"rel_permittivity",rel_permittivity);

				string material_tag;
				read_value_from_group<string>(Newmaterialsettings,"material_tag",material_tag);

				//create the material identifier for the new material
				MaterialId NewMaterialId;
				AddIsotropicMaterial(NewMaterialId,rel_permittivity,rel_permeability,electric_conductivity,magnetic_conductivity);
				//Add the new material to the material-collector object (even in check mode, since other objects refer to these materials)
				Materials.CreateMaterial(NewMaterialId,material_tag);
			}
		}
	}//if exists
}
