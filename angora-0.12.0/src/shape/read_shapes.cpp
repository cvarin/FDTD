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

#include "read_shapes.h"

#include "Cshapes.h"

extern bool check_mode;

extern int rank;

extern void MPI_exit(const int& exitcode);


void read_shapes(Cshapes &Shapes, const Config& fdtdconfig, const Config& validsettings)
{
	//Read shape settings
	string shape_setting_path = "Shapes";
	if (fdtdconfig.exists(shape_setting_path))
	{
		Setting& shape_settings = fdtdconfig.lookup(shape_setting_path);
		//check group for invalid settings
		CheckAngoraGroupSetting(shape_settings,validsettings);

		//Read planar layer settings
		string planar_layers_setting_path = "PlanarLayers";
		if (shape_settings.exists(planar_layers_setting_path))
		{
			Setting& planar_layer_collec_settings = read_list_from_group(shape_settings,planar_layers_setting_path);

			int num_of_planar_layers = planar_layer_collec_settings.getLength();
			for (int planar_layer_index=0; planar_layer_index<num_of_planar_layers; planar_layer_index++)
			{
				Setting& planar_layer_settings = planar_layer_collec_settings[planar_layer_index];	//go to the planar_layer_index'th planar-layer setting
				//check group for invalid settings
				CheckAngoraGroupSetting(planar_layer_settings,validsettings);

				if (SettingEnabledForGrid(planar_layer_settings))		//apply only if enabled for this grid
				{
					//read tag string
					string shape_tag;
					read_value_from_group<string>(planar_layer_settings,"shape_tag",shape_tag);
					//read the orientation of the layer
					string orientation;
					read_value_from_group<string>(planar_layer_settings,"orientation",orientation);
					//get the bounds of the layer
					//get starting position
					double min_coord;
					try{read_value_from_group<double>(planar_layer_settings,"min_coord",min_coord);}
					catch (AngoraInvalidSettingTypeException& exc)
					{//it could be string
						string min_coord_str;
						try{read_value_from_group<string>(planar_layer_settings,"min_coord",min_coord_str);}
						catch (AngoraInvalidSettingTypeException& exc)
						{//if it is neither numeric nor string, the exception is caught here
							throw AngoraInvalidSettingTypeException(planar_layer_settings["min_coord"],"should be a numeric value or the string \"-inf\"");
						}
						if (min_coord_str=="-inf")
						{
							/** FIXME: not the most ideal solution **/
							min_coord = -LIBSTD_DBL_MAX;  //negative infinity
						}
						else
						{
							throw AngoraInvalidSettingValueException(planar_layer_settings["min_coord"],"only allowable string value is \"-inf\"");
						}
					}
					//get end position
					double max_coord;
					try{read_value_from_group<double>(planar_layer_settings,"max_coord",max_coord);}
					catch (AngoraInvalidSettingTypeException& exc)
					{//it could be string
						string max_coord_str;
						try{read_value_from_group<string>(planar_layer_settings,"max_coord",max_coord_str);}
						catch (AngoraInvalidSettingTypeException& exc)
						{//if it is neither numeric nor string, the exception is caught here
							throw AngoraInvalidSettingTypeException(planar_layer_settings["max_coord"],"should be a numeric value or the string \"+inf\"");
						}
						if (max_coord_str=="+inf")
						{
							/** FIXME: not the most ideal solution **/
							max_coord = LIBSTD_DBL_MAX;  //positive infinity
						}
						else
						{
							throw AngoraInvalidSettingValueException(planar_layer_settings["max_coord"],"only allowable string value is \"+inf\"");
						}
					}
					//Add the planar layer to the shape collector object (even in check mode, since other objects refer to these shapes)
					Shapes.CreatePlanarLayer(min_coord,max_coord,orientation,shape_tag);
				}
			}
		}
		//Read rectangular box settings
//		else if (shape_setting_name=="RectangularBoxes")
//		{
//			Setting& rectangular_box_collec_settings = shape_setting;
		string rect_boxes_setting_path = "RectangularBoxes";
		if (shape_settings.exists(rect_boxes_setting_path))
		{
			Setting& rectangular_box_collec_settings = read_list_from_group(shape_settings,rect_boxes_setting_path);

			int num_of_rectangular_boxes = rectangular_box_collec_settings.getLength();
			for (int rectangular_box_index=0; rectangular_box_index<num_of_rectangular_boxes; rectangular_box_index++)
			{
				Setting& rectangular_box_settings = rectangular_box_collec_settings[rectangular_box_index];	//go to the rectangular_box_index'th rectangular box setting
				//check group for invalid settings
				CheckAngoraGroupSetting(rectangular_box_settings,validsettings);

				if (SettingEnabledForGrid(rectangular_box_settings))		//apply only if enabled for this grid
				{
					//read tag string
					string shape_tag;
					read_value_from_group<string>(rectangular_box_settings,"shape_tag",shape_tag);
					//get the bounds of the rectangular box
					//back_x
					double back_x;
					read_value_from_group<double>(rectangular_box_settings,"back_x",back_x);
					//front_x
					double front_x;
					read_value_from_group<double>(rectangular_box_settings,"front_x",front_x);
					//left_y
					double left_y;
					read_value_from_group<double>(rectangular_box_settings,"left_y",left_y);
					//right_y
					double right_y;
					read_value_from_group<double>(rectangular_box_settings,"right_y",right_y);
					//lower_z
					double lower_z;
					read_value_from_group<double>(rectangular_box_settings,"lower_z",lower_z);
					//upper_z
					double upper_z;
					read_value_from_group<double>(rectangular_box_settings,"upper_z",upper_z);
					//Add the rectangular box to the shape collector object (even in check mode, since other objects refer to these shapes)
					Shapes.CreateRectBox(back_x,front_x,left_y,right_y,lower_z,upper_z,shape_tag);
				}
			}
		}
		//Read sphere settings
		string spheres_setting_path = "Spheres";
		if (shape_settings.exists(spheres_setting_path))
		{
			Setting& sphere_collec_settings = read_list_from_group(shape_settings,spheres_setting_path);

			int num_of_spheres = sphere_collec_settings.getLength();
			for (int sphere_index=0; sphere_index<num_of_spheres; sphere_index++)
			{
				Setting& sphere_settings = sphere_collec_settings[sphere_index];
				//check group for invalid settings
				CheckAngoraGroupSetting(sphere_settings,validsettings);

				if (SettingEnabledForGrid(sphere_settings))		//apply only if enabled for this grid
				{
					//read tag string
					string shape_tag;
					read_value_from_group<string>(sphere_settings,"shape_tag",shape_tag);
					//get the center of the sphere
					//center_x
					double center_x;
					read_value_from_group<double>(sphere_settings,"center_x",center_x);
					//center_y
					double center_y;
					read_value_from_group<double>(sphere_settings,"center_y",center_y);
					//center_z
					double center_z;
					read_value_from_group<double>(sphere_settings,"center_z",center_z);
					//get the radius of the sphere
					double radius;
					read_value_from_group<double>(sphere_settings,"radius",radius);
					//Add the sphere to the shape collector object (even in check mode, since other objects refer to these shapes)
					Shapes.CreateSphere(center_x,center_y,center_z,radius,shape_tag);
				}
			}
		}
	}
}
