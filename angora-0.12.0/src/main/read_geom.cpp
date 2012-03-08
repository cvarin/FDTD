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

//Includes the routine that reads the geometry settings

#include "headers.h"

//base Angora exception class
#include "angora_excp.h"

#include "read_geom.h"

#include "shape/Cshapes.h"
#include "material/Cmats.h"

#include "place_obj.h"
#include "placegeom.h"

//for reading a material region from a file
#include "matfile.h"
//for creating a random material region
#include "random.h"

extern string config_filename;

extern bool check_mode;

extern string InputDir;

extern int NCELLS_Z,NPML;
extern int OriginX,OriginY,OriginZ;

extern int GridIndex;
extern int rank;

extern void add_slash_to_path(string& path);
extern bool is_absolute_path(string& path);

extern void MPI_exit(const int& exitcode);


void read_geom(const Config& fdtdconfig, const Config& validsettings, const Cshapes &Shapes, const Cmats& Materials)
{
	//Read the geometry settings
	string geometry_setting_path = "SimulationSpace";
	if (fdtdconfig.exists(geometry_setting_path))
	{
		Setting& Geometrysettings = fdtdconfig.lookup(geometry_setting_path);
		//check group for invalid settings
		CheckAngoraGroupSetting(Geometrysettings,validsettings);

		//treat the settings in order of appearance in the config file
        int num_of_geometry_settings = Geometrysettings.getLength();
        for (int geometry_setting_index=0; geometry_setting_index<num_of_geometry_settings; geometry_setting_index++)
        {
        	Setting& Geometrysetting = Geometrysettings[geometry_setting_index]; //get the current setting
            string geometry_setting_name = Geometrysetting.getName();       //get the name of the current setting

			//Read material objects
			if (geometry_setting_name=="Objects")
			{
				//redundant, but robust and shorthand
				Setting& Objectsettings = read_list_from_group(Geometrysettings,geometry_setting_name);

				//treat the settings in order of appearance in the config file
				int num_of_object_settings = Objectsettings.getLength();
				for (int object_setting_index=0; object_setting_index<num_of_object_settings; object_setting_index++)
				{
					Setting& Objectsetting = Objectsettings[object_setting_index]; //get the current setting
					//check group for invalid settings
					CheckAngoraGroupSetting(Objectsetting,validsettings);

					//read the shape
					string shape_tag;
					read_value_from_group<string>(Objectsetting,"shape_tag",shape_tag);
					const Cshape* shapeptr = Shapes[shape_tag];

					//read the material
					string material_tag;
					read_value_from_group<string>(Objectsetting,"material_tag",material_tag);
					MaterialId mat_id = Materials[material_tag];

					if (!check_mode)
					{
						//place the object into the grid
						place_obj(mat_id,shapeptr);
					}
				}
			}//if exists
/***************************************/
/** REMOVE THE FOLLOWING OPTION LATER **/
/***************************************/
			//Read dielectric slabs
			else if (geometry_setting_name=="MaterialSlabs")
			{
				//redundant, but robust and shorthand
				Setting& Slablistsettings = read_list_from_group(Geometrysettings,geometry_setting_name);

				int num_of_slabs = Slablistsettings.getLength();
				for (int slabindex=0; slabindex<num_of_slabs; slabindex++)
				{
					Setting& Slabsettings = Slablistsettings[slabindex];	//go to the slabindex'th slab setting
					//check group for invalid settings
					try{CheckAngoraGroupSetting(Slabsettings,validsettings);}
					/** check for deprecated options **/
					catch (AngoraInvalidSettingNameException& excp)
					{
						if (excp.getSettingName()=="start_position")
						{
							throw AngoraSettingDeprecatedException(Slabsettings["start_position"],"Use min_coord instead. (Conversion hint: if start_position is a, min_coord is a-1)");
						}
						else if (excp.getSettingName()=="end_position")
						{
							throw AngoraSettingDeprecatedException(Slabsettings["end_position"],"Use max_coord instead. (Conversion hint: if end_position is a, max_coord is also a)");
						}
						else
						{//rethrow the invalid setting exception
							throw;
						}
					}
					/** deprecated options checked **/

					if (SettingEnabledForGrid(Slabsettings))		//apply only if enabled for this grid
					{
						string tag;
						read_value_from_group<string>(Slabsettings,"tag",tag);

						MaterialId mat_id = Materials[tag];

						//z coordinates of the lower and upper interfaces w.r.t. the origin point (can be non-integer)
						double min_coord;
						try{
							read_value_from_group<double>(Slabsettings,"min_coord",min_coord);
							//the coordinate is w.r.t. the grid origin
							min_coord += OriginZ-1; //actual coordinate is 1 less than the grid index (which starts from 1)
						}
						catch (AngoraInvalidSettingTypeException& exc)
						{//it could be string
							string min_coord_str;
							try{read_value_from_group<string>(Slabsettings,"min_coord",min_coord_str);}
							catch (AngoraInvalidSettingTypeException& exc)
							{//if it is neither numeric nor string, the exception is caught here
								throw AngoraInvalidSettingTypeException(Slabsettings["min_coord"],"should be a numeric value or the string \"min\"");
							}
							if (min_coord_str=="min")
							{
								min_coord = 0;
							}
							else
							{
								throw AngoraInvalidSettingValueException(Slabsettings["min_coord"],"only allowable string value is \"min\"");
							}
						}
						double max_coord;
						try{
							read_value_from_group<double>(Slabsettings,"max_coord",max_coord);
							//the coordinate is w.r.t. the grid origin
							max_coord += OriginZ-1; //actual coordinate is 1 less than the grid index (which starts from 1)
						}
						catch (AngoraInvalidSettingTypeException& exc)
						{//it could be string
							string max_coord_str;
							try{read_value_from_group<string>(Slabsettings,"max_coord",max_coord_str);}
							catch (AngoraInvalidSettingTypeException& exc)
							{//if it is neither numeric nor string, the exception is caught here
								throw AngoraInvalidSettingTypeException(Slabsettings["max_coord"],"should be a numeric value or the string \"max\"");
							}
							if (max_coord_str=="max")
							{
								max_coord = NCELLS_Z+2*NPML;//not the cell index anymore: it's the z position
							}
							else
							{
								throw AngoraInvalidSettingValueException(Slabsettings["max_coord"],"only allowable string value is \"max\"");
							}
						}

						//If not in check mode, place the slab in grid
						if (!check_mode)
						{
							PlaceSlab(mat_id,min_coord,max_coord);
						}
					}
				}
			}
/***********************************/
/** REMOVE THE ABOVE OPTION LATER **/
/***********************************/
/***************************************/
/** REMOVE THE FOLLOWING OPTION LATER **/
/***************************************/
			else if (geometry_setting_name=="PECSlabs")
			{
				//redundant, but robust and shorthand
				Setting& PECSlablistsettings = read_list_from_group(Geometrysettings,geometry_setting_name);

				int num_of_PECslabs = PECSlablistsettings.getLength();
				for (int PECslabindex=0; PECslabindex<num_of_PECslabs; PECslabindex++)
				{
					Setting& PECSlabsettings = PECSlablistsettings[PECslabindex];	//go to the PECslabindex'th PEC slab setting
					//check group for invalid settings
					try{CheckAngoraGroupSetting(PECSlabsettings,validsettings);}
					/** check for deprecated options **/
					catch (AngoraInvalidSettingNameException& excp)
					{
						if (excp.getSettingName()=="start_position")
						{
							throw AngoraSettingDeprecatedException(PECSlabsettings["start_position"],"Use min_coord instead. (Conversion hint: if start_position is a, min_coord is a-1)");
						}
						else if (excp.getSettingName()=="end_position")
						{
							throw AngoraSettingDeprecatedException(PECSlabsettings["end_position"],"Use max_coord instead. (Conversion hint: if end_position is a, max_coord is also a)");
						}
						else
						{//rethrow the invalid setting exception
							throw;
						}
					}
					/** deprecated options checked **/

					if (SettingEnabledForGrid(PECSlabsettings))		//apply only if enabled for this grid
					{
						double min_coord;
						try{
							read_value_from_group<double>(PECSlabsettings,"min_coord",min_coord);
							//the coordinate is w.r.t. the grid origin
							min_coord += OriginZ-1; //actual coordinate is 1 less than the grid index (which starts from 1)
						}
						catch (AngoraInvalidSettingTypeException& exc)
						{//it could be string
							string min_coord_str;
							try{read_value_from_group<string>(PECSlabsettings,"min_coord",min_coord_str);}
							catch (AngoraInvalidSettingTypeException& exc)
							{//if it is neither numeric nor string, the exception is caught here
								throw AngoraInvalidSettingTypeException(PECSlabsettings["min_coord"],"should be a numeric value or the string \"min\"");
							}
							if (min_coord_str=="min")
							{
								min_coord = 0;
							}
							else
							{
								throw AngoraInvalidSettingValueException(PECSlabsettings["min_coord"],"only allowable string value is \"min\"");
							}
						}
						double max_coord;
						try{
							read_value_from_group<double>(PECSlabsettings,"max_coord",max_coord);
							//the coordinate is w.r.t. the grid origin
							max_coord += OriginZ-1; //actual coordinate is 1 less than the grid index (which starts from 1)
						}
						catch (AngoraInvalidSettingTypeException& exc)
						{//it could be string
							string max_coord_str;
							try{read_value_from_group<string>(PECSlabsettings,"max_coord",max_coord_str);}
							catch (AngoraInvalidSettingTypeException& exc)
							{//if it is neither numeric nor string, the exception is caught here
								throw AngoraInvalidSettingTypeException(PECSlabsettings["max_coord"],"should be a numeric value or the string \"max\"");
							}
							if (max_coord_str=="max")
							{
								max_coord = NCELLS_Z+2*NPML;//not the cell index anymore: it's the z position
							}
							else
							{
								throw AngoraInvalidSettingValueException(PECSlabsettings["max_coord"],"only allowable string value is \"max\"");
							}
						}

						//If not in check mode, place the PEC slab in grid
						if (!check_mode)
						{
							PlacePECSlab(min_coord,max_coord);
						}
					}
				}
			}
/***********************************/
/** REMOVE THE ABOVE OPTION LATER **/
/***********************************/
			else if (geometry_setting_name=="GroundPlanes")
			{
				//redundant, but robust and shorthand
				Setting& Groundplanelistsettings = read_list_from_group(Geometrysettings,geometry_setting_name);

				int num_of_groundplanes = Groundplanelistsettings.getLength();
				for (int groundplaneindex=0; groundplaneindex<num_of_groundplanes; groundplaneindex++)
				{
					Setting& Groundplanesettings = Groundplanelistsettings[groundplaneindex];	//go to the groundplaneindex'th ground-plane setting
					//check group for invalid settings
					try{CheckAngoraGroupSetting(Groundplanesettings,validsettings);}
					/** check for deprecated options **/
					catch (AngoraInvalidSettingNameException& excp)
					{
						if (excp.getSettingName()=="position_z")
						{
							throw AngoraSettingDeprecatedException(Groundplanesettings["position_z"],"Use \"coord\" instead.");
						}
						else
						{//rethrow the invalid setting exception
							throw;
						}
					}
					/** deprecated options checked **/

					if (SettingEnabledForGrid(Groundplanesettings))		//apply only if enabled for this grid
					{
						//get ground-plane position
						int coord;
						read_value_from_group<int>(Groundplanesettings,"coord",coord);
						coord += OriginZ-1; //actual coordinate is 1 less than the grid index (which starts from 1)

						//If not in check mode, place the slab in grid
						if (!check_mode)
						{
							PlaceGround(coord);
						}
					}
				}
			}
			else if (geometry_setting_name=="RandomMaterials")
			{
				Setting& Randommaterialsettings = Geometrysettings[geometry_setting_name];
				//check group for invalid settings
				CheckAngoraGroupSetting(Randommaterialsettings,validsettings);

				//treat the settings in order of appearance in the config file
				int num_of_random_material_settings = Randommaterialsettings.getLength();
				for (int random_material_setting_index=0; random_material_setting_index<num_of_random_material_settings; random_material_setting_index++)
				{
					Setting& Randommaterialsetting = Randommaterialsettings[random_material_setting_index]; //get the current setting
					string random_material_setting_name = Randommaterialsetting.getName();       //get the name of the current setting
					//Read Whittle-Matern-correlated random material settings
					if (random_material_setting_name=="WhittleMaternCorrelated")
					{
						//redundant, but robust and shorthand
						Setting& WMCollec = read_list_from_group(Randommaterialsettings,random_material_setting_name);

						int num_of_wms = WMCollec.getLength();
						for (int wmindex=0; wmindex<num_of_wms; wmindex++)
						{
							Setting& WMblocksettings = WMCollec[wmindex];	//go to the wmindex'th random block setting
							//check group for invalid settings
							CheckAngoraGroupSetting(WMblocksettings,validsettings);

							if (SettingEnabledForGrid(WMblocksettings))		//apply only if enabled for this grid
							{
								string constitutive_param_type;
								//read randomness type
								read_value_from_group<string>(WMblocksettings,"constitutive_param_type",constitutive_param_type);

								int random_seed;
								bool random_seed_exists = (WMblocksettings.exists("random_seed"));
								if (random_seed_exists)
								{
									try{
										read_value_from_group<int>(WMblocksettings,"random_seed",random_seed);
									}
									catch (AngoraInvalidSettingTypeException& exc)
									{//it could be string
										string random_seed_str;
										try{read_value_from_group<string>(WMblocksettings,"random_seed",random_seed_str);}
										catch (AngoraInvalidSettingTypeException& exc)
										{//if it is neither numeric nor string, the exception is caught here
											throw AngoraInvalidSettingTypeException(WMblocksettings["random_seed"],"should be a numeric value or the string \"run_index\"");
										}
										if (random_seed_str=="run_index")
										{
											random_seed = GridIndex;
										}
										else
										{
											throw AngoraInvalidSettingValueException(WMblocksettings["random_seed"],"only allowable string value is \"run_index\"");
										}
									}
								}

								double mean,std_dev,corr_len,m;
								read_value_from_group<double>(WMblocksettings,"mean",mean);
								read_value_from_group<double>(WMblocksettings,"std_dev",std_dev);
								read_value_from_group<double>(WMblocksettings,"corr_len",corr_len);
								read_value_from_group<double>(WMblocksettings,"m",m);

								double back_coord,front_coord,left_coord,right_coord,lower_coord,upper_coord;
								read_value_from_group<double>(WMblocksettings,"back_coord",back_coord);
								read_value_from_group<double>(WMblocksettings,"front_coord",front_coord);
								read_value_from_group<double>(WMblocksettings,"left_coord",left_coord);
								read_value_from_group<double>(WMblocksettings,"right_coord",right_coord);
								read_value_from_group<double>(WMblocksettings,"lower_coord",lower_coord);
								read_value_from_group<double>(WMblocksettings,"upper_coord",upper_coord);

								//these are actual coordinates, not cell indices
								back_coord += OriginX-1;
								front_coord += OriginX-1;
								left_coord += OriginY-1;
								right_coord += OriginY-1;
								lower_coord += OriginZ-1;
								upper_coord += OriginZ-1;

								//If not in check mode, add Whittle-Matern-correlated random block to grid
								if (!check_mode)
								{
									/** More flexible later? **/
									if (!random_seed_exists)
									{
										PlaceWMCorrRandomBlock<float>(constitutive_param_type,mean,std_dev,corr_len,m,back_coord,front_coord,left_coord,right_coord,lower_coord,upper_coord);
									}
									else
									{
										PlaceWMCorrRandomBlock<float>(constitutive_param_type,mean,std_dev,corr_len,m,back_coord,front_coord,left_coord,right_coord,lower_coord,upper_coord,random_seed);
									}
									/** More flexible later? **/
								}
							}
						}
					}
					//Read other random material settings
					// else if (random_material_setting_name==...)
				}
			}
			else if (geometry_setting_name=="MaterialsFromFiles")
			{
				//redundant, but robust and shorthand
				Setting& Matfilelistsettings = read_list_from_group(Geometrysettings,geometry_setting_name);

				int num_of_matfiles = Matfilelistsettings.getLength();
				for (int matfileindex=0; matfileindex<num_of_matfiles; matfileindex++)
				{
					Setting& Matfilesettings = Matfilelistsettings[matfileindex];	//go to the matfileindex'th material file setting
					//check group for invalid settings
					CheckAngoraGroupSetting(Matfilesettings,validsettings);

					if (SettingEnabledForGrid(Matfilesettings))		//apply only if enabled for this grid
					{
						string constitutive_param_type;
						//read constitutive parameter type
						read_value_from_group<string>(Matfilesettings,"constitutive_param_type",constitutive_param_type);

						//read material file name
						ostringstream matfilenamestream;
						string MatFileName,MatFileExtension,MatFullFileName;
						read_value_from_group<string>(Matfilesettings,"file_name",MatFileName);
						//if the path to the file is not absolute, prepend the base input path to get full path
						if (!is_absolute_path(MatFileName))
						{
							matfilenamestream << InputDir;
						}
						matfilenamestream << MatFileName;

						bool append_run_index_to_name;
						read_value_from_group<bool>(Matfilesettings,"append_run_index_to_name",append_run_index_to_name);
						if (append_run_index_to_name)
						{
							matfilenamestream << GridIndex;
						}
						if (!read_optional_value_from_group<string>(Matfilesettings,"file_extension",MatFileExtension))
						{
							MatFileExtension = "";
						}
						//append file extension
						if (!MatFileExtension.empty())
						{
							matfilenamestream << ".";
						}
						matfilenamestream << MatFileExtension;
						//get string from string stream
						MatFullFileName = matfilenamestream.str();

						int position_x,position_y,position_z;	// x, y and z positions of the material block
						read_value_from_group<int>(Matfilesettings,"position_x",position_x);
						read_value_from_group<int>(Matfilesettings,"position_y",position_y);
						read_value_from_group<int>(Matfilesettings,"position_z",position_z);

						//the position is relative to the grid origin
						position_x += OriginX;
						position_y += OriginY;
						position_z += OriginZ;

						//read the anchor (reference point for position_y,position_z) for the material block
						string anchor;
						if (!read_optional_value_from_group<string>(Matfilesettings,"anchor",anchor))
						{
							anchor = "center";
						}

						//read whether the material file is recorded in double or float type
						string datatype;
						read_value_from_group<string>(Matfilesettings,"datatype",datatype);
						if ((datatype!="double")&&(datatype!="float"))
						{

							throw AngoraInvalidSettingTypeException(Matfilesettings["datatype"],"should be \"double\" or \"float\"");
						}

						//If not in check mode, add material block to grid
						if (!check_mode)
						{
							int max_new_materials; //maximum number of new materials created to represent the different permittivities
							if (!read_optional_value_from_group<int>(Matfilesettings,"max_new_materials",max_new_materials))
							{//use the default number
								if (datatype=="double")
									PlaceMaterialRegionFromFile<double>(MatFullFileName,position_x,position_y,position_z,anchor,constitutive_param_type);
								else
									PlaceMaterialRegionFromFile<float>(MatFullFileName,position_x,position_y,position_z,anchor,constitutive_param_type);
							}
							else
							{//use the number given in the config file
								if (datatype=="double")
									PlaceMaterialRegionFromFile<double>(MatFullFileName,position_x,position_y,position_z,anchor,constitutive_param_type,max_new_materials);
								else
									PlaceMaterialRegionFromFile<float>(MatFullFileName,position_x,position_y,position_z,anchor,constitutive_param_type,max_new_materials);
							}
						}
					}
				}
			}
			else if (geometry_setting_name=="SurfaceProfilesFromFiles")
			{
				//redundant, but robust and shorthand
				Setting& Surffilelistsettings = read_list_from_group(Geometrysettings,geometry_setting_name);

				int num_of_surffiles = Surffilelistsettings.getLength();
				for (int surffileindex=0; surffileindex<num_of_surffiles; surffileindex++)
				{
					Setting& Surffilesettings = Surffilelistsettings[surffileindex];	//go to the surffileindex'th surface file setting
					//check group for invalid settings
					CheckAngoraGroupSetting(Surffilesettings,validsettings);

					if (SettingEnabledForGrid(Surffilesettings))		//apply only if enabled for this grid
					{
						//read material file name
						ostringstream surffilenamestream;
						string SurfFilePath,SurfFileName,SurfFileExtension,SurfFullFileName;
						//read path
						read_optional_value_from_group<string>(Surffilesettings,"filepath",SurfFilePath);
						//add slash to path if necessary
						add_slash_to_path(SurfFilePath);
						//if path is not absolute, prepend the base input path to get full path
						if (!is_absolute_path(SurfFilePath))
						{
							surffilenamestream << InputDir;
						}
						surffilenamestream << SurfFilePath;

						//read file name
						read_value_from_group<string>(Surffilesettings,"filename",SurfFileName);

						//add filename to string stream
						surffilenamestream << SurfFileName;
						bool append_run_index_to_name;
						read_value_from_group<bool>(Surffilesettings,"append_run_index_to_name",append_run_index_to_name);
						if (append_run_index_to_name)
						{
							surffilenamestream << GridIndex;
						}
						if (!read_optional_value_from_group<string>(Surffilesettings,"fileextension",SurfFileExtension))
						{
							SurfFileExtension = "surf";
						}
						//append file extension
						surffilenamestream << "." << SurfFileExtension;
						//get string from string stream
						SurfFullFileName = surffilenamestream.str();

						//read the filling material
						string material_tag;
						read_value_from_group<string>(Surffilesettings,"material_tag",material_tag);

						MaterialId mat_id = Materials[material_tag];

						int position_x,position_y,position_z;	// x, y and z positions of the dielectric slab
						read_value_from_group<int>(Surffilesettings,"position_x",position_x);
						read_value_from_group<int>(Surffilesettings,"position_y",position_y);
						read_value_from_group<int>(Surffilesettings,"position_z",position_z);

						//the position is relative to the grid origin
						position_x += OriginX;
						position_y += OriginY;
						position_z += OriginZ;

						//read the anchor (reference point for position_y,position_z) for the profiled material block
						string anchor;
						if (!read_optional_value_from_group<string>(Surffilesettings,"anchor",anchor))
						{
							anchor = "center";
						}

						//If not in check mode, add profiled material block to grid
						if (!check_mode)
						{
							PlaceSurfaceProfileFromFile(SurfFullFileName,mat_id,position_x,position_y,position_z,anchor);
						}
					}
				}
			}
			else if (geometry_setting_name=="SurfaceEngravingProfilesFromFiles")
			{
				//redundant, but robust and shorthand
				Setting& Surfengrfilelistsettings = read_list_from_group(Geometrysettings,geometry_setting_name);

				int num_of_surfengrfiles = Surfengrfilelistsettings.getLength();
				for (int surfengrfileindex=0; surfengrfileindex<num_of_surfengrfiles; surfengrfileindex++)
				{
					Setting& Surfengrfilesettings = Surfengrfilelistsettings[surfengrfileindex];	//go to the surfengrfileindex'th surface file setting
					//check group for invalid settings
					CheckAngoraGroupSetting(Surfengrfilesettings,validsettings);

					if (SettingEnabledForGrid(Surfengrfilesettings))	//apply only if enabled for this grid
					{
						//read surface-engraving file name
						ostringstream surfengrfilenamestream;
						string SurfEngrFilePath,SurfEngrFileName,SurfEngrFileExtension,SurfEngrFullFileName;
						//read path
						read_optional_value_from_group<string>(Surfengrfilesettings,"filepath",SurfEngrFilePath);
						//add slash to path if necessary
						add_slash_to_path(SurfEngrFilePath);
						//if path is not absolute, prepend the base input path to get full path
						if (!is_absolute_path(SurfEngrFilePath))
						{
							surfengrfilenamestream << InputDir;
						}
						surfengrfilenamestream << SurfEngrFilePath;

						//read file name
						read_value_from_group<string>(Surfengrfilesettings,"filename",SurfEngrFileName);
						//add filename to string stream
						surfengrfilenamestream << SurfEngrFileName;
						bool append_run_index_to_name;
						read_value_from_group<bool>(Surfengrfilesettings,"append_run_index_to_name",append_run_index_to_name);
						if (append_run_index_to_name)
						{
							surfengrfilenamestream << GridIndex;
						}
						if (!read_optional_value_from_group<string>(Surfengrfilesettings,"fileextension",SurfEngrFileExtension))
						{
							SurfEngrFileExtension = "surf";
						}
						//append file extension
						surfengrfilenamestream << "." << SurfEngrFileExtension;
						//get string from string stream
						SurfEngrFullFileName = surfengrfilenamestream.str();

						//read the filling material
						string material_tag;
						read_value_from_group<string>(Surfengrfilesettings,"material_tag",material_tag);

						MaterialId mat_id = Materials[material_tag];

						int position_x,position_y,position_z;	// x, y and z positions of the dielectric slab
						read_value_from_group<int>(Surfengrfilesettings,"position_x",position_x);
						read_value_from_group<int>(Surfengrfilesettings,"position_y",position_y);
						read_value_from_group<int>(Surfengrfilesettings,"position_z",position_z);

						//the position is relative to the grid origin
						position_x += OriginX;
						position_y += OriginY;
						position_z += OriginZ;

						//read the anchor (reference point for position_y,position_z) for the profiled material block
						string anchor;
						if (!read_optional_value_from_group<string>(Surfengrfilesettings,"anchor",anchor))
						{
							anchor = "center";
						}

						//If not in check mode, add profiled material block to grid
						if (!check_mode)
						{
							PlaceSurfaceEngravingProfileFromFile(SurfEngrFullFileName,mat_id,position_x,position_y,position_z,anchor);
						}
					}
				}
			}
			else if (geometry_setting_name=="PECMasksFromFiles")
			{
				//redundant, but robust and shorthand
				Setting& PECmaskfilelistsettings = read_list_from_group(Geometrysettings,geometry_setting_name);

				int num_of_PECmaskfiles = PECmaskfilelistsettings.getLength();
				for (int PECmaskfileindex=0; PECmaskfileindex<num_of_PECmaskfiles; PECmaskfileindex++)
				{
					Setting& PECmaskfilesettings = PECmaskfilelistsettings[PECmaskfileindex];	//go to the PECmaskfileindex'th mask file setting
					//check group for invalid settings
					CheckAngoraGroupSetting(PECmaskfilesettings,validsettings);

					if (SettingEnabledForGrid(PECmaskfilesettings))		//apply only if enabled for this grid
					{
						//read mask file name
						ostringstream PECmaskfilenamestream;
						string PECMaskFilePath,PECMaskFileName,PECMaskFileExtension,PECMaskFullFileName;
						//read path
						read_optional_value_from_group<string>(PECmaskfilesettings,"filepath",PECMaskFilePath);
						//add slash to path if necessary
						add_slash_to_path(PECMaskFilePath);
						//if path is not absolute, prepend the base input path to get full path
						if (!is_absolute_path(PECMaskFilePath))
						{
							PECmaskfilenamestream << InputDir;
						}
						PECmaskfilenamestream << PECMaskFilePath;

						//read file name
						read_value_from_group<string>(PECmaskfilesettings,"filename",PECMaskFileName);

						//add filename to string stream
						PECmaskfilenamestream << PECMaskFileName;
						bool append_run_index_to_name;
						read_value_from_group<bool>(PECmaskfilesettings,"append_run_index_to_name",append_run_index_to_name);
						if (append_run_index_to_name)
						{
							PECmaskfilenamestream << GridIndex;
						}
						if (!read_optional_value_from_group<string>(PECmaskfilesettings,"fileextension",PECMaskFileExtension))
						{
							PECMaskFileExtension = "geom";
						}
						//append file extension
						PECmaskfilenamestream << "." << PECMaskFileExtension;
						//get string from string stream
						PECMaskFullFileName = PECmaskfilenamestream.str();

						int position_x,position_y,position_z;	// x, y and z positions of the PEC mask
						read_value_from_group<int>(PECmaskfilesettings,"position_x",position_x);
						read_value_from_group<int>(PECmaskfilesettings,"position_y",position_y);
						read_value_from_group<int>(PECmaskfilesettings,"position_z",position_z);

						//the position is relative to the grid origin
						position_x += OriginX;
						position_y += OriginY;
						position_z += OriginZ;

						//read the anchor (reference point for position_y,position_z) for the PEC mask
						string anchor;
						if (!read_optional_value_from_group<string>(PECmaskfilesettings,"anchor",anchor))
						{
							anchor = "center";
						}

						//If not in check mode,add PEC mask to grid
						if (!check_mode)
						{
							PlacePECMaskFromFile(PECMaskFullFileName,position_x,position_y,position_z,anchor);
						}
					}
				}
			}
/***************************************/
/** REMOVE THE FOLLOWING OPTION LATER **/
/***************************************/
			else if (geometry_setting_name=="PECBlocks")
			{
				//redundant, but robust and shorthand
				Setting& PECblocklistsettings = read_list_from_group(Geometrysettings,geometry_setting_name);

				int num_of_PECblocks = PECblocklistsettings.getLength();
				for (int PECblockindex=0; PECblockindex<num_of_PECblocks; PECblockindex++)
				{
					Setting& PECblocksettings = PECblocklistsettings[PECblockindex];	//go to the PECblockindex'th PEC block setting
					//check group for invalid settings
					CheckAngoraGroupSetting(PECblocksettings,validsettings);

					if (SettingEnabledForGrid(PECblocksettings))		//apply only if enabled for this grid
					{
						int BlockBack,BlockFront,BlockLeft,BlockRight,BlockLower,BlockUpper;		// position variables for the PEC block
						read_value_from_group<int>(PECblocksettings,"BlockBack",BlockBack);
						read_value_from_group<int>(PECblocksettings,"BlockFront",BlockFront);
						read_value_from_group<int>(PECblocksettings,"BlockLeft",BlockLeft);
						read_value_from_group<int>(PECblocksettings,"BlockRight",BlockRight);
						read_value_from_group<int>(PECblocksettings,"BlockLower",BlockLower);
						read_value_from_group<int>(PECblocksettings,"BlockUpper",BlockUpper);

						//the position is relative to the grid origin
						BlockBack += OriginX;
						BlockFront += OriginX;
						BlockLeft += OriginY;
						BlockRight += OriginY;
						BlockLower += OriginZ;
						BlockUpper += OriginZ;

						//If not in check mode, place the PEC block in grid
						if (!check_mode)
						{
							PlacePECBlock(BlockBack,BlockFront,BlockLeft,BlockRight,BlockLower,BlockUpper);
						}
					}
				}
			}
/***********************************/
/** REMOVE THE ABOVE OPTION LATER **/
/***********************************/
        }
	}

/*******************************************************/
/******	BELOW PORTION SHOULD EVENTUALLY BE REMOVED *****/
/*******************************************************/
	//analyze layering structure
	analyze_layering();
	//find the minimum/maximum constitutive parameters in the grid
	find_extremal_constitutive_params();
/*******************************************************/
/******	ABOVE PORTION SHOULD EVENTUALLY BE REMOVED *****/
/*******************************************************/
}
