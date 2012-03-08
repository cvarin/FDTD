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

//Includes the routine that reads the recorder definitions

#include "headers.h"

#include "read_recorder.h"

//definition of Crecorder needed
#include "Crecorder.h"

//For file-directory  manipulations
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern bool check_mode;

extern string config_filename;
extern string OutputDir;
extern string RecorderOutputDir,MovieRecorderOutputDir,LineRecorderOutputDir,FieldValueRecorderOutputDir;
extern const string default_RecorderOutputDir;
extern const string default_MovieRecorderOutputDir;
extern const string default_LineRecorderOutputDir;
extern const string default_FieldValueRecorderOutputDir;
extern const string default_movie_filename;
extern const string default_movie_fileextension;
extern const string default_line_filename;
extern const string default_line_fileextension;
extern const string default_fieldvalue_filename;
extern const string default_fieldvalue_fileextension;

extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;

extern int GridIndex;
extern int rank;

extern void add_slash_to_path(string& path);
extern bool is_absolute_path(string& path);
extern int create_path(const string& path);

extern void MPI_exit(const int& exitcode);


void read_recorder(Crecorder &Recorder, const Config& fdtdconfig, const Config& validsettings)
{
	//read the recorder output directory name
	try{read_value_from_group<string>(fdtdconfig.getRoot(),"recorder_output_dir",RecorderOutputDir);}
	catch (AngoraSettingNotFoundException&)
	{
		RecorderOutputDir = default_RecorderOutputDir;
	}
	//add slash to path if necessary
	add_slash_to_path(RecorderOutputDir);
	//if path is not absolute, prepend the base path to get full path
	if (!is_absolute_path(RecorderOutputDir))
	{
		RecorderOutputDir = OutputDir + RecorderOutputDir;	//RecorderOutputDir is relative to the output directory
	}
	//create recorder output directory if it does not exist
	if (!check_mode)
	{
		if (create_path(RecorderOutputDir)<0)
		{
			/** throw exception **/
			if (rank==0) cout << "Could not create path " << RecorderOutputDir << endl;
		}
	}

	//Read recorder settings
	string recorder_setting_path = "Recorder";
	if (fdtdconfig.exists(recorder_setting_path))
	{
		Setting& Recordersettings = fdtdconfig.lookup(recorder_setting_path);
		//check group for invalid settings
		CheckAngoraGroupSetting(Recordersettings,validsettings);

		//Read movie-recorder output directory name
		try{read_value_from_group<string>(Recordersettings,"movie_recorder_output_dir",MovieRecorderOutputDir);}
		catch (AngoraSettingNotFoundException&)
		{
			MovieRecorderOutputDir = default_MovieRecorderOutputDir;
		}
		//add slash to path if necessary
		add_slash_to_path(MovieRecorderOutputDir);
		//if path is not absolute, prepend the base path to get full path
		if (!is_absolute_path(MovieRecorderOutputDir))
		{
			MovieRecorderOutputDir = RecorderOutputDir + MovieRecorderOutputDir;	//MovieRecorderOutputDir is relative to the recorder output directory
		}
		//create movie-recorder output directory if it does not exist
		if (!check_mode)
		{
			if (create_path(MovieRecorderOutputDir)<0)
			{
				/** throw exception **/
				if (rank==0) cout << "Could not create path " << MovieRecorderOutputDir << endl;
			}
		}

		//Read line-recorder output directory name
		try{read_value_from_group<string>(Recordersettings,"line_recorder_output_dir",LineRecorderOutputDir);}
		catch (AngoraSettingNotFoundException&)
		{
			LineRecorderOutputDir = default_LineRecorderOutputDir;
		}
		//add slash to path if necessary
		add_slash_to_path(LineRecorderOutputDir);
		//if path is not absolute, prepend the base path to get full path
		if (!is_absolute_path(LineRecorderOutputDir))
		{
			LineRecorderOutputDir = RecorderOutputDir + LineRecorderOutputDir;	//LineRecorderOutputDir is relative to the recorder output directory
		}
		//create line-recorder output directory if it does not exist
		if (!check_mode)
		{
			if (create_path(LineRecorderOutputDir)<0)
			{
				/** throw exception **/
				if (rank==0) cout << "Could not create path " << LineRecorderOutputDir << endl;
			}
		}

		//Read field-value-recorder output directory name
		try{read_value_from_group<string>(Recordersettings,"field_value_recorder_output_dir",FieldValueRecorderOutputDir);}
		catch (AngoraSettingNotFoundException&)
		{
			FieldValueRecorderOutputDir = default_FieldValueRecorderOutputDir;
		}
		//add slash to path if necessary
		add_slash_to_path(FieldValueRecorderOutputDir);
		//if path is not absolute, prepend the base path to get full path
		if (!is_absolute_path(FieldValueRecorderOutputDir))
		{
			FieldValueRecorderOutputDir = RecorderOutputDir + FieldValueRecorderOutputDir;	//FieldValueRecorderOutputDir is relative to the recorder output directory
		}
		//create field-value-recorder output directory if it does not exist
		if (!check_mode)
		{
			if (create_path(FieldValueRecorderOutputDir)<0)
			{
				/** throw exception **/
				if (rank==0) cout << "Could not create path " << FieldValueRecorderOutputDir << endl;
			}
		}

		//Read movie-recorder settings
		string movie_setting_path = "MovieRecorders";
		if (Recordersettings.exists(movie_setting_path))
		{
			Setting& MovieCollec = read_list_from_group(Recordersettings,movie_setting_path);

			int num_of_movies = MovieCollec.getLength();
			for (int movieindex=0; movieindex<num_of_movies; movieindex++)
			{
				Setting& Moviesettings = MovieCollec[movieindex];	//go to the movieindex'th movie-recorder setting
				//check group for invalid settings
				CheckAngoraGroupSetting(Moviesettings,validsettings);

				if (SettingEnabledForGrid(Moviesettings))		//apply only if enabled for this grid
				{
					//read recorded section
					string recorded_section;
					read_value_from_group<string>(Moviesettings,"recorded_section",recorded_section);
					//read recording position
					/** FIXME : This should be double !! **/
					int recorded_position;
					read_value_from_group<int>(Moviesettings,"recorded_position",recorded_position);
					//read recording position
					string recorded_component;
					read_value_from_group<string>(Moviesettings,"recorded_component",recorded_component);
					//read recording scale
					string recording_scale;
					read_value_from_group<string>(Moviesettings,"recording_scale",recording_scale);
					//read recording type
					string recording_type;
					read_value_from_group<string>(Moviesettings,"recording_type",recording_type);

					//read movie-file name
					ostringstream moviefilenamestream;
					string MovieFilePath;
					//read path
					try{read_value_from_group<string>(Moviesettings,"movie_dir",MovieFilePath);}
					catch (AngoraSettingNotFoundException&)
					{//do nothing if does not exist, but type is checked in the try block
					}
					//add slash to path if necessary
					add_slash_to_path(MovieFilePath);
					//if path is not absolute, prepend the base path to get full path
					if (!is_absolute_path(MovieFilePath))
					{
						MovieFilePath = MovieRecorderOutputDir + MovieFilePath;	//MovieFilePath is relative to the recorder output directory
					}
					//create directory if it does not exist
					if (!check_mode)
					{
						if (create_path(MovieFilePath)<0)
						{
							/** throw exception **/
							if (rank==0) cout << "Could not create path " << MovieFilePath << endl;
						}
					}

					string MovieFileName;
					//read file name
					if (!read_optional_value_from_group<string>(Moviesettings,"movie_file_name",MovieFileName))
					{
						MovieFileName = default_movie_filename;
					}
					string MovieFileExtension;
					//read file extension
					if (!read_optional_value_from_group<string>(Moviesettings,"movie_file_extension",MovieFileExtension))
					{
						MovieFileExtension = default_movie_fileextension;
					}
					//construct full filename
					moviefilenamestream << MovieFilePath << MovieFileName << "_" << recorded_component << "_"
						<< recorded_section << "_" << GridIndex << "_" << movieindex;
					if (MovieFileExtension!="")
					{
						moviefilenamestream << "." << MovieFileExtension;
					}
					//get string from string stream
					string MovieFullFileName;
					MovieFullFileName = moviefilenamestream.str();

					//does the movie recorder ONLY record the geometry properties (NOT the field values at each time step)?
					bool OnlyRecordsGeometry;
					try{read_value_from_group<bool>(Moviesettings,"only_records_material_info",OnlyRecordsGeometry);}
					catch (AngoraSettingNotFoundException&)
					{
						OnlyRecordsGeometry = false;
					}

					//If not in check mode, add movie recorder to the Recorder object
					if (!check_mode)
					{
						Recorder.AddMovieRecorder(recorded_section,recorded_position,recorded_component,recording_scale,recording_type,MovieFullFileName,OnlyRecordsGeometry);
					}
				}
			}
		}

		//Read line-recorder settings
		string line_setting_path = "LineRecorders";
		if (Recordersettings.exists(line_setting_path))
		{
			Setting& LineCollec = read_list_from_group(Recordersettings,line_setting_path);

			int num_of_lines = LineCollec.getLength();
			for (int lineindex=0; lineindex<num_of_lines; lineindex++)
			{
				Setting& Linesettings = LineCollec[lineindex];	//go to the lineindex'th line-recorder setting
				//check group for invalid settings
				CheckAngoraGroupSetting(Linesettings,validsettings);

				if (SettingEnabledForGrid(Linesettings))		//apply only if enabled for this grid
				{
					//read line orientation
					string line_orientation;
					read_value_from_group<string>(Linesettings,"line_orientation",line_orientation);
					//read line positions
					int line_position_x1,line_position_x2;
					read_value_from_group<int>(Linesettings,"line_position_x1",line_position_x1);
					read_value_from_group<int>(Linesettings,"line_position_x2",line_position_x2);
					//read recorded component
					string recorded_component;
					read_value_from_group<string>(Linesettings,"recorded_component",recorded_component);
					//read recording scale
					string recording_scale;
					read_value_from_group<string>(Linesettings,"recording_scale",recording_scale);

					//read path
					string LineFilePath;
					try{read_value_from_group<string>(Linesettings,"line_dir",LineFilePath);}
					catch (AngoraSettingNotFoundException&)
					{//do nothing if does not exist, but type is checked in the try block
					}
					//add slash to path if necessary
					add_slash_to_path(LineFilePath);
					//if path is not absolute, prepend the base path to get full path
					if (!is_absolute_path(LineFilePath))
					{
						LineFilePath = LineRecorderOutputDir + LineFilePath;	//LineFilePath is relative to the recorder output directory
					}
					//create directory if it does not exist
					if (!check_mode)
					{
						if (create_path(LineFilePath)<0)
						{
							/** throw exception **/
							if (rank==0) cout << "Could not create path " << LineFilePath << endl;
						}
					}

					//read file name
					ostringstream linefilenamestream;
					string LineFileName;
					try{read_value_from_group<string>(Linesettings,"line_file_name",LineFileName);}
					catch (AngoraSettingNotFoundException&)
					{
						LineFileName = default_line_filename;
					}
					string LineFileExtension;
					//read file extension
					if (!read_optional_value_from_group<string>(Linesettings,"line_file_extension",LineFileExtension))
					{
						LineFileExtension = default_line_fileextension;
					}

					//construct full filename
					string LineFullFileName;
					if (line_orientation=="x_directed")
					{
						linefilenamestream << LineFilePath << LineFileName << "_" << recorded_component << "_X_" << GridIndex << "_" << lineindex;
					}
					else if (line_orientation=="y_directed")
					{
						linefilenamestream << LineFilePath << LineFileName << "_" << recorded_component << "_Y_" << GridIndex << "_" << lineindex;
					}
					else if (line_orientation=="z_directed")
					{
						linefilenamestream << LineFilePath << LineFileName << "_" << recorded_component << "_Z_" << GridIndex << "_" << lineindex;
					}

					if (LineFileExtension!="")
					{
						linefilenamestream << "." << LineFileExtension;
					}
					LineFullFileName = linefilenamestream.str();

					//If not in check mode, add line recorder to the Recorder object
					if (!check_mode)
					{
						//Add line recorder to the Recorder object
						Recorder.AddLineRecorder(line_orientation,line_position_x1,line_position_x2,recorded_component,recording_scale,LineFullFileName);
					}
				}
			}
		}

		//Read field-value-recorder settings
		string fieldvalue_setting_path = "FieldValueRecorders";
		if (Recordersettings.exists(fieldvalue_setting_path))
		{
			Setting& FieldValueCollec = read_list_from_group(Recordersettings,fieldvalue_setting_path);

			int num_of_fieldvalues = FieldValueCollec.getLength();
			for (int fieldvalueindex=0; fieldvalueindex<num_of_fieldvalues; fieldvalueindex++)
			{
				Setting& Fieldvaluesettings = FieldValueCollec[fieldvalueindex];	//go to the fieldvalueindex'th field-value-recorder setting
				//check group for invalid settings
				CheckAngoraGroupSetting(Fieldvaluesettings,validsettings);

				if (SettingEnabledForGrid(Fieldvaluesettings))		//apply only if enabled for this grid
				{
					//read the recording positions
					int position_x,position_y,position_z;
					read_value_from_group<int>(Fieldvaluesettings,"position_x",position_x);
					read_value_from_group<int>(Fieldvaluesettings,"position_y",position_y);
					read_value_from_group<int>(Fieldvaluesettings,"position_z",position_z);

					//read recorded component
					string recorded_component;
					read_value_from_group<string>(Fieldvaluesettings,"recorded_component",recorded_component);
					//read recording scale
					string recording_scale;
					read_value_from_group<string>(Fieldvaluesettings,"recording_scale",recording_scale);

					//read path
					string FieldValueFilePath;
					try{read_value_from_group<string>(Fieldvaluesettings,"field_value_dir",FieldValueFilePath);}
					catch (AngoraSettingNotFoundException&)
					{//do nothing if does not exist, but type is checked in the try block
					}
					//add slash to path if necessary
					add_slash_to_path(FieldValueFilePath);
					//if path is not absolute, prepend the base path to get full path
					if (!is_absolute_path(FieldValueFilePath))
					{
						FieldValueFilePath = FieldValueRecorderOutputDir + FieldValueFilePath;	//FieldValueFilePath is relative to the recorder output directory
					}
					//create directory if it does not exist
					if (!check_mode)
					{
						if (create_path(FieldValueFilePath)<0)
						{
							/** throw exception **/
							if (rank==0) cout << "Could not create path " << FieldValueFilePath << endl;
						}
					}

					//read the filename
					ostringstream fieldvaluefilenamestream;
					string FieldValueFileName;
					try{read_value_from_group<string>(Fieldvaluesettings,"field_value_file_name",FieldValueFileName);}
					catch (AngoraSettingNotFoundException&)
					{
						FieldValueFileName = default_fieldvalue_filename;
					}
					string FieldFileExtension;
					//read file extension
					if (!read_optional_value_from_group<string>(Fieldvaluesettings,"field_value_file_extension",FieldFileExtension))
					{
						FieldFileExtension = default_fieldvalue_fileextension;
					}

					//construct full filename
					string FieldValueFullFileName;
					fieldvaluefilenamestream << FieldValueFilePath << FieldValueFileName << "_" << recorded_component << "_" << GridIndex << "_" << fieldvalueindex;
					if (FieldFileExtension!="")
					{
						fieldvaluefilenamestream << "." << FieldFileExtension;
					}
					FieldValueFullFileName = fieldvaluefilenamestream.str();

					//If not in check mode, add field-value recorder to the Recorder object
					if (!check_mode)
					{
						//Add field-value recorder to the Recorder object
						Recorder.AddFieldValueRecorder(position_x,position_y,position_z,recorded_component,recording_scale,FieldValueFullFileName);
					}
				}
			}
		}
	}
}
