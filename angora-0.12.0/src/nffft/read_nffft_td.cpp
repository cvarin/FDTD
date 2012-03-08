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

//Includes the routine that reads the TIME_DOMAIN near-field-to-far-field transform (NFFFT) definitions

#include "headers.h"

#include "read_nffft_td.h"

//definition of Cnffft_td needed
#include "td/Cnffft_td.h"

//definition of TrDataType_td needed
#include "td/Ctr_td.h"

//For file-directory  manipulations
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern bool check_mode;

extern string config_filename;
extern string OutputDir;
extern string TimeDomainNFFFTOutputDir;
extern const string default_TimeDomainNFFFTOutputDir;
extern const string default_TimeDomainNFFFT_filename;

extern int OriginX,OriginY,OriginZ;

extern int GridIndex;
extern int rank;

extern void add_slash_to_path(string& path);
extern bool is_absolute_path(string& path);
extern int create_path(const string& path);


void read_nffft_td(Cnffft_td &NFFFT_td, const Config& fdtdconfig, const Config& validsettings, const Cpointsources &PointSources)
{
	//Read NFFFT output directory name
	try{read_value_from_group<string>(fdtdconfig.getRoot(),"TimeDomainNFFFTOutputDir",TimeDomainNFFFTOutputDir);}
	catch (AngoraSettingNotFoundException&)
	{
		TimeDomainNFFFTOutputDir = default_TimeDomainNFFFTOutputDir;
	}
	//add slash to path if necessary
	add_slash_to_path(TimeDomainNFFFTOutputDir);
	//if path is not absolute, prepend the base path to get full path
	if (!is_absolute_path(TimeDomainNFFFTOutputDir))
	{
		TimeDomainNFFFTOutputDir = OutputDir + TimeDomainNFFFTOutputDir;	//TimeDomainNFFFTOutputDir is relative to the output directory
	}
	//create NFFFT output directory if it does not exist
	if (!check_mode)
	{
		if (create_path(TimeDomainNFFFTOutputDir)<0)
		{
			/** throw exception **/
			if (rank==0) cout << "Could not create path " << TimeDomainNFFFTOutputDir << endl;
		}
	}

	//Read time-domain NFFFT settings
	string nffft_td_setting_path = "TimeDomainNFFFT";
	if (fdtdconfig.exists(nffft_td_setting_path))
	{
		Setting& TimeDomainNFFFTsettings = read_list_from_group(fdtdconfig.getRoot(),nffft_td_setting_path);

		//read transformer data
		int num_of_td_transformers = TimeDomainNFFFTsettings.getLength();
		for (int tdtrindex=0; tdtrindex<num_of_td_transformers; tdtrindex++)
		{
			Setting& TimeDomainTransformersettings = TimeDomainNFFFTsettings[tdtrindex];	//go to the tdtrindex'th time-domain transformer setting
			//check group for invalid settings
			CheckAngoraGroupSetting(TimeDomainTransformersettings,validsettings);

			if (SettingEnabledForGrid(TimeDomainTransformersettings))		//apply only if enabled for this grid
			{
				TrDataType_td TimeDomainTrData;

				read_value_from_group<double>(TimeDomainTransformersettings,"THETA",TimeDomainTrData.THETA);
				TimeDomainTrData.THETA *= M_PI/180;	//convert to radians

				read_value_from_group<double>(TimeDomainTransformersettings,"PHI",TimeDomainTrData.PHI);
				TimeDomainTrData.PHI *= M_PI/180;	//convert to radians

				try{
					read_value_from_group<int>(TimeDomainTransformersettings,"NFFFTMarginBackX",TimeDomainTrData.NFFFTMarginBackX);
					read_value_from_group<int>(TimeDomainTransformersettings,"NFFFTMarginFrontX",TimeDomainTrData.NFFFTMarginFrontX);
					read_value_from_group<int>(TimeDomainTransformersettings,"NFFFTMarginLeftY",TimeDomainTrData.NFFFTMarginLeftY);
					read_value_from_group<int>(TimeDomainTransformersettings,"NFFFTMarginRightY",TimeDomainTrData.NFFFTMarginRightY);
					read_value_from_group<int>(TimeDomainTransformersettings,"NFFFTMarginLowerZ",TimeDomainTrData.NFFFTMarginLowerZ);
					read_value_from_group<int>(TimeDomainTransformersettings,"NFFFTMarginUpperZ",TimeDomainTrData.NFFFTMarginUpperZ);
				}
				catch (AngoraSettingNotFoundException&)
				{
					//it's OK if these don't exist, but type exceptions are not caught here
				}

				try{
					read_value_from_group<double>(TimeDomainTransformersettings,"NFFFTOriginX",TimeDomainTrData.NFFFTOriginX);
					//always with respect to the grid origin
					TimeDomainTrData.NFFFTOriginX += OriginX;
				}
				catch (AngoraSettingNotFoundException&)
				{//do nothing if it does not exist
				}
				try{
					read_value_from_group<double>(TimeDomainTransformersettings,"NFFFTOriginY",TimeDomainTrData.NFFFTOriginY);
					//always with respect to the grid origin
					TimeDomainTrData.NFFFTOriginY += OriginY;
				}
				catch (AngoraSettingNotFoundException&)
				{//do nothing if it does not exist
				}
				try{
					read_value_from_group<double>(TimeDomainTransformersettings,"NFFFTOriginZ",TimeDomainTrData.NFFFTOriginZ);
					//always with respect to the grid origin
					TimeDomainTrData.NFFFTOriginZ += OriginZ;
				}
				catch (AngoraSettingNotFoundException&)
				{//do nothing if it does not exist
				}

				//read far-field file name
				ostringstream farfieldfilenamestream;
				string FarFieldFilePath,FarFieldFileName,FarFieldFullFileName;
				//read path
				try{read_value_from_group<string>(TimeDomainTransformersettings,"FarFieldFilePath",FarFieldFilePath);}
				catch (AngoraSettingNotFoundException&)
				{//do nothing if does not exist, but type is checked in the try block
				}
				//add slash to path if necessary
				add_slash_to_path(FarFieldFilePath);
				//if path is not absolute, prepend the base path to get full path
				if (!is_absolute_path(FarFieldFilePath))
				{
					FarFieldFilePath = TimeDomainNFFFTOutputDir + FarFieldFilePath;	//FarFieldFilePath is relative to the recorder output directory
				}
				//create directory if it does not exist
				if (!check_mode)
				{
					if (create_path(FarFieldFilePath)<0)
					{
						/** throw exception **/
						if (rank==0) cout << "Could not create path " << FarFieldFilePath << endl;
					}
				}

				//read file name
				try{read_value_from_group<string>(TimeDomainTransformersettings,"FarFieldFileName",FarFieldFileName);}
				catch (AngoraSettingNotFoundException&)
				{
					FarFieldFileName = default_TimeDomainNFFFT_filename;
				}

				//construct full filename
				farfieldfilenamestream << FarFieldFilePath << FarFieldFileName << "_" << GridIndex << "_" << tdtrindex << ".ff";
				//get string from string stream
				FarFieldFullFileName = farfieldfilenamestream.str();

				//If not in check mode, add the transformer to the Cnffft_td object
				if (!check_mode)
				{
					TimeDomainTrData.PointSourcesPtr = &PointSources;//pointer to PointSources object is used in theoretical FF calculations
					//finally, add the transformer to the Cnffft_td object
					NFFFT_td.AddTransformer(TimeDomainTrData,FarFieldFullFileName);
				}
			}
		}
	}
}
