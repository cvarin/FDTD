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

//Includes the routine that reads the PHASOR-DOMAIN near-field-to-far-field transform (NFFFT_pd) definitions

#include "headers.h"

#include "read_nffft_pd.h"

//definition of Cnffft_pd needed
#include "pd/Cnffft_pd.h"

//definition of TrDataType_pd needed
#include "pd/Ctr_pd.h"

//For file-directory  manipulations
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern bool check_mode;

extern string config_filename;
extern string OutputDir;
extern string PhasorDomainNFFFTOutputDir;
extern const string default_PhasorDomainNFFFTOutputDir;
extern const string default_PhasorDomainNFFFT_filename;

extern int GridIndex;
extern int rank;

extern void add_slash_to_path(string& path);
extern bool is_absolute_path(string& path);
extern int create_path(const string& path);


void read_nffft_pd(Cnffft_pd &NFFFT_pd, const Config& fdtdconfig, const Config& validsettings, const Cpointsources &PointSources, const Ctfsf &TFSF)
{
	//Read NFFFT_pd output directory name
	try{read_value_from_group<string>(fdtdconfig.getRoot(),"PhasorDomainNFFFTOutputDir",PhasorDomainNFFFTOutputDir);}
	catch (AngoraSettingNotFoundException&)
	{
		PhasorDomainNFFFTOutputDir = default_PhasorDomainNFFFTOutputDir;
	}
	//add slash to path if necessary
	add_slash_to_path(PhasorDomainNFFFTOutputDir);
	//if path is not absolute, prepend the base path to get full path
	if (!is_absolute_path(PhasorDomainNFFFTOutputDir))
	{
		PhasorDomainNFFFTOutputDir = OutputDir + PhasorDomainNFFFTOutputDir;	//PhasorDomainNFFFTOutputDir is relative to the output directory
	}
	//create NFFFT output directory if it does not exist
	if (!check_mode)
	{
		if (create_path(PhasorDomainNFFFTOutputDir)<0)
		{
			/** throw exception **/
		}
	}

	//Read phasor-domain NFFFT settings
	string nffft_pd_setting_path = "PhasorDomainNFFFT";
	if (fdtdconfig.exists(nffft_pd_setting_path))
	{
		Setting& PhasorDomainNFFFTsettings = read_list_from_group(fdtdconfig.getRoot(),nffft_pd_setting_path);

		//read phasor-domain transformer data
		int num_of_pd_transformers = PhasorDomainNFFFTsettings.getLength();
		for (int pdtrindex=0; pdtrindex<num_of_pd_transformers; pdtrindex++)
		{
			Setting& PhasorDomainTransformersettings = PhasorDomainNFFFTsettings[pdtrindex];	//go to the pdtrindex'th phasor-domain transformer setting
			//check group for invalid settings
			CheckAngoraGroupSetting(PhasorDomainTransformersettings,validsettings);

			if (SettingEnabledForGrid(PhasorDomainTransformersettings))		//apply only if enabled for this grid
			{
				//Add transformer(s) to the NFFFT_pd object
				int num_of_lambdas;
				read_value_from_group<int>(PhasorDomainTransformersettings,"num_of_lambdas",num_of_lambdas);
				if (num_of_lambdas<=0)
				{
					throw AngoraInvalidSettingValueException(PhasorDomainTransformersettings["num_of_lambdas"], "should be positive");
				}
				Array<double,1> lambda(num_of_lambdas);
				double lambda_min,lambda_max;
				read_value_from_group<double>(PhasorDomainTransformersettings,"lambda_min",lambda_min);
				read_value_from_group<double>(PhasorDomainTransformersettings,"lambda_max",lambda_max);
				if (lambda_min>=lambda_max)
				{
					throw AngoraInvalidSettingValueException(PhasorDomainTransformersettings["lambda_max"],"should be greater than \"lambda_min\"");
				}

				string lambda_spacing_type;
				if (num_of_lambdas>1)
				{
					read_value_from_group<string>(PhasorDomainTransformersettings,"lambda_spacing_type",lambda_spacing_type);
					if ((lambda_spacing_type!="log")&&(lambda_spacing_type!="lambda-linear")&&(lambda_spacing_type!="k-linear"))
					{
						throw AngoraInvalidSettingValueException(PhasorDomainTransformersettings["lambda_spacing_type"], "should be \"log\", \"lambda-linear\", or \"k-linear\"");
					}
				}
				bool do_not_include_first_lambda,do_not_include_last_lambda;
				try{read_value_from_group<bool>(PhasorDomainTransformersettings,"do_not_include_first_lambda",do_not_include_first_lambda);}
				catch (AngoraSettingNotFoundException&)
				{
					do_not_include_first_lambda = false;
				}
				try{read_value_from_group<bool>(PhasorDomainTransformersettings,"do_not_include_last_lambda",do_not_include_last_lambda);}
				catch (AngoraSettingNotFoundException&)
				{
					  do_not_include_last_lambda = false;
				}
				if (lambda.size()==1)
				{
					lambda = lambda_min;
				}
				else if (lambda_spacing_type=="lambda-linear")
				{//place lambda values linearly between lambda_min and lambda_max
					double lambda_start,lambda_step;
					if (do_not_include_first_lambda)
					{
						if (do_not_include_last_lambda)
						{//NEITHER the first nor the last lambda is included
							lambda_step = (lambda_max-lambda_min)/lambda.size();
							lambda_start = lambda_min + lambda_step/2;
						}
						else
						{//first lambda is NOT included, last lambda is included
							lambda_step = (lambda_max-lambda_min)/lambda.size();
							lambda_start = lambda_min + lambda_step;
						}
					}
					else
					{
						if (do_not_include_last_lambda)
						{//first lambda is included, last lambda is NOT included
							lambda_step = (lambda_max-lambda_min)/lambda.size();
							lambda_start = lambda_min;
						}
						else
						{//BOTH the first and the last lambda are included
							lambda_step = (lambda_max-lambda_min)/(lambda.size()-1);
							lambda_start = lambda_min;
						}
					}
					for (int i=0;i<lambda.size();i++)
					{
						//linearly-spaced lambda values
						lambda(i) = lambda_start + i*lambda_step;
					}
				}
				else if (lambda_spacing_type=="k-linear")
				{//place k values linearly between 2pi/lambda_max and 2pi/lambda_min
					double k,k_start,k_step;
					if (do_not_include_first_lambda)
					{
						if (do_not_include_last_lambda)
						{//NEITHER the first nor the last lambda is included
							k_step = (2*M_PI/lambda_min-2*M_PI/lambda_max)/lambda.size();
							k_start = 2*M_PI/lambda_max + k_step/2;
						}
						else
						{//first lambda is NOT included, last lambda is included
							k_step = (2*M_PI/lambda_min-2*M_PI/lambda_max)/lambda.size();
							k_start = 2*M_PI/lambda_max;
						}
					}
					else
					{
						if (do_not_include_last_lambda)
						{//first lambda is included, last lambda is NOT included
							k_step = (2*M_PI/lambda_min-2*M_PI/lambda_max)/lambda.size();
							k_start = 2*M_PI/lambda_max + k_step;
						}
						else
						{//BOTH the first and the last lambda are included
							k_step = (2*M_PI/lambda_min-2*M_PI/lambda_max)/(lambda.size()-1);
							k_start = 2*M_PI/lambda_max;
						}
					}
					for (int i=0;i<lambda.size();i++)
					{
						//linearly-spaced k values
						k = k_start + i*k_step;
						//convert to lambda
						lambda(i) = 2*M_PI/k;
					}
				}
				else if (lambda_spacing_type=="log")
				{//place k values logarithmically between 2pi/lambda_max and 2pi/lambda_min
					double klog,klog_start,klog_step;
					if (do_not_include_first_lambda)
					{
						if (do_not_include_last_lambda)
						{//NEITHER the first nor the last lambda is included
							klog_step = (log10(2*M_PI/lambda_min)-log10(2*M_PI/lambda_max))/lambda.size();
							klog_start = log10(2*M_PI/lambda_max) + klog_step/2;
						}
						else
						{//first lambda is NOT included, last lambda is included
							klog_step = (log10(2*M_PI/lambda_min)-log10(2*M_PI/lambda_max))/lambda.size();
							klog_start = log10(2*M_PI/lambda_max);
						}
					}
					else
					{
						if (do_not_include_last_lambda)
						{//first lambda is included, last lambda is NOT included
							klog_step = (log10(2*M_PI/lambda_min)-log10(2*M_PI/lambda_max))/lambda.size();
							klog_start = log10(2*M_PI/lambda_max) + klog_step;
						}
						else
						{//BOTH the first and the last lambda are included
							klog_step = (log10(2*M_PI/lambda_min)-log10(2*M_PI/lambda_max))/(lambda.size()-1);
							klog_start = log10(2*M_PI/lambda_max);
						}
					}
					for (int i=0;i<lambda.size();i++)
					{
						//logarithmically-spaced k values
						klog = klog_start + i*klog_step;
						//convert to lambda
						lambda(i) = 2*M_PI/pow(10.0,klog);
					}
				}

				//read the observation direction parameters
				int num_of_dirs_1,num_of_dirs_2;
				read_value_from_group<int>(PhasorDomainTransformersettings,"num_of_dirs_1",num_of_dirs_1);
				if (num_of_dirs_1<=0)
				{
					throw AngoraInvalidSettingValueException(PhasorDomainTransformersettings["num_of_dirs_1"], "should be positive");
				}
				read_value_from_group<int>(PhasorDomainTransformersettings,"num_of_dirs_2",num_of_dirs_2);
				if (num_of_dirs_2<=0)
				{
					throw AngoraInvalidSettingValueException(PhasorDomainTransformersettings["num_of_dirs_2"], "should be positive");
				}
				Array<double,1> dir1(num_of_dirs_1);
				Array<double,1> dir2(num_of_dirs_2);
				double dir1_min,dir1_max,dir2_min,dir2_max;
				read_value_from_group<double>(PhasorDomainTransformersettings,"dir1_min",dir1_min);
				read_value_from_group<double>(PhasorDomainTransformersettings,"dir1_max",dir1_max);
				read_value_from_group<double>(PhasorDomainTransformersettings,"dir2_min",dir2_min);
				read_value_from_group<double>(PhasorDomainTransformersettings,"dir2_max",dir2_max);
				double max_s;
				try{
					read_value_from_group<double>(PhasorDomainTransformersettings,"limit_to_s",max_s);
					if ((max_s<0)||(max_s>1))
					{
						throw AngoraInvalidSettingValueException(PhasorDomainTransformersettings["limit_to_s"], "should be between 0 and 1");
					}
				}
				catch (AngoraSettingNotFoundException&)
				{
					max_s = 1.0;
				}

				//read the direction specification
				string direction_spec;
				read_value_from_group<string>(PhasorDomainTransformersettings,"direction_spec",direction_spec);
//				if (!PhasorDomainTransformersettings.lookupValue("direction_spec",direction_spec))
//				{
//					cout << "Error: Direction specification (direction_spec) for phasor-domain near-field-to-far-field transformer (NFFFT) " << pdtrindex << " not found or invalid in config file " << config_filename << "." << endl << endl;
//					exit(-1);
//				}
				if ((direction_spec!="theta-phi")&&(direction_spec!="dircosx-dircosy-upper")&&(direction_spec!="dircosx-dircosy-lower"))
				{
					throw AngoraInvalidSettingValueException(PhasorDomainTransformersettings["direction_spec"], "should be \"theta-phi\". \"dircosx-dircosy-upper\", or \"dircosx-dircosy-lower\"");
				}

				if (dir1.size()==1)
				{
					dir1 = dir1_min;
				}
				else
				{
					for (int i=0;i<dir1.size();i++)
					{
						dir1(i) = dir1_min+i*(dir1_max-dir1_min)/(dir1.size()-1);
					}
				}
				if (dir2.size()==1)
				{
					dir2 = dir2_min;
				}
				else
				{
					for (int i=0;i<dir2.size();i++)
					{
						dir2(i) = dir2_min+i*(dir2_max-dir2_min)/(dir2.size()-1);
					}
				}

				//define the transformer data object
				TrDataType_pd PhasorDomainTrData(lambda,dir1,dir2,max_s,direction_spec);	//call the constructor for each new transformer
				try{
					read_value_from_group<int>(PhasorDomainTransformersettings,"NFFFTMarginBackX",PhasorDomainTrData.NFFFTMarginBackX);
					read_value_from_group<int>(PhasorDomainTransformersettings,"NFFFTMarginFrontX",PhasorDomainTrData.NFFFTMarginFrontX);
					read_value_from_group<int>(PhasorDomainTransformersettings,"NFFFTMarginLeftY",PhasorDomainTrData.NFFFTMarginLeftY);
					read_value_from_group<int>(PhasorDomainTransformersettings,"NFFFTMarginRightY",PhasorDomainTrData.NFFFTMarginRightY);
					read_value_from_group<int>(PhasorDomainTransformersettings,"NFFFTMarginLowerZ",PhasorDomainTrData.NFFFTMarginLowerZ);
					read_value_from_group<int>(PhasorDomainTransformersettings,"NFFFTMarginUpperZ",PhasorDomainTrData.NFFFTMarginUpperZ);
				}
				catch (AngoraSettingNotFoundException&)
				{
					//it's OK if these don't exist, but type exceptions are not caught here
				}

				try{read_value_from_group<bool>(PhasorDomainTransformersettings,"scale_with_wavelength",PhasorDomainTrData.scale_with_wavelength);}
				catch (AngoraSettingNotFoundException&)
				{
				    PhasorDomainTrData.scale_with_wavelength = false;
				}

				try{
					read_value_from_group<double>(PhasorDomainTransformersettings,"NFFFTOriginX",PhasorDomainTrData.NFFFTOriginX);
					//always with respect to the grid origin
					PhasorDomainTrData.NFFFTOriginX += OriginX;
				}
				catch (AngoraSettingNotFoundException&)
				{//do nothing if it does not exist
				}
				try{
					read_value_from_group<double>(PhasorDomainTransformersettings,"NFFFTOriginY",PhasorDomainTrData.NFFFTOriginY);
					//always with respect to the grid origin
					PhasorDomainTrData.NFFFTOriginY += OriginY;
				}
				catch (AngoraSettingNotFoundException&)
				{//do nothing if it does not exist
				}
				try{
					read_value_from_group<double>(PhasorDomainTransformersettings,"NFFFTOriginZ",PhasorDomainTrData.NFFFTOriginZ);
					//always with respect to the grid origin
					PhasorDomainTrData.NFFFTOriginZ += OriginZ;
				}
				catch (AngoraSettingNotFoundException&)
				{//do nothing if it does not exist
				}


				//read far-field file name
				ostringstream farfieldfilenamestream;
				string FarFieldFilePath,FarFieldFileName,FarFieldFullFileName;
				//read path
				try{read_value_from_group<string>(PhasorDomainTransformersettings,"FarFieldFilePath",FarFieldFilePath);}
				catch (AngoraSettingNotFoundException&)
				{//do nothing if does not exist, but type is checked in the try block
				}
				//add slash to path if necessary
				add_slash_to_path(FarFieldFilePath);
				//if path is not absolute, prepend the base path to get full path
				if (!is_absolute_path(FarFieldFilePath))
				{
					FarFieldFilePath = PhasorDomainNFFFTOutputDir + FarFieldFilePath;	//FarFieldFilePath is relative to the recorder output directory
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
				try{read_value_from_group<string>(PhasorDomainTransformersettings,"FarFieldFileName",FarFieldFileName);}
				catch (AngoraSettingNotFoundException&)
				{
					FarFieldFileName = default_PhasorDomainNFFFT_filename;
				}

				//construct full filename
				farfieldfilenamestream << FarFieldFilePath << FarFieldFileName << "_" << GridIndex << "_" << pdtrindex << ".hd5";
				//get string from string stream
				FarFieldFullFileName = farfieldfilenamestream.str();

				//If not in check mode, add the transformer to the Cnffft_pd object
				if (!check_mode)
				{
					PhasorDomainTrData.PointSourcesPtr = &PointSources;//pointer to PointSources object is used in theoretical FF calculations
					PhasorDomainTrData.TFSFPtr = &TFSF;//pointer to TFSF object is used for getting the singular plane-wave components in the far field
					//finally, add the transformer to the Cnffft_pd object
					NFFFT_pd.AddTransformer(PhasorDomainTrData,FarFieldFullFileName);
				}
			}
		}
	}
}
