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

//Initializes the main grid, executes the main simulation.

//Global variables
#include "globals.h"

//include the HDF5 header, if not disabled
#ifndef HDF5_DISABLE
#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif
#endif

//Routines for reading data from the config file
#include "read_global.h"
#include "read_mat.h"
#include "read_geom.h"
#include "read_basedirs.h"
#include "read_multgridinfo.h"
#include "read_shapes.h"
#include "read_wfs.h"
#include "read_tfsf.h"
#include "read_nffft_td.h"
#include "read_nffft_pd.h"
#include "read_imaging.h"
#include "read_recorder.h"
#include "read_pointsources.h"
#include "read_logging_data.h"

Config fdtdconfig; // configuration object
string config_filename; //configuration filename
string ConfigOutputDir; //configuration file output path
const string default_config_path = ""; //default configuration path
const string default_config_filename = "angora.cfg"; //default configuration filename
const string default_ConfigOutputDir = "cfg/"; //default configuration file output path
bool max_field_value_set_in_configfile=false;		//is the max. field value fixed in the config file?
bool dB_accuracy_set_in_configfile=false;		//is the dB accuracy fixed in the config file?

//Global routines
#include "parallel.h"
#include "init.h"
#include "initgeom.h"
#include "update.h"
#include "time_axis.h"

//Object declarations
#include "Cabsorb.h"
#include "Cshapes.h"
#include "Cmats.h"
#include "Cwfs.h"
#include "Cpointsources.h"
#include "Cnffft_td.h"
#include "Cnffft_pd.h"
#include "Cimgs.h"
#include "Cpw.h"
#include "Cfb.h"
#include "Ctfsf.h"
#include "Crecorder.h"
#include "Cestimator.h"

#include <cstring>
//Timing variables
#include <ctime>
time_t simulation_starting_time, simulation_finishing_time;	//timing variables that mark the beginning and end of the whole simulation

//For file-directory  manipulations
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

//For the "basename" function (Unix only)
#include <libgen.h>

//Command-line parsing library argp (GNU only)
#include <argp.h>

//number of non-option command-line arguments
const int num_of_arguments = 1; //only the config file at the moment
//is the check-mode enabled? (default:no)
bool check_mode = false;

//Get these from config.h (GNU only)
const char *argp_program_version = PACKAGE_VERSION;
const char *argp_program_bug_address = PACKAGE_BUGREPORT;

static char doc[] = "Run three-dimensional electromagnetic finite-difference time-domain simulation using the parameters in CONFIGFILE.";

static char args_doc[] = "CONFIGFILE";

// List of options
// http://www.gnu.org/s/libc/manual/html_node/Argp-Option-Vectors.html#Argp-Option-Vectors
static struct argp_option options[] = {
       {"check", 'c', 0, 0, "Check mode (do not run simulation, only check CONFIGFILE for warnings and errors)" },
       {0, 0, "CONFIGFILE", OPTION_ARG_OPTIONAL, 0},
       { 0 }
     };

// Parse a single option
// http://www.gnu.org/s/libc/manual/html_node/Argp-Parser-Functions.html#Argp-Parser-Functions
static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
	ostringstream config_filenamestream;

	switch (key)
		{
		case 'c':
		check_mode = true;
		break;

		case ARGP_KEY_ARG:
		if (state->arg_num >= num_of_arguments)
			/* Too many arguments. */
			argp_usage (state);
		//the first argument is the configuration file name
		config_filename = arg;
		break;

		case ARGP_KEY_END:
		if (state->arg_num < num_of_arguments)
		{/* Not enough arguments. */
// 			argp_usage (state);
			//assign default config filename
			config_filenamestream << currentworkdir << default_config_path << default_config_filename;
			config_filename = config_filenamestream.str();
			if (rank==0)
			{
			    cout << "No config file specified. Using default config file " << config_filename << endl << endl;
			}
		}
		break;

		default:
		return ARGP_ERR_UNKNOWN;
		}
	return 0;
}

// The argp parser
// http://www.gnu.org/s/libc/manual/html_node/Argp-Parsers.html#Argp-Parsers
static struct argp argp = { options, parse_opt, args_doc, doc };

extern void add_slash_to_path(string& path);
extern int create_path(const string& path);

//Config object that holds all the valid settings
extern const Config& valid_angora_settings();


//*******************************************************************
//	MAIN PROGRAM
//*******************************************************************
int main(int argc, char *argv[])
{
#ifndef MPI_DISABLE
 	MPI_Init(&argc,&argv);
 	atexit((void (*)())MPI_Finalize);
#endif
//	//Initialize MPI by declaring the MPI object
//	//Important: This must be declared before everything, especially before those using parallel I/O
//	CMPI GlobalMPI(argc,argv);	//constructor does MPI initialization, destructor does finalization
//								//not declared global since MPI_Finalize() may not be allowed after main() is finished

	try{

	//Assign the group number, rank number and group size
	init_parallel();

	//get the current working directory
	currentworkdir = getcwd(NULL,0);	//GNU specific working directory resolution
	//add slash to path if necessary
	add_slash_to_path(currentworkdir);

	//Parse the command line using the struct argp
	//http://www.gnu.org/s/libc/manual/html_node/Argp.html#Argp
	argp_parse (&argp, argc, argv, 0, 0, NULL);

	//get the package version from the "config.h" file (stored in string format by autoconf in preprocessor macro PACKAGE_VERSION)
#ifndef PACKAGE_VERSION
	cout << "Warning: Package version variable ""PACKAGE_VERSION"" undefined. Assuming version ""0.0.0""." << endl << endl;
	angora_version_major = 0;
	angora_version_minor = 0;
	angora_version_revision = 0;
#else
	//extract the major,minor and revision numbers from the string
	char separator = '.';
	string package_version_string = PACKAGE_VERSION;
	//major version number
	size_t first_sep_pos = package_version_string.find(separator);
	string major_ver = package_version_string.substr(0,first_sep_pos);
	istringstream major_ver_str(major_ver);
	if (!(major_ver_str >> angora_version_major))
	{
		throw AngoraDeveloperException("Invalid major version number (" + major_ver + ") for the Angora package");
	}
	//major version number
	size_t second_sep_pos = package_version_string.find(separator,first_sep_pos+1);
	string minor_ver = package_version_string.substr(first_sep_pos+1,second_sep_pos-first_sep_pos-1);
	istringstream minor_ver_str(minor_ver);
	if (!(minor_ver_str >> angora_version_minor))
	{
		throw AngoraDeveloperException("Invalid minor version number (" + minor_ver + ") for the Angora package");
	}
	//revision number
	string revision_num = package_version_string.substr(second_sep_pos+1);
	istringstream revision_str(revision_num);
	if (!(revision_str >> angora_version_revision))
	{
		throw AngoraDeveloperException("Invalid revision number (" + revision_num + ") for the Angora package");
	}
#endif

	//print message if check mode is enabled
	if ((check_mode)&&(rank==0))
	{
		cout << "Check mode enabled." << endl << endl;
	}

	//read the Config object that holds the valid Angora settings
	const Config& validsettings = valid_angora_settings();

	//Open the configuration file
	if (rank==0){ cout << "Reading from config file " << config_filename << "..." << endl << endl; }

	//read the configuration from the config file
	fdtdconfig.readFile(config_filename.c_str());

	//check the root setting (which is a group) for invalid settings
	try{
	CheckAngoraGroupSetting(fdtdconfig.getRoot(),validsettings);
	}
	/** check for deprecated options **/
	catch (AngoraInvalidSettingNameException& excp)
	{
		if (excp.getSettingName()=="Geometry")
		{
			throw AngoraSettingDeprecatedException(fdtdconfig.getRoot()["Geometry"],"Use \"SimulationSpace\" instead.");
		}
		else
		{//rethrow the invalid setting exception
			throw;
		}
	}
	/** deprecated options checked **/

	//do type conversions automatically
	fdtdconfig.setAutoConvert(true);
	if (rank==0){ cout << "Config file " << config_filename << " read successfully." << endl << endl; }

	//read the base directories
	read_basedirs(fdtdconfig,validsettings);

	//write out the config file, if not disabled
	if (rank==0)
	{
		bool auto_save_cfg;
		if (!fdtdconfig.lookupValue("auto_save_cfg",auto_save_cfg))
		{
			auto_save_cfg = false;
		}
		if ((auto_save_cfg)&&(!check_mode)) //don't write out log file in check mode
		{
			//Read the config output path
			if (!fdtdconfig.lookupValue("cfg_output_dir",ConfigOutputDir))
			{
				ConfigOutputDir = default_ConfigOutputDir;
			}
			//add slash to path if necessary
			add_slash_to_path(ConfigOutputDir);
			ConfigOutputDir = OutputDir + ConfigOutputDir;	//ConfigOutputDir is relative to the output directory
			//create output directory if it does not exist
			if (!check_mode)
			{
				if (create_path(ConfigOutputDir)<0)
				{
					/** throw exception **/
				}
			}

			//the starting time will be appended to the stored log file name, so...
			//Just read the friggin time
			time_t just_read_the_friggin_time = time(NULL);
			tm *struct_just_read_the_friggin_time = localtime(&just_read_the_friggin_time);
			const int time_size=40;
			char *str_just_read_the_friggin_time = new char[time_size];	//time & date in string format
			strftime(str_just_read_the_friggin_time,time_size,"_%m-%d-%y_%I-%M-%S%p",struct_just_read_the_friggin_time);
			//get the basename from config_filename (currently Unix-only):
			//first, copy "config_filename" to a non-const char array
			char* char_config_filename = new char [config_filename.size()+1];
			strcpy(char_config_filename,config_filename.c_str());
			const char* char_config_basename = basename(char_config_filename);
			string config_basename = char_config_basename;	//cast from char* to string
			//append the time to the log file name
			string ConfigOutputFileName = ConfigOutputDir + char_config_basename + str_just_read_the_friggin_time;	//use the default config file name
			//write the configuration to file
			try{fdtdconfig.writeFile(ConfigOutputFileName.c_str());}
			catch (FileIOException& e)
			{
				cerr << "Warning: Could not write configuration to file " << ConfigOutputFileName << "." << endl << endl;
			}
		}
	}

	//read the multiple-grid information (how many grids, which of them are disabled, etc.)
	read_multgridinfo(fdtdconfig,validsettings);  //only for executable version

	//read some global variables
	read_global(fdtdconfig,validsettings);

	//Initialize the grid
	init_grid();  //later: make constructor of global object?

	//begin multiple-grid loop
	for (GridIndex=0;GridIndex<number_of_runs;GridIndex++)
	{
		if ((rank==0)&&(GridIndex==0))
		{
			cout << "----------------------------------------------------------" << endl;
			cout << "Number of runs: " << number_of_runs << endl;
			cout << "----------------------------------------------------------" << endl;
		}
		if (!grid_is_enabled(GridIndex))
		{
			if (rank==0)
			{
				cout << "Run " << GridIndex << " is disabled. Skipping... " << endl << endl;
			}
			continue;	//if grid is not enabled, continue to the next grid
		}

		bool ShowNodeDist;
		if (!fdtdconfig.lookupValue("ShowNodeDist",ShowNodeDist))
		{
			ShowNodeDist = true;
		}
#ifndef MPI_DISABLE
		if ((rank==0)&&ShowNodeDist)
		{
			cout << "Cartesian node distribution in run " << GridIndex << ": "
				<< "[" << nodes_x << "," << nodes_y << "," << nodes_z << "]" << endl;
		}
#endif

		//initialize the geometry
		init_geom();

		///////*************************************************************************************/
		///////*************************************************************************************/
		///////************************        MATERIAL DEFINITIONS          ***********************/
		///////*************************************************************************************/
		///////*************************************************************************************/

		//read the material definitions from the config file
		Cmats Materials;
		read_mat(Materials,fdtdconfig,validsettings);

		///////*********************************************************************************************/
		///////*********************************************************************************************/
		///////************************        GEOMETRIC SHAPE  DEFINITIONS          ***********************/
		///////*********************************************************************************************/
		///////*********************************************************************************************/

		//read the geometric-shape information from the config file
		Cshapes Shapes;
		read_shapes(Shapes,fdtdconfig,validsettings);

		///////************************************************************************************/
		///////************************************************************************************/
		///////************************        OBJECT  DEFINITIONS          ***********************/
		///////************************************************************************************/
		///////************************************************************************************/

		//read geometry into the grid from the config file
		read_geom(fdtdconfig,validsettings,Shapes,Materials);

//		//read objects into the grid from the config file
//		read_obj(fdtdconfig,validsettings,Shapes,Materials);

		///////*********************************************************************************************/
		///////*********************************************************************************************/
		///////************************       ABSORBING-BOUNDARY DEFINITIONS        ************************/
		///////*********************************************************************************************/
		///////*********************************************************************************************/

		//Define the absorbing-boundary object
		Cabsorb Absorb;
		/** More flexible later? **/
		//adds CPMLs on all faces of the grid
		if (!check_mode)
		{
			//first, read the feature size (in grid cells) of the scattering or radiating object
			// (needed for choosing the alpha parameter in the CPML)
			double CPML_feature_size;
			if (!fdtdconfig.lookupValue("CPML_feature_size",CPML_feature_size))
			{
				//choose the maximum of NCELLS_X,NCELLS_Y,NCELLS_Z
				CPML_feature_size = max(max(NCELLS_X,NCELLS_Y),NCELLS_Z);
			}
			Absorb.AddCPMLs(NPML,CPML_feature_size);
		}
		/** More flexible later? **/

		///////*********************************************************************************************/
		///////*********************************************************************************************/
		///////************************          TIME-WAVEFORM DEFINITIONS           ***********************/
		///////*********************************************************************************************/
		///////*********************************************************************************************/

		//read the time-waveform information from the config file
		Cwfs Waveforms;
		read_wfs(Waveforms,fdtdconfig,validsettings);


		///////*********************************************************************************************/
		///////*********************************************************************************************/
		///////************************          POINT SOURCE DEFINITIONS           ************************/
		///////*********************************************************************************************/
		///////*********************************************************************************************/

		//Define the PointSources object (representing a collection of point sources)
		Cpointsources PointSources;
		//read the point-source information from the config file
		read_pointsources(PointSources,fdtdconfig,validsettings,Waveforms);	//reads info from the Waveforms object


		///////*********************************************************************************************/
		///////*********************************************************************************************/
		///////***********************      TF/SF BOUNDARY DEFINITIONS      ********************************/
		///////*********************************************************************************************/
		///////*********************************************************************************************/

		//Define the TFSF object (representing a collection of incident wave sources)
		Ctfsf TFSF;
		//read the TFSF information from the config file
		read_tfsf(TFSF,fdtdconfig,validsettings,Waveforms);	//reads info from the Waveforms object


		///////*********************************************************************************************/
		///////*********************************************************************************************/
		///////**********      NEAR-FIELD-TO-FAR-FIELD TRANSFORMER DEFINITIONS      ************************/
		///////*********************************************************************************************/
		///////*********************************************************************************************/
		//Define the NFFFT_td object (representing a collection of time-domain transformers)
		Cnffft_td NFFFT_td;
		read_nffft_td(NFFFT_td,fdtdconfig,validsettings,PointSources);

		//Define the NFFFT_pd object (representing a collection of phasor-domain transformers)
		Cnffft_pd NFFFT_pd;
		read_nffft_pd(NFFFT_pd,fdtdconfig,validsettings,PointSources,TFSF);		//reads info from the PointSources and TFSF objects

		///////*********************************************************************************************/
		///////*********************************************************************************************/
		///////***************************      OPTICAL-IMAGE DEFINITIONS      *****************************/
		///////*********************************************************************************************/
		///////*********************************************************************************************/
		//Define the Cimgs object (representing a collection of optical images)
		Cimgs OpticalImages;
		read_imaging(OpticalImages,fdtdconfig,validsettings,TFSF);	//reads info from the TFSF object

		///////*********************************************************************************************/
		///////*********************************************************************************************/
		///////***********************           RECORDER DEFINITIONS       ********************************/
		///////*********************************************************************************************/
		///////*********************************************************************************************/
		//Define the Recorder object (representing a collection of movie or field value recorders)
		Crecorder Recorder;
		read_recorder(Recorder,fdtdconfig,validsettings);

		///////*********************************************************************************************/
		///////*********************************************************************************************/
		///////***********************           MAIN TIME LOOP             ********************************/
		///////*********************************************************************************************/
		///////*********************************************************************************************/

		bool print_initial_time_value=false;
		fdtdconfig.lookupValue("print_initial_time_value",print_initial_time_value);
		if ((rank==0)&&(print_initial_time_value))
		{
			cout << "Initial time value: " << get_initial_time_value() << " sec (" << (int)round(get_initial_time_value()/dt) << " time steps)" << endl;
		}

		if (rank==0)
		{
//			cout << endl << "Accuracy in grid " << GridIndex << " is " << accuracy << " dB." << endl;
			cout << endl << "Run " << GridIndex << ":" << endl;
//			cout << endl << "In run " << GridIndex << ":" << endl
//				<< " Number of planar material layers : " << number_of_layers << endl
//				<< " Number of electric dipoles : " << PointSources.NumberOfElectricDipoles() << endl
//				<< " Number of plane waves : " << TFSF.NumberOfPlaneWaves() << endl
//				<< " Number of focused beams : " << TFSF.NumberOfFocusedBeams() << endl
//				<< " Number of time-domain transformers : " << NFFFT_td.NumberOfTransformers() << endl
//				<< " Number of phasor-domain transformers : " << NFFFT_pd.NumberOfTransformers() << endl
//				<< " Number of optical images : " << OpticalImages.NumberOfOpticalImages() << endl << endl;
		}

		if (!check_mode)//don't run simulation or do post processing in check_mode
		{
			//instantiate the estimator object
			Cestimator Estimator;
			bool logging_enabled;
			if (!read_optional_value_from_group(fdtdconfig.getRoot(),"enable_logging",logging_enabled))
			{
				logging_enabled = true;
			}
			if (logging_enabled)
			{
				read_logging_data(Estimator,fdtdconfig);
			}

			for (int n=0;n<NSTEPS;n++)
			{
				//use the current values, not the updated ones
				Recorder.Record();	//record everything that needs to be recorded at this time step
				NFFFT_td.UpdateFarField(n);
				NFFFT_pd.UpdateFarField(n);
				OpticalImages.UpdateFarField(n);

				updateE(n);		//Update E-field in the main grid
				Absorb.UpdateE(); //Update the E-field at the boundaries
				TFSF.CorrectE(n);	//Apply TF/SF corrections to the E-field

				PointSources.ApplySources(n);

				updateH(n);		//Update H-field in the main grid
				Absorb.UpdateH(); //Update the H-field at the boundaries
				TFSF.CorrectH(n);	//Apply TF/SF corrections to the H-field

				exchangeH();	//Exchange necessary H-field information before next E-field update

				if (rank==0)
				{
					cout << "Main simulation : " << ceil(100.0*(n+1)/NSTEPS) << "% completed ("
							<< n+1 << " of " << NSTEPS << ")\r" ;

					if (logging_enabled)
					{
						Estimator.MakeEstimate(n);	//if (sample_size) samples are taken, make estimate and output to the log file
					}
				}
			}

			if (rank==0)
			{
				cout << endl << "Finished main simulation in run " << GridIndex << endl;
			}

			///////*********************************************************************************************/
			///////*********************************************************************************************/
			///////***********************            POST-PROCESSING           ********************************/
			///////*********************************************************************************************/
			///////*********************************************************************************************/

			//construct and write out the far-field information
			//time-domain NFFFT
			NFFFT_td.ConstructFarField();
			//phasor-domain NFFFT
			NFFFT_pd.ConstructFarField();

			//construct the optical images
			OpticalImages.ConstructImages();

			//get the actual finishing time and elapsed time, and write them to log file
			if (rank==0)
			{
				time_t actual_finishing_time,actual_elapsed_time;
				if (logging_enabled)
				{
					Estimator.WriteActualFinishingTime(actual_finishing_time,actual_elapsed_time);
				}
				tm *struct_actual_finishing_time = localtime(&actual_finishing_time);
				const int time_size=40;
				char *str_actual_finishing_time = new char[time_size];	//time & date in string format
				strftime(str_actual_finishing_time,time_size,"%x %I:%M:%S%p",struct_actual_finishing_time);
				// Obtain coordinated universal time:
				tm *struct_actual_elapsed_time = gmtime(&actual_elapsed_time);
				cout << "Finished on " <<  str_actual_finishing_time << endl
					<< "Time elapsed: " << struct_actual_elapsed_time->tm_hour << " hours, "
					<< struct_actual_elapsed_time->tm_min << " minutes and "
					<< struct_actual_elapsed_time->tm_sec << " seconds" << endl ;
				cout << "----------------------------------------------------------" << endl << endl;
			}

#ifndef MPI_DISABLE
			MPI_Barrier(MPI_SubComm);	//wait until all nodes in the subcommunicator are finished
#endif
		}//check_mode
	} //GridIndex loop

	} //try block

	/** catch syntax and IO exceptions while reading the config file **/
	catch (ParseException& e)
	{
		if (rank==0)
		{
			cerr << config_filename << ": line " << e.getLine() << ": parsing error: " << e.getError() << endl;
		}
		MPI_exit(-1);
	}
	catch (FileIOException& e)
	{
		if (rank==0)
		{
			cerr << "IO error: Could not read config file " << config_filename << endl;
		}
		MPI_exit(-1);
	}
	/** catch any config or internal exceptions **/
	catch (AngoraInvalidArgumentException& exc)
	{
		if (rank==0)
		{
			cerr << exc.getError() << endl;
		}
		MPI_exit(-1);
	}
	catch (AngoraSettingException& exc)
	{
		if (rank==0)
		{
			cerr << exc.getError() << endl;
		}
		MPI_exit(-1);
	}
	/** catch exceptions raised by HDF5 operations **/
#ifndef HDF5_DISABLE
	catch( FileIException& error )
	{
		if (rank==0)
		{
			error.printError();
		}
		MPI_exit(-1);
	}
	catch( DataSetIException& error )
	{
		if (rank==0)
		{
			error.printError();
		}
		MPI_exit(-1);
	}
	catch( DataSpaceIException& error )
	{
		if (rank==0)
		{
			error.printError();
		}
		MPI_exit(-1);
	}
	catch( DataTypeIException& error )
	{
		if (rank==0)
		{
			error.printError();
		}
		MPI_exit(-1);
	}
#endif
	/** catch the exceptions due to bugs **/
	catch( AngoraDeveloperException& error )
	{
		if (rank==0)
		{
			cerr << error.getError() << endl;
		}
//		MPI_exit(-1);
		abort();
	}
	// Isn't what() usually printed to cerr for unhandled exceptions?
	/** report unhandled exceptions and rethrow **/
	catch (AngoraException& exc)
	{
		if (rank==0)
		{
			cerr << exc.getError() << endl;
		}
		MPI_exit(-1);
//		throw exc;
	}
	return EXIT_SUCCESS;
}
