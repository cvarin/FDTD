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

/** README:
Due to a HDF5 script bug (ver. 1.8.4), each time this file is changed, you should do "rm createcfg.o" to recompile.
**/

//Reads the reference config file and creates the C++ source file for the built-in reference config object

#include "cfg_utils.h"

//GNU autoconf config header
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>

////Config object that holds all the valid settings
//extern const Config& valid_angora_settings();

using namespace std;

namespace{
string reference_config_filename;
};

//Command-line parsing library argp (GNU only)
#include <argp.h>

const int num_of_arguments = 1;

//Get these from config.h (GNU only)
const char *argp_program_version = PACKAGE_VERSION;
const char *argp_program_bug_address = PACKAGE_BUGREPORT;

static char doc[] = "Reads the reference config file REF_CONFIGFILE and creates the C++ source file for the built-in reference config object.";

static char args_doc[] = "REF_CONFIGFILE";

//// List of options
//// http://www.gnu.org/s/libc/manual/html_node/Argp-Option-Vectors.html#Argp-Option-Vectors
//static struct argp_option options[] = {
//       { 0 }
//     };

// Parse a single option
// http://www.gnu.org/s/libc/manual/html_node/Argp-Parser-Functions.html#Argp-Parser-Functions
static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
	switch (key)
		{
		case ARGP_KEY_ARG:
		if (state->arg_num >= num_of_arguments)
			/* Too many arguments. */
			argp_usage (state);
		reference_config_filename = arg;
		break;

		case ARGP_KEY_END:
		if (state->arg_num < num_of_arguments)
			/* Not enough arguments. */
             argp_usage (state);
		break;

		default:
		return ARGP_ERR_UNKNOWN;
		}
	return 0;
}

// The argp parser
// http://www.gnu.org/s/libc/manual/html_node/Argp-Parsers.html#Argp-Parsers
//static struct argp argp = { options, parse_opt, args_doc, doc };
static struct argp argp = { NULL, parse_opt, args_doc, doc };


void write_valid_settings_source_file(const Setting& group_setting, const string& group_setting_name, ofstream& source_file)
{//reads the settings in group_setting, writes them into source_file
	for (int setting_index = 0; setting_index < group_setting.getLength(); setting_index++)
	{
		Setting& sub_setting = group_setting[setting_index];
		if (sub_setting.getType()==Setting::TypeInt)
		{
			source_file << group_setting_name << ".add(\"" << sub_setting.getName() << "\",Setting::TypeInt);" << endl;
		}
		else if (sub_setting.getType()==Setting::TypeInt64)
		{
			source_file << group_setting_name << ".add(\"" << sub_setting.getName() << "\",Setting::TypeInt64);" << endl;
		}
		else if (sub_setting.getType()==Setting::TypeFloat)
		{
			source_file << group_setting_name << ".add(\"" << sub_setting.getName() << "\",Setting::TypeFloat);" << endl;
		}
		else if (sub_setting.getType()==Setting::TypeBoolean)
		{
			source_file << group_setting_name << ".add(\"" << sub_setting.getName() << "\",Setting::TypeBoolean);" << endl;
		}
		else if (sub_setting.getType()==Setting::TypeString)
		{
			source_file << group_setting_name << ".add(\"" << sub_setting.getName() << "\",Setting::TypeString);" << endl;
		}
		else if (sub_setting.getType()==Setting::TypeArray)
		{
			source_file << group_setting_name << ".add(\"" << sub_setting.getName() << "\",Setting::TypeArray);" << endl;
		}
		else if (sub_setting.getType()==Setting::TypeList)
		{
			string list_setting_name = sub_setting.getName();
			source_file << "Setting& " << list_setting_name << "=" << group_setting_name << ".add(\"" << list_setting_name << "\",Setting::TypeList);" << endl;
			//recurse into the subsettings
			for (int sub_setting_index = 0; sub_setting_index < sub_setting.getLength(); sub_setting_index++)
			{
				Setting& sub_sub_setting = sub_setting[sub_setting_index];
				if (!sub_sub_setting.isGroup())
				{
					throw AngoraDeveloperException(reference_config_filename + ": " + sub_sub_setting.getPath() + " can only be an unnamed group.");
				}
				else
				{
					stringstream new_group_setting_name_stream;
					new_group_setting_name_stream << list_setting_name << "_" << sub_setting_index;
					string new_group_setting_name = new_group_setting_name_stream.str();
					source_file << "Setting& " << new_group_setting_name << "=" << list_setting_name << ".add(Setting::TypeGroup);" << endl;
					write_valid_settings_source_file(sub_sub_setting, new_group_setting_name_stream.str(), source_file);
				}
			}
		}
		else if (sub_setting.getType()==Setting::TypeGroup)
		{
			string sub_group_setting_name = sub_setting.getName();
			source_file << "Setting& " << sub_group_setting_name << "=" << group_setting_name << ".add(\"" << sub_group_setting_name << "\",Setting::TypeGroup);" << endl;
			write_valid_settings_source_file(sub_setting, sub_group_setting_name, source_file);
		}
	}
}

int main(int argc, char *argv[])
{
	//Parse the command line using the struct argp
	//http://www.gnu.org/s/libc/manual/html_node/Argp.html#Argp
	argp_parse (&argp, argc, argv, 0, 0, NULL);

	//read the reference config file listing the valid Angora settings
	Config reference_config;
	try{reference_config.readFile(reference_config_filename.c_str());}
	catch (ParseException& e)
	{
		stringstream errstream;
		errstream << reference_config_filename << ": line " << e.getLine() << ": parsing error in reference config file: "  << e.getError();
		throw AngoraDeveloperException(errstream.str());
	}
	catch (FileIOException& e)
	{
		throw AngoraDeveloperException("the reference config file " + reference_config_filename + " not found");
	}

	//create the .cc source file
	ofstream built_in_config_source_file("valid_settings.cc");
	built_in_config_source_file << "Setting& root_setting = validsettings.getRoot();" << endl;
	write_valid_settings_source_file(reference_config.getRoot(),"root_setting",built_in_config_source_file);
	built_in_config_source_file << "return validsettings;" << endl;
	built_in_config_source_file.close();

	return EXIT_SUCCESS;
}
