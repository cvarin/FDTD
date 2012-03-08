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

#ifndef CFG_UTILS_H
#define CFG_UTILS_H

//defines some exceptions raised by configuration-reading functions

#include "angora_excp.h"

//use the libconfig library
#include <libconfig.h++>
using namespace libconfig;

#if !(((LIBCONFIG_VER_MAJOR == 1) && (LIBCONFIG_VER_MINOR >= 4)) \
               || (LIBCONFIG_VER_MAJOR > 1))
            /* library version is < 1.4, getSourceFile() does not exist */
extern string config_filename;
#endif



class AngoraSettingException: public AngoraException
{// base class for all config-setting exceptions in Angora
public:
  AngoraSettingException(const unsigned int& line_number, const string& cfg_filename) : _line_number(line_number), _cfg_filename(cfg_filename) {};
  virtual ~AngoraSettingException() throw() {};

  unsigned int getLineNumber() const {return _line_number;}
  string getFileName() const {return _cfg_filename;}

protected:
 const unsigned int _line_number;
 const string _cfg_filename;
};

class AngoraSettingEmptyException: public AngoraSettingException
{// exception for a an empty list setting
public:
#if (((LIBCONFIG_VER_MAJOR == 1) && (LIBCONFIG_VER_MINOR >= 4)) \
               || (LIBCONFIG_VER_MAJOR > 1))
            /* use the getSourceFile() feature in libconfig 1.4 and later */
  AngoraSettingEmptyException(const Setting& empty_setting) :_setting_path(empty_setting.getPath()), AngoraSettingException(empty_setting.getSourceLine(),empty_setting.getSourceFile()) {};
#else
  AngoraSettingEmptyException(const Setting& empty_setting) :_setting_path(empty_setting.getPath()), AngoraSettingException(empty_setting.getSourceLine(),config_filename) {};
#endif
  virtual ~AngoraSettingEmptyException() throw() {};

  virtual const string getError() const
  {
	ostringstream _msgstr;
	_msgstr << getFileName() << ": line " << getLineNumber() << ": error: " << _setting_path << " cannot be empty";
	return _msgstr.str();
  };

protected:
 const string _setting_path;
};

class AngoraSettingDeprecatedException: public AngoraSettingException
{// exception for a a deprecated list setting
public:
#if (((LIBCONFIG_VER_MAJOR == 1) && (LIBCONFIG_VER_MINOR >= 4)) \
               || (LIBCONFIG_VER_MAJOR > 1))
            /* use the getSourceFile() feature in libconfig 1.4 and later */
  AngoraSettingDeprecatedException(const Setting& deprecated_setting, const string& help_string) :_setting_path(deprecated_setting.getPath()), _help_string(help_string), AngoraSettingException(deprecated_setting.getSourceLine(),deprecated_setting.getSourceFile()) {};
#else
  AngoraSettingDeprecatedException(const Setting& deprecated_setting, const string& help_string) :_setting_path(deprecated_setting.getPath()), _help_string(help_string), AngoraSettingException(deprecated_setting.getSourceLine(),config_filename) {};
#endif
  virtual ~AngoraSettingDeprecatedException() throw() {};

  virtual const string getError() const
  {
	ostringstream _msgstr;
	if (_help_string=="")
  	{
  		_msgstr << getFileName() << ": line " << getLineNumber() << ": error: " << _setting_path << " is deprecated.";
  	}
	else
  	{
  		_msgstr << getFileName() << ": line " << getLineNumber() << ": error: " << _setting_path << " is a deprecated setting name. " << _help_string;
  	}
	return _msgstr.str();
  };

protected:
 const string _setting_path;
 const string _help_string;
};

class AngoraSettingNotFoundException: public AngoraSettingException
{// exception for a non-existent setting
public:
#if (((LIBCONFIG_VER_MAJOR == 1) && (LIBCONFIG_VER_MINOR >= 4)) \
               || (LIBCONFIG_VER_MAJOR > 1))
            /* use the getSourceFile() feature in libconfig 1.4 and later */
  AngoraSettingNotFoundException(const string& setting_name, const Setting& parent_setting) :_setting_name(setting_name), AngoraSettingException(parent_setting.getSourceLine(),parent_setting.getSourceFile()) {};
#else
  AngoraSettingNotFoundException(const string& setting_name, const Setting& parent_setting) :_setting_name(setting_name), AngoraSettingException(parent_setting.getSourceLine(),config_filename) {};
#endif
  virtual ~AngoraSettingNotFoundException() throw() {};

  virtual const string getError() const
  {
	ostringstream _msgstr;
	_msgstr << getFileName() << ": line " << getLineNumber() << ": error: " << _setting_name << " not found";
	return _msgstr.str();
  };

protected:
 const string _setting_name;
};

class AngoraInvalidSettingNameException: public AngoraSettingException
{// exception for an unknown setting
public:
#if (((LIBCONFIG_VER_MAJOR == 1) && (LIBCONFIG_VER_MINOR >= 4)) \
               || (LIBCONFIG_VER_MAJOR > 1))
            /* use the getSourceFile() feature in libconfig 1.4 and later */
  AngoraInvalidSettingNameException(const Setting& invalid_setting) :_setting_path(invalid_setting.getPath()), _setting_name(invalid_setting.getName()), AngoraSettingException(invalid_setting.getSourceLine(),invalid_setting.getSourceFile()) {};
#else
  AngoraInvalidSettingNameException(const Setting& invalid_setting) :_setting_path(invalid_setting.getPath()), _setting_name(invalid_setting.getName()), AngoraSettingException(invalid_setting.getSourceLine(),config_filename) {};
#endif
  virtual ~AngoraInvalidSettingNameException() throw() {};

  virtual const string getError() const
  {//error message
	ostringstream _msgstr;
	_msgstr << getFileName() << ": line " << getLineNumber() << ": error: " << _setting_path << " is an invalid setting name";
	return _msgstr.str();
  };

  const string getSettingName() const
  {//name of the (named) setting (this exception wouldn't be thrown had there not been a name)
	return _setting_name;
  };

protected:
 const string _setting_path;
 const string _setting_name;
};

class AngoraInvalidSettingTypeException: public AngoraSettingException
{// exception for an invalid data type for a setting
public:
#if (((LIBCONFIG_VER_MAJOR == 1) && (LIBCONFIG_VER_MINOR >= 4)) \
               || (LIBCONFIG_VER_MAJOR > 1))
            /* use the getSourceFile() feature in libconfig 1.4 and later */
  AngoraInvalidSettingTypeException(const Setting& invalid_setting, const string& allowed_types) :_setting_path(invalid_setting.getPath()), _allowed_types(allowed_types), AngoraSettingException(invalid_setting.getSourceLine(),invalid_setting.getSourceFile()) {};
#else
  AngoraInvalidSettingTypeException(const Setting& invalid_setting, const string& allowed_types) :_setting_path(invalid_setting.getPath()), _allowed_types(allowed_types), AngoraSettingException(invalid_setting.getSourceLine(),config_filename) {};
#endif
  virtual ~AngoraInvalidSettingTypeException() throw() {};

  virtual const string getError() const
  {//error message
	ostringstream _msgstr;
	if (_allowed_types=="")
  	{
  		_msgstr << getFileName() << ": line " << getLineNumber() << ": error: " << _setting_path << " has invalid type";
  	}
	else
  	{
  		_msgstr << getFileName() << ": line " << getLineNumber() << ": error: " << _setting_path << " has invalid type (" << _allowed_types << ")";
  	}
	return _msgstr.str();
  };

protected:
 const string _setting_path;
 const string _allowed_types;
};

class AngoraInvalidSettingValueException: public AngoraSettingException
{// exception for an invalid value for a setting
public:
#if (((LIBCONFIG_VER_MAJOR == 1) && (LIBCONFIG_VER_MINOR >= 4)) \
               || (LIBCONFIG_VER_MAJOR > 1))
            /* use the getSourceFile() feature in libconfig 1.4 and later */
  AngoraInvalidSettingValueException(const Setting& invalid_setting, const string& allowed_values) :_setting_path(invalid_setting.getPath()), _allowed_values(allowed_values), AngoraSettingException(invalid_setting.getSourceLine(),invalid_setting.getSourceFile()) {};
#else
  AngoraInvalidSettingValueException(const Setting& invalid_setting, const string& allowed_values) :_setting_path(invalid_setting.getPath()), _allowed_values(allowed_values), AngoraSettingException(invalid_setting.getSourceLine(),config_filename) {};
#endif
  virtual ~AngoraInvalidSettingValueException() throw() {};

  virtual const string getError() const
  {
	ostringstream _msgstr;
	if (_allowed_values=="")
  	{
  		_msgstr << getFileName() << ": line " << getLineNumber() << ": error: " << _setting_path << " has an invalid value";
  	}
	else
  	{
  		_msgstr << getFileName() << ": line " << getLineNumber() << ": error: " << _setting_path << " has an invalid value (" << _allowed_values << ")";
  	}
	return _msgstr.str();
  };

protected:
 const string _setting_path;
 const string _allowed_values;
};

template<typename SettingType>
void read_value_from_group(const Setting& parent_group, const string& setting_name, SettingType& setting_variable)
{
	try{setting_variable = parent_group[setting_name];}
	catch (SettingNotFoundException& exc)
	{
		throw AngoraSettingNotFoundException(setting_name,parent_group);
	}
	catch (SettingTypeException& exc)
	{
		if (typeid(SettingType)==typeid(double))
		{
			throw AngoraInvalidSettingTypeException(parent_group[setting_name],"should be double");
		}
		if (typeid(SettingType)==typeid(float))
		{
			throw AngoraInvalidSettingTypeException(parent_group[setting_name],"should be float");
		}
		if (typeid(SettingType)==typeid(int))
		{
			throw AngoraInvalidSettingTypeException(parent_group[setting_name],"should be integer");
		}
		if (typeid(SettingType)==typeid(unsigned int))
		{
			throw AngoraInvalidSettingTypeException(parent_group[setting_name],"should be unsigned integer");
		}
		if (typeid(SettingType)==typeid(short int))
		{
			throw AngoraInvalidSettingTypeException(parent_group[setting_name],"should be short integer");
		}
		if (typeid(SettingType)==typeid(bool))
		{
			throw AngoraInvalidSettingTypeException(parent_group[setting_name],"should be boolean");
		}
		if (typeid(SettingType)==typeid(char))
		{
			throw AngoraInvalidSettingTypeException(parent_group[setting_name],"should be char");
		}
	}
}

template<>
void read_value_from_group(const Setting& parent_group, const string& setting_name, string& setting_variable);

template<typename SettingType>
bool read_optional_value_from_group(const Setting& parent_group, const string& setting_name, SettingType& setting_variable)
{
	try{
		read_value_from_group(parent_group,setting_name,setting_variable);
		//No exceptions: return true to signal successful attempt
		return true;
	}
	catch (AngoraSettingNotFoundException& exc)
	{
		//Nothing is done, since the setting is optional
		//However, return false to signal unsuccessful reading attempt.
		return false;
		//Type mismatch exceptions (AngoraInvalidSettingTypeException) are not caught here!
	}
	//if type mismatch exception (AngoraInvalidSettingTypeException) is thrown, the return value is undefined
}

Setting& read_list_from_group(const Setting& parent_group, const string& list_name);
//
////Setting& read_optional_list_from_group(const Setting& parent_group, const string& list_name);
//
Setting& read_array_from_group(const Setting& parent_group, const string& array_name);
//
////Setting& read_optional_array_from_group(const Setting& parent_group, const string& array_name);
//
Setting& read_group_from_group(const Setting& parent_group, const string& group_name);
//
////Setting& read_optional_group_from_group(const Setting& parent_group, const string& group_name);

bool SettingEnabledForGrid(const Setting& mySetting); //determines if the setting "mySetting" is enabled for the current grid

void CheckAngoraGroupSetting(const Setting& mygroup, const Config& validsettings); //check the group "mygroup" for invalid settings

#endif // CFG_UTILS_H
