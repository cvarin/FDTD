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

#ifndef ANGORA_EXCP_H
#define ANGORA_EXCP_H

//Common exception classes in Angora

#include "headers.h"


class AngoraException: public exception
{// base class for all exceptions in Angora
public:
  virtual ~AngoraException() throw() {};

  virtual const string getError() const=0;

  virtual const char* what() const throw()
  {
  	return getError().c_str();
  }
};

class AngoraDeveloperException: public AngoraException
{// exception class for an error caused by a bug in the code
public:
  AngoraDeveloperException(const string& error_msg) :_error_msg(error_msg){};
  virtual ~AngoraDeveloperException() throw() {};

  virtual const string getError() const
  {
  	return "Development error: " + _error_msg  + ": To help improve " + PACKAGE_NAME + ", please send a bug report to " + PACKAGE_BUGREPORT + ".";
  };

  virtual const char* what() const throw()
  {
  	return getError().c_str();
  }
protected:
 const string _error_msg;
};

class AngoraInvalidArgumentException: public AngoraException
{// exception for an invalid argument to a function (for catching exceptions outside functions)
// Holds the invalid argument in string type: The specific conversion is done by the derived class AngoraInvalidArgumentExceptionWithType
public:
  AngoraInvalidArgumentException(const string& func_name, const string& valid_args = "") :_func_name(func_name), _valid_args(valid_args) {};
  virtual ~AngoraInvalidArgumentException() throw() {};

  virtual const string getError() const
  {//error message
  	ostringstream _msgstr;
	if (_func_name=="")
		_msgstr << "error: invalid argument " << getArgument() << " " << _valid_args;
	else
		_msgstr << "error: invalid argument " << getArgument() << " to " << _func_name << " " << _valid_args;
  	return _msgstr.str();
  };

  virtual const char* what() const throw()
  {
  	return getError().c_str();
  }

  virtual const string getArgument() const = 0;  //returns the invalid argument in string type

  const string getFunctionName() const
  {//name of the function from which exception is raised
    return _func_name;
  }
protected:
 const string _func_name;
 const string _valid_args;
};

template<typename T>
class AngoraInvalidArgumentExceptionWithType: public AngoraInvalidArgumentException
{// derived exception class for invalid arguments of different types (for throwing exceptions inside functions)
public:
  AngoraInvalidArgumentExceptionWithType(const string& func_name, const T& arg, const string& valid_args = "") :AngoraInvalidArgumentException(func_name,valid_args), _arg(arg)
  {};
  virtual ~AngoraInvalidArgumentExceptionWithType() throw() {};

  const string getArgument() const
  {//the string argument that is invalid (converted to string using C++ string streams)
  	ostringstream _argstr;
  	_argstr << _arg;
  	return _argstr.str();
  }
private:
 const T _arg;
};

#endif // ANGORA_EXCP_H
