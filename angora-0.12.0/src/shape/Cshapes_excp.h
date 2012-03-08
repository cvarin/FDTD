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

#ifndef CSHAPES_EXCP_H
#define CSHAPES_EXCP_H

//the exceptions thrown by Cshapes and derived classes

//base Angora exception class
#include "angora_excp.h"


class NamedShapeExistsException: public AngoraException
{// exception raised when the shape with a given name already exists
public:
  NamedShapeExistsException(const string& shape_tag) :_shape_tag(shape_tag){};
  virtual ~NamedShapeExistsException() throw() {};

  virtual const string getError() const
  {//error message
  	ostringstream _msgstr;
	_msgstr << "error: a shape with tag \"" << _shape_tag << "\" already exists";
  	return _msgstr.str();
  }
protected:
 const string _shape_tag;
};

class NamedShapeNotFoundException: public AngoraException
{// exception raised when the shape with a given name does not exist
public:
  NamedShapeNotFoundException(const string& shape_tag) :_shape_tag(shape_tag){};
  virtual ~NamedShapeNotFoundException() throw() {};

  virtual const string getError() const
  {//error message
  	ostringstream _msgstr;
	_msgstr << "error: shape with tag \"" << _shape_tag << "\" not found";
  	return _msgstr.str();
  }
protected:
 const string _shape_tag;
};

#endif // CSHAPES_EXCP_H
