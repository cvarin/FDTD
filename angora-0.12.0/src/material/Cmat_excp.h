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

#ifndef CMAT_EXCP_H
#define CMAT_EXCP_H

//the exceptions thrown by Cmats

//base Angora exception class
#include "angora_excp.h"


class NamedMaterialExistsException: public AngoraException
{// exception raised when the material with a given name already exists
public:
  NamedMaterialExistsException(const string& material_tag) :_material_tag(material_tag){};
  virtual ~NamedMaterialExistsException() throw() {};

  virtual const string getError() const
  {//error message
  	ostringstream _msgstr;
	_msgstr << "error: a material with tag \"" << _material_tag << "\" already exists";
  	return _msgstr.str();
  }
protected:
 const string _material_tag;
};

class NamedMaterialNotFoundException: public AngoraException
{// exception raised when the material with a given name does not exist
public:
  NamedMaterialNotFoundException(const string& material_tag) :_material_tag(material_tag){};
  virtual ~NamedMaterialNotFoundException() throw() {};

  virtual const string getError() const
  {//error message
  	ostringstream _msgstr;
	_msgstr << "error: material with tag \"" << _material_tag << "\" not found";
  	return _msgstr.str();
  }
protected:
 const string _material_tag;
};

#endif // CMAT_EXCP_H
