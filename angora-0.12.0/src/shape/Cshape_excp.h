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

#ifndef CSHAPE_EXCP_H
#define CSHAPE_EXCP_H

//the exceptions thrown by Cshape and derived classes

//base Angora exception class
#include "angora_excp.h"


class AngoraInvalidBoundsException: public AngoraException
{// exception for invalid coordinates for a planar layer
public:
  AngoraInvalidBoundsException(const string& shape, const double& min_coord, const double& max_coord) :_shape(shape), _min_coord(min_coord),_max_coord(max_coord){};
  virtual ~AngoraInvalidBoundsException() throw() {};

  virtual const string getError() const
  {//error message
  	ostringstream _msgstr;
	_msgstr << "error: invalid bounds for " << _shape << ": " << _min_coord << ">" << _max_coord;
  	return _msgstr.str();
  }
protected:
 const string _shape;
 const double _min_coord, _max_coord;
};

#endif // CSHAPE_EXCP_H
