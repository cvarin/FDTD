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

#ifndef CWFS_EXCP_H
#define CWFS_EXCP_H

//base Angora exception class
#include "angora_excp.h"


class NamedWaveformExistsException: public AngoraException
{// exception raised when the waveform with a given name already exists
public:
  NamedWaveformExistsException(const string& waveform_tag) :_waveform_tag(waveform_tag){};
  virtual ~NamedWaveformExistsException() throw() {};

  virtual const string getError() const
  {//error message
  	ostringstream _msgstr;
	_msgstr << "error: a waveform with tag \"" << _waveform_tag << "\" already exists";
  	return _msgstr.str();
  }
protected:
 const string _waveform_tag;
};

class NamedWaveformNotFoundException: public AngoraException
{// exception raised when the waveform with a given name does not exist
public:
  NamedWaveformNotFoundException(const string& waveform_tag) :_waveform_tag(waveform_tag){};
  virtual ~NamedWaveformNotFoundException() throw() {};

  virtual const string getError() const
  {//error message
  	ostringstream _msgstr;
	_msgstr << "error: waveform with tag \"" << _waveform_tag << "\" not found";
  	return _msgstr.str();
  }
protected:
 const string _waveform_tag;
};

#endif // CWFS_EXCP_H
